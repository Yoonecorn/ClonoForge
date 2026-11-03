//! Deterministic external sort of full FASTQ records. Scratch runs are versioned,
//! sorted, length-bounded and SHA-256 authenticated; they are never completion markers.
use crate::{
    consensus::Config,
    io::{input, input_timed, Reader},
    model::{invalid, ReadRecord},
    publication::{inject, Failure},
};
use sha2::{Digest, Sha256};
use std::{
    cmp::Ordering,
    fs::{self, File, OpenOptions},
    io::{BufReader, BufWriter, Read, Write},
    path::{Path, PathBuf},
    sync::{
        atomic::{AtomicU64, Ordering as AtomicOrdering},
        Arc,
    },
    time::Instant,
};

#[derive(Clone, Debug, Default)]
pub struct InputTimings {
    pub detail_enabled: bool,
    pub record_read_parse_inclusive_ns: u64,
    pub raw_file_read_ns: u64,
    pub decoded_read_inclusive_ns: u64,
    pub raw_bytes: u64,
    pub decoded_bytes: u64,
    pub compressed: bool,
    pub spill: SpillBreakdown,
    pub read_parse_header_spill_ns: u64,
    pub header_finalize_verify_ns: u64,
    pub duplicate_scan_barcode_spill_ns: u64,
    pub barcode_finalize_verify_ns: u64,
}
#[derive(Clone, Debug, Default)]
pub struct OutputTimings {
    pub detail_enabled: bool,
    pub consensus_group_elapsed_sum_ns: u64,
    pub consensus_groups: u64,
    pub spill: SpillBreakdown,
    pub group_pipeline_ns: u64,
    pub output_finalize_merge_ns: u64,
}
#[derive(Clone, Debug, Default)]
pub struct SpillBreakdown {
    pub initial_run_build_verify_ns: u64,
    pub merge_pass_inclusive_ns: u64,
    pub flush_ns: u64,
    pub sync_ns: u64,
    pub initial_runs: u64,
    pub merge_passes: u64,
}
#[derive(Default)]
struct SpillCounters {
    initial_run_build_verify_ns: AtomicU64,
    merge_pass_inclusive_ns: AtomicU64,
    flush_ns: AtomicU64,
    sync_ns: AtomicU64,
    initial_runs: AtomicU64,
    merge_passes: AtomicU64,
}
fn count_add(a: &AtomicU64, v: u64) {
    let _ = a.fetch_update(AtomicOrdering::Relaxed, AtomicOrdering::Relaxed, |old| {
        Some(old.saturating_add(v))
    });
}
impl SpillCounters {
    fn snapshot(&self) -> SpillBreakdown {
        let load = |v: &AtomicU64| v.load(AtomicOrdering::Relaxed);
        SpillBreakdown {
            initial_run_build_verify_ns: load(&self.initial_run_build_verify_ns),
            merge_pass_inclusive_ns: load(&self.merge_pass_inclusive_ns),
            flush_ns: load(&self.flush_ns),
            sync_ns: load(&self.sync_ns),
            initial_runs: load(&self.initial_runs),
            merge_passes: load(&self.merge_passes),
        }
    }
}
fn timed_ns(t: Instant) -> u64 {
    u64::try_from(t.elapsed().as_nanos()).unwrap_or(u64::MAX)
}
fn paired_profile_enabled() -> bool {
    std::env::var_os("CLONOFORGE_PROFILE").is_some_and(|value| value == "1")
}
fn paired_sort_geometry(stage: &str, records: u64, metrics: &Option<Arc<SpillCounters>>) {
    if let Some(metrics) = metrics {
        let m = metrics.snapshot();
        let _ = writeln!(
            std::io::stderr().lock(),
            "{}",
            serde_json::json!({
                "phase":"paired_sort_geometry", "stage":stage, "records":records,
                "initial_nonempty_runs":m.initial_runs, "merge_generations":m.merge_passes,
                "empty_run_written":records == 0,
                "initial_build_verify_ns":m.initial_run_build_verify_ns,
                "merge_inclusive_ns":m.merge_pass_inclusive_ns,
                "flush_ns":m.flush_ns,"sync_ns":m.sync_ns
            })
        );
    }
}
const MAGIC: &[u8; 8] = b"CFRUN002";
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Order {
    Header,
    Barcode,
    OutputPass,
    OutputFail,
}
impl Order {
    fn byte(self) -> u8 {
        match self {
            Self::Header => 0,
            Self::Barcode => 1,
            Self::OutputPass => 2,
            Self::OutputFail => 3,
        }
    }
}
#[derive(Clone, Debug)]
pub struct Limits {
    pub chunk_bytes: usize,
    pub max_record_bytes: usize,
    pub fan_in: usize,
}
impl Limits {
    pub fn validate(&self) -> std::io::Result<()> {
        if self.fan_in < 2
            || self.chunk_bytes == 0
            || self.max_record_bytes == 0
            || self.max_record_bytes > u32::MAX as usize
        {
            return Err(invalid("external: invalid limits"));
        }
        Ok(())
    }
}
#[derive(Clone, Debug)]
struct Item {
    record: ReadRecord,
    header: String,
    barcode: String,
    first: u64,
}
impl Item {
    fn new(record: ReadRecord, cfg: &Config) -> std::io::Result<Self> {
        if record.phred.is_none() || !record.bases.iter().all(|b| b"ACGTN".contains(b)) {
            return Err(invalid("input: only uppercase ACGTN FASTQ supported"));
        }
        let barcode = record
            .annotations
            .iter()
            .find(|(k, _)| k == &cfg.barcode_field)
            .map(|(_, v)| v.clone())
            .ok_or_else(|| invalid("grouping: missing selected barcode"))?;
        if barcode.is_empty()
            || !barcode
                .bytes()
                .all(|b| b.is_ascii_alphanumeric() || b == b'_')
        {
            return Err(invalid("grouping: unsafe or empty barcode"));
        }
        let header = record.header();
        let first = record.ordinal;
        Ok(Self {
            record,
            header,
            barcode,
            first,
        })
    }
    fn weight(&self) -> std::io::Result<usize> {
        // Includes spare space for keys, owned strings, and Vec/sort slots.
        let n = self
            .header
            .len()
            .checked_mul(40)
            .and_then(|x| x.checked_add(self.barcode.len().checked_mul(2)?))
            .and_then(|x| x.checked_add(self.record.bases.len().checked_mul(4)?))
            .and_then(|x| x.checked_add(4096))
            .ok_or_else(|| invalid("external: record size overflow"))?;
        Ok(n)
    }
}
fn cmp(a: &Item, b: &Item, order: Order) -> Ordering {
    match order {
        Order::Header => a
            .header
            .cmp(&b.header)
            .then(a.record.ordinal.cmp(&b.record.ordinal)),
        Order::Barcode => a
            .barcode
            .cmp(&b.barcode)
            .then(a.record.ordinal.cmp(&b.record.ordinal)),
        Order::OutputPass | Order::OutputFail => a
            .first
            .cmp(&b.first)
            .then(a.record.ordinal.cmp(&b.record.ordinal)),
    }
}
fn name(dir: &Path, order: Order, pass: u64, index: u64) -> PathBuf {
    dir.join(format!("run-{}-{pass:016x}-{index:016x}.bin", order.byte()))
}
struct RunWriter {
    file: BufWriter<File>,
    hash: Sha256,
    expected: u64,
    written: u64,
    failure: Option<Failure>,
    metrics: Option<Arc<SpillCounters>>,
}
impl RunWriter {
    fn new(
        path: &Path,
        order: Order,
        count: u64,
        failure: Option<Failure>,
        metrics: Option<Arc<SpillCounters>>,
    ) -> std::io::Result<Self> {
        let mut w = Self {
            file: BufWriter::new(OpenOptions::new().write(true).create_new(true).open(path)?),
            hash: Sha256::new(),
            expected: count,
            written: 0,
            failure,
            metrics,
        };
        w.bytes(MAGIC)?;
        w.bytes(&[order.byte()])?;
        w.bytes(&count.to_le_bytes())?;
        Ok(w)
    }
    fn bytes(&mut self, b: &[u8]) -> std::io::Result<()> {
        self.file.write_all(b)?;
        self.hash.update(b);
        Ok(())
    }
    fn item(&mut self, x: &Item) -> std::io::Result<()> {
        let h = u32::try_from(x.header.len()).map_err(|_| invalid("external: header too long"))?;
        let n = u32::try_from(x.record.bases.len())
            .map_err(|_| invalid("external: sequence too long"))?;
        self.bytes(&h.to_le_bytes())?;
        inject(self.failure, Failure::SpillWrite)?;
        self.bytes(&n.to_le_bytes())?;
        self.bytes(&x.record.ordinal.to_le_bytes())?;
        self.bytes(&x.first.to_le_bytes())?;
        self.bytes(x.header.as_bytes())?;
        self.bytes(&x.record.bases)?;
        self.bytes(
            x.record
                .phred
                .as_ref()
                .ok_or_else(|| invalid("external: missing quality"))?,
        )?;
        self.written = self
            .written
            .checked_add(1)
            .ok_or_else(|| invalid("external: count overflow"))?;
        Ok(())
    }
    fn finish(mut self) -> std::io::Result<()> {
        if self.written != self.expected {
            return Err(invalid("external: run count mismatch"));
        }
        let digest = self.hash.finalize();
        self.file.write_all(&digest)?;
        inject(self.failure, Failure::SpillFlush)?;
        let start = self.metrics.as_ref().map(|_| Instant::now());
        self.file.flush()?;
        if let (Some(m), Some(t)) = (&self.metrics, start) {
            count_add(&m.flush_ns, timed_ns(t));
        }
        inject(self.failure, Failure::SpillSync)?;
        let start = self.metrics.as_ref().map(|_| Instant::now());
        self.file.get_ref().sync_all()?;
        if let (Some(m), Some(t)) = (&self.metrics, start) {
            count_add(&m.sync_ns, timed_ns(t));
        }
        Ok(())
    }
}
struct RunReader {
    file: BufReader<File>,
    hash: Sha256,
    left: u64,
    expected: u64,
    order: Order,
    previous: Option<(String, u64, u64)>,
    max_record: usize,
    cfg: Config,
    done: bool,
}
impl RunReader {
    fn open(path: &Path, order: Order, max_record: usize, cfg: &Config) -> std::io::Result<Self> {
        let mut s = Self {
            file: BufReader::new(File::open(path)?),
            hash: Sha256::new(),
            left: 0,
            expected: 0,
            order,
            previous: None,
            max_record,
            cfg: cfg.clone(),
            done: false,
        };
        let mut magic = [0; 8];
        s.bytes(&mut magic)?;
        let mut o = [0];
        s.bytes(&mut o)?;
        let mut c = [0; 8];
        s.bytes(&mut c)?;
        if &magic != MAGIC || o[0] != order.byte() {
            return Err(invalid("external: wrong run format or order"));
        }
        s.left = u64::from_le_bytes(c);
        s.expected = s.left;
        Ok(s)
    }
    fn bytes(&mut self, b: &mut [u8]) -> std::io::Result<()> {
        self.file.read_exact(b)?;
        self.hash.update(&*b);
        Ok(())
    }
    fn next(&mut self) -> std::io::Result<Option<Item>> {
        if self.done {
            return Ok(None);
        }
        if self.left == 0 {
            let mut stored = [0; 32];
            self.file.read_exact(&mut stored)?;
            let actual = self.hash.clone().finalize();
            if stored.as_slice() != actual.as_slice() {
                return Err(invalid("external: checksum mismatch"));
            }
            let mut extra = [0];
            if self.file.read(&mut extra)? != 0 {
                return Err(invalid("external: trailing bytes"));
            }
            self.done = true;
            return Ok(None);
        }
        let mut h = [0; 4];
        let mut n = [0; 4];
        let mut o = [0; 8];
        let mut first = [0; 8];
        self.bytes(&mut h)?;
        self.bytes(&mut n)?;
        self.bytes(&mut o)?;
        self.bytes(&mut first)?;
        let h = u32::from_le_bytes(h) as usize;
        let n = u32::from_le_bytes(n) as usize;
        let sz = h
            .checked_add(
                n.checked_mul(2)
                    .ok_or_else(|| invalid("external: length overflow"))?,
            )
            .ok_or_else(|| invalid("external: length overflow"))?;
        if sz > self.max_record || h == 0 || n == 0 {
            return Err(invalid("external: record exceeds limit"));
        }
        let mut hb = vec![0; h];
        let mut bases = vec![0; n];
        let mut q = vec![0; n];
        self.bytes(&mut hb)?;
        self.bytes(&mut bases)?;
        self.bytes(&mut q)?;
        let header = String::from_utf8(hb).map_err(|_| invalid("external: non UTF8 header"))?;
        let r = ReadRecord::from_header(&header, bases, Some(q), u64::from_le_bytes(o))?;
        let mut item = if matches!(self.order, Order::OutputPass | Order::OutputFail) {
            let header = r.header();
            Item {
                first: u64::from_le_bytes(first),
                barcode: String::new(),
                header,
                record: r,
            }
        } else {
            Item::new(r, &self.cfg)?
        };
        item.first = u64::from_le_bytes(first);
        let key = match self.order {
            Order::Header => (item.header.clone(), 0, item.record.ordinal),
            Order::Barcode => (item.barcode.clone(), 0, item.record.ordinal),
            Order::OutputPass | Order::OutputFail => {
                (String::new(), item.first, item.record.ordinal)
            }
        };
        if self.previous.as_ref().is_some_and(|prev| prev >= &key) {
            return Err(invalid("external: unordered or duplicated run records"));
        }
        self.previous = Some(key);
        self.left -= 1;
        Ok(Some(item))
    }
}
fn verify(path: &Path, order: Order, lim: &Limits, cfg: &Config) -> std::io::Result<u64> {
    let mut r = RunReader::open(path, order, lim.max_record_bytes, cfg)?;
    let count = r.expected;
    while r.next()?.is_some() {}
    Ok(count)
}
fn write_run(
    path: &Path,
    order: Order,
    batch: &mut [Item],
    failure: Option<Failure>,
    metrics: Option<Arc<SpillCounters>>,
) -> std::io::Result<()> {
    batch.sort_unstable_by(|a, b| cmp(a, b, order));
    let mut w = RunWriter::new(
        path,
        order,
        u64::try_from(batch.len()).map_err(|_| invalid("external: count overflow"))?,
        failure,
        metrics,
    )?;
    for x in batch {
        w.item(x)?
    }
    w.finish()
}
#[allow(clippy::too_many_arguments)]
fn merge(
    dir: &Path,
    order: Order,
    pass: u64,
    n: u64,
    lim: &Limits,
    cfg: &Config,
    failure: Option<Failure>,
    metrics: Option<Arc<SpillCounters>>,
) -> std::io::Result<(u64, u64)> {
    let next_pass = pass
        .checked_add(1)
        .ok_or_else(|| invalid("external: merge pass overflow"))?;
    let mut i = 0u64;
    let mut out = 0u64;
    while i < n {
        let end = i
            .checked_add(lim.fan_in as u64)
            .ok_or_else(|| invalid("external: merge index overflow"))?
            .min(n);
        let mut readers = Vec::new();
        let mut total = 0u64;
        for j in i..end {
            let p = name(dir, order, pass, j);
            total = total
                .checked_add(verify(&p, order, lim, cfg)?)
                .ok_or_else(|| invalid("external: merge count overflow"))?;
            readers.push(RunReader::open(&p, order, lim.max_record_bytes, cfg)?)
        }
        let mut heads = Vec::new();
        for r in &mut readers {
            heads.push(r.next()?)
        }
        let np = name(dir, order, next_pass, out);
        let mut w = RunWriter::new(&np, order, total, failure, metrics.clone())?;
        loop {
            let best = heads
                .iter()
                .enumerate()
                .filter_map(|(i, x)| x.as_ref().map(|x| (i, x)))
                .min_by(|a, b| cmp(a.1, b.1, order))
                .map(|x| x.0);
            let Some(k) = best else { break };
            w.item(heads[k].as_ref().unwrap())?;
            inject(failure, Failure::MergePartial)?;
            heads[k] = readers[k].next()?;
        }
        w.finish()?;
        verify(&np, order, lim, cfg)?;
        for j in i..end {
            fs::remove_file(name(dir, order, pass, j))?
        }
        i = end;
        out = out
            .checked_add(1)
            .ok_or_else(|| invalid("external: run count overflow"))?;
    }
    Ok((next_pass, out))
}
struct Builder<'a> {
    dir: &'a Path,
    order: Order,
    lim: &'a Limits,
    cfg: &'a Config,
    batch: Vec<Item>,
    weight: usize,
    runs: u64,
    failure: Option<Failure>,
    metrics: Option<Arc<SpillCounters>>,
}
impl<'a> Builder<'a> {
    fn new(
        dir: &'a Path,
        order: Order,
        lim: &'a Limits,
        cfg: &'a Config,
        failure: Option<Failure>,
        metrics: Option<Arc<SpillCounters>>,
    ) -> Self {
        Self {
            dir,
            order,
            lim,
            cfg,
            batch: Vec::new(),
            weight: 0,
            runs: 0,
            failure,
            metrics,
        }
    }
    fn push(&mut self, x: Item) -> std::io::Result<()> {
        let len = x
            .header
            .len()
            .checked_add(
                x.record
                    .bases
                    .len()
                    .checked_mul(2)
                    .ok_or_else(|| invalid("external: length overflow"))?,
            )
            .ok_or_else(|| invalid("external: length overflow"))?;
        let weight = x.weight()?;
        if len > self.lim.max_record_bytes || weight > self.lim.chunk_bytes {
            return Err(invalid("external: record exceeds memory admission"));
        }
        if self
            .weight
            .checked_add(weight)
            .ok_or_else(|| invalid("external: memory accounting overflow"))?
            > self.lim.chunk_bytes
        {
            self.flush()?
        }
        self.weight = self
            .weight
            .checked_add(weight)
            .ok_or_else(|| invalid("external: memory accounting overflow"))?;
        self.batch.push(x);
        Ok(())
    }
    fn flush(&mut self) -> std::io::Result<()> {
        if self.batch.is_empty() {
            return Ok(());
        }
        let p = name(self.dir, self.order, 0, self.runs);
        inject(self.failure, Failure::SpillCreate)?;
        let start = self.metrics.as_ref().map(|_| Instant::now());
        write_run(
            &p,
            self.order,
            &mut self.batch,
            self.failure,
            self.metrics.clone(),
        )?;
        inject(self.failure, Failure::SpillValidate)?;
        verify(&p, self.order, self.lim, self.cfg)?;
        if let (Some(m), Some(t)) = (&self.metrics, start) {
            count_add(&m.initial_run_build_verify_ns, timed_ns(t));
            count_add(&m.initial_runs, 1);
        }
        self.runs = self
            .runs
            .checked_add(1)
            .ok_or_else(|| invalid("external: run count overflow"))?;
        self.batch.clear();
        self.weight = 0;
        Ok(())
    }
    fn finish(mut self) -> std::io::Result<PathBuf> {
        self.flush()?;
        if self.runs == 0 {
            let p = name(self.dir, self.order, 0, 0);
            write_run(&p, self.order, &mut [], self.failure, self.metrics.clone())?;
            return Ok(p);
        }
        let (mut pass, mut n) = (0, self.runs);
        while n > 1 {
            inject(self.failure, Failure::Merge)?;
            let start = self.metrics.as_ref().map(|_| Instant::now());
            (pass, n) = merge(
                self.dir,
                self.order,
                pass,
                n,
                self.lim,
                self.cfg,
                self.failure,
                self.metrics.clone(),
            )?;
            if let (Some(m), Some(t)) = (&self.metrics, start) {
                count_add(&m.merge_pass_inclusive_ns, timed_ns(t));
                count_add(&m.merge_passes, 1);
            }
        }
        Ok(name(self.dir, self.order, pass, 0))
    }
}
/// Scratch directory must already exist and be owned exclusively by this invocation.
/// On error all incomplete state is retained for diagnosis; it is never a success marker.
pub fn sort_input(path: &Path, cfg: &Config, dir: &Path, lim: &Limits) -> std::io::Result<PathBuf> {
    sort_input_with_failure(path, cfg, dir, lim, None)
}
/// Sort a private generated stream directly by barcode. The caller must have
/// independently checked complete-header uniqueness (including excluded reads).
/// This deliberately does not replace the duplicate-checking public input path.
pub fn sort_generated_by_barcode(
    path: &Path,
    cfg: &Config,
    dir: &Path,
    lim: &Limits,
) -> std::io::Result<PathBuf> {
    cfg.validate()?;
    lim.validate()?;
    let metrics = paired_profile_enabled().then(|| Arc::new(SpillCounters::default()));
    let mut count = 0u64;
    let mut src = Reader::fastq_limited(input(path)?, lim.max_record_bytes)?;
    let mut builder = Builder::new(dir, Order::Barcode, lim, cfg, None, metrics.clone());
    while let Some(record) = src.next_record()? {
        builder.push(Item::new(record, cfg)?)?;
        count = count
            .checked_add(1)
            .ok_or_else(|| invalid("external: count overflow"))?;
    }
    let result = builder.finish()?;
    verify(&result, Order::Barcode, lim, cfg)?;
    paired_sort_geometry("generated_barcode", count, &metrics);
    Ok(result)
}

/// Check complete-header uniqueness with one bounded header sort. This is used
/// by the paired identity index, whose constant barcode requires no second sort.
pub fn check_unique_headers(
    path: &Path,
    cfg: &Config,
    dir: &Path,
    lim: &Limits,
) -> std::io::Result<()> {
    cfg.validate()?;
    lim.validate()?;
    let metrics = paired_profile_enabled().then(|| Arc::new(SpillCounters::default()));
    let mut count = 0u64;
    let mut src = Reader::fastq_limited(input(path)?, lim.max_record_bytes)?;
    let mut builder = Builder::new(dir, Order::Header, lim, cfg, None, metrics.clone());
    while let Some(record) = src.next_record()? {
        builder.push(Item::new(record, cfg)?)?;
        count = count
            .checked_add(1)
            .ok_or_else(|| invalid("external: count overflow"))?;
    }
    let result = builder.finish()?;
    verify(&result, Order::Header, lim, cfg)?;
    let mut reader = RunReader::open(&result, Order::Header, lim.max_record_bytes, cfg)?;
    let mut previous = None::<String>;
    while let Some(item) = reader.next()? {
        if previous.as_ref() == Some(&item.header) {
            return Err(invalid("grouping: duplicate complete header"));
        }
        previous = Some(item.header);
    }
    paired_sort_geometry("identity_header", count, &metrics);
    Ok(())
}

/// Retain eligible paired groups from an authenticated barcode run without
/// sorting or FASTQ materialization. At most `max_group_count` records are kept;
/// larger and unequal-length groups are discarded as complete groups. Counts
/// must match the caller's independently completed preparation disposition.
pub fn filter_paired_groups(
    input: &Path,
    cfg: &Config,
    output: &Path,
    lim: &Limits,
    max_group_count: usize,
    expected_qualified: u64,
) -> std::io::Result<PathBuf> {
    filter_paired_groups_impl(
        input,
        cfg,
        output,
        lim,
        max_group_count,
        expected_qualified,
        false,
    )
}
#[cfg(feature = "experimental-large")]
pub fn filter_paired_groups_experimental(
    input: &Path,
    cfg: &Config,
    output: &Path,
    lim: &Limits,
    max_group_count: usize,
    expected_qualified: u64,
) -> std::io::Result<PathBuf> {
    if ![100, 250, 500].contains(&max_group_count) {
        return Err(invalid("paired: unsupported experimental group limit"));
    }
    filter_paired_groups_impl(
        input,
        cfg,
        output,
        lim,
        max_group_count,
        expected_qualified,
        true,
    )
}
fn filter_paired_groups_impl(
    input: &Path,
    cfg: &Config,
    output: &Path,
    lim: &Limits,
    max_group_count: usize,
    expected_qualified: u64,
    experimental: bool,
) -> std::io::Result<PathBuf> {
    cfg.validate()?;
    lim.validate()?;
    if max_group_count == 0 || max_group_count > if experimental { 500 } else { 50 } {
        return Err(invalid("paired: invalid group limit"));
    }
    let input_records = verify(input, Order::Barcode, lim, cfg)?;
    let mut reader = RunReader::open(input, Order::Barcode, lim.max_record_bytes, cfg)?;
    let mut writer = RunWriter::new(output, Order::Barcode, expected_qualified, None, None)?;
    let mut group = Vec::<Item>::new();
    let mut barcode = String::new();
    let mut count = 0u64;
    let mut length = 0usize;
    let mut unequal = false;
    while let Some(item) = reader.next()? {
        if count > 0 && item.barcode != barcode {
            if count <= max_group_count as u64 && !unequal {
                for record in &group {
                    writer.item(record)?;
                }
            }
            group.clear();
            count = 0;
            unequal = false;
        }
        if count == 0 {
            barcode = item.barcode.clone();
            length = item.record.bases.len();
        }
        count = count
            .checked_add(1)
            .ok_or_else(|| invalid("paired: count overflow"))?;
        unequal |= item.record.bases.len() != length;
        if count > max_group_count as u64 || unequal {
            group.clear();
        } else {
            group.push(item);
        }
    }
    if count <= max_group_count as u64 && !unequal {
        for record in &group {
            writer.item(record)?;
        }
    }
    writer.finish()?;
    verify(output, Order::Barcode, lim, cfg)?;
    if paired_profile_enabled() {
        let _ = writeln!(
            std::io::stderr().lock(),
            "{}",
            serde_json::json!({
                "phase":"paired_filter_geometry", "input_records":input_records,
                "output_records":expected_qualified,"input_full_reads":2,
                "output_full_verifications":1,"maximum_buffered_records":max_group_count
            })
        );
    }
    Ok(output.to_owned())
}
pub fn sort_input_with_failure(
    path: &Path,
    cfg: &Config,
    dir: &Path,
    lim: &Limits,
    failure: Option<Failure>,
) -> std::io::Result<PathBuf> {
    sort_input_profiled(path, cfg, dir, lim, failure, None)
}
pub fn sort_input_profiled(
    path: &Path,
    cfg: &Config,
    dir: &Path,
    lim: &Limits,
    failure: Option<Failure>,
    mut timing: Option<&mut InputTimings>,
) -> std::io::Result<PathBuf> {
    cfg.validate()?;
    lim.validate()?;
    inject(failure, Failure::Parse)?;
    let stage = Instant::now();
    let detailed = timing.as_ref().is_some_and(|t| t.detail_enabled);
    let spill_metrics = detailed.then(|| Arc::new(SpillCounters::default()));
    let (source, read_metrics, compressed) = if detailed {
        let (source, metrics, compressed) = input_timed(path)?;
        (source, Some(metrics), compressed)
    } else {
        (input(path)?, None, false)
    };
    let mut src = Reader::fastq_limited(source, lim.max_record_bytes)?;
    let mut b = Builder::new(dir, Order::Header, lim, cfg, failure, spill_metrics.clone());
    let mut parse_ns = 0u64;
    loop {
        let start = detailed.then(Instant::now);
        let next = src.next_record()?;
        if let Some(start) = start {
            parse_ns = parse_ns.saturating_add(timed_ns(start));
        }
        let Some(record) = next else { break };
        b.push(Item::new(record, cfg)?)?;
    }
    if let (Some(t), Some(m)) = (timing.as_deref_mut(), read_metrics) {
        t.record_read_parse_inclusive_ns = parse_ns;
        t.raw_file_read_ns = m.raw_ns.load(AtomicOrdering::Relaxed);
        t.decoded_read_inclusive_ns = m.decoded_ns.load(AtomicOrdering::Relaxed);
        t.raw_bytes = m.raw_bytes.load(AtomicOrdering::Relaxed);
        t.decoded_bytes = m.decoded_bytes.load(AtomicOrdering::Relaxed);
        t.compressed = compressed;
    }
    if let Some(t) = timing.as_deref_mut() {
        t.read_parse_header_spill_ns = timed_ns(stage);
    }
    let stage = Instant::now();
    let header = b.finish()?;
    verify(&header, Order::Header, lim, cfg)?;
    if let Some(t) = timing.as_deref_mut() {
        t.header_finalize_verify_ns = timed_ns(stage);
    }
    let stage = Instant::now();
    let mut reader = RunReader::open(&header, Order::Header, lim.max_record_bytes, cfg)?;
    let mut b = Builder::new(
        dir,
        Order::Barcode,
        lim,
        cfg,
        failure,
        spill_metrics.clone(),
    );
    inject(failure, Failure::DuplicateDetection)?;
    let mut last = None::<String>;
    while let Some(x) = reader.next()? {
        if last.as_ref() == Some(&x.header) {
            return Err(invalid("grouping: duplicate complete header"));
        }
        last = Some(x.header.clone());
        b.push(x)?
    }
    if let Some(t) = timing.as_deref_mut() {
        t.duplicate_scan_barcode_spill_ns = timed_ns(stage);
    }
    let stage = Instant::now();
    let result = b.finish()?;
    verify(&result, Order::Barcode, lim, cfg)?;
    if let Some(t) = timing {
        t.barcode_finalize_verify_ns = timed_ns(stage);
        if let Some(metrics) = spill_metrics {
            t.spill = metrics.snapshot();
        }
    }
    Ok(result)
}
/// Read a completely validated barcode-sorted run. Consumers must still handle
/// iterator errors (e.g. external mutation after validation).
pub fn for_each_record(
    path: &Path,
    cfg: &Config,
    lim: &Limits,
    mut f: impl FnMut(ReadRecord) -> std::io::Result<()>,
) -> std::io::Result<()> {
    verify(path, Order::Barcode, lim, cfg)?;
    let mut r = RunReader::open(path, Order::Barcode, lim.max_record_bytes, cfg)?;
    while let Some(x) = r.next()? {
        f(x.record)?
    }
    Ok(())
}

/// Produce verified, externally ordered pass/fail runs from a barcode-sorted run.
/// The maximum live group weight is checked before retaining each member.
/// Output runs are sorted by the group's first input ordinal, then read ordinal.
pub fn sort_outputs(
    input: &Path,
    cfg: &Config,
    dir: &Path,
    lim: &Limits,
    max_group_bytes: usize,
    workers: usize,
    window: usize,
) -> std::io::Result<(PathBuf, PathBuf)> {
    sort_outputs_configured(
        input,
        cfg,
        dir,
        lim,
        GroupExecution {
            max_group_bytes,
            workers,
            window,
            max_group_count: 50,
            failure: None,
        },
    )
}
#[derive(Clone, Copy)]
pub struct GroupExecution {
    pub max_group_bytes: usize,
    pub workers: usize,
    pub window: usize,
    pub max_group_count: usize,
    pub failure: Option<Failure>,
}
pub fn sort_outputs_configured(
    input: &Path,
    cfg: &Config,
    dir: &Path,
    lim: &Limits,
    exec: GroupExecution,
) -> std::io::Result<(PathBuf, PathBuf)> {
    sort_outputs_profiled(input, cfg, dir, lim, exec, None)
}
pub fn sort_outputs_profiled(
    input: &Path,
    cfg: &Config,
    dir: &Path,
    lim: &Limits,
    exec: GroupExecution,
    timing: Option<&mut OutputTimings>,
) -> std::io::Result<(PathBuf, PathBuf)> {
    sort_outputs_impl(input, cfg, dir, lim, exec, timing, false)
}
#[cfg(feature = "experimental-large")]
pub fn sort_outputs_experimental(
    input: &Path,
    cfg: &Config,
    dir: &Path,
    lim: &Limits,
    exec: GroupExecution,
    timing: Option<&mut OutputTimings>,
) -> std::io::Result<(PathBuf, PathBuf)> {
    if ![100, 250, 500].contains(&exec.max_group_count) {
        return Err(invalid("grouping: unsupported experimental limit"));
    }
    sort_outputs_impl(input, cfg, dir, lim, exec, timing, true)
}
fn sort_outputs_impl(
    input: &Path,
    cfg: &Config,
    dir: &Path,
    lim: &Limits,
    exec: GroupExecution,
    mut timing: Option<&mut OutputTimings>,
    experimental: bool,
) -> std::io::Result<(PathBuf, PathBuf)> {
    let stage = Instant::now();
    let detailed = timing.as_ref().is_some_and(|t| t.detail_enabled);
    let spill_metrics = detailed.then(|| Arc::new(SpillCounters::default()));
    let mut kernel_sum_ns = 0u64;
    let mut kernel_groups = 0u64;
    let GroupExecution {
        max_group_bytes,
        workers,
        window,
        max_group_count,
        failure,
    } = exec;
    if max_group_count == 0 || max_group_count > if experimental { 500 } else { 50 } {
        return Err(invalid("grouping: invalid configured read limit"));
    }
    use crate::{consensus::consensus, ordered_executor::ordered_map_inner};
    struct Groups {
        reader: RunReader,
        pending: Option<Item>,
        max_bytes: usize,
        max_count: usize,
        done: bool,
    }
    impl Iterator for Groups {
        type Item = std::io::Result<Vec<Item>>;
        fn next(&mut self) -> Option<Self::Item> {
            if self.done {
                return None;
            }
            let start = match self.pending.take() {
                Some(x) => x,
                None => match self.reader.next() {
                    Ok(Some(x)) => x,
                    Ok(None) => {
                        self.done = true;
                        return None;
                    }
                    Err(e) => {
                        self.done = true;
                        return Some(Err(e));
                    }
                },
            };
            let barcode = start.barcode.clone();
            let mut result = Vec::new();
            let mut weight = 0usize;
            let mut next = Some(start);
            loop {
                let item = match next.take() {
                    Some(i) => i,
                    None => match self.reader.next() {
                        Ok(Some(x)) => x,
                        Ok(None) => {
                            self.done = true;
                            break;
                        }
                        Err(e) => {
                            self.done = true;
                            return Some(Err(e));
                        }
                    },
                };
                if item.barcode != barcode {
                    self.pending = Some(item);
                    break;
                }
                if result.len() >= self.max_count {
                    self.done = true;
                    return Some(Err(invalid(
                        "grouping: group exceeds configured read limit",
                    )));
                }
                if result
                    .first()
                    .is_some_and(|r: &Item| r.record.bases.len() != item.record.bases.len())
                {
                    self.done = true;
                    return Some(Err(invalid("grouping: unequal read lengths")));
                }
                weight = match item.weight().and_then(|n| {
                    weight
                        .checked_add(n)
                        .ok_or_else(|| invalid("external: group weight overflow"))
                }) {
                    Ok(n) => n,
                    Err(e) => {
                        self.done = true;
                        return Some(Err(e));
                    }
                };
                if weight > self.max_bytes {
                    self.done = true;
                    return Some(Err(invalid("external: group exceeds admitted weight")));
                }
                result.push(item);
            }
            Some(Ok(result))
        }
    }
    cfg.validate()?;
    lim.validate()?;
    verify(input, Order::Barcode, lim, cfg)?;
    let mut output_limits = lim.clone();
    output_limits.max_record_bytes = output_limits
        .max_record_bytes
        .checked_add(16)
        .ok_or_else(|| invalid("external: output record limit overflow"))?;
    output_limits.validate()?;
    let reader = RunReader::open(input, Order::Barcode, lim.max_record_bytes, cfg)?;
    let groups = Groups {
        reader,
        pending: None,
        max_bytes: max_group_bytes,
        max_count: max_group_count,
        done: false,
    };
    let mut pass = Builder::new(
        dir,
        Order::OutputPass,
        &output_limits,
        cfg,
        failure,
        spill_metrics.clone(),
    );
    let mut fail = Builder::new(
        dir,
        Order::OutputFail,
        &output_limits,
        cfg,
        failure,
        spill_metrics.clone(),
    );
    inject(failure, Failure::Group)?;
    #[cfg(feature = "experimental-large")]
    let worker_load = std::env::var_os("CLONOFORGE_WORKER_LOAD")
        .filter(|v| v == "1")
        .map(|_| {
            std::sync::Mutex::new(std::collections::BTreeMap::<String, (u64, u64, u64)>::new())
        });
    ordered_map_inner(
        groups,
        workers,
        window,
        |mut group| {
            inject(failure, Failure::Worker)?;
            let first = group[0].record.ordinal;
            if group.len() >= cfg.min_count {
                let barcode = group[0].barcode.clone();
                let records: Vec<_> = group.iter().map(|x| x.record.clone()).collect();
                let kernel_start = detailed.then(Instant::now);
                let record = if max_group_count > 50 {
                    #[cfg(feature = "experimental-large")]
                    {
                        crate::consensus::experimental_consensus(&records, &barcode, cfg)?
                    }
                    #[cfg(not(feature = "experimental-large"))]
                    {
                        return Err(invalid("experimental group support disabled"));
                    }
                } else {
                    consensus(&records, &barcode, cfg)?
                };
                let kernel_ns = kernel_start.map(timed_ns).unwrap_or(0);
                #[cfg(feature = "experimental-large")]
                if let Some(load) = &worker_load {
                    let mut load = load
                        .lock()
                        .map_err(|_| invalid("worker load lock poisoned"))?;
                    let entry = load
                        .entry(format!("{:?}", std::thread::current().id()))
                        .or_default();
                    entry.0 += 1;
                    entry.1 += records.len() as u64;
                    entry.2 = entry.2.max(records.len() as u64);
                }
                Ok((
                    true,
                    vec![Item {
                        header: record.header(),
                        record,
                        barcode: String::new(),
                        first,
                    }],
                    kernel_ns,
                ))
            } else if cfg.failed {
                for item in &mut group {
                    item.first = first;
                }
                Ok((false, group, 0))
            } else {
                Ok((false, Vec::new(), 0))
            }
        },
        |(is_pass, items, kernel_ns)| {
            if detailed && is_pass {
                kernel_sum_ns = kernel_sum_ns.saturating_add(kernel_ns);
                kernel_groups = kernel_groups.saturating_add(1);
            }
            let target = if is_pass { &mut pass } else { &mut fail };
            for item in items {
                target.push(item)?;
            }
            Ok(())
        },
        if failure == Some(Failure::WorkerEarlyExit) {
            2
        } else {
            0
        },
    )?;
    #[cfg(feature = "experimental-large")]
    if let Some(load) = worker_load {
        let load = load
            .into_inner()
            .map_err(|_| invalid("worker load lock poisoned"))?;
        eprintln!(
            "{}",
            serde_json::json!({"phase":"experimental_worker_load", "workers":load})
        );
    }
    if let Some(t) = timing.as_deref_mut() {
        t.group_pipeline_ns = timed_ns(stage);
        t.consensus_group_elapsed_sum_ns = kernel_sum_ns;
        t.consensus_groups = kernel_groups;
    }
    let stage = Instant::now();
    inject(failure, Failure::OutputOrder)?;
    let result = (pass.finish()?, fail.finish()?);
    if let Some(t) = timing {
        t.output_finalize_merge_ns = timed_ns(stage);
        if let Some(metrics) = spill_metrics {
            t.spill = metrics.snapshot();
        }
    }
    Ok(result)
}

/// Validate a complete result run before invoking any sink callback.
pub fn for_each_output(
    path: &Path,
    pass: bool,
    cfg: &Config,
    lim: &Limits,
    mut f: impl FnMut(ReadRecord) -> std::io::Result<()>,
) -> std::io::Result<()> {
    let order = if pass {
        Order::OutputPass
    } else {
        Order::OutputFail
    };
    let mut output_limits = lim.clone();
    output_limits.max_record_bytes = output_limits
        .max_record_bytes
        .checked_add(16)
        .ok_or_else(|| invalid("external: output record limit overflow"))?;
    output_limits.validate()?;
    verify(path, order, &output_limits, cfg)?;
    let mut reader = RunReader::open(path, order, output_limits.max_record_bytes, cfg)?;
    while let Some(item) = reader.next()? {
        f(item.record)?;
    }
    Ok(())
}
/// Resume an already open parser into external sorting. Timings here exclude the initial in-memory trial.
#[allow(clippy::too_many_arguments)]
pub fn sort_reader_with_seed<R: std::io::BufRead>(
    mut source: Reader<R>,
    groups: Vec<crate::group::Group>,
    pending: ReadRecord,
    cfg: &Config,
    dir: &Path,
    lim: &Limits,
    failure: Option<Failure>,
    mut timing: Option<&mut InputTimings>,
) -> std::io::Result<PathBuf> {
    cfg.validate()?;
    lim.validate()?;
    let detailed = timing.as_ref().is_some_and(|t| t.detail_enabled);
    let metrics = detailed.then(|| Arc::new(SpillCounters::default()));
    let stage = Instant::now();
    let mut b = Builder::new(dir, Order::Header, lim, cfg, failure, metrics.clone());
    for group in groups {
        for record in group.records {
            b.push(Item::new(record, cfg)?)?;
        }
    }
    b.push(Item::new(pending, cfg)?)?;
    let mut parse_ns = 0u64;
    loop {
        let start = detailed.then(Instant::now);
        let next = source.next_record()?;
        if let Some(start) = start {
            parse_ns = parse_ns.saturating_add(timed_ns(start));
        }
        let Some(record) = next else { break };
        b.push(Item::new(record, cfg)?)?;
    }
    if let Some(t) = timing.as_deref_mut() {
        t.record_read_parse_inclusive_ns =
            t.record_read_parse_inclusive_ns.saturating_add(parse_ns);
        t.read_parse_header_spill_ns = t.read_parse_header_spill_ns.saturating_add(timed_ns(stage));
    }
    let stage = Instant::now();
    let header = b.finish()?;
    verify(&header, Order::Header, lim, cfg)?;
    if let Some(t) = timing.as_deref_mut() {
        t.header_finalize_verify_ns = timed_ns(stage);
    }
    let stage = Instant::now();
    let mut reader = RunReader::open(&header, Order::Header, lim.max_record_bytes, cfg)?;
    let mut b = Builder::new(dir, Order::Barcode, lim, cfg, failure, metrics.clone());
    inject(failure, Failure::DuplicateDetection)?;
    let mut last = None::<String>;
    while let Some(x) = reader.next()? {
        if last.as_ref() == Some(&x.header) {
            return Err(invalid("grouping: duplicate complete header"));
        }
        last = Some(x.header.clone());
        b.push(x)?;
    }
    if let Some(t) = timing.as_deref_mut() {
        t.duplicate_scan_barcode_spill_ns = timed_ns(stage);
    }
    let stage = Instant::now();
    let result = b.finish()?;
    verify(&result, Order::Barcode, lim, cfg)?;
    if let Some(t) = timing {
        t.barcode_finalize_verify_ns = timed_ns(stage);
        if let Some(m) = metrics {
            t.spill = m.snapshot();
        }
    }
    Ok(result)
}

#[cfg(test)]
mod adversarial_runs {
    use super::*;
    use std::sync::atomic::{AtomicUsize, Ordering as AtomicOrdering};
    static NEXT: AtomicUsize = AtomicUsize::new(0);
    struct Scratch(PathBuf);
    impl Scratch {
        fn new() -> Self {
            let p = std::env::temp_dir().join(format!(
                "clonoforge-adversarial-{}-{}",
                std::process::id(),
                NEXT.fetch_add(1, AtomicOrdering::Relaxed)
            ));
            fs::create_dir(&p).unwrap();
            Self(p)
        }
    }
    impl Drop for Scratch {
        fn drop(&mut self) {
            let _ = fs::remove_dir_all(&self.0);
        }
    }
    fn limits() -> Limits {
        Limits {
            chunk_bytes: 10_000,
            max_record_bytes: 1024,
            fan_in: 2,
        }
    }
    fn item(id: &str, ordinal: u64) -> Item {
        let r = ReadRecord::from_header(
            &format!("{id}|BARCODE=A"),
            b"A".to_vec(),
            Some(vec![20]),
            ordinal,
        )
        .unwrap();
        Item::new(r, &Config::default()).unwrap()
    }
    fn prepare(dir: &Path) {
        for (i, id) in ["a", "b"].iter().enumerate() {
            let mut rows = [item(id, i as u64)];
            write_run(
                &name(dir, Order::Header, 0, i as u64),
                Order::Header,
                &mut rows,
                None,
                None,
            )
            .unwrap();
        }
    }
    #[test]
    fn missing_run_and_duplicated_run_are_rejected_without_deleting_inputs() {
        let t = Scratch::new();
        let cfg = Config::default();
        let lim = limits();
        prepare(&t.0);
        let p0 = name(&t.0, Order::Header, 0, 0);
        let p1 = name(&t.0, Order::Header, 0, 1);
        fs::remove_file(&p1).unwrap();
        assert!(merge(&t.0, Order::Header, 0, 2, &lim, &cfg, None, None).is_err());
        assert!(p0.exists());
        fs::copy(&p0, &p1).unwrap();
        assert!(merge(&t.0, Order::Header, 0, 2, &lim, &cfg, None, None).is_err());
        assert!(p0.exists() && p1.exists());
        let result = name(&t.0, Order::Header, 1, 0);
        assert!(verify(&result, Order::Header, &lim, &cfg).is_err());
    }
    #[test]
    fn valid_digest_cannot_hide_reordered_records() {
        let t = Scratch::new();
        let cfg = Config::default();
        let lim = limits();
        let p = name(&t.0, Order::Header, 0, 0);
        let mut w = RunWriter::new(&p, Order::Header, 2, None, None).unwrap();
        w.item(&item("z", 1)).unwrap();
        w.item(&item("a", 0)).unwrap();
        w.finish().unwrap();
        assert!(verify(&p, Order::Header, &lim, &cfg)
            .unwrap_err()
            .to_string()
            .contains("unordered"));
    }
    #[test]
    fn interrupted_merge_keeps_source_runs_and_incomplete_output() {
        let t = Scratch::new();
        let cfg = Config::default();
        let lim = limits();
        prepare(&t.0);
        let sentinel = t.0.join("foreign");
        fs::write(&sentinel, b"keep").unwrap();
        assert!(merge(
            &t.0,
            Order::Header,
            0,
            2,
            &lim,
            &cfg,
            Some(Failure::MergePartial),
            None
        )
        .is_err());
        assert!(
            name(&t.0, Order::Header, 0, 0).exists() && name(&t.0, Order::Header, 0, 1).exists()
        );
        let incomplete = name(&t.0, Order::Header, 1, 0);
        assert!(verify(&incomplete, Order::Header, &lim, &cfg).is_err());
        assert_eq!(fs::read(&sentinel).unwrap(), b"keep");
    }
}
