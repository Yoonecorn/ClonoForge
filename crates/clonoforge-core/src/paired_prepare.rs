//! Preliminary target-only paired preprocessing. Caller owns the scratch
//! directory and must publish only after M3 completion and full verification.
use crate::{
    consensus::Config,
    external::{self, Limits},
    filter::{Disposition, Filter},
    io::{input, write_record, Format},
    model::{invalid, ReadRecord},
    pair::{matched_id, HeaderMode},
    primer::trim_exact_5prime,
    raw_fastq::RawReader,
};
use std::{
    fs::OpenOptions,
    io::{BufWriter, Result, Write},
    path::Path,
};

#[derive(Clone, Debug)]
pub struct Options {
    pub header_mode: HeaderMode,
    pub key_start: usize,
    pub key_len: usize,
    pub target_is_first: bool,
    pub remove_target_prefix: usize,
    pub primer: Option<Vec<u8>>,
    pub filter: Filter,
    pub max_raw_record: usize,
    pub max_group_count: usize,
}
/// Opt-in scientific extensions; the historical `Options` and `prepare` retain
/// the exact M4 contract when this is default.
#[derive(Clone, Debug, Default)]
pub struct Extensions {
    pub assembly: Option<crate::assembly::AssemblyConfig>,
    pub primer_iupac: Option<Vec<u8>>,
    pub primer_anchor: Option<crate::primer_ext::Anchor>,
}
#[derive(Default, Debug)]
pub struct Counts {
    pub pairs: u64,
    pub assembly_no_overlap: u64,
    pub assembly_ambiguous: u64,
    pub assembled: u64,
    pub invalid_key: u64,
    pub short_target: u64,
    pub primer_miss: u64,
    pub too_short: u64,
    pub too_many_n: u64,
    pub low_quality: u64,
    pub retained: u64,
    pub group_excluded: u64,
    pub group_excluded_groups: u64,
    pub oversize_records: u64,
    pub oversize_groups: u64,
    pub unequal_records: u64,
    pub unequal_groups: u64,
    pub qualified: u64,
    pub qualified_groups: u64,
}
fn inject(point: &str) -> Result<()> {
    if std::env::var("CLONOFORGE_M4_FAIL").ok().as_deref() == Some(point) {
        return Err(invalid(&format!("paired: injected {point}")));
    }
    Ok(())
}
fn inc(n: &mut u64) -> Result<()> {
    *n = n
        .checked_add(1)
        .ok_or_else(|| invalid("paired: count overflow"))?;
    Ok(())
}
fn hex(s: &str) -> String {
    use std::fmt::Write as _;
    let mut out = String::with_capacity(s.len().saturating_mul(2) + 1);
    out.push('R');
    for b in s.bytes() {
        write!(&mut out, "{b:02X}").unwrap();
    }
    out
}
/// Scratch must be a newly-created empty directory. Both outputs and sorted
/// runs are private intermediates. No completion marker is produced here.
pub fn prepare(
    first: &Path,
    second: &Path,
    scratch: &Path,
    opt: &Options,
    limits: &Limits,
) -> Result<Counts> {
    prepare_extended(first, second, scratch, opt, limits, &Extensions::default())
}

/// Extended bounded preprocessing. Assembly and expanded primers are explicitly
/// selected and never change the historical target-only defaults.
pub fn prepare_extended(
    first: &Path,
    second: &Path,
    scratch: &Path,
    opt: &Options,
    limits: &Limits,
    extensions: &Extensions,
) -> Result<Counts> {
    prepare_impl(first, second, scratch, opt, limits, extensions, false)
}
#[cfg(feature = "experimental-large")]
pub fn prepare_experimental(
    first: &Path,
    second: &Path,
    scratch: &Path,
    opt: &Options,
    limits: &Limits,
    extensions: &Extensions,
) -> Result<Counts> {
    if ![100, 250, 500].contains(&opt.max_group_count) {
        return Err(invalid("paired: unsupported experimental group limit"));
    }
    prepare_impl(first, second, scratch, opt, limits, extensions, true)
}
fn prepare_impl(
    first: &Path,
    second: &Path,
    scratch: &Path,
    opt: &Options,
    limits: &Limits,
    extensions: &Extensions,
    experimental: bool,
) -> Result<Counts> {
    if opt.key_len == 0 || opt.key_start.checked_add(opt.key_len).is_none() {
        return Err(invalid("paired: invalid key interval"));
    }
    if opt.max_group_count == 0 || opt.max_group_count > if experimental { 500 } else { 50 } {
        return Err(invalid("paired: invalid group limit"));
    }
    limits.validate()?;
    opt.filter.validate()?;
    if let Some(cfg) = &extensions.assembly {
        cfg.validate()?;
        if !opt.target_is_first {
            return Err(invalid("paired: assembly requires target first"));
        }
    }
    if let Some(primer) = &extensions.primer_iupac {
        crate::primer_ext::validate_iupac(primer)?;
        if primer.len() > opt.max_raw_record {
            return Err(invalid("paired: primer exceeds record bound"));
        }
        if opt.primer.is_some() {
            return Err(invalid(
                "paired: exact and IUPAC primers are mutually exclusive",
            ));
        }
    }
    if opt.max_raw_record == 0
        || opt.max_raw_record > 16 * 1024 * 1024
        || limits.max_record_bytes > 16 * 1024 * 1024
    {
        return Err(invalid("paired: invalid record bounds"));
    }
    if opt
        .primer
        .as_ref()
        .is_some_and(|p| p.is_empty() || !p.iter().all(|b| b"ACGT".contains(b)))
    {
        return Err(invalid("paired: invalid primer"));
    }
    let all_start = std::time::Instant::now();
    let mut a = RawReader::new(input(first)?, opt.max_raw_record)?;
    let mut b = RawReader::new(input(second)?, opt.max_raw_record)?;
    let mut identities = BufWriter::new(
        OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(scratch.join("identities.fastq"))?,
    );
    let mut prepared = BufWriter::new(
        OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(scratch.join("prepared.fastq"))?,
    );
    let mut counts = Counts::default();
    loop {
        let x = a.next_record()?;
        let y = b.next_record()?;
        let (x, y) = match (x, y) {
            (None, None) => break,
            (Some(x), Some(y)) => (x, y),
            _ => return Err(invalid("pair: input counts differ")),
        };
        let id = hex(matched_id(&x.header, &y.header, opt.header_mode)?);
        inc(&mut counts.pairs)?;
        inject("read")?;
        if id
            .len()
            .checked_add(12)
            .is_none_or(|n| n > limits.max_record_bytes)
        {
            return Err(invalid("paired: identity exceeds admitted record limit"));
        }
        let identity = ReadRecord::from_header(
            &format!("{id}|BARCODE=X"),
            b"A".to_vec(),
            Some(vec![0]),
            counts.pairs - 1,
        )?;
        write_record(&mut identities, &identity, Format::Fastq)?;
        let end = opt.key_start + opt.key_len;
        if x.bases.len() < end
            || !x.bases[opt.key_start..end]
                .iter()
                .all(|b| b"ACGT".contains(b))
        {
            inc(&mut counts.invalid_key)?;
            continue;
        }
        let key = std::str::from_utf8(&x.bases[opt.key_start..end])
            .map_err(|_| invalid("paired: key encoding"))?;
        let assembly_mate = extensions
            .assembly
            .as_ref()
            .map(|_| (y.bases.clone(), y.phred.clone()));
        let (mut bases, mut quality) = if opt.target_is_first {
            (x.bases.clone(), x.phred.clone())
        } else {
            (y.bases, y.phred)
        };
        if opt.remove_target_prefix >= bases.len() {
            inc(&mut counts.short_target)?;
            continue;
        }
        bases.drain(..opt.remove_target_prefix);
        quality.drain(..opt.remove_target_prefix);
        if let (Some(cfg), Some((mate_bases, mate_quality))) = (&extensions.assembly, assembly_mate)
        {
            use crate::assembly::{assemble_suffix_prefix, AssemblyDisposition};
            match assemble_suffix_prefix(&bases, &quality, &mate_bases, &mate_quality, cfg)? {
                AssemblyDisposition::Assembled(result) => {
                    bases = result.bases;
                    quality = result.phred;
                    inc(&mut counts.assembled)?;
                }
                AssemblyDisposition::NoOverlap => {
                    inc(&mut counts.assembly_no_overlap)?;
                    continue;
                }
                AssemblyDisposition::AmbiguousOverlap => {
                    inc(&mut counts.assembly_ambiguous)?;
                    continue;
                }
            }
        }
        if let Some(primer) = &extensions.primer_iupac {
            if bases.len() == primer.len() {
                // A matching full-read primer consumes the entire target.
                let mut trial_bases = bases.clone();
                let mut trial_quality = quality.clone();
                if extensions.primer_anchor == Some(crate::primer_ext::Anchor::ThreePrime) {
                    trial_bases.insert(0, b'A');
                    trial_quality.insert(0, 0);
                } else {
                    trial_bases.push(b'A');
                    trial_quality.push(0);
                }
                if crate::primer_ext::trim_iupac_anchored(
                    &mut trial_bases,
                    &mut trial_quality,
                    primer,
                    extensions
                        .primer_anchor
                        .unwrap_or(crate::primer_ext::Anchor::FivePrime),
                )? {
                    inc(&mut counts.short_target)?;
                    continue;
                }
            }
            if !crate::primer_ext::trim_iupac_anchored(
                &mut bases,
                &mut quality,
                primer,
                extensions
                    .primer_anchor
                    .unwrap_or(crate::primer_ext::Anchor::FivePrime),
            )? {
                inc(&mut counts.primer_miss)?;
                continue;
            }
        }
        if let Some(primer) = &opt.primer {
            if bases == *primer {
                inc(&mut counts.short_target)?;
                continue;
            }
            if !trim_exact_5prime(&mut bases, &mut quality, primer)? {
                inc(&mut counts.primer_miss)?;
                continue;
            }
        }
        match opt.filter.apply(&bases, &quality)? {
            Disposition::Pass => (),
            Disposition::TooShort => {
                inc(&mut counts.too_short)?;
                continue;
            }
            Disposition::TooManyNs => {
                inc(&mut counts.too_many_n)?;
                continue;
            }
            Disposition::LowMeanQuality => {
                inc(&mut counts.low_quality)?;
                continue;
            }
        }
        let generated_size = id
            .len()
            .checked_add(64)
            .and_then(|n| n.checked_add(key.len()))
            .and_then(|n| n.checked_add(bases.len().checked_mul(2)?))
            .ok_or_else(|| invalid("paired: generated record size overflow"))?;
        if generated_size > limits.max_record_bytes {
            return Err(invalid("paired: generated record exceeds admitted limit"));
        }
        let record = ReadRecord::from_header(
            &format!("{id}|BARCODE={key}"),
            bases,
            Some(quality),
            counts.pairs - 1,
        )?;
        write_record(&mut prepared, &record, Format::Fastq)?;
        inc(&mut counts.retained)?;
    }
    identities.flush()?;
    identities.get_ref().sync_all()?;
    prepared.flush()?;
    prepared.get_ref().sync_all()?;
    let read_filter_ns = all_start.elapsed().as_nanos();
    inject("identity-sort")?;
    let index_start = std::time::Instant::now();
    let index = scratch.join("index");
    std::fs::create_dir(&index)?;
    external::check_unique_headers(
        &scratch.join("identities.fastq"),
        &Config::default(),
        &index,
        limits,
    )?;
    let identity_sort_ns = index_start.elapsed().as_nanos();
    inject("group-sort")?;
    let group_sort_start = std::time::Instant::now();
    let group_dir = scratch.join("group-index");
    std::fs::create_dir(&group_dir)?;
    let run = external::sort_generated_by_barcode(
        &scratch.join("prepared.fastq"),
        &Config::default(),
        &group_dir,
        limits,
    )?;
    let group_sort_ns = group_sort_start.elapsed().as_nanos();
    let disposition_start = std::time::Instant::now();
    let mut qualified = BufWriter::new(
        OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(scratch.join("order-input.fastq"))?,
    );
    let mut key = String::new();
    let mut group = Vec::<ReadRecord>::new();
    let mut group_size = 0u64;
    let mut group_length = None;
    let mut unequal = false;
    fn emit(
        out: &mut BufWriter<std::fs::File>,
        group: &mut Vec<ReadRecord>,
        size: u64,
        unequal: bool,
        max_group: usize,
        counts: &mut Counts,
    ) -> Result<()> {
        if size == 0 {
            return Ok(());
        }
        if size > max_group as u64 || unequal {
            counts.group_excluded = counts
                .group_excluded
                .checked_add(size)
                .ok_or_else(|| invalid("paired: count overflow"))?;
            inc(&mut counts.group_excluded_groups)?;
            if size > max_group as u64 {
                counts.oversize_records = counts
                    .oversize_records
                    .checked_add(size)
                    .ok_or_else(|| invalid("paired: count overflow"))?;
                inc(&mut counts.oversize_groups)?;
            } else {
                counts.unequal_records = counts
                    .unequal_records
                    .checked_add(size)
                    .ok_or_else(|| invalid("paired: count overflow"))?;
                inc(&mut counts.unequal_groups)?;
            }
        } else {
            for record in group.iter() {
                let mut ordered = record.clone();
                let original = ordered.annotations[0].1.clone();
                ordered.annotations[0].1 = format!("{:020}", record.ordinal);
                ordered.annotations.push(("ORIGINAL".into(), original));
                write_record(out, &ordered, Format::Fastq)?;
            }
            counts.qualified = counts
                .qualified
                .checked_add(size)
                .ok_or_else(|| invalid("paired: count overflow"))?;
            inc(&mut counts.qualified_groups)?;
        }
        group.clear();
        Ok(())
    }
    external::for_each_record(&run, &Config::default(), limits, |r| {
        let next = r
            .annotations
            .iter()
            .find(|(k, _)| k == "BARCODE")
            .ok_or_else(|| invalid("paired: missing group key"))?
            .1
            .clone();
        if group_size > 0 && next != key {
            emit(
                &mut qualified,
                &mut group,
                group_size,
                unequal,
                opt.max_group_count,
                &mut counts,
            )?;
            group_size = 0;
            group_length = None;
            unequal = false;
        }
        key = next;
        group_size = group_size
            .checked_add(1)
            .ok_or_else(|| invalid("paired: group count overflow"))?;
        if group_length.is_some_and(|n| n != r.bases.len()) {
            unequal = true;
        }
        group_length = Some(r.bases.len());
        if unequal || group_size > opt.max_group_count as u64 {
            group.clear();
        } else {
            group.push(r);
        }
        Ok(())
    })?;
    emit(
        &mut qualified,
        &mut group,
        group_size,
        unequal,
        opt.max_group_count,
        &mut counts,
    )?;
    qualified.flush()?;
    qualified.get_ref().sync_all()?;
    let disposition_ns = disposition_start.elapsed().as_nanos();
    inject("order-sort")?;
    let order_sort_start = std::time::Instant::now();
    let order_dir = scratch.join("order-index");
    std::fs::create_dir(&order_dir)?;
    let ordered = external::sort_generated_by_barcode(
        &scratch.join("order-input.fastq"),
        &Config::default(),
        &order_dir,
        limits,
    )?;
    let mut final_out = BufWriter::new(
        OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(scratch.join("qualified.fastq"))?,
    );
    external::for_each_record(&ordered, &Config::default(), limits, |mut record| {
        let original = record
            .annotations
            .pop()
            .ok_or_else(|| invalid("paired: missing original barcode"))?;
        if original.0 != "ORIGINAL" {
            return Err(invalid("paired: invalid original barcode"));
        }
        record.annotations[0].1 = original.1;
        write_record(&mut final_out, &record, Format::Fastq)
    })?;
    final_out.flush()?;
    final_out.get_ref().sync_all()?;
    inject("qualified")?;
    let order_sort_and_write_ns = order_sort_start.elapsed().as_nanos();
    if counts.qualified.checked_add(counts.group_excluded) != Some(counts.retained) {
        return Err(invalid("paired: group count mismatch"));
    }
    let excluded = [
        counts.invalid_key,
        counts.assembly_no_overlap,
        counts.assembly_ambiguous,
        counts.short_target,
        counts.primer_miss,
        counts.too_short,
        counts.too_many_n,
        counts.low_quality,
        counts.retained,
    ]
    .into_iter()
    .try_fold(0u64, |a, b| a.checked_add(b));
    if excluded != Some(counts.pairs) {
        return Err(invalid("paired: disposition count mismatch"));
    }
    let reuse_start = std::time::Instant::now();
    let filter = if experimental {
        #[cfg(feature = "experimental-large")]
        {
            external::filter_paired_groups_experimental
        }
        #[cfg(not(feature = "experimental-large"))]
        {
            external::filter_paired_groups
        }
    } else {
        external::filter_paired_groups
    };
    filter(
        &run,
        &Config::default(),
        &scratch.join("group-qualified.run"),
        limits,
        opt.max_group_count,
        counts.qualified,
    )?;
    let group_reuse_ns = reuse_start.elapsed().as_nanos();
    if std::env::var_os("CLONOFORGE_PROFILE").is_some_and(|v| v == "1") {
        let _ = writeln!(
            std::io::stderr().lock(),
            "{}",
            serde_json::json!({"phase":"paired_prepare", "group_reuse_ns":group_reuse_ns, "read_filter_ns":read_filter_ns, "identity_sort_ns":identity_sort_ns, "group_sort_ns":group_sort_ns, "group_disposition_ns":disposition_ns, "order_sort_and_write_ns":order_sort_and_write_ns, "total_ns":all_start.elapsed().as_nanos()})
        );
    }
    Ok(counts)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn groups_and_duplicate_pairs() {
        let root =
            std::env::temp_dir().join(format!("clonoforge-paired-test-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&root);
        std::fs::create_dir(&root).unwrap();
        let first = root.join("r1");
        let second = root.join("r2");
        std::fs::write(
            &first,
            b"@a\nAAGG\n+\nIIII\n@b\nAACC\n+\nIIII\n@c\nTTGG\n+\nIIII\n@d\nNNGG\n+\nIIII\n",
        )
        .unwrap();
        std::fs::write(
            &second,
            b"@a\nCCCC\n+\nIIII\n@b\nCCCC\n+\nIIII\n@c\nCCCC\n+\nIIII\n@d\nCCCC\n+\nIIII\n",
        )
        .unwrap();
        let opt = Options {
            header_mode: HeaderMode::Identical,
            key_start: 0,
            key_len: 2,
            target_is_first: false,
            remove_target_prefix: 0,
            primer: None,
            filter: Filter {
                min_length: 1,
                max_n: 4,
                min_mean_phred: 0,
            },
            max_raw_record: 1024,
            max_group_count: 1,
        };
        let lim = Limits {
            chunk_bytes: 65536,
            max_record_bytes: 1024,
            fan_in: 2,
        };
        let scratch = root.join("one");
        std::fs::create_dir(&scratch).unwrap();
        let counts = prepare(&first, &second, &scratch, &opt, &lim).unwrap();
        assert_eq!(
            (
                counts.pairs,
                counts.invalid_key,
                counts.retained,
                counts.group_excluded,
                counts.qualified
            ),
            (4, 1, 3, 2, 1)
        );
        let data = std::fs::read_to_string(scratch.join("qualified.fastq")).unwrap();
        assert!(data.contains("BARCODE=TT"));
        assert!(!data.contains("BARCODE=AA"));
        std::fs::write(&first, b"@a\nAAGG\n+\nIIII\n@a\nNNGG\n+\nIIII\n").unwrap();
        std::fs::write(&second, b"@a\nCCCC\n+\nIIII\n@a\nCCCC\n+\nIIII\n").unwrap();
        let scratch = root.join("two");
        std::fs::create_dir(&scratch).unwrap();
        assert!(prepare(&first, &second, &scratch, &opt, &lim).is_err());
        let mut first_many = String::new();
        let mut second_many = String::new();
        for i in 0..5000 {
            use std::fmt::Write as _;
            writeln!(first_many, "@r{i}\nAAGG\n+\nIIII").unwrap();
            writeln!(second_many, "@r{i}\nCCCC\n+\nIIII").unwrap();
        }
        std::fs::write(&first, first_many).unwrap();
        std::fs::write(&second, second_many).unwrap();
        let scratch = root.join("three");
        std::fs::create_dir(&scratch).unwrap();
        let stress_limits = Limits {
            chunk_bytes: 1024 * 1024,
            ..lim.clone()
        };
        let many = prepare(&first, &second, &scratch, &opt, &stress_limits).unwrap();
        assert_eq!(
            (
                many.pairs,
                many.oversize_records,
                many.oversize_groups,
                many.qualified
            ),
            (5000, 5000, 1, 0)
        );
        std::fs::remove_dir_all(root).unwrap();
    }
}
