//! Streaming, record-bounded sequence I/O and same-directory output staging.
use crate::model::{invalid, ReadRecord};
use flate2::read::MultiGzDecoder;
use std::fs::{File, OpenOptions};
use std::io::{BufRead, BufReader, Read, Result, Write};
use std::path::{Path, PathBuf};
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Format {
    Fasta,
    Fastq,
}
pub fn input(path: &Path) -> Result<Box<dyn BufRead>> {
    let mut file = BufReader::new(File::open(path)?);
    let gzip = file.fill_buf()?.starts_with(&[0x1f, 0x8b]);
    if gzip {
        Ok(Box::new(BufReader::new(MultiGzDecoder::new(file))))
    } else {
        Ok(Box::new(file))
    }
}
/// Opt-in nested input read instrumentation. Decoder durations include raw reads.
#[derive(Default)]
pub(crate) struct ReadMetrics {
    pub raw_ns: std::sync::atomic::AtomicU64,
    pub decoded_ns: std::sync::atomic::AtomicU64,
    pub raw_bytes: std::sync::atomic::AtomicU64,
    pub decoded_bytes: std::sync::atomic::AtomicU64,
}
struct TimedRead<R> {
    inner: R,
    metrics: std::sync::Arc<ReadMetrics>,
    raw: bool,
}
impl<R: Read> Read for TimedRead<R> {
    fn read(&mut self, out: &mut [u8]) -> Result<usize> {
        use std::sync::atomic::Ordering;
        let start = std::time::Instant::now();
        let result = self.inner.read(out);
        let elapsed = u64::try_from(start.elapsed().as_nanos()).unwrap_or(u64::MAX);
        let (time, bytes) = if self.raw {
            (&self.metrics.raw_ns, &self.metrics.raw_bytes)
        } else {
            (&self.metrics.decoded_ns, &self.metrics.decoded_bytes)
        };
        let _ = time.fetch_update(Ordering::Relaxed, Ordering::Relaxed, |v| {
            Some(v.saturating_add(elapsed))
        });
        if let Ok(n) = result {
            let _ = bytes.fetch_update(Ordering::Relaxed, Ordering::Relaxed, |v| {
                Some(v.saturating_add(n as u64))
            });
        }
        result
    }
}
/// Raw duration is nested in decoded duration for gzip; the latter is not pure
/// decompression CPU. `compressed` describes file magic, not its extension.
pub(crate) fn input_timed(
    path: &Path,
) -> Result<(Box<dyn BufRead>, std::sync::Arc<ReadMetrics>, bool)> {
    let metrics = std::sync::Arc::new(ReadMetrics::default());
    let raw = TimedRead {
        inner: File::open(path)?,
        metrics: metrics.clone(),
        raw: true,
    };
    let mut file = BufReader::new(raw);
    let compressed = file.fill_buf()?.starts_with(&[0x1f, 0x8b]);
    let decoder: Box<dyn Read> = if compressed {
        Box::new(MultiGzDecoder::new(file))
    } else {
        Box::new(file)
    };
    let decoded = TimedRead {
        inner: decoder,
        metrics: metrics.clone(),
        raw: false,
    };
    Ok((Box::new(BufReader::new(decoded)), metrics, compressed))
}

pub struct Reader<R> {
    source: R,
    format: Format,
    pending: Option<String>,
    ordinal: u64,
    limit: usize,
    max_fastq_record: Option<usize>,
}
const MAX_LINE: usize = 16 * 1024 * 1024;
impl<R: BufRead> Reader<R> {
    pub fn new(source: R, format: Format) -> Self {
        Self {
            source,
            format,
            pending: None,
            ordinal: 0,
            limit: MAX_LINE,
            max_fastq_record: None,
        }
    }
    /// Read a staged FASTQ with room for the bounded generated CONSCOUNT header.
    pub(crate) fn staged_fastq(source: R) -> Self {
        Self {
            source,
            format: Format::Fastq,
            pending: None,
            ordinal: 0,
            limit: MAX_LINE + 16,
            max_fastq_record: None,
        }
    }
    /// Bounded FASTQ parser. Limits apply before header annotation expansion
    /// and quality decoding. Lines count their terminating newline as in `new`.
    pub fn fastq_limited(source: R, max_record: usize) -> Result<Self> {
        if max_record == 0 || max_record > MAX_LINE {
            return Err(invalid("invalid FASTQ record limit"));
        }
        let limit = max_record.checked_add(3).unwrap_or(MAX_LINE).min(MAX_LINE);
        Ok(Self {
            source,
            format: Format::Fastq,
            pending: None,
            ordinal: 0,
            limit,
            max_fastq_record: Some(max_record),
        })
    }
    fn line(&mut self) -> Result<Option<String>> {
        let mut bytes = Vec::new();
        let n = self
            .source
            .by_ref()
            .take((self.limit + 1) as u64)
            .read_until(b'\n', &mut bytes)?;
        if n == 0 {
            return Ok(None);
        }
        if n > self.limit {
            return Err(invalid(if self.max_fastq_record.is_some() {
                "line exceeds configured limit"
            } else {
                "line exceeds 16 MiB limit"
            }));
        }
        if bytes.last() == Some(&b'\n') {
            bytes.pop();
            if bytes.last() == Some(&b'\r') {
                bytes.pop();
            }
        }
        String::from_utf8(bytes)
            .map(Some)
            .map_err(|_| invalid("non-UTF8 input"))
    }
    pub fn next_record(&mut self) -> Result<Option<ReadRecord>> {
        let header = match self
            .pending
            .take()
            .map(Some)
            .map(Ok)
            .unwrap_or_else(|| self.line())?
        {
            Some(h) => h,
            None => return Ok(None),
        };
        let prefix = if self.format == Format::Fastq {
            '@'
        } else {
            '>'
        };
        let header = header
            .strip_prefix(prefix)
            .ok_or_else(|| invalid("invalid record header"))?;
        let mut bases = Vec::new();
        let phred = if self.format == Format::Fastq {
            bases = self
                .line()?
                .ok_or_else(|| invalid("truncated FASTQ sequence"))?
                .into_bytes();
            if let Some(maximum) = self.max_fastq_record {
                let size = bases
                    .len()
                    .checked_mul(2)
                    .and_then(|n| n.checked_add(header.len()))
                    .ok_or_else(|| invalid("FASTQ record length overflow"))?;
                if size > maximum {
                    return Err(invalid("FASTQ record exceeds configured limit"));
                }
            }
            let plus = self
                .line()?
                .ok_or_else(|| invalid("truncated FASTQ separator"))?;
            if plus != "+" && plus.strip_prefix('+') != Some(header) {
                return Err(invalid("invalid FASTQ separator"));
            }
            let q = self
                .line()?
                .ok_or_else(|| invalid("truncated FASTQ qualities"))?;
            if !q.bytes().all(|b| (33..=126).contains(&b)) {
                return Err(invalid("invalid Phred+33 encoding"));
            }
            if self.max_fastq_record.is_some() && q.len() != bases.len() {
                return Err(invalid("invalid quality length or Phred value"));
            }
            Some(q.bytes().map(|b| b - 33).collect())
        } else {
            while let Some(line) = self.line()? {
                if line.starts_with('>') {
                    self.pending = Some(line);
                    break;
                }
                if line.is_empty() {
                    return Err(invalid("blank FASTA sequence line"));
                }
                if bases
                    .len()
                    .checked_add(line.len())
                    .is_none_or(|n| n > self.limit)
                {
                    return Err(invalid("record exceeds 16 MiB limit"));
                }
                bases.extend_from_slice(line.as_bytes());
            }
            None
        };
        let record = ReadRecord::from_header(header, bases, phred, self.ordinal)?;
        self.ordinal = self
            .ordinal
            .checked_add(1)
            .ok_or_else(|| invalid("ordinal overflow"))?;
        Ok(Some(record))
    }
}
pub fn write_record<W: Write>(out: &mut W, record: &ReadRecord, format: Format) -> Result<()> {
    record.validate()?;
    // Reparse public fields to reject manually constructed malformed headers.
    let checked = ReadRecord::from_header(
        &record.header(),
        record.bases.clone(),
        record.phred.clone(),
        record.ordinal,
    )?;
    if checked.id != record.id || checked.annotations != record.annotations {
        return Err(invalid("noncanonical header fields"));
    }
    match format {
        Format::Fasta => {
            if record.phred.is_some() {
                return Err(invalid("FASTA conversion would discard qualities"));
            }
            writeln!(out, ">{}", record.header())?;
            out.write_all(&record.bases)?;
            writeln!(out)?;
        }
        Format::Fastq => {
            let q = record
                .phred
                .as_ref()
                .ok_or_else(|| invalid("FASTQ requires qualities"))?;
            writeln!(out, "@{}", record.header())?;
            out.write_all(&record.bases)?;
            writeln!(out, "\n+")?;
            for value in q {
                out.write_all(&[value + 33])?;
            }
            writeln!(out)?;
        }
    }
    Ok(())
}
pub struct StagedOutput {
    file: File,
    temporary: PathBuf,
    destination: PathBuf,
}
impl StagedOutput {
    pub fn new(destination: &Path) -> Result<Self> {
        // create_new and hard_link publication forbid clobbering existing paths.
        let mut temporary = destination.as_os_str().to_owned();
        temporary.push(format!(".partial.{}", std::process::id()));
        let temporary = PathBuf::from(temporary);
        let file = OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(&temporary)?;
        Ok(Self {
            file,
            temporary,
            destination: destination.to_owned(),
        })
    }
    pub fn commit(mut self) -> Result<()> {
        self.file.flush()?;
        self.file.sync_all()?;
        std::fs::hard_link(&self.temporary, &self.destination)?;
        std::fs::remove_file(&self.temporary)?;
        Ok(())
    }
}
impl Write for StagedOutput {
    fn write(&mut self, bytes: &[u8]) -> Result<usize> {
        self.file.write(bytes)
    }
    fn flush(&mut self) -> Result<()> {
        self.file.flush()
    }
}
impl Drop for StagedOutput {
    fn drop(&mut self) {
        let _ = std::fs::remove_file(&self.temporary);
    }
}

#[cfg(test)]
mod bounds_tests {
    use super::*;
    #[test]
    fn ordinal_overflow_is_error() {
        let mut r = Reader::new(std::io::Cursor::new(b">r\nA\n"), Format::Fasta);
        r.ordinal = u64::MAX;
        assert!(r.next_record().is_err());
    }
    #[test]
    fn aggregate_fasta_bound() {
        let mut r = Reader::new(std::io::Cursor::new(b">r\nAAAA\nAAAA\n"), Format::Fasta);
        r.limit = 6;
        assert!(r.next_record().is_err());
    }
    #[test]
    fn staged_fastq_accepts_generated_header_above_legacy_line_limit() {
        let header = format!("@{}|CONSCOUNT=1\nA\n+\nI\n", "A".repeat(MAX_LINE - 4));
        let mut legacy = Reader::new(std::io::Cursor::new(header.as_bytes()), Format::Fastq);
        assert!(legacy.next_record().is_err());
        let mut staged = Reader::staged_fastq(std::io::Cursor::new(header.as_bytes()));
        assert!(staged.next_record().unwrap().is_some());
    }
}
