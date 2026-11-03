//! Bounded raw four-line FASTQ parsing, including space-delimited instrument headers.
use crate::model::invalid;
use std::io::{BufRead, Read, Result};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RawRecord {
    pub header: String,
    pub bases: Vec<u8>,
    pub phred: Vec<u8>,
}

pub struct RawReader<R> {
    source: R,
    max: usize,
}
impl<R: BufRead> RawReader<R> {
    pub fn new(source: R, max: usize) -> Result<Self> {
        if max == 0 || max > 16 * 1024 * 1024 {
            return Err(invalid("raw FASTQ: invalid record limit"));
        }
        Ok(Self { source, max })
    }
    fn line(&mut self) -> Result<Option<Vec<u8>>> {
        let mut bytes = Vec::new();
        let n = self
            .source
            .by_ref()
            .take((self.max + 3) as u64)
            .read_until(b'\n', &mut bytes)?;
        if n == 0 {
            return Ok(None);
        }
        if n > self.max + 2 {
            return Err(invalid("raw FASTQ: line too long"));
        }
        if bytes.last() == Some(&b'\n') {
            bytes.pop();
            if bytes.last() == Some(&b'\r') {
                bytes.pop();
            }
        }
        if bytes.len() > self.max {
            return Err(invalid("raw FASTQ: line too long"));
        }
        Ok(Some(bytes))
    }
    pub fn next_record(&mut self) -> Result<Option<RawRecord>> {
        let Some(header) = self.line()? else {
            return Ok(None);
        };
        let header = header
            .strip_prefix(b"@")
            .ok_or_else(|| invalid("raw FASTQ: expected header"))?;
        if header.is_empty()
            || header.first() == Some(&b' ')
            || header.last() == Some(&b' ')
            || header.windows(2).any(|w| w == b"  ")
            || !header.iter().all(|b| b.is_ascii_graphic() || *b == b' ')
        {
            return Err(invalid("raw FASTQ: invalid header"));
        }
        let bases = self
            .line()?
            .ok_or_else(|| invalid("raw FASTQ: truncated sequence"))?;
        let sep = self
            .line()?
            .ok_or_else(|| invalid("raw FASTQ: truncated separator"))?;
        let quality = self
            .line()?
            .ok_or_else(|| invalid("raw FASTQ: truncated quality"))?;
        if sep != b"+" && !(sep.first() == Some(&b'+') && &sep[1..] == header) {
            return Err(invalid("raw FASTQ: invalid separator"));
        }
        if bases.is_empty() || !bases.iter().all(|b| b"ACGTN".contains(b)) {
            return Err(invalid("raw FASTQ: unsupported sequence"));
        }
        if quality.len() != bases.len() || !quality.iter().all(|b| (33..=126).contains(b)) {
            return Err(invalid("raw FASTQ: invalid qualities"));
        }
        if header
            .len()
            .checked_add(bases.len())
            .and_then(|n| n.checked_add(quality.len()))
            .is_none_or(|n| n > self.max)
        {
            return Err(invalid("raw FASTQ: record exceeds limit"));
        }
        Ok(Some(RawRecord {
            header: String::from_utf8(header.to_vec())
                .map_err(|_| invalid("raw FASTQ: invalid header encoding"))?,
            bases,
            phred: quality.into_iter().map(|b| b - 33).collect(),
        }))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn bounded_raw_parser() {
        let mut r = RawReader::new(&b"@x 1:N:0:AC\nACN\n+\n!I~\n"[..], 30).unwrap();
        assert!(RawReader::new(&b"@ x\nA\n+\n!\n"[..], 30)
            .unwrap()
            .next_record()
            .is_err());
        assert!(RawReader::new(&b"@x\nA\n+\n!"[..], 30)
            .unwrap()
            .next_record()
            .unwrap()
            .is_some());
        let x = r.next_record().unwrap().unwrap();
        assert_eq!(x.header, "x 1:N:0:AC");
        assert_eq!(x.phred, vec![0, 40, 93]);
        assert!(r.next_record().unwrap().is_none());
        for input in [&b"@x\nAC\n+\n!\n"[..], &b"@x\nAc\n+\n!!\n"[..]] {
            assert!(RawReader::new(input, 100).unwrap().next_record().is_err());
        }
        assert!(RawReader::new(&b"@x\nAC\n+\n!!\n"[..], 4)
            .unwrap()
            .next_record()
            .is_err());
    }
}
