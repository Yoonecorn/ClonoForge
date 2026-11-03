//! Original record model; ordered fields never depend on hash iteration.
use std::io::{Error, ErrorKind, Result};
pub fn invalid(message: &str) -> Error {
    Error::new(ErrorKind::InvalidData, message)
}
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ReadRecord {
    pub id: String,
    pub annotations: Vec<(String, String)>,
    pub bases: Vec<u8>,
    pub phred: Option<Vec<u8>>,
    pub ordinal: u64,
}
impl ReadRecord {
    pub fn from_header(
        header: &str,
        bases: Vec<u8>,
        phred: Option<Vec<u8>>,
        ordinal: u64,
    ) -> Result<Self> {
        let mut fields = header.split('|');
        let id = fields.next().unwrap_or_default();
        if id.is_empty() || !id.bytes().all(|b| b.is_ascii_graphic() && b != b'=') {
            return Err(invalid("invalid identifier"));
        }
        let mut annotations = Vec::new();
        for field in fields {
            let (key, value) = field
                .split_once('=')
                .ok_or_else(|| invalid("annotation requires KEY=VALUE"))?;
            if key.is_empty()
                || !key
                    .bytes()
                    .all(|b| b.is_ascii_uppercase() || b.is_ascii_digit() || b == b'_')
                || value.is_empty()
                || !value.bytes().all(|b| b.is_ascii_graphic() && b != b'=')
            {
                return Err(invalid("invalid annotation"));
            }
            if annotations.iter().any(|(k, _)| k == key) {
                return Err(invalid("duplicate annotation key"));
            }
            annotations.push((key.to_owned(), value.to_owned()));
        }
        let record = Self {
            id: id.to_owned(),
            annotations,
            bases,
            phred,
            ordinal,
        };
        record.validate()?;
        Ok(record)
    }
    pub fn header(&self) -> String {
        let mut header = self.id.clone();
        for (key, value) in &self.annotations {
            header.push('|');
            header.push_str(key);
            header.push('=');
            header.push_str(value);
        }
        header
    }
    pub fn validate(&self) -> Result<()> {
        if self.bases.is_empty()
            || !self
                .bases
                .iter()
                .all(|b| b"ACGTURYSWKMBDHVNacgturyswkmbdhvn-.".contains(b))
        {
            return Err(invalid("empty sequence or unsupported sequence symbol"));
        }
        if let Some(q) = &self.phred {
            if q.len() != self.bases.len() || q.iter().any(|q| *q > 93) {
                return Err(invalid("invalid quality length or Phred value"));
            }
        }
        Ok(())
    }
}
