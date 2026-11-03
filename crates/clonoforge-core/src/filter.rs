//! Simple deterministic per-read filters; no probabilistic quality inference.
use crate::model::invalid;
use std::io::Result;

#[derive(Clone, Copy, Debug)]
pub struct Filter {
    pub min_length: usize,
    pub max_n: usize,
    /// Minimum arithmetic mean Phred, compared using integer arithmetic.
    pub min_mean_phred: u8,
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Disposition {
    Pass,
    TooShort,
    TooManyNs,
    LowMeanQuality,
}
impl Filter {
    pub fn validate(self) -> Result<()> {
        if self.min_length == 0 || self.min_mean_phred > 93 {
            return Err(invalid("filter: invalid configuration"));
        }
        Ok(())
    }
    pub fn apply(self, bases: &[u8], qualities: &[u8]) -> Result<Disposition> {
        self.validate()?;
        if bases.len() != qualities.len()
            || bases.is_empty()
            || !bases.iter().all(|b| b"ACGTN".contains(b))
            || qualities.iter().any(|q| *q > 93)
        {
            return Err(invalid("filter: invalid configuration or input"));
        }
        if bases.len() < self.min_length {
            return Ok(Disposition::TooShort);
        }
        if bases.iter().filter(|b| **b == b'N').count() > self.max_n {
            return Ok(Disposition::TooManyNs);
        }
        let sum = qualities
            .iter()
            .try_fold(0u64, |s, q| s.checked_add(u64::from(*q)))
            .ok_or_else(|| invalid("filter: quality sum overflow"))?;
        let threshold = (bases.len() as u64)
            .checked_mul(u64::from(self.min_mean_phred))
            .ok_or_else(|| invalid("filter: quality threshold overflow"))?;
        if sum < threshold {
            return Ok(Disposition::LowMeanQuality);
        }
        Ok(Disposition::Pass)
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn ordered_dispositions_and_mean_boundary() {
        let f = Filter {
            min_length: 2,
            max_n: 0,
            min_mean_phred: 10,
        };
        assert_eq!(f.apply(b"N", &[0]).unwrap(), Disposition::TooShort);
        assert_eq!(f.apply(b"NA", &[0, 0]).unwrap(), Disposition::TooManyNs);
        assert_eq!(
            f.apply(b"AA", &[9, 10]).unwrap(),
            Disposition::LowMeanQuality
        );
        assert_eq!(f.apply(b"AA", &[9, 11]).unwrap(), Disposition::Pass);
        assert!(f.apply(b"AA", &[0]).is_err());
    }
}
