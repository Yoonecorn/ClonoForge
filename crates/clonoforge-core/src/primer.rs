//! Explicit exact, anchored primer handling; no mismatch or indel inference.
use crate::model::invalid;
use std::io::Result;

/// Remove an exact 5-prime primer and its corresponding qualities. Returns
/// false on a mismatch without changing either input.
pub fn trim_exact_5prime(
    bases: &mut Vec<u8>,
    qualities: &mut Vec<u8>,
    primer: &[u8],
) -> Result<bool> {
    if primer.is_empty()
        || !primer.iter().all(|b| b"ACGT".contains(b))
        || bases.len() != qualities.len()
    {
        return Err(invalid("primer: invalid sequence or quality"));
    }
    if !bases.starts_with(primer) {
        return Ok(false);
    }
    if bases.len() == primer.len() {
        return Ok(false);
    }
    bases.drain(..primer.len());
    qualities.drain(..primer.len());
    Ok(true)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn exact_only_and_quality_co_trim() {
        let mut b = b"ACGT".to_vec();
        let mut q = vec![1, 2, 3, 4];
        assert!(!trim_exact_5prime(&mut b, &mut q, b"AG").unwrap());
        assert_eq!(b, b"ACGT");
        assert_eq!(q, vec![1, 2, 3, 4]);
        assert!(trim_exact_5prime(&mut b, &mut q, b"AC").unwrap());
        assert_eq!(b, b"GT");
        assert_eq!(q, vec![3, 4]);
        assert!(!trim_exact_5prime(&mut b, &mut q, b"GT").unwrap());
        assert!(trim_exact_5prime(&mut b, &mut q, b"N").is_err());
    }
}
