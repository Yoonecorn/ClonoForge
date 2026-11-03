//! Anchored IUPAC primer matching without mismatches or indels.
use crate::model::invalid;
use std::io::Result;
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Anchor {
    FivePrime,
    ThreePrime,
}
fn mask(b: u8) -> Option<u8> {
    Some(match b {
        b'A' => 1,
        b'C' => 2,
        b'G' => 4,
        b'T' => 8,
        b'R' => 5,
        b'Y' => 10,
        b'S' => 6,
        b'W' => 9,
        b'K' => 12,
        b'M' => 3,
        b'B' => 14,
        b'D' => 13,
        b'H' => 11,
        b'V' => 7,
        b'N' => 15,
        _ => return None,
    })
}
pub fn validate_iupac(primer: &[u8]) -> Result<()> {
    if primer.is_empty() || primer.len() > 4096 || primer.iter().any(|b| mask(*b).is_none()) {
        return Err(invalid(
            "primer IUPAC: invalid primer or length exceeds 4096",
        ));
    }
    Ok(())
}
/// A read N never satisfies a primer position, including primer N. The primer
/// symbols denote allowed canonical bases, not inferred read calls. A full-read
/// match returns false because trimming would leave an empty target.
pub fn trim_iupac_anchored(
    bases: &mut Vec<u8>,
    qualities: &mut Vec<u8>,
    primer: &[u8],
    anchor: Anchor,
) -> Result<bool> {
    validate_iupac(primer)?;
    if bases.len() != qualities.len()
        || !bases.iter().all(|b| b"ACGTN".contains(b))
        || qualities.iter().any(|q| *q > 93)
    {
        return Err(invalid("primer IUPAC: invalid read or qualities"));
    }
    if primer.len() >= bases.len() {
        return Ok(false);
    }
    let start = match anchor {
        Anchor::FivePrime => 0,
        Anchor::ThreePrime => bases.len() - primer.len(),
    };
    if !bases[start..start + primer.len()]
        .iter()
        .zip(primer)
        .all(|(b, p)| *b != b'N' && mask(*b).unwrap() & mask(*p).unwrap() != 0)
    {
        return Ok(false);
    }
    match anchor {
        Anchor::FivePrime => {
            bases.drain(..primer.len());
            qualities.drain(..primer.len());
        }
        Anchor::ThreePrime => {
            bases.truncate(start);
            qualities.truncate(start);
        }
    }
    Ok(true)
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn exhaustive_iupac_table_and_quality_trim() {
        for (symbol, allowed) in [
            (b'A', "A"),
            (b'C', "C"),
            (b'G', "G"),
            (b'T', "T"),
            (b'R', "AG"),
            (b'Y', "CT"),
            (b'S', "CG"),
            (b'W', "AT"),
            (b'K', "GT"),
            (b'M', "AC"),
            (b'B', "CGT"),
            (b'D', "AGT"),
            (b'H', "ACT"),
            (b'V', "ACG"),
            (b'N', "ACGT"),
        ] {
            for b in b"ACGTN" {
                let mut seq = vec![*b, b'T'];
                let mut q = vec![1, 9];
                let hit =
                    trim_iupac_anchored(&mut seq, &mut q, &[symbol], Anchor::FivePrime).unwrap();
                assert_eq!(hit, allowed.as_bytes().contains(b));
                assert_eq!(q, if hit { vec![9] } else { vec![1, 9] });
            }
        }
    }
    #[test]
    fn terminal_misses_empty_and_invalid() {
        let mut b = b"AACG".to_vec();
        let mut q = vec![1, 2, 3, 4];
        assert!(trim_iupac_anchored(&mut b, &mut q, b"YS", Anchor::ThreePrime).unwrap());
        assert_eq!(b, b"AA");
        assert_eq!(q, vec![1, 2]);
        assert!(!trim_iupac_anchored(&mut b, &mut q, b"NN", Anchor::FivePrime).unwrap());
        assert!(trim_iupac_anchored(&mut b, &mut q, b"ac", Anchor::FivePrime).is_err());
        assert!(validate_iupac(b"").is_err());
        assert!(validate_iupac(&vec![b'A'; 4097]).is_err());
    }
}
