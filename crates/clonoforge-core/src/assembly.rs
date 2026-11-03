//! Conservative, bounded ungapped overlap assembly; independent of M3 consensus.
use crate::model::invalid;
use std::io::Result;

/// All eligible suffix/prefix overlaps are considered. Exactly one is required.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct AssemblyConfig {
    pub min_overlap: usize,
    pub max_mismatches: usize,
    pub max_read_length: usize,
    pub quality_cap: u8,
}
impl Default for AssemblyConfig {
    fn default() -> Self {
        Self {
            min_overlap: 20,
            max_mismatches: 0,
            max_read_length: 4096,
            quality_cap: 60,
        }
    }
}
impl AssemblyConfig {
    pub fn validate(&self) -> Result<()> {
        if self.min_overlap == 0
            || self.max_read_length == 0
            || self.max_read_length > 4096
            || self.min_overlap > self.max_read_length
            || self.max_mismatches >= self.min_overlap
            || self.quality_cap > 93
        {
            return Err(invalid("assembly: invalid bounded configuration"));
        }
        Ok(())
    }
}
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Assembly {
    pub bases: Vec<u8>,
    pub phred: Vec<u8>,
    pub overlap: usize,
    pub mismatches: usize,
}
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum AssemblyDisposition {
    Assembled(Assembly),
    NoOverlap,
    AmbiguousOverlap,
}
fn validate_read(b: &[u8], q: &[u8], max: usize) -> Result<()> {
    if b.is_empty()
        || b.len() > max
        || b.len() != q.len()
        || !b.iter().all(|x| b"ACGTN".contains(x))
        || q.iter().any(|x| *x > 93)
    {
        return Err(invalid(
            "assembly: invalid read or configured length exceeded",
        ));
    }
    Ok(())
}
fn complement(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        _ => b'N',
    }
}
/// Likelihood of observing a call given each possible true nucleotide.
/// Q0, Q1 and N carry no information; errors are independent symmetric substitutions.
fn likelihood(call: u8, q: u8) -> [f64; 4] {
    if call == b'N' {
        return [0.25; 4];
    }
    let e = 10f64.powf(-f64::from(q) / 10.0).min(0.75);
    let mut p = [e / 3.0; 4];
    p[b"ACGT".iter().position(|b| *b == call).unwrap()] = 1.0 - e;
    p
}
/// Uniform-prior posterior under the stated independent substitution model.
/// Exact posterior ties become N/Q0. Quality is floor(-10 log10(error)), capped.
fn reconcile(a: u8, qa: u8, b: u8, qb: u8, cap: u8) -> (u8, u8) {
    if a == b'N' && b == b'N' {
        return (b'N', 0);
    }
    // A single informative observation must retain its exact input quality.
    if a == b'N' || qa <= 1 {
        return if b == b'N' || qb <= 1 {
            (b'N', 0)
        } else {
            (b, qb.min(cap))
        };
    }
    if b == b'N' || qb <= 1 {
        return (a, qa.min(cap));
    }
    let x = likelihood(a, qa);
    let y = likelihood(b, qb);
    let w = std::array::from_fn::<_, 4, _>(|i| x[i] * y[i]);
    let mut best = 0;
    for i in 1..4 {
        if w[i] > w[best] {
            best = i;
        }
    }
    if (0..4).any(|i| i != best && w[i] == w[best]) {
        return (b'N', 0);
    }
    let total: f64 = w.iter().sum();
    // Sum alternate mass directly, avoiding cancellation at high Q.
    let error: f64 = w
        .iter()
        .enumerate()
        .filter(|(i, _)| *i != best)
        .map(|(_, p)| *p)
        .sum::<f64>()
        / total;
    let q = (-10.0 * error.log10()).floor().clamp(0.0, f64::from(cap)) as u8;
    (b"ACGT"[best], q)
}
/// Assemble forward R1 suffix to reverse-complemented raw R2 prefix.
/// N is a mismatch for eligibility, even when both calls are N. No indels,
/// containment beyond a full shorter-read overlap, adapters or sequence fill.
/// Work <= max_read_length squared; sequence/quality vectors <= 5*max_read_length bytes (excluding allocator overhead).
pub fn assemble_suffix_prefix(
    left_bases: &[u8],
    left_phred: &[u8],
    right_bases: &[u8],
    right_phred: &[u8],
    config: &AssemblyConfig,
) -> Result<AssemblyDisposition> {
    config.validate()?;
    validate_read(left_bases, left_phred, config.max_read_length)?;
    validate_read(right_bases, right_phred, config.max_read_length)?;
    let right: Vec<u8> = right_bases.iter().rev().map(|b| complement(*b)).collect();
    let mut candidate = None;
    for n in config.min_overlap..=left_bases.len().min(right.len()) {
        let mut mismatches = 0;
        for (a, b) in left_bases[left_bases.len() - n..].iter().zip(&right[..n]) {
            if a != b || *a == b'N' {
                mismatches += 1;
            }
            if mismatches > config.max_mismatches {
                break;
            }
        }
        if mismatches <= config.max_mismatches {
            if candidate.is_some() {
                return Ok(AssemblyDisposition::AmbiguousOverlap);
            }
            candidate = Some((n, mismatches));
        }
    }
    let Some((overlap, mismatches)) = candidate else {
        return Ok(AssemblyDisposition::NoOverlap);
    };
    let prefix = left_bases.len() - overlap;
    let mut bases = Vec::with_capacity(left_bases.len() + right.len() - overlap);
    let mut phred = Vec::with_capacity(bases.capacity());
    bases.extend_from_slice(&left_bases[..prefix]);
    phred.extend_from_slice(&left_phred[..prefix]);
    for i in 0..overlap {
        let (b, q) = reconcile(
            left_bases[prefix + i],
            left_phred[prefix + i],
            right[i],
            right_phred[right.len() - 1 - i],
            config.quality_cap,
        );
        bases.push(b);
        phred.push(q);
    }
    bases.extend_from_slice(&right[overlap..]);
    phred.extend(right_phred[..right.len() - overlap].iter().rev());
    Ok(AssemblyDisposition::Assembled(Assembly {
        bases,
        phred,
        overlap,
        mismatches,
    }))
}

#[cfg(test)]
mod tests {
    use super::*;
    fn c(min: usize, mismatches: usize) -> AssemblyConfig {
        AssemblyConfig {
            min_overlap: min,
            max_mismatches: mismatches,
            ..Default::default()
        }
    }
    #[test]
    fn geometry_orientation_and_quality() {
        let AssemblyDisposition::Assembled(x) = assemble_suffix_prefix(
            b"AACCGT",
            &[10, 11, 20, 20, 20, 20],
            b"TTACGG",
            &[31, 30, 20, 20, 20, 20],
            &c(4, 0),
        )
        .unwrap() else {
            panic!()
        };
        assert_eq!(x.bases, b"AACCGTAA");
        assert_eq!(x.overlap, 4);
        assert_eq!(x.mismatches, 0);
        assert_eq!(x.phred, vec![10, 11, 44, 44, 44, 44, 30, 31]);
    }
    #[test]
    fn dispositions_and_bounds() {
        assert_eq!(
            assemble_suffix_prefix(b"AAAA", &[20; 4], b"TTTT", &[20; 4], &c(2, 0)).unwrap(),
            AssemblyDisposition::AmbiguousOverlap
        );
        assert_eq!(
            assemble_suffix_prefix(b"ACGT", &[20; 4], b"AAAA", &[20; 4], &c(3, 0)).unwrap(),
            AssemblyDisposition::NoOverlap
        );
        assert_eq!(
            assemble_suffix_prefix(b"AC", &[20; 2], b"GT", &[20; 2], &c(3, 0)).unwrap(),
            AssemblyDisposition::NoOverlap
        );
        assert!(assemble_suffix_prefix(b"AC", &[20], b"GT", &[20; 2], &c(2, 0)).is_err());
        assert!(assemble_suffix_prefix(b"AX", &[20; 2], b"GT", &[20; 2], &c(2, 0)).is_err());
        assert!(assemble_suffix_prefix(b"AC", &[94; 2], b"GT", &[20; 2], &c(2, 0)).is_err());
        assert!(c(0, 0).validate().is_err());
        assert!(c(2, 2).validate().is_err());
        assert!(AssemblyConfig {
            max_read_length: 4097,
            ..Default::default()
        }
        .validate()
        .is_err());
    }
    #[test]
    fn posterior_cases() {
        assert_eq!(reconcile(b'A', 20, b'A', 20, 60), (b'A', 44));
        assert_eq!(reconcile(b'A', 20, b'C', 20, 60), (b'N', 0));
        assert_eq!(reconcile(b'A', 30, b'C', 20, 60), (b'A', 10));
        assert_eq!(reconcile(b'N', 93, b'C', 20, 60), (b'C', 20));
        assert_eq!(reconcile(b'N', 93, b'N', 20, 60), (b'N', 0));
        assert_eq!(reconcile(b'A', 0, b'C', 0, 60), (b'N', 0));
        assert_eq!(reconcile(b'A', 93, b'A', 93, 60), (b'A', 60));
        for qa in 0..=93 {
            for qb in 0..=93 {
                for a in b"ACGTN" {
                    for b in b"ACGTN" {
                        assert_eq!(reconcile(*a, qa, *b, qb, 60), reconcile(*b, qb, *a, qa, 60));
                    }
                }
            }
        }
    }
}
