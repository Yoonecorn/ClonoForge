//! Pure bounded consensus kernel derived from qualified M2A observations.
use crate::model::{invalid, ReadRecord};
use std::io::{Error, ErrorKind, Result};

#[derive(Clone, Debug, PartialEq)]
pub struct Config {
    pub barcode_field: String,
    pub min_count: usize,
    pub frequency: f64,
    pub min_quality: u8,
    pub failed: bool,
}
impl Default for Config {
    fn default() -> Self {
        Self {
            barcode_field: "BARCODE".into(),
            min_count: 1,
            frequency: 0.6,
            min_quality: 0,
            failed: false,
        }
    }
}
impl Config {
    pub fn validate(&self) -> Result<()> {
        if self.barcode_field.is_empty()
            || !self
                .barcode_field
                .bytes()
                .all(|b| b.is_ascii_uppercase() || b.is_ascii_digit() || b == b'_')
        {
            return Err(Error::new(
                ErrorKind::InvalidInput,
                "configuration: barcode field must be a nonempty uppercase annotation key",
            ));
        }
        if self.min_count == 0
            || !self.frequency.is_finite()
            || !(0.0..=1.0).contains(&self.frequency)
            || self.min_quality > 93
        {
            return Err(Error::new(
                ErrorKind::InvalidInput,
                "configuration: require n>0, finite frequency in [0,1], quality in [0,93]",
            ));
        }
        Ok(())
    }
}

/// Compare against the correctly rounded binary64 count fraction used by the
/// qualified reference. Counts <=50 convert to binary64 exactly; one division
/// produces the required fraction with no intermediate decimal rounding or
/// multiplication. Comparing the *unrounded rational* with the threshold would
/// differ at a threshold equal to the rounded fraction, e.g. 2/3.
pub fn frequency_passes(count: u32, eligible: u32, threshold: f64) -> Result<bool> {
    frequency_bounded(count, eligible, threshold, 50)
}
fn frequency_bounded(count: u32, eligible: u32, threshold: f64, limit: usize) -> Result<bool> {
    if eligible == 0
        || eligible > limit as u32
        || count > eligible
        || !threshold.is_finite()
        || !(0.0..=1.0).contains(&threshold)
    {
        return Err(invalid("scientific: invalid bounded frequency comparison"));
    }
    Ok(f64::from(count) / f64::from(eligible) >= threshold)
}

/// One independent position, including the separately qualified singleton path.
/// Observations must contain one through 50 uppercase ACGTN bases and Q0..93.
pub fn position(observations: &[(u8, u8)], frequency: f64, min_quality: u8) -> Result<(u8, u8)> {
    position_bounded(observations, frequency, min_quality, 50)
}
fn position_bounded(
    observations: &[(u8, u8)],
    frequency: f64,
    min_quality: u8,
    limit: usize,
) -> Result<(u8, u8)> {
    if observations.is_empty()
        || observations.len() > limit
        || !frequency.is_finite()
        || !(0.0..=1.0).contains(&frequency)
        || min_quality > 93
    {
        return Err(invalid("scientific: position outside bounded contract"));
    }
    if observations
        .iter()
        .any(|(b, q)| !b"ACGTN".contains(b) || *q > 93)
    {
        return Err(invalid("scientific: unsupported base or quality"));
    }
    if observations.len() == 1 {
        let (base, quality) = observations[0];
        return Ok((if quality < min_quality { b'N' } else { base }, quality));
    }
    let mut sums = [0u32; 4];
    let mut counts = [0u32; 4];
    let mut total = 0u32;
    let mut eligible = 0u32;
    for &(base, quality) in observations {
        total = total
            .checked_add(u32::from(quality))
            .ok_or_else(|| invalid("scientific: total-quality overflow"))?;
        if let Some(i) = b"ACGT".iter().position(|b| *b == base) {
            eligible = eligible
                .checked_add(1)
                .ok_or_else(|| invalid("scientific: count overflow"))?;
            counts[i] = counts[i]
                .checked_add(1)
                .ok_or_else(|| invalid("scientific: count overflow"))?;
            sums[i] = sums[i]
                .checked_add(u32::from(quality))
                .ok_or_else(|| invalid("scientific: support-quality overflow"))?;
        }
    }
    let mut scores = [0u32; 4];
    if total != 0 {
        for i in 0..4 {
            scores[i] = sums[i]
                .checked_mul(sums[i])
                .ok_or_else(|| invalid("scientific: squared-quality overflow"))?
                .checked_div(total)
                .ok_or_else(|| invalid("scientific: quality division failure"))?;
        }
    }
    let mut winner = None;
    for i in 0..4 {
        if counts[i] > 0 && winner.is_none_or(|previous| scores[i] > scores[previous]) {
            winner = Some(i);
        }
    }
    let Some(i) = winner else {
        return Ok((b'N', 0));
    };
    let quality = scores[i].min(90) as u8;
    if !frequency_bounded(counts[i], eligible, frequency, limit)? || quality < min_quality {
        Ok((b'N', 0))
    } else {
        Ok((b"ACGT"[i], quality))
    }
}

/// Transform a validated eligible group without performing file I/O.
pub fn consensus(records: &[ReadRecord], barcode: &str, config: &Config) -> Result<ReadRecord> {
    consensus_bounded(records, barcode, config, 50)
}
#[cfg(feature = "experimental-large")]
pub fn experimental_consensus(
    records: &[ReadRecord],
    barcode: &str,
    config: &Config,
) -> Result<ReadRecord> {
    consensus_bounded(records, barcode, config, 500)
}
fn consensus_bounded(
    records: &[ReadRecord],
    barcode: &str,
    config: &Config,
    limit: usize,
) -> Result<ReadRecord> {
    config.validate()?;
    if records.is_empty() || records.len() > limit || records.len() < config.min_count {
        return Err(invalid(
            "scientific: group is empty, ineligible or exceeds 50 reads",
        ));
    }
    if barcode.is_empty()
        || !barcode
            .bytes()
            .all(|b| b.is_ascii_alphanumeric() || b == b'_')
    {
        return Err(invalid("scientific: unsafe barcode"));
    }
    let length = records[0].bases.len();
    for record in records {
        record.validate()?;
        if record.bases.len() != length
            || !record.bases.iter().all(|b| b"ACGTN".contains(b))
            || record.phred.is_none()
        {
            return Err(invalid(
                "scientific: group requires equal-length uppercase ACGTN FASTQ records",
            ));
        }
    }
    let mut bases = Vec::with_capacity(length);
    let mut qualities = Vec::with_capacity(length);
    let mut observations = Vec::with_capacity(records.len());
    for i in 0..length {
        observations.clear();
        for record in records {
            let q = record
                .phred
                .as_ref()
                .ok_or_else(|| invalid("scientific: missing qualities"))?;
            observations.push((record.bases[i], q[i]));
        }
        let (base, quality) =
            position_bounded(&observations, config.frequency, config.min_quality, limit)?;
        bases.push(base);
        qualities.push(quality);
    }
    ReadRecord::from_header(
        &format!("{barcode}|CONSCOUNT={}", records.len()),
        bases,
        Some(qualities),
        records[0].ordinal,
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    fn call(bases: &[u8], qualities: &[u8], frequency: f64, cutoff: u8) -> (u8, u8) {
        position(
            &bases
                .iter()
                .copied()
                .zip(qualities.iter().copied())
                .collect::<Vec<_>>(),
            frequency,
            cutoff,
        )
        .unwrap()
    }
    #[test]
    fn unanimous_and_all_n() {
        for base in b"ACGT" {
            assert_eq!(call(&[*base; 3], &[20; 3], 0.6, 0), (*base, 60));
        }
        assert_eq!(call(b"NN", &[93, 93], 0.0, 0), (b'N', 0));
    }
    #[test]
    fn candidate_selection_precedes_cap_and_uses_floored_scores() {
        assert_eq!(call(b"TTCGTT", &[2, 9, 34, 8, 15, 9], 0.0, 0), (b'C', 15));
        assert_eq!(call(b"TTCGTT", &[2, 9, 34, 8, 15, 9], 0.6, 0), (b'N', 0));
        assert_eq!(call(b"CAT", &[1, 0, 1], 0.0, 0), (b'A', 0));
        assert_eq!(
            call(b"TGTGGT", &[92, 86, 88, 85, 84, 87], 0.0, 0),
            (b'T', 90)
        );
    }
    #[test]
    fn sums_reject_count_and_max_individual() {
        assert_eq!(call(b"AAC", &[1, 1, 3], 0.0, 0), (b'C', 1));
        assert_eq!(call(b"AAC", &[15, 15, 25], 0.0, 0), (b'A', 16));
        assert_eq!(call(b"AAC", &[0, 0, 40], 0.6, 0), (b'N', 0));
    }
    #[test]
    fn all_canonical_tie_pairs_and_reversals() {
        let alphabet = b"ACGT";
        for i in 0..4 {
            for j in (i + 1)..4 {
                for bases in [[alphabet[i], alphabet[j]], [alphabet[j], alphabet[i]]] {
                    assert_eq!(call(&bases, &[20, 20], 0.5, 0), (alphabet[i], 10));
                    assert_eq!(call(&bases, &[0, 0], 0.0, 0), (alphabet[i], 0));
                }
            }
        }
        for p in [b"ACGT", b"TGCA", b"TCGA", b"GACT"] {
            assert_eq!(call(p, &[30; 4], 0.0, 0), (b'A', 7));
        }
        assert_eq!(call(b"GG", &[0, 0], 0.6, 0), (b'G', 0));
    }
    #[test]
    fn missing_candidates_frequency_and_total() {
        assert_eq!(call(b"ANN", &[20, 20, 20], 1.0, 0), (b'A', 6));
        assert_eq!(call(b"AAN", &[10, 10, 30], 1.0, 0), (b'A', 8));
        assert_eq!(call(b"NNN", &[0, 20, 93], 0.0, 0), (b'N', 0));
    }
    #[test]
    fn adjacent_binary64_frequency_boundaries() {
        for (count, total) in [(1, 2), (3, 5), (2, 3), (1, 1)] {
            let f = f64::from(count) / f64::from(total);
            assert!(frequency_passes(count, total, f).unwrap());
            assert!(frequency_passes(count, total, f64::from_bits(f.to_bits() - 1)).unwrap());
            if f < 1.0 {
                assert!(!frequency_passes(count, total, f64::from_bits(f.to_bits() + 1)).unwrap());
            }
        }
        assert_eq!(call(b"AAACC", &[30; 5], 0.6, 0), (b'A', 54));
        assert_eq!(
            call(b"AAACC", &[30; 5], f64::from_bits(0.6f64.to_bits() + 1), 0),
            (b'N', 0)
        );
    }
    #[test]
    fn floors_clipping_and_quality_thresholds() {
        assert_eq!(call(b"AC", &[15, 15], 0.0, 0), (b'A', 7));
        assert_eq!(call(b"AACGT", &[49, 50, 88, 88, 88], 0.0, 0), (b'A', 27));
        for (qs, q) in [([44, 45], 89), ([45, 45], 90), ([45, 46], 90)] {
            assert_eq!(call(b"AA", &qs, 0.6, 0), (b'A', q));
        }
        assert_eq!(call(b"AA", &[10, 10], 0.6, 20), (b'A', 20));
        assert_eq!(call(b"AA", &[10, 10], 0.6, 21), (b'N', 0));
        assert_eq!(call(b"AC", &[0, 0], 0.0, 0), (b'A', 0));
    }
    #[test]
    fn singleton_path_preserves_quality() {
        for q in [0, 20, 93] {
            for base in b"ACGTN" {
                assert_eq!(call(&[*base], &[q], 1.0, 0), (*base, q));
                if q < 93 {
                    assert_eq!(call(&[*base], &[q], 0.0, q + 1), (b'N', q));
                }
            }
        }
    }
    fn record(id: &str, bases: &[u8], q: &[u8], ordinal: u64) -> ReadRecord {
        ReadRecord::from_header(id, bases.to_vec(), Some(q.to_vec()), ordinal).unwrap()
    }
    #[test]
    fn positions_and_output_header_independent() {
        let rs = [
            record("a", b"ACGT", &[1, 10, 30, 93], 8),
            record("b", b"ACGT", &[2, 20, 40, 0], 9),
            record("c", b"ACGT", &[3, 30, 50, 1], 10),
        ];
        let out = consensus(&rs, "Family_1", &Config::default()).unwrap();
        assert_eq!(out.header(), "Family_1|CONSCOUNT=3");
        assert_eq!(out.bases, b"ACGT");
        assert_eq!(out.phred, Some(vec![6, 60, 90, 90]));
        assert_eq!(out.ordinal, 8);
    }
    #[test]
    fn rejects_unsupported_and_invalid_configuration() {
        for frequency in [-1.0, 1.1, f64::NAN, f64::INFINITY] {
            assert!(Config {
                frequency,
                ..Config::default()
            }
            .validate()
            .is_err());
        }
        assert!(Config {
            min_count: 0,
            ..Config::default()
        }
        .validate()
        .is_err());
        assert!(Config {
            min_quality: 94,
            ..Config::default()
        }
        .validate()
        .is_err());
        assert!(Config {
            barcode_field: "barcode".into(),
            ..Config::default()
        }
        .validate()
        .is_err());
        assert!(position(&[], 0.6, 0).is_err());
        assert_eq!(position(&[(b'A', 20); 11], 0.6, 0).unwrap(), (b'A', 90));
        assert!(position(&[(b'A', 20); 51], 0.6, 0).is_err());
        for o in [(b'a', 20), (b'R', 20), (b'-', 20), (b'A', 94)] {
            assert!(position(&[o], 0.6, 0).is_err());
        }
        let rs = [
            record("a", b"A", &[20], 0),
            record("b", b"AA", &[20, 20], 1),
        ];
        assert!(consensus(&rs, "X", &Config::default()).is_err());
        assert!(frequency_passes(1, 0, 0.6).is_err());
        assert!(frequency_passes(2, 1, 0.6).is_err());
    }
}
