//! Strict lockstep pairing primitives. Pairing never searches ahead or silently
//! discards unmatched reads. Identity uniqueness across a stream is the caller's
//! responsibility and requires a bounded external index.
use crate::model::invalid;
use std::io::Result;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum HeaderMode {
    /// Both complete raw headers must be equal (e.g. derived barcode/target files).
    Identical,
    /// First whitespace-delimited token ends in /1 and /2 respectively;
    /// optional trailing descriptions must match exactly.
    SlashSuffix,
    /// SRA paired headers with a shared accession token and exact `1/1`, `1/2` descriptions.
    SraNumbered,
}

pub fn matched_id<'a>(first: &'a str, second: &str, mode: HeaderMode) -> Result<&'a str> {
    if first.is_empty()
        || second.is_empty()
        || first.starts_with(' ')
        || second.starts_with(' ')
        || first.ends_with(' ')
        || second.ends_with(' ')
        || first.contains("  ")
        || second.contains("  ")
        || first.bytes().any(|b| !b.is_ascii_graphic() && b != b' ')
        || second.bytes().any(|b| !b.is_ascii_graphic() && b != b' ')
    {
        return Err(invalid("pair: invalid header"));
    }
    match mode {
        HeaderMode::Identical => {
            if first != second {
                return Err(invalid("pair: headers differ"));
            }
            Ok(first.split(' ').next().unwrap())
        }
        HeaderMode::SraNumbered => {
            let (a, ad) = first
                .split_once(' ')
                .ok_or_else(|| invalid("pair: missing SRA tag"))?;
            let (b, bd) = second
                .split_once(' ')
                .ok_or_else(|| invalid("pair: missing SRA tag"))?;
            let spot = a.rsplit_once('.').map(|(_, spot)| spot).unwrap_or("");
            let valid_spot = !spot.is_empty()
                && !spot.starts_with('0')
                && spot.bytes().all(|b| b.is_ascii_digit());
            if a.is_empty()
                || a != b
                || !valid_spot
                || ad.strip_suffix("/1") != Some(spot)
                || bd.strip_suffix("/2") != Some(spot)
            {
                return Err(invalid("pair: mismatched SRA headers"));
            }
            Ok(a)
        }
        HeaderMode::SlashSuffix => {
            let (a, ad) = first.split_once(' ').unwrap_or((first, ""));
            let (b, bd) = second.split_once(' ').unwrap_or((second, ""));
            let a = a
                .strip_suffix("/1")
                .ok_or_else(|| invalid("pair: first header lacks /1"))?;
            let b = b
                .strip_suffix("/2")
                .ok_or_else(|| invalid("pair: second header lacks /2"))?;
            if a.is_empty() || a != b || ad != bd || a.ends_with("/1") || a.ends_with("/2") {
                return Err(invalid("pair: mismatched or ambiguous headers"));
            }
            Ok(a)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn strict_identity() {
        assert_eq!(
            matched_id("SRR1.2 2/1", "SRR1.2 2/2", HeaderMode::SraNumbered).unwrap(),
            "SRR1.2"
        );
        assert_eq!(
            matched_id("SRR1.1 1/1", "SRR1.1 1/2", HeaderMode::SraNumbered).unwrap(),
            "SRR1.1"
        );
        assert!(matched_id("SRR1.1 1/1", "SRR1.2 1/2", HeaderMode::SraNumbered).is_err());
        assert!(matched_id("  ", "  ", HeaderMode::Identical).is_err());
        assert!(matched_id(" x", " x", HeaderMode::Identical).is_err());
        assert_eq!(
            matched_id("x/1 foo", "x/2 foo", HeaderMode::SlashSuffix).unwrap(),
            "x"
        );
        assert!(matched_id("x/1", "y/2", HeaderMode::SlashSuffix).is_err());
        assert!(matched_id("x/1/1", "x/1/2", HeaderMode::SlashSuffix).is_err());
        assert!(matched_id("x/1", "x/1", HeaderMode::SlashSuffix).is_err());
        assert_eq!(
            matched_id("x 1:N:0:A", "x 1:N:0:A", HeaderMode::Identical).unwrap(),
            "x"
        );
        assert!(matched_id("x 1:N:0:A", "x 2:N:0:A", HeaderMode::Identical).is_err());
    }
}
