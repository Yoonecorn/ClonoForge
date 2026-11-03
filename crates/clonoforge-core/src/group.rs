//! Bounded consensus admission; stable vector order with ordered lookup only.
use crate::{
    consensus::Config,
    io::{input, Format, Reader},
    model::{invalid, ReadRecord},
};
use std::{
    collections::{BTreeMap, BTreeSet},
    io::BufRead,
    path::Path,
};
#[derive(Debug)]
pub struct Group {
    pub barcode: String,
    pub records: Vec<ReadRecord>,
}
pub fn load(path: &Path, config: &Config) -> std::io::Result<Vec<Group>> {
    config.validate()?;
    admit(Reader::new(input(path)?, Format::Fastq), config)
}
pub fn admit<R: BufRead>(mut reader: Reader<R>, config: &Config) -> std::io::Result<Vec<Group>> {
    config.validate()?;
    let mut groups: Vec<Group> = Vec::new();
    let mut lookup = BTreeMap::new();
    let mut headers = BTreeSet::new();
    while let Some(record) = reader
        .next_record()
        .map_err(|e| invalid(&format!("input: {e}")))?
    {
        if !record.bases.iter().all(|b| b"ACGTN".contains(b)) {
            return Err(invalid("input: only uppercase ACGTN supported"));
        }
        if record.phred.is_none() {
            return Err(invalid("input: FASTQ qualities required"));
        }
        let value = record
            .annotations
            .iter()
            .find(|(k, _)| k == &config.barcode_field)
            .map(|(_, v)| v)
            .ok_or_else(|| invalid("grouping: missing selected barcode"))?;
        if value.is_empty()
            || !value
                .bytes()
                .all(|b| b.is_ascii_alphanumeric() || b == b'_')
        {
            return Err(invalid("grouping: unsafe or empty barcode"));
        }
        if !headers.insert(record.header()) {
            return Err(invalid("grouping: duplicate complete header"));
        }
        let index = if let Some(i) = lookup.get(value) {
            *i
        } else {
            let i = groups.len();
            lookup.insert(value.clone(), i);
            groups.push(Group {
                barcode: value.clone(),
                records: Vec::new(),
            });
            i
        };
        let group = &mut groups[index];
        if group.records.len() >= 50 {
            return Err(invalid("grouping: group exceeds 50 reads"));
        }
        if group
            .records
            .first()
            .is_some_and(|r| r.bases.len() != record.bases.len())
        {
            return Err(invalid("grouping: unequal read lengths"));
        }
        group.records.push(record);
    }
    Ok(groups)
}

/// A single-pass bounded grouping attempt. A spill transition retains the
/// original parser and input ordinal, so the input is never reopened.
pub enum GroupAttempt<R> {
    Complete(Vec<Group>),
    Spill {
        reader: Reader<R>,
        groups: Vec<Group>,
        pending: ReadRecord,
    },
}
pub fn consume_capped<R: BufRead>(
    mut reader: Reader<R>,
    config: &Config,
    cap: usize,
    max_group_count: usize,
    max_group_bytes: usize,
    failure: Option<crate::publication::Failure>,
    mut parse_ns: Option<&mut u64>,
) -> std::io::Result<GroupAttempt<R>> {
    config.validate()?;
    if cap == 0 || max_group_count == 0 || max_group_count > 50 {
        return Err(invalid("grouping: invalid in-memory limits"));
    }
    let mut groups: Vec<Group> = Vec::new();
    let mut lookup: BTreeMap<String, usize> = BTreeMap::new();
    let mut headers = BTreeSet::new();
    let mut used = 0usize;
    let mut group_weights: Vec<usize> = Vec::new();
    loop {
        let timer = parse_ns.as_ref().map(|_| std::time::Instant::now());
        let next = reader
            .next_record()
            .map_err(|e| invalid(&format!("input: {e}")))?;
        if let (Some(ns), Some(timer)) = (parse_ns.as_deref_mut(), timer) {
            *ns = ns.saturating_add(u64::try_from(timer.elapsed().as_nanos()).unwrap_or(u64::MAX));
        }
        let Some(record) = next else { break };
        let value = record
            .annotations
            .iter()
            .find(|(k, _)| k == &config.barcode_field)
            .map(|(_, v)| v);
        let bad = !record.bases.iter().all(|b| b"ACGTN".contains(b))
            || record.phred.is_none()
            || value.is_none_or(|v| {
                v.is_empty() || !v.bytes().all(|b| b.is_ascii_alphanumeric() || b == b'_')
            });
        if bad {
            return Ok(GroupAttempt::Spill {
                reader,
                groups,
                pending: record,
            });
        }
        let value = value.unwrap();
        let header = record.header();
        let weight = header
            .len()
            .checked_mul(40)
            .and_then(|n| n.checked_add(value.len().checked_mul(2)?))
            .and_then(|n| n.checked_add(record.bases.len().checked_mul(4)?))
            .and_then(|n| n.checked_add(4096))
            .ok_or_else(|| invalid("grouping: memory arithmetic overflow"))?;
        let next = used
            .checked_add(weight)
            .ok_or_else(|| invalid("grouping: memory arithmetic overflow"))?;
        // Defer all group and duplicate errors to the external path, which
        // checks complete input framing first, then duplicates, then groups.
        let existing = lookup.get(value).copied();
        let group_invalid = existing.is_some_and(|i| {
            groups[i].records.len() >= max_group_count
                || groups[i]
                    .records
                    .first()
                    .is_some_and(|r| r.bases.len() != record.bases.len())
                || group_weights[i]
                    .checked_add(weight)
                    .is_none_or(|n| n > max_group_bytes)
        });
        crate::publication::inject(failure, crate::publication::Failure::DuplicateDetection)?;
        crate::publication::inject(failure, crate::publication::Failure::Group)?;
        if next > cap
            || headers.contains(&header)
            || group_invalid
            || (existing.is_none() && weight > max_group_bytes)
        {
            return Ok(GroupAttempt::Spill {
                reader,
                groups,
                pending: record,
            });
        }
        used = next;
        headers.insert(header);
        let i = if let Some(i) = existing {
            i
        } else {
            let i = groups.len();
            lookup.insert(value.clone(), i);
            groups.push(Group {
                barcode: value.clone(),
                records: Vec::new(),
            });
            group_weights.push(0);
            i
        };
        group_weights[i] = group_weights[i]
            .checked_add(weight)
            .ok_or_else(|| invalid("grouping: group weight overflow"))?;
        groups[i].records.push(record);
    }
    Ok(GroupAttempt::Complete(groups))
}
