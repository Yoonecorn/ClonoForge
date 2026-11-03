//! Structural comparison; record identifiers must be unique.
use crate::io::{Format, Reader};
use crate::model::ReadRecord;
use std::collections::BTreeMap;
use std::io::BufRead;
#[derive(Debug, PartialEq, Eq)]
pub struct Difference {
    pub category: &'static str,
    pub id: String,
}
pub fn compare<A: BufRead, B: BufRead>(a: A, b: B, format: Format) -> Vec<Difference> {
    fn load<R: BufRead>(r: R, f: Format, differences: &mut Vec<Difference>) -> Vec<ReadRecord> {
        let mut reader = Reader::new(r, f);
        let mut records = Vec::new();
        loop {
            match reader.next_record() {
                Ok(Some(r)) => records.push(r),
                Ok(None) => break,
                Err(e) => {
                    differences.push(Difference {
                        category: "parse",
                        id: e.to_string(),
                    });
                    break;
                }
            }
        }
        records
    }
    let mut differences = Vec::new();
    let a = load(a, format, &mut differences);
    let b = load(b, format, &mut differences);
    if a.iter().map(|r| &r.id).collect::<Vec<_>>() != b.iter().map(|r| &r.id).collect::<Vec<_>>() {
        differences.push(Difference {
            category: "record-order",
            id: String::new(),
        });
    }
    fn index<'a>(
        records: &'a [ReadRecord],
        differences: &mut Vec<Difference>,
    ) -> BTreeMap<&'a str, &'a ReadRecord> {
        let mut map = BTreeMap::new();
        for r in records {
            if map.insert(r.id.as_str(), r).is_some() {
                differences.push(Difference {
                    category: "duplicate-record",
                    id: r.id.clone(),
                });
            }
        }
        map
    }
    let left = index(&a, &mut differences);
    let right = index(&b, &mut differences);
    for id in left
        .keys()
        .chain(right.keys())
        .copied()
        .collect::<std::collections::BTreeSet<_>>()
    {
        match (left.get(id), right.get(id)) {
            (Some(a), Some(b)) => {
                for (different, category) in [
                    (a.header() != b.header(), "header"),
                    (a.annotations != b.annotations, "annotation"),
                    (a.bases != b.bases, "sequence"),
                    (a.phred != b.phred, "quality"),
                ] {
                    if different {
                        differences.push(Difference {
                            category,
                            id: id.to_owned(),
                        });
                    }
                }
            }
            _ => differences.push(Difference {
                category: "missing-record",
                id: id.to_owned(),
            }),
        }
    }
    differences
}
