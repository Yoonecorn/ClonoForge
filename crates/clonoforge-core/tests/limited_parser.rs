use clonoforge_core::io::{Format, Reader};
use std::io::Cursor;
#[test]
fn exact_record_limit_and_legacy_default() {
    let data = b"@r|BARCODE=A\nAC\n+\nII\n";
    let mut unrestricted = Reader::new(Cursor::new(data), Format::Fastq);
    let expected = unrestricted.next_record().unwrap().unwrap();
    let len = expected.header().len() + 2 * expected.bases.len();
    let mut bounded = Reader::fastq_limited(Cursor::new(data), len).unwrap();
    assert_eq!(bounded.next_record().unwrap().unwrap(), expected);
    assert!(bounded.next_record().unwrap().is_none());
    let mut too_small = Reader::fastq_limited(Cursor::new(data), len - 1).unwrap();
    assert!(too_small.next_record().is_err());
}
#[test]
fn lines_and_record_bound_before_header_expansion() {
    let giant = format!(
        "@r|BARCODE=A|{}\nA\n+\nI\n",
        (0..10000)
            .map(|i| format!("K{i}=X"))
            .collect::<Vec<_>>()
            .join("|")
    );
    let mut bounded = Reader::fastq_limited(Cursor::new(giant.as_bytes()), 64).unwrap();
    assert!(bounded
        .next_record()
        .unwrap_err()
        .to_string()
        .contains("line exceeds"));
    assert!(Reader::fastq_limited(Cursor::new(b""), 0).is_err());
    assert!(Reader::fastq_limited(Cursor::new(b""), 16 * 1024 * 1024 + 1).is_err());
}
#[test]
fn crlf_and_bad_quality_handling() {
    let data = b"@r|BARCODE=A\r\nAC\r\n+\r\nII\r\n";
    let mut r = Reader::fastq_limited(Cursor::new(data), 15).unwrap();
    assert_eq!(r.next_record().unwrap().unwrap().bases, b"AC");
    let mut bad = Reader::fastq_limited(Cursor::new(b"@r|BARCODE=A\nAC\n+\nI\n"), 100).unwrap();
    assert!(bad.next_record().is_err());
}
