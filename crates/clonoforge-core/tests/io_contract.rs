use clonoforge_core::{compare::compare, io::*, model::ReadRecord};
use std::io::{self, Cursor, Write};
fn records(bytes: &[u8], format: Format) -> io::Result<Vec<ReadRecord>> {
    let mut reader = Reader::new(Cursor::new(bytes), format);
    let mut out = Vec::new();
    while let Some(r) = reader.next_record()? {
        out.push(r);
    }
    Ok(out)
}
#[test]
fn generated_roundtrips() {
    for n in 1..128 {
        for f in [Format::Fasta, Format::Fastq] {
            let bases: Vec<u8> = (0..n).map(|i| b"ACGTRYn.-"[i % 9]).collect();
            let q = if f == Format::Fastq {
                Some((0..n).map(|i| (i % 94) as u8).collect())
            } else {
                None
            };
            let r = ReadRecord::from_header("read|BARCODE=ACGT|COUNT=2", bases, q, 0).unwrap();
            let mut bytes = Vec::new();
            write_record(&mut bytes, &r, f).unwrap();
            assert_eq!(records(&bytes, f).unwrap(), [r]);
        }
    }
}
#[test]
fn annotation_errors_and_order() {
    for h in [
        "",
        "a b",
        "r|A=",
        "r|=x",
        "r|A=x|A=y",
        "r|a=x",
        "r|A=x=y",
        "r|A",
        "r|",
    ] {
        assert!(
            ReadRecord::from_header(h, b"A".to_vec(), None, 0).is_err(),
            "{h}"
        );
    }
    let r = ReadRecord::from_header("r|Z=x|A=y", b"N".to_vec(), None, 0).unwrap();
    assert_eq!(r.header(), "r|Z=x|A=y");
    assert_eq!(r.annotations.iter().find(|(k, _)| k == "MISSING"), None);
}
#[test]
fn malformed() {
    for b in [
        "@r\nA\n+\n",
        "@r\nA\n",
        "@r\nA\nwrong\n!\n",
        "@r\nAA\n+\n!\n",
        "@r\nA\n+\n \n",
        "@r\nX\n+\n!\n",
        "\n",
        "@\nA\n+\n!\n",
    ] {
        assert!(records(b.as_bytes(), Format::Fastq).is_err(), "{b:?}");
    }
    for b in [">r\n", ">r\n\n", "A\n", ">r\nXX\n"] {
        assert!(records(b.as_bytes(), Format::Fasta).is_err());
    }
    assert!(records(&vec![b'A'; 16 * 1024 * 1024 + 1], Format::Fasta).is_err());
}
#[test]
fn multiline_crlf_ordinals() {
    let r = records(b">a\r\nAC\r\nGT\r\n>b\nN", Format::Fasta).unwrap();
    assert_eq!(r[0].bases, b"ACGT");
    assert_eq!(r[1].ordinal, 1);
    assert_eq!(
        records(b"@a\r\nA\r\n+a\r\n!", Format::Fastq).unwrap()[0].phred,
        Some(vec![0])
    );
}
#[test]
fn gzip_and_corruption() {
    let bytes = b"@r\nAC\n+\n!~\n";
    let mut gzip = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
    gzip.write_all(bytes).unwrap();
    let compressed = gzip.finish().unwrap();
    let path = std::env::temp_dir().join(format!("clonoforge-gzip-{}", std::process::id()));
    let mut two = compressed.clone();
    two.extend_from_slice(&compressed);
    std::fs::write(&path, two).unwrap();
    let mut reader = Reader::new(input(&path).unwrap(), Format::Fastq);
    assert_eq!(reader.next_record().unwrap().unwrap().bases, b"AC");
    assert_eq!(reader.next_record().unwrap().unwrap().ordinal, 1);
    assert!(reader.next_record().unwrap().is_none());
    std::fs::write(&path, &compressed[..compressed.len() - 4]).unwrap();
    let mut reader = Reader::new(input(&path).unwrap(), Format::Fastq);
    let mut failed = false;
    for _ in 0..3 {
        if reader.next_record().is_err() {
            failed = true;
            break;
        }
    }
    assert!(failed);
    std::fs::remove_file(path).unwrap();
}
#[test]
fn comparator_categories() {
    let a = b"@a|A=x\nAC\n+\n!!\n@b\nA\n+\n!\n";
    assert!(compare(Cursor::new(a), Cursor::new(a), Format::Fastq).is_empty());
    for (b, category) in [
        ("@b\nA\n+\n!\n@a|A=x\nAC\n+\n!!\n", "record-order"),
        ("@a|A=y\nAC\n+\n!!\n", "annotation"),
        ("@a|A=x\nAT\n+\n!!\n", "sequence"),
        ("@a|A=x\nAC\n+\n!#\n", "quality"),
        ("@a|A=x\nAC\n+\n!!\n", "missing-record"),
        ("@a\nA\n+\n!\n@a\nA\n+\n!\n", "duplicate-record"),
        ("bad", "parse"),
    ] {
        assert!(compare(Cursor::new(a), Cursor::new(b), Format::Fastq)
            .iter()
            .any(|d| d.category == category));
    }
}
struct Failing;
impl Write for Failing {
    fn write(&mut self, _: &[u8]) -> io::Result<usize> {
        Err(io::Error::other("injected"))
    }
    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}
#[test]
fn io_faults() {
    let r = ReadRecord::from_header("r", b"A".to_vec(), Some(vec![0]), 0).unwrap();
    assert!(write_record(&mut Failing, &r, Format::Fastq).is_err());
    struct BadRead;
    impl io::Read for BadRead {
        fn read(&mut self, _: &mut [u8]) -> io::Result<usize> {
            Err(io::Error::other("injected"))
        }
    }
    assert!(Reader::new(io::BufReader::new(BadRead), Format::Fastq)
        .next_record()
        .is_err());
}
