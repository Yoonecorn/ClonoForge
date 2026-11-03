//! Independently generated Python Decimal finite reference, not production-derived.
use clonoforge_core::assembly::{assemble_suffix_prefix, AssemblyConfig, AssemblyDisposition};
#[test]
fn independent_decimal_reference() {
    let reference: serde_json::Value = serde_json::from_str(include_str!(
        "../../../tests/fixtures/science/finite-reference.json"
    ))
    .unwrap();
    let mut assembled = 0;
    for (index, f) in reference["fixtures"].as_array().unwrap().iter().enumerate() {
        let bytes = |k: &str| f[k].as_str().unwrap().as_bytes();
        let q = |k: &str| {
            f[k].as_array()
                .unwrap()
                .iter()
                .map(|x| x.as_u64().unwrap() as u8)
                .collect::<Vec<_>>()
        };
        let cfg = AssemblyConfig {
            min_overlap: f["min_overlap"].as_u64().unwrap() as usize,
            max_mismatches: f["max_mismatches"].as_u64().unwrap() as usize,
            ..Default::default()
        };
        let actual =
            assemble_suffix_prefix(bytes("left"), &q("lq"), bytes("right"), &q("rq"), &cfg)
                .unwrap();
        let e = &f["expected"];
        match e["disposition"].as_str().unwrap() {
            "none" => assert_eq!(actual, AssemblyDisposition::NoOverlap, "case {index}"),
            "ambiguous" => assert_eq!(
                actual,
                AssemblyDisposition::AmbiguousOverlap,
                "case {index}"
            ),
            "assembled" => {
                assembled += 1;
                let AssemblyDisposition::Assembled(a) = actual else {
                    panic!("case {index}: missing assembly")
                };
                assert_eq!(
                    a.bases,
                    e["bases"].as_str().unwrap().as_bytes(),
                    "case {index}"
                );
                let eq = e["phred"]
                    .as_array()
                    .unwrap()
                    .iter()
                    .map(|x| x.as_u64().unwrap() as u8)
                    .collect::<Vec<_>>();
                assert_eq!(a.phred, eq, "quality case {index}");
                assert_eq!(
                    a.overlap,
                    e["overlap"].as_u64().unwrap() as usize,
                    "case {index}"
                );
                assert_eq!(
                    a.mismatches,
                    e["mismatches"].as_u64().unwrap() as usize,
                    "case {index}"
                );
            }
            _ => panic!(),
        }
    }
    assert!(assembled > 2000);
}
