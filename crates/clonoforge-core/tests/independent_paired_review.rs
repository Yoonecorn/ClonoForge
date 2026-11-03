//! Independent review cases prepared from molecule geometry and hand calculations.
use clonoforge_core::{
    assembly::{assemble_suffix_prefix, AssemblyConfig, AssemblyDisposition},
    primer_ext::{trim_iupac_anchored, Anchor},
};
#[test]
fn independently_calculated_molecule_and_likelihoods() {
    // Forward molecule GACGTT: R1 GACG, raw R2 AACGT -> oriented ACGTT.
    // Unique 3-base overlap. Q10/Q10 gives Q23; Q40/Q10 gives Q54;
    // Q0 is uninformative, leaving the second observation's Q30.
    let cfg = AssemblyConfig {
        min_overlap: 3,
        max_mismatches: 0,
        max_read_length: 8,
        quality_cap: 60,
    };
    let AssemblyDisposition::Assembled(value) = assemble_suffix_prefix(
        b"GACG",
        &[7, 10, 40, 0],
        b"AACGT",
        &[13, 12, 30, 10, 10],
        &cfg,
    )
    .unwrap() else {
        panic!("unique overlap rejected")
    };
    assert_eq!(value.bases, b"GACGTT");
    assert_eq!(value.phred, [7, 23, 54, 30, 12, 13]);
    assert_eq!((value.overlap, value.mismatches), (3, 0));
    let mut bases = value.bases;
    let mut quality = value.phred;
    // YY means two pyrimidines, here TT; only terminal qualities removed.
    assert!(trim_iupac_anchored(&mut bases, &mut quality, b"YY", Anchor::ThreePrime).unwrap());
    assert_eq!(bases, b"GACG");
    assert_eq!(quality, [7, 23, 54, 30]);
}
#[test]
fn exhaustive_short_geometry_matches_independent_hamming_enumerator() {
    fn dna(mut code: usize, n: usize) -> Vec<u8> {
        (0..n)
            .map(|_| {
                let b = b"ACGT"[code % 4];
                code /= 4;
                b
            })
            .collect()
    }
    for a in 0..256 {
        for b in 0..256 {
            let left = dna(a, 4);
            let oriented = dna(b, 4);
            let raw: Vec<u8> = oriented
                .iter()
                .rev()
                .map(|x| match x {
                    b'A' => b'T',
                    b'C' => b'G',
                    b'G' => b'C',
                    _ => b'A',
                })
                .collect();
            let offsets: Vec<usize> = (0..=1)
                .filter(|shift| left[*shift..] == oriented[..4 - *shift])
                .collect();
            let got = assemble_suffix_prefix(
                &left,
                &[30; 4],
                &raw,
                &[30; 4],
                &AssemblyConfig {
                    min_overlap: 3,
                    max_mismatches: 0,
                    max_read_length: 4,
                    quality_cap: 60,
                },
            )
            .unwrap();
            match offsets.len() {
                0 => assert_eq!(got, AssemblyDisposition::NoOverlap),
                1 => {
                    let AssemblyDisposition::Assembled(x) = got else {
                        panic!("missed unique placement")
                    };
                    let mut expected = left[..offsets[0]].to_vec();
                    expected.extend_from_slice(&oriented);
                    assert_eq!(x.bases, expected);
                    assert_eq!(x.overlap, 4 - offsets[0]);
                }
                _ => assert_eq!(got, AssemblyDisposition::AmbiguousOverlap),
            }
        }
    }
}
