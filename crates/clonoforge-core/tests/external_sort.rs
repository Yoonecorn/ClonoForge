use clonoforge_core::{
    consensus::Config,
    external::{self, Limits},
};
use std::{
    fs,
    io::{Seek, SeekFrom, Write},
    path::PathBuf,
    sync::atomic::{AtomicUsize, Ordering},
};
static NEXT: AtomicUsize = AtomicUsize::new(0);
struct Temp(PathBuf);
impl Temp {
    fn new() -> Self {
        let p = std::env::temp_dir().join(format!(
            "cf-ext-{}-{}",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        fs::create_dir(&p).unwrap();
        Self(p)
    }
}
impl Drop for Temp {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}
fn lim() -> Limits {
    Limits {
        chunk_bytes: 10000,
        max_record_bytes: 100,
        fan_in: 2,
    }
}
fn text(n: usize) -> String {
    (0..n)
        .map(|i| {
            format!(
                "@r{i}|BARCODE={}\nAC\n+\nII\n",
                if i % 2 == 0 { "B" } else { "A" }
            )
        })
        .collect()
}
#[test]
fn many_runs_and_group_order() {
    let t = Temp::new();
    let input = t.0.join("input");
    fs::write(&input, text(25)).unwrap();
    let run = external::sort_input(&input, &Config::default(), &t.0, &lim()).unwrap();
    let mut items = Vec::new();
    external::for_each_record(&run, &Config::default(), &lim(), |r| {
        items.push((r.annotations[0].1.clone(), r.ordinal));
        Ok(())
    })
    .unwrap();
    let mut expected = (0..25)
        .map(|i| (if i % 2 == 0 { "B" } else { "A" }.to_string(), i))
        .collect::<Vec<_>>();
    expected.sort();
    assert_eq!(items, expected);
}
#[test]
fn global_duplicates_across_runs_are_rejected() {
    let t = Temp::new();
    let input = t.0.join("input");
    let mut s = text(25);
    s.push_str("@r0|BARCODE=B\nAC\n+\nII\n");
    fs::write(&input, s).unwrap();
    assert!(
        external::sort_input(&input, &Config::default(), &t.0, &lim())
            .unwrap_err()
            .to_string()
            .contains("duplicate complete header")
    );
}
#[test]
fn corruption_and_truncation_are_rejected_before_callback() {
    for truncate in [true, false] {
        let t = Temp::new();
        let input = t.0.join("input");
        fs::write(&input, text(5)).unwrap();
        let run = external::sort_input(&input, &Config::default(), &t.0, &lim()).unwrap();
        if truncate {
            let f = fs::OpenOptions::new().write(true).open(&run).unwrap();
            f.set_len(22).unwrap()
        } else {
            let mut f = fs::OpenOptions::new().write(true).open(&run).unwrap();
            f.seek(SeekFrom::Start(18)).unwrap();
            f.write_all(&[9]).unwrap()
        }
        let mut seen = 0;
        assert!(
            external::for_each_record(&run, &Config::default(), &lim(), |_| {
                seen += 1;
                Ok(())
            })
            .is_err()
        );
        assert_eq!(seen, 0);
    }
}
#[test]
fn rejects_oversized_record() {
    let t = Temp::new();
    let input = t.0.join("input");
    fs::write(&input, text(1)).unwrap();
    let mut l = lim();
    l.chunk_bytes = 1;
    assert!(external::sort_input(&input, &Config::default(), &t.0, &l).is_err());
}

#[test]
fn input_record_limit_allows_consensus_header_expansion() {
    let t = Temp::new();
    let input = t.0.join("input");
    fs::write(&input, b"@r|BARCODE=A\nA\n+\nI\n").unwrap();
    let limits = Limits {
        chunk_bytes: 10000,
        max_record_bytes: 13,
        fan_in: 2,
    };
    let out = t.0.join("out");
    clonoforge_core::bounded_workflow::run(
        &input,
        &Config::default(),
        &out,
        "p",
        &t.0,
        (&limits, 10000, 1, 1),
        None,
    )
    .unwrap();
    assert_eq!(
        fs::read(out.join("p_consensus-pass.fastq")).unwrap(),
        b"@A|CONSCOUNT=1\nA\n+\nI\n"
    );
}

#[test]
fn spill_cleanup_failure_prevents_marker_and_allows_retry() {
    use clonoforge_core::{bounded_workflow, publication::Failure};
    let t = Temp::new();
    let input = t.0.join("input");
    let out = t.0.join("out");
    fs::create_dir(&out).unwrap();
    fs::write(out.join("user"), b"keep").unwrap();
    fs::write(&input, b"@r|BARCODE=A\nA\n+\nI\n").unwrap();
    let mut timing = bounded_workflow::StageTimings::default();
    let err = bounded_workflow::run_configured_profiled(
        &input,
        &Config::default(),
        &out,
        "p",
        &t.0,
        (&lim(), 10000, 1, 1, 50),
        Some(Failure::Cleanup),
        Some(&mut timing),
    )
    .unwrap_err();
    assert!(err.to_string().contains("injected Cleanup"));
    assert!(!out.join("p_completion.json").exists());
    assert!(!out.join("p_consensus-pass.fastq").exists());
    assert_eq!(fs::read(out.join("user")).unwrap(), b"keep");
    assert!(fs::read_dir(&t.0).unwrap().any(|e| e
        .unwrap()
        .file_name()
        .to_string_lossy()
        .starts_with(".clonoforge-p-spill-incomplete")));
    bounded_workflow::run_configured(
        &input,
        &Config::default(),
        &out,
        "p",
        &t.0,
        (&lim(), 10000, 1, 1, 50),
        None,
    )
    .unwrap();
    assert!(out.join("p_completion.json").exists());
}

#[test]
fn corrected_candidate_reference_through_external_workflow() {
    let t = Temp::new();
    let root = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../tests/fixtures/candidate");
    let cfg = Config {
        frequency: 0.0,
        failed: true,
        ..Config::default()
    };
    let limits = Limits {
        chunk_bytes: 2 * 1024 * 1024,
        max_record_bytes: 16 * 1024,
        fan_in: 8,
    };
    let out = t.0.join("published");
    clonoforge_core::bounded_workflow::run(
        &root.join("input.fastq"),
        &cfg,
        &out,
        "probe",
        &t.0,
        (&limits, 1024 * 1024, 2, 4),
        None,
    )
    .unwrap();
    assert_eq!(
        fs::read(out.join("probe_consensus-pass.fastq")).unwrap(),
        fs::read(root.join("ref/p_consensus-pass.fastq")).unwrap()
    );
    assert!(!out.join("probe_consensus-fail.fastq").exists());
}
#[test]
fn workflow_outputs_match_m2b_publisher_in_first_encounter_order() {
    use clonoforge_core::{bounded_workflow, group, publication};
    let t = Temp::new();
    let input = t.0.join("input.fastq");
    let data = (0..70)
        .map(|i| {
            format!(
                "@r{i}|BARCODE={}\nAC\n+\nII\n",
                match i % 7 {
                    0 => "Z",
                    1 => "A",
                    2 => "Y",
                    3 => "B",
                    4 => "X",
                    5 => "C",
                    _ => "W",
                }
            )
        })
        .collect::<String>();
    fs::write(&input, data).unwrap();
    let cfg = Config {
        min_count: 1,
        failed: true,
        ..Config::default()
    };
    let mem = group::load(&input, &cfg).unwrap();
    let expected = t.0.join("expected");
    let actual = t.0.join("actual");
    publication::publish(&mem, &cfg, &expected, "probe", None).unwrap();
    bounded_workflow::run(
        &input,
        &cfg,
        &actual,
        "probe",
        &t.0,
        (&lim(), 100_000, 1, 2),
        None,
    )
    .unwrap();
    for suffix in ["consensus-pass.fastq", "completion.json"] {
        assert_eq!(
            fs::read(expected.join(format!("probe_{suffix}"))).unwrap(),
            fs::read(actual.join(format!("probe_{suffix}"))).unwrap()
        );
    }
}
#[test]
fn workflow_failures_keep_existing_paths_and_emit_no_manifest() {
    use clonoforge_core::{bounded_workflow, publication::Failure};
    let t = Temp::new();
    let input = t.0.join("input");
    fs::write(&input, text(5)).unwrap();
    let cfg = Config::default();
    let out = t.0.join("out");
    assert!(bounded_workflow::run(
        &input,
        &cfg,
        &out,
        "probe",
        &t.0,
        (&lim(), 100_000, 1, 2),
        Some(Failure::ManifestPublish)
    )
    .is_err());
    assert!(!out.join("probe_completion.json").exists());
    assert!(!out.join("probe_consensus-pass.fastq").exists());
}
#[test]
fn mixed_fail_order_and_lazy_categories_match_reference_publisher() {
    use clonoforge_core::{bounded_workflow, group, publication};
    let t = Temp::new();
    let input = t.0.join("input");
    fs::write(&input,"@r0|BARCODE=Z\nAC\n+\nII\n@r1|BARCODE=B\nAC\n+\nII\n@r2|BARCODE=Y\nGT\n+\nII\n@r3|BARCODE=Z\nAC\n+\nII\n@r4|BARCODE=A\nGT\n+\nII\n@r5|BARCODE=Y\nGT\n+\nII\n").unwrap();
    for min in [1, 2, 3] {
        let cfg = Config {
            min_count: min,
            failed: true,
            ..Config::default()
        };
        let groups = group::load(&input, &cfg).unwrap();
        let exp = t.0.join(format!("exp{min}"));
        let got = t.0.join(format!("got{min}"));
        publication::publish(&groups, &cfg, &exp, "probe", None).unwrap();
        bounded_workflow::run(
            &input,
            &cfg,
            &got,
            "probe",
            &t.0,
            (&lim(), 100_000, 1, 2),
            None,
        )
        .unwrap();
        for suffix in [
            "consensus-pass.fastq",
            "consensus-fail.fastq",
            "completion.json",
        ] {
            let a = exp.join(format!("probe_{suffix}"));
            let b = got.join(format!("probe_{suffix}"));
            assert_eq!(a.exists(), b.exists());
            if a.exists() {
                assert_eq!(fs::read(a).unwrap(), fs::read(b).unwrap());
            }
        }
    }
}
#[test]
fn empty_workflow_has_only_completion_manifest() {
    use clonoforge_core::bounded_workflow;
    let t = Temp::new();
    let input = t.0.join("empty");
    fs::write(&input, b"").unwrap();
    let out = t.0.join("out");
    bounded_workflow::run(
        &input,
        &Config::default(),
        &out,
        "probe",
        &t.0,
        (&lim(), 100_000, 1, 2),
        None,
    )
    .unwrap();
    assert!(out.join("probe_completion.json").exists());
    assert!(!out.join("probe_consensus-pass.fastq").exists());
    assert!(!out.join("probe_consensus-fail.fastq").exists());
}
#[test]
fn bounded_parallel_outputs_match_for_1_2_8_workers() {
    use clonoforge_core::bounded_workflow;
    let t = Temp::new();
    let input = t.0.join("input");
    let data = (0..300)
        .map(|i| format!("@r{i}|BARCODE=G{}\nACGT\n+\nIIII\n", (i * 19) % 60))
        .collect::<String>();
    fs::write(&input, data).unwrap();
    let cfg = Config {
        min_count: 2,
        failed: true,
        ..Config::default()
    };
    let mut reference = None;
    for workers in [1, 2, 8] {
        for rep in 0..2 {
            let out = t.0.join(format!("out{workers}_{rep}"));
            bounded_workflow::run(
                &input,
                &cfg,
                &out,
                "probe",
                &t.0,
                (&lim(), 100_000, workers, 8),
                None,
            )
            .unwrap();
            let bytes = ["consensus-pass.fastq", "completion.json"]
                .map(|suffix| fs::read(out.join(format!("probe_{suffix}"))).unwrap());
            if let Some(r) = &reference {
                assert_eq!(&bytes, r);
            } else {
                reference = Some(bytes);
            }
        }
    }
}
#[test]
fn parallel_group_failure_never_publishes_success() {
    use clonoforge_core::bounded_workflow;
    let t = Temp::new();
    let input = t.0.join("input");
    let mut data = (0..100)
        .map(|i| format!("@r{i}|BARCODE=G{}\nA\n+\nI\n", i % 20))
        .collect::<String>();
    data.push_str("@bad|BARCODE=G5\nAA\n+\nII\n");
    fs::write(&input, data).unwrap();
    for workers in [1, 2, 8] {
        let out = t.0.join(format!("err{workers}"));
        let err = bounded_workflow::run(
            &input,
            &Config::default(),
            &out,
            "probe",
            &t.0,
            (&lim(), 100_000, workers, 8),
            None,
        )
        .unwrap_err();
        assert!(err.to_string().contains("unequal read lengths"));
        assert!(!out.join("probe_completion.json").exists());
    }
}

#[test]
fn injected_external_failures_preserve_destinations_and_allow_retry() {
    use clonoforge_core::{bounded_workflow, publication::Failure};
    let cases = [
        Failure::Admission,
        Failure::Parse,
        Failure::SpillCreate,
        Failure::SpillFlush,
        Failure::SpillValidate,
        Failure::Merge,
        Failure::DuplicateDetection,
        Failure::Group,
        Failure::Worker,
        Failure::OutputOrder,
        Failure::Stage,
        Failure::WritePass,
        Failure::WriteFail,
        Failure::FlushPass,
        Failure::FlushFail,
        Failure::ValidatePass,
        Failure::ValidateFail,
        Failure::ManifestWrite,
        Failure::ManifestFlush,
        Failure::ManifestValidate,
        Failure::PublishPass,
        Failure::PublishFail,
        Failure::DirectorySync,
        Failure::ManifestPublish,
    ];
    for (i, failure) in cases.into_iter().enumerate() {
        let t = Temp::new();
        let input = t.0.join("input");
        let out = t.0.join("out");
        fs::create_dir(&out).unwrap();
        let mut data = String::new();
        for x in 0..12 {
            data.push_str(&format!(
                "@r{x}|BARCODE={}\nA\n+\nI\n",
                if x == 0 { "B" } else { "A" }
            ));
        }
        fs::write(&input, data).unwrap();
        let sentinel = out.join("user");
        fs::write(&sentinel, b"sentinel").unwrap();
        let cfg = Config {
            min_count: 2,
            failed: true,
            ..Config::default()
        };
        let result = bounded_workflow::run_configured(
            &input,
            &cfg,
            &out,
            "probe",
            &t.0,
            (&lim(), 100_000, 2, 4, 50),
            Some(failure),
        );
        assert!(result.is_err(), "{i}: {failure:?}");
        assert!(!out.join("probe_completion.json").exists());
        assert!(!out.join("probe_consensus-pass.fastq").exists());
        assert_eq!(fs::read(&sentinel).unwrap(), b"sentinel");
        bounded_workflow::run_configured(
            &input,
            &cfg,
            &out,
            "probe",
            &t.0,
            (&lim(), 100_000, 2, 4, 50),
            None,
        )
        .unwrap();
        assert!(out.join("probe_completion.json").exists());
    }
}
#[test]
fn configured_group_count_is_enforced_before_success() {
    use clonoforge_core::bounded_workflow;
    let t = Temp::new();
    let input = t.0.join("input");
    fs::write(&input, text(5)).unwrap();
    let out = t.0.join("out");
    assert!(bounded_workflow::run_configured(
        &input,
        &Config::default(),
        &out,
        "probe",
        &t.0,
        (&lim(), 100_000, 2, 4, 2),
        None
    )
    .unwrap_err()
    .to_string()
    .contains("group exceeds configured"));
    assert!(!out.join("probe_completion.json").exists());
}

#[test]
fn actual_parser_and_decompressor_failures_do_not_publish() {
    use clonoforge_core::bounded_workflow;
    use flate2::{write::GzEncoder, Compression};
    let t = Temp::new();
    let mut gz = GzEncoder::new(Vec::new(), Compression::default());
    gz.write_all(b"@r|BARCODE=A\nA\n+\nI\n").unwrap();
    let mut bad_gz = gz.finish().unwrap();
    bad_gz.truncate(bad_gz.len() - 4);
    for (i, bytes) in [b"@r|BARCODE=A\nA\n+\n".to_vec(), bad_gz]
        .into_iter()
        .enumerate()
    {
        let input = t.0.join(format!("bad{i}"));
        let out = t.0.join(format!("out{i}"));
        fs::write(&input, bytes).unwrap();
        assert!(bounded_workflow::run_configured(
            &input,
            &Config::default(),
            &out,
            "probe",
            &t.0,
            (&lim(), 100_000, 2, 4, 50),
            None
        )
        .is_err());
        assert!(!out.join("probe_completion.json").exists());
    }
}

#[test]
fn partial_spill_write_flush_sync_and_merge_failures_are_not_success() {
    use clonoforge_core::{bounded_workflow, publication::Failure};
    use sha2::{Digest, Sha256};
    for failure in [
        Failure::SpillWrite,
        Failure::SpillFlush,
        Failure::SpillSync,
        Failure::MergePartial,
    ] {
        let t = Temp::new();
        let input = t.0.join("in");
        let out = t.0.join("out");
        fs::create_dir(&out).unwrap();
        fs::write(&input, text(12)).unwrap();
        let sentinel = t.0.join("unrelated");
        fs::write(&sentinel, b"keep").unwrap();
        let cfg = Config::default();
        let err = bounded_workflow::run_configured(
            &input,
            &cfg,
            &out,
            "probe",
            &t.0,
            (&lim(), 100_000, 2, 4, 50),
            Some(failure),
        )
        .unwrap_err();
        assert!(err.to_string().contains("spill path"));
        assert!(!out.join("probe_completion.json").exists());
        assert!(!out.join("probe_consensus-pass.fastq").exists());
        assert_eq!(fs::read(&sentinel).unwrap(), b"keep");
        let scratch = fs::read_dir(&t.0)
            .unwrap()
            .map(|e| e.unwrap().path())
            .find(|p| {
                p.file_name()
                    .unwrap()
                    .to_string_lossy()
                    .starts_with(".clonoforge-probe-spill-incomplete-")
            })
            .unwrap();
        let runs: Vec<_> = fs::read_dir(&scratch)
            .unwrap()
            .map(|e| e.unwrap().path())
            .filter(|p| p.extension().is_some_and(|e| e == "bin"))
            .collect();
        assert!(!runs.is_empty());
        if matches!(failure, Failure::SpillWrite | Failure::MergePartial) {
            let partial = runs.iter().any(|p| {
                let data = fs::read(p).unwrap();
                data.len() < 32
                    || data[data.len() - 32..] != Sha256::digest(&data[..data.len() - 32])[..]
            });
            assert!(partial, "no incomplete run detected for {failure:?}");
        }
        bounded_workflow::run_configured(
            &input,
            &cfg,
            &out,
            "probe",
            &t.0,
            (&lim(), 100_000, 2, 4, 50),
            None,
        )
        .unwrap();
        assert!(out.join("probe_completion.json").exists());
    }
}
