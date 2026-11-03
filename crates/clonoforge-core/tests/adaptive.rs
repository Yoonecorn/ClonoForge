use clonoforge_core::{bounded_workflow, consensus::Config, external::Limits};
use std::{
    fs,
    path::PathBuf,
    sync::atomic::{AtomicUsize, Ordering},
};
static NEXT: AtomicUsize = AtomicUsize::new(0);
struct Temp(PathBuf);
impl Temp {
    fn new() -> Self {
        let p = std::env::temp_dir().join(format!(
            "cf-adaptive-{}-{}",
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
fn limits(chunk: usize) -> Limits {
    Limits {
        chunk_bytes: chunk,
        max_record_bytes: 256,
        fan_in: 2,
    }
}
#[test]
fn no_spill_and_transition_match_external_for_pass_fail_and_workers() {
    let t = Temp::new();
    let input = t.0.join("input");
    let data = (0..42)
        .map(|i| {
            format!(
                "@r{i}|BARCODE={}\nAC\n+\nII\n",
                if i < 40 {
                    format!("G{}", i % 10)
                } else {
                    format!("S{i}")
                }
            )
        })
        .collect::<String>();
    fs::write(&input, data).unwrap();
    let cfg = Config {
        min_count: 2,
        failed: true,
        ..Config::default()
    };
    for workers in [1, 2, 8] {
        let exp = t.0.join(format!("exp{workers}"));
        let mem = t.0.join(format!("mem{workers}"));
        let spill = t.0.join(format!("spill{workers}"));
        let big = limits(4_000_000);
        let small = limits(10_000);
        bounded_workflow::run_configured(
            &input,
            &cfg,
            &exp,
            "probe",
            &t.0,
            (&big, 100_000, workers, 8, 50),
            None,
        )
        .unwrap();
        let (_, spilled) = bounded_workflow::run_adaptive(
            &input,
            &cfg,
            &mem,
            "probe",
            &t.0,
            (&big, 100_000, workers, 8, 50),
        )
        .unwrap();
        assert!(!spilled);
        let (_, spilled) = bounded_workflow::run_adaptive(
            &input,
            &cfg,
            &spill,
            "probe",
            &t.0,
            (&small, 100_000, workers, 8, 50),
        )
        .unwrap();
        assert!(spilled);
        for suffix in [
            "consensus-pass.fastq",
            "consensus-fail.fastq",
            "completion.json",
        ] {
            let p = format!("probe_{suffix}");
            assert_eq!(
                fs::read(exp.join(&p)).unwrap(),
                fs::read(mem.join(&p)).unwrap()
            );
            assert_eq!(
                fs::read(exp.join(&p)).unwrap(),
                fs::read(spill.join(&p)).unwrap()
            );
        }
    }
}
#[test]
fn single_pass_transition_preserves_external_error_precedence() {
    let t = Temp::new();
    let cfg = Config::default();
    let big = limits(4_000_000);
    for (i, (data, kind)) in [
        (
            "@r|BARCODE=A\nA\n+\nI\n@r|BARCODE=A\nA\n+\nI\n@s|BARCODE=B\nA\n+\n",
            "truncated",
        ),
        (
            "@r|BARCODE=A\nA\n+\nI\n@s|BARCODE=A\nAA\n+\nII\n@r|BARCODE=A\nA\n+\nI\n",
            "duplicate complete header",
        ),
    ]
    .into_iter()
    .enumerate()
    {
        let input = t.0.join(format!("in{i}"));
        fs::write(&input, data).unwrap();
        let out = t.0.join(format!("out{i}"));
        let err = bounded_workflow::run_adaptive(
            &input,
            &cfg,
            &out,
            "probe",
            &t.0,
            (&big, 100_000, 2, 4, 50),
        )
        .unwrap_err();
        assert!(err.to_string().contains(kind), "{err}");
        assert!(!out.join("probe_completion.json").exists());
    }
}
#[test]
fn no_spill_transaction_failures_cancel_and_retry() {
    use clonoforge_core::publication::Failure;
    let failures = [
        Failure::Admission,
        Failure::Parse,
        Failure::DuplicateDetection,
        Failure::Group,
        Failure::Worker,
        Failure::WorkerEarlyExit,
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
    for failure in failures {
        let t = Temp::new();
        let input = t.0.join("in");
        let out = t.0.join("out");
        fs::create_dir(&out).unwrap();
        fs::write(out.join("user"), b"keep").unwrap();
        let mut data = String::new();
        for i in 0..7 {
            data.push_str(&format!(
                "@r{i}|BARCODE={}\nA\n+\nI\n",
                if i == 0 {
                    "S"
                } else if i % 2 == 0 {
                    "A"
                } else {
                    "B"
                }
            ));
        }
        fs::write(&input, data).unwrap();
        let lim = limits(1_000_000);
        let cfg = Config {
            min_count: 2,
            failed: true,
            ..Config::default()
        };
        let resource = (&lim, 100_000, 2, 4, 50);
        let mut timing = bounded_workflow::StageTimings::default();
        let result = bounded_workflow::run_adaptive_profiled(
            &input,
            &cfg,
            &out,
            "probe",
            &t.0,
            resource,
            Some(failure),
            Some(&mut timing),
        );
        assert!(result.is_err(), "{failure:?}");
        assert!(!out.join("probe_completion.json").exists());
        assert!(!out.join("probe_consensus-pass.fastq").exists());
        assert!(!out.join("probe_consensus-fail.fastq").exists());
        assert_eq!(fs::read(out.join("user")).unwrap(), b"keep");
        let (_, spilled) =
            bounded_workflow::run_adaptive(&input, &cfg, &out, "probe", &t.0, resource).unwrap();
        assert!(!spilled);
        assert!(out.join("probe_completion.json").exists());
    }
}
#[test]
fn no_spill_stage_cleanup_failure_prevents_marker_and_allows_retry() {
    use clonoforge_core::publication::Failure;
    let t = Temp::new();
    let input = t.0.join("in");
    let out = t.0.join("out");
    fs::create_dir(&out).unwrap();
    fs::write(out.join("user"), b"keep").unwrap();
    fs::write(&input, b"@r|BARCODE=A\nA\n+\nI\n").unwrap();
    let lim = limits(1_000_000);
    let cfg = Config::default();
    let mut timing = bounded_workflow::StageTimings::default();
    let err = bounded_workflow::run_adaptive_profiled(
        &input,
        &cfg,
        &out,
        "probe",
        &t.0,
        (&lim, 100_000, 1, 2, 50),
        Some(Failure::StageCleanup),
        Some(&mut timing),
    )
    .unwrap_err();
    assert!(err.to_string().contains("injected stage cleanup failure"));
    assert!(!out.join("probe_completion.json").exists());
    assert!(!out.join("probe_consensus-pass.fastq").exists());
    assert_eq!(fs::read(out.join("user")).unwrap(), b"keep");
    let (_, spilled) = bounded_workflow::run_adaptive(
        &input,
        &cfg,
        &out,
        "probe",
        &t.0,
        (&lim, 100_000, 1, 2, 50),
    )
    .unwrap();
    assert!(!spilled);
    assert!(out.join("probe_completion.json").exists());
}
#[test]
fn seeded_spill_failure_hooks_do_not_publish_and_allow_retry() {
    use clonoforge_core::publication::Failure;
    for failure in [
        Failure::SpillCreate,
        Failure::SpillWrite,
        Failure::SpillFlush,
        Failure::SpillSync,
        Failure::SpillValidate,
        Failure::Merge,
        Failure::MergePartial,
    ] {
        let t = Temp::new();
        let input = t.0.join("in");
        let out = t.0.join("out");
        let data = (0..18)
            .map(|i| format!("@r{i}|BARCODE=G{}\nA\n+\nI\n", i % 6))
            .collect::<String>();
        fs::write(&input, data).unwrap();
        let lim = limits(10_000);
        let cfg = Config::default();
        let r = (&lim, 100_000, 2, 4, 50);
        let mut timing = bounded_workflow::StageTimings::default();
        assert!(
            bounded_workflow::run_adaptive_profiled(
                &input,
                &cfg,
                &out,
                "probe",
                &t.0,
                r,
                Some(failure),
                Some(&mut timing)
            )
            .is_err(),
            "{failure:?}"
        );
        assert!(!out.join("probe_completion.json").exists());
        let (_, spilled) =
            bounded_workflow::run_adaptive(&input, &cfg, &out, "probe", &t.0, r).unwrap();
        assert!(spilled);
        assert!(out.join("probe_completion.json").exists());
    }
}
#[test]
fn seeded_spill_cleanup_failure_prevents_marker_and_allows_retry() {
    use clonoforge_core::publication::Failure;
    let t = Temp::new();
    let input = t.0.join("in");
    let out = t.0.join("out");
    fs::create_dir(&out).unwrap();
    fs::write(out.join("user"), b"keep").unwrap();
    let data = (0..18)
        .map(|i| format!("@r{i}|BARCODE=G{}\nA\n+\nI\n", i % 6))
        .collect::<String>();
    fs::write(&input, data).unwrap();
    let lim = limits(10_000);
    let cfg = Config::default();
    let resource = (&lim, 100_000, 2, 4, 50);
    let mut timing = bounded_workflow::StageTimings::default();
    let e = bounded_workflow::run_adaptive_profiled(
        &input,
        &cfg,
        &out,
        "probe",
        &t.0,
        resource,
        Some(Failure::Cleanup),
        Some(&mut timing),
    )
    .unwrap_err();
    assert!(e.to_string().contains("injected Cleanup"));
    assert!(!out.join("probe_completion.json").exists());
    assert!(!out.join("probe_consensus-pass.fastq").exists());
    assert_eq!(fs::read(out.join("user")).unwrap(), b"keep");
    let (_, spilled) =
        bounded_workflow::run_adaptive(&input, &cfg, &out, "probe", &t.0, resource).unwrap();
    assert!(spilled);
    assert!(out.join("probe_completion.json").exists());
}
