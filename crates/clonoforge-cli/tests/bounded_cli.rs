use std::{fs, process::Command};
#[test]
fn thread_parity_and_admission_rejection() {
    let root = std::env::temp_dir().join(format!("clonoforge-bounded-cli-{}", std::process::id()));
    fs::create_dir(&root).unwrap();
    let input = root.join("input.fastq");
    let content = (0..120)
        .map(|i| format!("@r{i}|BARCODE=B{}\nACGT\n+\nIIII\n", i % 17))
        .collect::<String>();
    fs::write(&input, content).unwrap();
    let mut expected = None;
    for (index, workers) in [1, 2, 8, 2].iter().enumerate() {
        let out = root.join(format!("out{index}"));
        let result = Command::new(env!("CARGO_BIN_EXE_clonoforge"))
            .arg("build-consensus")
            .arg(&input)
            .arg("--outdir")
            .arg(&out)
            .args([
                "--outname",
                "probe",
                "--memory-mib",
                "2048",
                "--threads",
                &workers.to_string(),
            ])
            .output()
            .unwrap();
        assert!(
            result.status.success(),
            "{}",
            String::from_utf8_lossy(&result.stderr)
        );
        assert!(String::from_utf8_lossy(&result.stderr).contains(&format!("workers {workers}")));
        let got = (
            fs::read(out.join("probe_consensus-pass.fastq")).unwrap(),
            fs::read(out.join("probe_completion.json")).unwrap(),
        );
        if let Some(ref prior) = expected {
            assert_eq!(prior, &got);
        } else {
            expected = Some(got);
        }
    }
    let out = root.join("rejected");
    let result = Command::new(env!("CARGO_BIN_EXE_clonoforge"))
        .arg("build-consensus")
        .arg(&input)
        .arg("--outdir")
        .arg(&out)
        .args(["--memory-mib", "1"])
        .output()
        .unwrap();
    assert_eq!(result.status.code(), Some(2));
    assert!(!out.join("input_completion.json").exists());
    let limited = root.join("group-limited");
    let result = Command::new(env!("CARGO_BIN_EXE_clonoforge"))
        .arg("build-consensus")
        .arg(&input)
        .arg("--outdir")
        .arg(&limited)
        .args(["--max-group-count", "2"])
        .output()
        .unwrap();
    assert_eq!(result.status.code(), Some(2));
    assert!(!limited.join("input_completion.json").exists());
    let nested = root.join("new-parent").join("nested-output");
    let result = Command::new(env!("CARGO_BIN_EXE_clonoforge"))
        .arg("build-consensus")
        .arg(&input)
        .arg("--outdir")
        .arg(&nested)
        .output()
        .unwrap();
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    assert!(nested.join("input_completion.json").exists());
    fs::remove_dir_all(root).unwrap();
}

#[test]
fn plain_gzip_equivalence_and_truncated_stream() {
    let root = std::env::temp_dir().join(format!("clonoforge-gzip-cli-{}", std::process::id()));
    fs::create_dir(&root).unwrap();
    let fixture = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/fixtures/io");
    let mut expected = None;
    for (i, name) in ["valid.fastq", "valid.fastq.gz"].iter().enumerate() {
        let out = root.join(format!("out{i}"));
        let result = Command::new(env!("CARGO_BIN_EXE_clonoforge"))
            .arg("build-consensus")
            .arg(fixture.join(name))
            .arg("--outdir")
            .arg(&out)
            .args(["--outname", "probe"])
            .output()
            .unwrap();
        assert!(
            result.status.success(),
            "{}",
            String::from_utf8_lossy(&result.stderr)
        );
        let receipt: serde_json::Value = serde_json::from_slice(&result.stdout).unwrap();
        assert_eq!(receipt["input_mode"], if i == 0 { "plain" } else { "gzip" });
        assert_eq!(receipt["success"], true);
        let got = (
            fs::read(out.join("probe_consensus-pass.fastq")).unwrap(),
            fs::read(out.join("probe_completion.json")).unwrap(),
        );
        if let Some(ref prior) = expected {
            assert_eq!(prior, &got);
        } else {
            expected = Some(got);
        }
    }
    let bytes = fs::read(fixture.join("valid.fastq.gz")).unwrap();
    let truncated = root.join("truncated.fastq.gz");
    fs::write(&truncated, &bytes[..bytes.len() - 4]).unwrap();
    let out = root.join("bad");
    let result = Command::new(env!("CARGO_BIN_EXE_clonoforge"))
        .arg("build-consensus")
        .arg(&truncated)
        .arg("--outdir")
        .arg(&out)
        .args(["--outname", "probe"])
        .output()
        .unwrap();
    assert_eq!(result.status.code(), Some(2));
    assert!(!out.join("probe_completion.json").exists());
    fs::remove_dir_all(root).unwrap();
}

#[test]
fn profiling_preserves_scientific_and_manifest_bytes() {
    let root = std::env::temp_dir().join(format!("clonoforge-profile-cli-{}", std::process::id()));
    fs::create_dir(&root).unwrap();
    let input = root.join("input.fastq");
    fs::write(&input, b"@r|BARCODE=A\nA\n+\nI\n").unwrap();
    let mut expected = None;
    for profile in [false, true] {
        let out = root.join(if profile { "profile" } else { "plain" });
        let mut command = Command::new(env!("CARGO_BIN_EXE_clonoforge"));
        command
            .arg("build-consensus")
            .arg(&input)
            .arg("--outdir")
            .arg(&out)
            .args(["--outname", "probe"]);
        if profile {
            command.env("CLONOFORGE_PROFILE", "1");
        } else {
            command.env_remove("CLONOFORGE_PROFILE");
        }
        let result = command.output().unwrap();
        assert!(
            result.status.success(),
            "{}",
            String::from_utf8_lossy(&result.stderr)
        );
        let receipt: serde_json::Value = serde_json::from_slice(&result.stdout).unwrap();
        assert_eq!(receipt["profile"].is_object(), profile);
        assert_eq!(receipt["spilled"], false);
        if profile {
            assert_eq!(receipt["profile"]["execution_path"], "no-spill");
            assert_eq!(
                receipt["profile"]["raw_bytes"],
                fs::metadata(&input).unwrap().len()
            );
            assert_eq!(
                receipt["profile"]["decoded_bytes"],
                fs::metadata(&input).unwrap().len()
            );
            assert_eq!(receipt["profile"]["consensus_groups"], 1);
            assert!(
                receipt["profile"]["record_read_parse_inclusive_ns"]
                    .as_u64()
                    .unwrap()
                    > 0
            );
            assert!(
                receipt["profile"]["consensus_group_elapsed_sum_ns"]
                    .as_u64()
                    .unwrap()
                    > 0
            );
            assert!(receipt["profile"]["group_pipeline_ns"].is_null());
        }
        let got = (
            fs::read(out.join("probe_consensus-pass.fastq")).unwrap(),
            fs::read(out.join("probe_completion.json")).unwrap(),
        );
        if let Some(ref prior) = expected {
            assert_eq!(prior, &got);
        } else {
            expected = Some(got);
        }
    }
    fs::remove_dir_all(root).unwrap();
}

#[test]
fn closed_stdout_does_not_turn_published_success_into_failure() {
    let root = std::env::temp_dir().join(format!("clonoforge-pipe-cli-{}", std::process::id()));
    fs::create_dir(&root).unwrap();
    let input = root.join("input.fastq");
    let out = root.join("out");
    fs::write(&input, b"@r|BARCODE=A\nA\n+\nI\n").unwrap();
    let mut child = Command::new(env!("CARGO_BIN_EXE_clonoforge"))
        .arg("build-consensus")
        .arg(&input)
        .arg("--outdir")
        .arg(&out)
        .stdout(std::process::Stdio::piped())
        .spawn()
        .unwrap();
    drop(child.stdout.take());
    assert!(child.wait().unwrap().success());
    assert!(out.join("input_completion.json").exists());
    fs::remove_dir_all(root).unwrap();
}

#[test]
fn profiling_and_force_spill_report_actual_paths() {
    let root = std::env::temp_dir().join(format!("cf-path-profile-{}", std::process::id()));
    fs::create_dir(&root).unwrap();
    let input = root.join("in");
    let data = (0..400)
        .map(|i| format!("@r{i}|BARCODE=G{}\nA\n+\nI\n", i % 20))
        .collect::<String>();
    fs::write(&input, &data).unwrap();
    for (name, force) in [("seeded", false), ("forced", true)] {
        let out = root.join(name);
        let mut c = Command::new(env!("CARGO_BIN_EXE_clonoforge"));
        c.env("CLONOFORGE_PROFILE", "1")
            .arg("build-consensus")
            .arg(&input)
            .arg("--outdir")
            .arg(&out)
            .args(["--outname", "probe", "--memory-mib", "512"]);
        if force {
            c.arg("--force-spill");
        }
        let r = c.output().unwrap();
        assert!(r.status.success(), "{}", String::from_utf8_lossy(&r.stderr));
        let v: serde_json::Value = serde_json::from_slice(&r.stdout).unwrap();
        assert_eq!(v["spilled"], true);
        assert_eq!(v["forced_spill"], force);
        assert_eq!(
            v["profile"]["execution_path"],
            if force {
                "forced-spill"
            } else {
                "seeded-spill"
            }
        );
        assert_eq!(v["profile"]["raw_bytes"], data.len());
        assert_eq!(v["profile"]["decoded_bytes"], data.len());
        assert!(
            v["profile"]["record_read_parse_inclusive_ns"]
                .as_u64()
                .unwrap()
                > 0
        );
        assert!(v["profile"]["input_initial_runs"].as_u64().unwrap() > 0);
    }
    fs::remove_dir_all(root).unwrap();
}
