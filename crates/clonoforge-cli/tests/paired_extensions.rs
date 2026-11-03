use std::{fs, path::Path, process::Command};
fn command(root: &Path, out: &str, workers: &str) -> Command {
    let mut c = Command::new(env!("CARGO_BIN_EXE_clonoforge"));
    c.env("CLONOFORGE_PROFILE", "1").args([
        "preprocess-paired",
        root.join("r1").to_str().unwrap(),
        root.join("r2").to_str().unwrap(),
        "--output",
        root.join(out).to_str().unwrap(),
        "--header",
        "identical",
        "--key-start",
        "0",
        "--key-len",
        "2",
        "--target",
        "first",
        "--trim",
        "2",
        "--assembly-min-overlap",
        "4",
        "--threads",
        workers,
        "--memory-mib",
        "2048",
        "--primer-iupac",
        "ACR",
    ]);
    c
}
#[test]
fn opt_in_assembly_dispositions_parity_and_publication() {
    let root = std::env::temp_dir().join(format!("clonoforge-extensions-{}", std::process::id()));
    let _ = fs::remove_dir_all(&root);
    fs::create_dir(&root).unwrap();
    fs::write(root.join("r1"), b"@a\nAAACGTACTG\n+\nIIIIIIIIII\n@b\nCCACGTACTG\n+\nIIIIIIIIII\n@c\nTTAAAAAAAA\n+\nIIIIIIIIII\n").unwrap();
    fs::write(
        root.join("r2"),
        b"@a\nTTAACAGT\n+\nIIIIIIII\n@b\nCCCCCCCC\n+\nIIIIIIII\n@c\nTTTTTTTT\n+\nIIIIIIII\n",
    )
    .unwrap();
    for n in ["1", "2", "8"] {
        let result = command(&root, n, n).output().unwrap();
        assert!(
            result.status.success(),
            "{}",
            String::from_utf8_lossy(&result.stderr)
        );
        let stderr = String::from_utf8(result.stderr).unwrap();
        let profile: serde_json::Value =
            serde_json::from_str(stderr.lines().last().unwrap()).unwrap();
        assert_eq!(profile["actual_workers"].as_u64(), Some(n.parse().unwrap()));
    }
    for file in [
        "qualified.fastq",
        "m3_consensus-pass.fastq",
        "m3_completion.json",
        "completion.json",
    ] {
        let one = fs::read(root.join("1").join(file)).unwrap();
        for n in ["2", "8"] {
            assert_eq!(one, fs::read(root.join(n).join(file)).unwrap(), "{file}");
        }
    }
    let q = fs::read_to_string(root.join("1/qualified.fastq")).unwrap();
    assert_eq!(q.lines().nth(1), Some("TACTGTTAA"));
    let m: serde_json::Value =
        serde_json::from_slice(&fs::read(root.join("1/completion.json")).unwrap()).unwrap();
    assert_eq!(m["schema"], "clonoforge-paired-v2");
    assert_eq!(m["target_only"], false);
    for (key, value) in [
        ("pairs", 3),
        ("assembled", 1),
        ("assembly_no_overlap", 1),
        ("assembly_ambiguous", 1),
        ("qualified", 1),
    ] {
        assert_eq!(m["counts"][key], value);
    }
    // A 3-prime degenerate primer is applied to the assembled R1-oriented target.
    let template = command(&root, "three-prime", "1");
    let mut args: Vec<String> = template
        .get_args()
        .map(|s| s.to_str().unwrap().to_owned())
        .collect();
    let primer_index = args.iter().position(|a| a == "--primer-iupac").unwrap() + 1;
    args[primer_index] = "TAR".to_owned();
    let result = Command::new(env!("CARGO_BIN_EXE_clonoforge"))
        .args(args)
        .args(["--primer-anchor", "3"])
        .output()
        .unwrap();
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    let three = fs::read_to_string(root.join("three-prime/qualified.fastq")).unwrap();
    assert_eq!(three.lines().nth(1), Some("ACGTACTGT"));
    for flag in ["--primer", "--assembly-max-read-length"] {
        let mut c = command(&root, flag, "1");
        c.args([flag, if flag == "--primer" { "ACG" } else { "0" }]);
        assert!(!c.output().unwrap().status.success());
        assert!(!root.join(flag).exists());
    }
    assert!(!command(&root, "1", "1").output().unwrap().status.success());
    let before = fs::read(root.join("1/completion.json")).unwrap();
    for phase in ["read", "group-sort", "consensus", "precommit"] {
        let fail = command(&root, phase, "1")
            .env("CLONOFORGE_M4_FAIL", phase)
            .output()
            .unwrap();
        assert!(!fail.status.success());
        assert!(!root.join(phase).exists());
        assert!(command(&root, phase, "1")
            .output()
            .unwrap()
            .status
            .success());
    }
    assert_eq!(before, fs::read(root.join("1/completion.json")).unwrap());
    fs::remove_dir_all(root).unwrap();
}
