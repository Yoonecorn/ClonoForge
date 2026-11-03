use std::{fs, path::Path, process::Command};
fn run(a: &Path, b: &Path, out: &Path, workers: &str) -> std::process::Output {
    Command::new(env!("CARGO_BIN_EXE_clonoforge"))
        .env("CLONOFORGE_PROFILE", "1")
        .args([
            "preprocess-paired",
            a.to_str().unwrap(),
            b.to_str().unwrap(),
            "--output",
            out.to_str().unwrap(),
            "--header",
            "identical",
            "--key-start",
            "0",
            "--key-len",
            "2",
            "--target",
            "second",
            "--trim",
            "0",
            "--max-group",
            "2",
            "--threads",
            workers,
            "--memory-mib",
            "2048",
        ])
        .output()
        .unwrap()
}
#[test]
fn paired_order_worker_identity_and_failures() {
    let root = std::env::temp_dir().join(format!("clonoforge-paired-cli-{}", std::process::id()));
    let _ = fs::remove_dir_all(&root);
    fs::create_dir(&root).unwrap();
    let a = root.join("a");
    let b = root.join("b");
    fs::write(&a, b"@x\nTT\n+\nII\n@y\nAA\n+\nII\n@z\nTT\n+\nII\n").unwrap();
    fs::write(&b, b"@x\nAC\n+\nII\n@y\nTG\n+\nII\n@z\nAC\n+\nII\n").unwrap();
    for n in ["1", "2", "8"] {
        let out = root.join(n);
        let result = run(&a, &b, &out, n);
        assert!(
            result.status.success(),
            "{}",
            String::from_utf8_lossy(&result.stderr)
        );
        let stderr = String::from_utf8(result.stderr).unwrap();
        let report: serde_json::Value =
            serde_json::from_str(stderr.lines().last().unwrap()).unwrap();
        assert_eq!(report["actual_workers"].as_u64(), Some(n.parse().unwrap()));
    }
    for name in [
        "qualified.fastq",
        "m3_consensus-pass.fastq",
        "m3_completion.json",
        "completion.json",
    ] {
        let one = fs::read(root.join("1").join(name)).unwrap();
        for n in ["2", "8"] {
            assert_eq!(one, fs::read(root.join(n).join(name)).unwrap(), "{name}");
        }
    }
    let q = fs::read_to_string(root.join("1/qualified.fastq")).unwrap();
    let keys: Vec<_> = q
        .lines()
        .filter(|l| l.starts_with('@'))
        .map(|l| l.split("BARCODE=").nth(1).unwrap())
        .collect();
    assert_eq!(keys, ["TT", "AA", "TT"]);
    let smoke = Command::new(env!("CARGO_BIN_EXE_clonoforge"))
        .args([
            "preprocess-paired",
            a.to_str().unwrap(),
            b.to_str().unwrap(),
            "--output",
            root.join("default").to_str().unwrap(),
            "--header",
            "identical",
            "--key-start",
            "0",
            "--key-len",
            "2",
            "--target",
            "second",
            "--trim",
            "0",
        ])
        .output()
        .unwrap();
    assert!(
        smoke.status.success(),
        "{}",
        String::from_utf8_lossy(&smoke.stderr)
    );
    for phase in [
        "stage",
        "read",
        "identity-sort",
        "group-sort",
        "order-sort",
        "qualified",
        "prepare",
        "consensus",
        "precommit",
    ] {
        let path = root.join(format!("injected-{phase}"));
        let fail = Command::new(env!("CARGO_BIN_EXE_clonoforge"))
            .env("CLONOFORGE_M4_FAIL", phase)
            .args([
                "preprocess-paired",
                a.to_str().unwrap(),
                b.to_str().unwrap(),
                "--output",
                path.to_str().unwrap(),
                "--header",
                "identical",
                "--key-start",
                "0",
                "--key-len",
                "2",
                "--target",
                "second",
                "--trim",
                "0",
                "--max-group",
                "2",
                "--memory-mib",
                "2048",
            ])
            .output()
            .unwrap();
        assert!(!fail.status.success());
        assert!(!path.exists());
        assert!(run(&a, &b, &path, "1").status.success());
    }
    let output_closed = root.join("closed-stdout");
    let mut child = Command::new(env!("CARGO_BIN_EXE_clonoforge"))
        .args([
            "preprocess-paired",
            a.to_str().unwrap(),
            b.to_str().unwrap(),
            "--output",
            output_closed.to_str().unwrap(),
            "--header",
            "identical",
            "--key-start",
            "0",
            "--key-len",
            "2",
            "--target",
            "second",
            "--trim",
            "0",
            "--max-group",
            "2",
            "--memory-mib",
            "2048",
        ])
        .stdout(std::process::Stdio::piped())
        .spawn()
        .unwrap();
    drop(child.stdout.take());
    assert!(child.wait().unwrap().success());
    assert!(output_closed.join("completion.json").exists());
    let collision = run(&a, &b, &root.join("1"), "1");
    assert!(!collision.status.success());
    fs::write(&b, [0x1f, 0x8b, 0x08, 0, 0, 0, 0, 0, 0, 3, 0xff]).unwrap();
    assert!(!run(&a, &b, &root.join("bad-gzip"), "1").status.success());
    assert!(!root.join("bad-gzip").exists());
    fs::write(&b, b"@x\nAC\n+\nII\n").unwrap();
    assert!(!run(&a, &b, &root.join("bad"), "1").status.success());
    assert!(!root.join("bad").exists());
    fs::remove_dir_all(root).unwrap();
}
