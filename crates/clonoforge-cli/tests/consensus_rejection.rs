use std::{fs, process::Command};
#[test]
fn rejects_invalid_or_unsupported_before_completion() {
    let root = std::env::temp_dir().join(format!("clonoforge-cli-reject-{}", std::process::id()));
    fs::create_dir(&root).unwrap();
    let input = root.join("input.fastq");
    let out = root.join("out");
    let binary = env!("CARGO_BIN_EXE_clonoforge");
    let run = |extra: &[&str]| {
        Command::new(binary)
            .arg("build-consensus")
            .arg(&input)
            .arg("--outdir")
            .arg(&out)
            .args(extra)
            .output()
            .unwrap()
    };
    fs::write(&input, b"@r|BARCODE=A\nA\n+\n!\n").unwrap();
    for args in [
        vec!["--freq", "nan"],
        vec!["--freq", "inf"],
        vec!["--freq", "1.1"],
        vec!["--freq", "-0.1"],
        vec!["-n", "0"],
        vec!["-n", "-1"],
        vec!["-q", "94"],
        vec!["-q", "-1"],
        vec!["--bf", "barcode"],
        vec!["--outname", "../escape"],
        vec!["--dep"],
        vec!["--maxgap", "0.5"],
        vec!["--nproc", "2"],
        vec!["--gzip-output"],
        vec!["--freq"],
    ] {
        let result = run(&args);
        assert_eq!(result.status.code(), Some(2), "{args:?}");
        assert!(!out.join("input_completion.json").exists());
    }
    for bytes in [
        b"\x1f\x8b".as_slice(),
        b">r|BARCODE=A\nA\n",
        b"@r|BARCODE=A\na\n+\n!\n",
        b"@r|BARCODE=A\nR\n+\n!\n",
        b"@r|BARCODE=A\n-\n+\n!\n",
        b"@r|BARCODE=A\nA\n+\n",
        b"@r|BARCODE=\nA\n+\n!\n",
        b"@r|BARCODE=A|BARCODE=B\nA\n+\n!\n",
    ] {
        fs::write(&input, bytes).unwrap();
        assert_eq!(run(&[]).status.code(), Some(2));
        assert!(!out.join("input_completion.json").exists());
    }
    assert!(Command::new(binary)
        .args(["build-consensus", "--help"])
        .output()
        .unwrap()
        .status
        .success());
    fs::remove_dir_all(root).unwrap();
}
