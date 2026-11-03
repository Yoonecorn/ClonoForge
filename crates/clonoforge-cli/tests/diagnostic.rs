use std::{fs, process::Command};
#[test]
fn cli_publication_and_determinism() {
    let dir = std::env::temp_dir().join(format!("clonoforge-cli-{}", std::process::id()));
    fs::create_dir_all(&dir).unwrap();
    let input = dir.join("input.fq");
    let output = dir.join("output.fq");
    let run = |cmd: &str, out: Option<&std::path::Path>| {
        let mut c = Command::new(env!("CARGO_BIN_EXE_clonoforge"));
        c.arg(cmd).arg("fastq").arg(&input);
        if let Some(p) = out {
            c.arg(p);
        }
        c.output().unwrap()
    };
    fs::write(&input, b"@r\nAC\n+\n!~\n").unwrap();
    assert!(run("validate", None).status.success());
    assert_eq!(run("inspect", None).stdout, run("inspect", None).stdout);
    assert!(run("normalize", Some(&output)).status.success());
    assert_eq!(fs::read(&input).unwrap(), fs::read(&output).unwrap());
    assert_eq!(run("normalize", Some(&output)).status.code(), Some(2));
    fs::remove_file(&output).unwrap();
    fs::write(&input, b"@r\nA\n+\n!\n@bad\n").unwrap();
    assert_eq!(run("normalize", Some(&output)).status.code(), Some(2));
    assert!(!output.exists());
    assert_eq!(fs::read_dir(&dir).unwrap().count(), 1);
    assert!(Command::new(env!("CARGO_BIN_EXE_clonoforge"))
        .arg("--version")
        .output()
        .unwrap()
        .status
        .success());
    fs::remove_dir_all(dir).unwrap();
}
