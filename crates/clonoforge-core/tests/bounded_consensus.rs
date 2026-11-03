use clonoforge_core::{
    consensus::Config,
    group,
    io::{Format, Reader},
    publication::{self, Failure},
};
use std::{
    fs,
    io::{Cursor, Write},
    path::{Path, PathBuf},
    sync::atomic::{AtomicUsize, Ordering},
};
static NEXT: AtomicUsize = AtomicUsize::new(0);
struct Scratch(PathBuf);
impl Scratch {
    fn new() -> Self {
        let p = std::env::temp_dir().join(format!(
            "clonoforge-m2b-{}-{}",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        fs::create_dir(&p).unwrap();
        Self(p)
    }
}
#[test]
fn gzip_content_detection_and_corruption() {
    use flate2::{write::GzEncoder, Compression};
    let scratch = Scratch::new();
    let plain = scratch.0.join("plain.fastq");
    let compressed = scratch.0.join("compressed.fastq");
    fs::write(&plain, mixed()).unwrap();
    let mut encoder = GzEncoder::new(Vec::new(), Compression::default());
    encoder.write_all(mixed().as_bytes()).unwrap();
    let bytes = encoder.finish().unwrap();
    fs::write(&compressed, &bytes).unwrap();
    let config = Config::default();
    let a = group::load(&plain, &config).unwrap();
    let b = group::load(&compressed, &config).unwrap();
    assert_eq!(a.len(), b.len());
    for (left, right) in a.iter().zip(&b) {
        assert_eq!(left.barcode, right.barcode);
        assert_eq!(left.records, right.records);
    }
    fs::write(&compressed, &bytes[..bytes.len() - 4]).unwrap();
    assert!(group::load(&compressed, &config).is_err());
}
impl Drop for Scratch {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}
fn admit(text: &str, c: &Config) -> std::io::Result<Vec<group::Group>> {
    group::admit(Reader::new(Cursor::new(text), Format::Fastq), c)
}
fn mixed() -> String {
    "@same|BARCODE=B|X=z\nAC\n+\nII\n@same|BARCODE=A\nGT\n+\n!~\n@same|BARCODE=B|X=w\nAC\n+\nII\n"
        .into()
}
#[test]
fn grouping_contract() {
    let c = Config::default();
    let groups = admit(&mixed(), &c).unwrap();
    assert_eq!(
        groups
            .iter()
            .map(|g| g.barcode.as_str())
            .collect::<Vec<_>>(),
        ["B", "A"]
    );
    assert_eq!(
        groups[0]
            .records
            .iter()
            .map(|r| r.ordinal)
            .collect::<Vec<_>>(),
        [0, 2]
    );
    let mut custom = c.clone();
    custom.barcode_field = "TAG".into();
    assert!(admit("@r|TAG=a\nA\n+\n!\n", &custom).is_ok());
    let x = admit("@r|BARCODE=a\nA\n+\n!\n@s|BARCODE=A\nA\n+\n!\n", &c).unwrap();
    assert_eq!(x.len(), 2);
    for t in [
        "@r|BARCODE=A\nA\n+\n!\n@r|BARCODE=A\nA\n+\n!\n",
        "@r\nA\n+\n!\n",
        "@r|BARCODE=\nA\n+\n!\n",
        "@r|BARCODE=A|BARCODE=B\nA\n+\n!\n",
        "@r|BARCODE=A-B\nA\n+\n!\n",
        "@r|BARCODE=A\nA\n+\n!\n@s|BARCODE=A\nAA\n+\n!!\n",
        "@r|BARCODE=A\na\n+\n!\n",
        "@r|BARCODE=A\nR\n+\n!\n",
    ] {
        assert!(admit(t, &c).is_err(), "{t}");
    }
    let fifty_one = (0..51)
        .map(|i| format!("@r{i}|BARCODE=A\nA\n+\n!\n"))
        .collect::<String>();
    assert!(admit(&fifty_one, &c).is_err());
}
fn verify_manifest(out: &Path, prefix: &str) -> serde_json::Value {
    let path = out.join(format!("{prefix}_completion.json"));
    let v: serde_json::Value = serde_json::from_slice(&fs::read(path).unwrap()).unwrap();
    assert_eq!(v["success"], true);
    assert_eq!(v["schema_version"], 1);
    assert_eq!(v["command"], "clonoforge build-consensus");
    assert_eq!(v["configuration"]["prefix"], prefix);
    assert_eq!(v["configuration"]["logical_processes"], 1);
    for entry in v["files"].as_array().unwrap() {
        let path = out.join(entry["name"].as_str().unwrap());
        let (size, hash) = publication::digest(&path).unwrap();
        assert_eq!(entry["bytes"], size);
        assert_eq!(entry["sha256"], hash);
        let text = fs::read_to_string(&path).unwrap();
        assert_eq!(entry["records"], text.lines().count() / 4);
    }
    for n in v["absent_categories"].as_array().unwrap() {
        assert!(!out.join(n.as_str().unwrap()).exists());
    }
    v
}
#[test]
fn output_shapes_and_determinism() {
    let scratch = Scratch::new();
    for (name, text, n, failed, expected) in [
        ("pass", mixed(), 1, true, 1),
        ("fail", mixed(), 3, true, 1),
        ("suppressed", mixed(), 3, false, 0),
        ("mixed", mixed(), 2, true, 2),
        ("empty", String::new(), 1, true, 0),
    ] {
        let c = Config {
            min_count: n,
            failed,
            ..Config::default()
        };
        let groups = admit(&text, &c).unwrap();
        let out = scratch.0.join(name);
        publication::publish(&groups, &c, &out, "probe", None).unwrap();
        let m = verify_manifest(&out, "probe");
        assert_eq!(m["files"].as_array().unwrap().len(), expected);
        let repeated = scratch.0.join(format!("{name}_repeat"));
        publication::publish(&groups, &c, &repeated, "probe", None).unwrap();
        assert_eq!(
            fs::read(out.join("probe_completion.json")).unwrap(),
            fs::read(repeated.join("probe_completion.json")).unwrap()
        );
        if name == "mixed" {
            assert_eq!(
                fs::read_to_string(out.join("probe_consensus-fail.fastq")).unwrap(),
                "@same|BARCODE=A\nGT\n+\n!~\n"
            );
        }
        let before = fs::read(out.join("probe_completion.json")).unwrap();
        assert!(publication::publish(&groups, &c, &out, "probe", None).is_err());
        assert_eq!(before, fs::read(out.join("probe_completion.json")).unwrap());
    }
}
#[test]
fn injected_failure_at_every_boundary() {
    let scratch = Scratch::new();
    let c = Config {
        min_count: 2,
        failed: true,
        ..Config::default()
    };
    let groups = admit(&mixed(), &c).unwrap();
    for (i, point) in [
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
    ]
    .into_iter()
    .enumerate()
    {
        let out = scratch.0.join(format!("out{i}"));
        assert!(publication::publish(&groups, &c, &out, "probe", Some(point)).is_err());
        assert!(!out.join("probe_completion.json").exists());
        assert!(!out.join("probe_consensus-pass.fastq").exists());
        assert!(!out.join("probe_consensus-fail.fastq").exists());
        assert!(fs::read_dir(&scratch.0).unwrap().any(|x| x
            .unwrap()
            .file_name()
            .to_string_lossy()
            .contains("incomplete")));
        publication::publish(&groups, &c, &out, "probe", None).unwrap();
        verify_manifest(&out, "probe");
    }
}
#[test]
fn lazy_no_clobber_and_broken_symlink() {
    let scratch = Scratch::new();
    let c = Config::default();
    let groups = admit(&mixed(), &c).unwrap();
    for name in [
        "probe_consensus-pass.fastq",
        "probe_consensus-fail.fastq",
        "probe_completion.json",
    ] {
        let out = scratch.0.join(name);
        fs::create_dir(&out).unwrap();
        fs::write(out.join(name), b"user-data").unwrap();
        assert!(publication::publish(&groups, &c, &out, "probe", None).is_err());
        assert_eq!(fs::read(out.join(name)).unwrap(), b"user-data");
    }
    #[cfg(unix)]
    {
        let out = scratch.0.join("symlink");
        fs::create_dir(&out).unwrap();
        std::os::unix::fs::symlink("missing", out.join("probe_completion.json")).unwrap();
        assert!(publication::publish(&groups, &c, &out, "probe", None).is_err());
        assert!(fs::symlink_metadata(out.join("probe_completion.json"))
            .unwrap()
            .file_type()
            .is_symlink());
    }
}

#[test]
fn lock_and_publication_collisions() {
    let scratch = Scratch::new();
    let c = Config {
        min_count: 2,
        failed: true,
        ..Config::default()
    };
    let groups = admit(&mixed(), &c).unwrap();
    for (i, (point, name)) in [
        (Failure::CollisionPass, "probe_consensus-pass.fastq"),
        (Failure::CollisionFail, "probe_consensus-fail.fastq"),
        (Failure::CollisionManifest, "probe_completion.json"),
    ]
    .into_iter()
    .enumerate()
    {
        let out = scratch.0.join(format!("collision{i}"));
        assert!(publication::publish(&groups, &c, &out, "probe", Some(point)).is_err());
        assert_eq!(fs::read(out.join(name)).unwrap(), b"external-user-file");
        if name != "probe_completion.json" {
            assert!(!out.join("probe_completion.json").exists());
        } else {
            assert!(serde_json::from_slice::<serde_json::Value>(
                &fs::read(out.join(name)).unwrap()
            )
            .is_err());
        }
    }
    let out = scratch.0.join("lock");
    fs::create_dir(&out).unwrap();
    let lock = fs::OpenOptions::new()
        .read(true)
        .write(true)
        .create_new(true)
        .open(out.join(".probe.transaction.lock"))
        .unwrap();
    lock.try_lock().unwrap();
    assert!(publication::publish(&groups, &c, &out, "probe", None).is_err());
    drop(lock);
    publication::publish(&groups, &c, &out, "probe", None).unwrap();
    verify_manifest(&out, "probe");
}
#[test]
fn late_publish_failure_retains_distinct_pending_marker_and_rolls_back_owned_categories() {
    use clonoforge_core::publication::Failure;
    let scratch = Scratch::new();
    let out = scratch.0.join("out");
    fs::create_dir(&out).unwrap();
    fs::write(out.join("user"), b"keep").unwrap();
    let cfg = Config::default();
    let groups = admit("@r|BARCODE=A\nA\n+\nI\n", &cfg).unwrap();
    let err = publication::publish(&groups, &cfg, &out, "probe", Some(Failure::ManifestPublish))
        .unwrap_err();
    assert!(err
        .to_string()
        .contains("incomplete pending marker retained"));
    assert!(!out.join("probe_completion.json").exists());
    assert!(!out.join("probe_consensus-pass.fastq").exists());
    assert_eq!(fs::read(out.join("user")).unwrap(), b"keep");
    let pending: Vec<_> = fs::read_dir(&out)
        .unwrap()
        .map(|x| x.unwrap().path())
        .filter(|p| {
            p.file_name()
                .unwrap()
                .to_string_lossy()
                .starts_with(".clonoforge-probe.pending-manifest.")
        })
        .collect();
    assert_eq!(pending.len(), 1);
    let value: serde_json::Value = serde_json::from_slice(&fs::read(&pending[0]).unwrap()).unwrap();
    assert_eq!(value["success"], true);
    assert_ne!(pending[0].file_name().unwrap(), "probe_completion.json");
    publication::publish(&groups, &cfg, &out, "probe", None).unwrap();
    assert!(out.join("probe_completion.json").exists());
    assert!(pending[0].exists());
}
