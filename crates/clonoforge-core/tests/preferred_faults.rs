use clonoforge_core::{
    admission::Request, bounded_workflow, consensus::Config, external::Limits, publication::Failure,
};
use std::{
    fs,
    path::PathBuf,
    sync::atomic::{AtomicUsize, Ordering},
};
static NEXT: AtomicUsize = AtomicUsize::new(0);
struct Scratch(PathBuf);
impl Scratch {
    fn new() -> Self {
        let p = std::env::temp_dir().join(format!(
            "cf-preferred-fault-{}-{}",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        fs::create_dir(&p).unwrap();
        Self(p)
    }
}
impl Drop for Scratch {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}
#[test]
fn preferred_chunk_faults_never_publish_and_retry() {
    let (r, p) = Request::select(1024 * 1024 * 1024, 2, 16384, 50).unwrap();
    assert_eq!(r.workers, 2);
    assert_eq!(r.sort_chunk_bytes, 8 * 1024 * 1024);
    assert_eq!(r.merge_fan_in, 4);
    let lim = Limits {
        chunk_bytes: r.sort_chunk_bytes as usize,
        max_record_bytes: r.max_record_bytes as usize,
        fan_in: r.merge_fan_in as usize,
    };
    for failure in [
        Failure::SpillWrite,
        Failure::SpillSync,
        Failure::MergePartial,
        Failure::Worker,
        Failure::ManifestPublish,
    ] {
        let t = Scratch::new();
        let input = t.0.join("in");
        let out = t.0.join("out");
        fs::create_dir(&out).unwrap();
        let mut text = String::new();
        for i in 0..2200 {
            text.push_str(&format!("@r{i}|BARCODE=G{}\nA\n+\nI\n", i % 100));
        }
        fs::write(&input, text).unwrap();
        let sentinel = out.join("user");
        fs::write(&sentinel, b"keep").unwrap();
        let cfg = Config::default();
        let resource = (
            &lim,
            p.group_weight as usize,
            r.workers as usize,
            r.window as usize,
            r.max_group_records as usize,
        );
        let result = bounded_workflow::run_configured(
            &input,
            &cfg,
            &out,
            "probe",
            &t.0,
            resource,
            Some(failure),
        );
        assert!(result.is_err(), "{failure:?}");
        assert!(!out.join("probe_completion.json").exists());
        assert_eq!(fs::read(&sentinel).unwrap(), b"keep");
        bounded_workflow::run_configured(&input, &cfg, &out, "probe", &t.0, resource, None)
            .unwrap();
        assert!(out.join("probe_completion.json").exists());
    }
}
