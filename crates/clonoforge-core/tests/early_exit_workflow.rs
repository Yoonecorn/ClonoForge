use clonoforge_core::{
    bounded_workflow, consensus::Config, external::Limits, publication::Failure,
};
use std::{fs, sync::mpsc, time::Duration};
#[test]
fn early_worker_exit_cancels_pipeline_and_releases_lock_for_retry() {
    let (tx, rx) = mpsc::channel();
    std::thread::spawn(move || {
        let root = std::env::temp_dir().join(format!("cf-early-exit-{}", std::process::id()));
        fs::create_dir(&root).unwrap();
        let input = root.join("input");
        let out = root.join("out");
        fs::create_dir(&out).unwrap();
        fs::write(out.join("user"), b"keep").unwrap();
        let text = (0..100)
            .map(|i| format!("@r{i}|BARCODE=G{}\nA\n+\nI\n", i % 20))
            .collect::<String>();
        fs::write(&input, text).unwrap();
        let lim = Limits {
            chunk_bytes: 10000,
            max_record_bytes: 100,
            fan_in: 2,
        };
        let cfg = Config::default();
        let resource = (&lim, 100_000, 4, 8, 50);
        let failure = bounded_workflow::run_configured(
            &input,
            &cfg,
            &out,
            "probe",
            &root,
            resource,
            Some(Failure::WorkerEarlyExit),
        )
        .unwrap_err()
        .to_string();
        assert!(failure.contains("worker exited unexpectedly"), "{failure}");
        assert!(!out.join("probe_completion.json").exists());
        assert!(!out.join("probe_consensus-pass.fastq").exists());
        assert_eq!(fs::read(out.join("user")).unwrap(), b"keep");
        bounded_workflow::run_configured(&input, &cfg, &out, "probe", &root, resource, None)
            .unwrap();
        assert!(out.join("probe_completion.json").exists());
        fs::remove_dir_all(&root).unwrap();
        tx.send(()).unwrap();
    });
    rx.recv_timeout(Duration::from_secs(30))
        .expect("pipeline did not terminate and retry");
}
