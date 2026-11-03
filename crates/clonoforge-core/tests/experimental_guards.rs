#![cfg(feature = "experimental-large")]
use clonoforge_core::{
    admission::Request,
    bounded_workflow, consensus,
    external::{self, GroupExecution, Limits},
    filter::Filter,
    pair::HeaderMode,
    paired_prepare::{self, Extensions, Options},
};
use std::path::Path;
#[test]
fn old_public_apis_retain_fifty_and_experimental_options_are_explicit() {
    let path = Path::new("/nonexistent-clonoforge-experimental-guard");
    let lim = Limits {
        chunk_bytes: 1024 * 1024,
        max_record_bytes: 256,
        fan_in: 2,
    };
    let mut opt = Options {
        header_mode: HeaderMode::Identical,
        key_start: 0,
        key_len: 8,
        target_is_first: false,
        remove_target_prefix: 0,
        primer: None,
        filter: Filter {
            min_length: 1,
            max_n: usize::MAX,
            min_mean_phred: 0,
        },
        max_raw_record: 256,
        max_group_count: 100,
    };
    assert!(
        paired_prepare::prepare_extended(path, path, path, &opt, &lim, &Extensions::default())
            .unwrap_err()
            .to_string()
            .contains("group limit")
    );
    opt.max_group_count = 51;
    assert!(paired_prepare::prepare_experimental(
        path,
        path,
        path,
        &opt,
        &lim,
        &Extensions::default()
    )
    .unwrap_err()
    .to_string()
    .contains("unsupported"));
    let cfg = consensus::Config::default();
    assert!(consensus::position(&[(b'A', 20); 51], 0.6, 0).is_err());
    assert!(
        external::filter_paired_groups(path, &cfg, path, &lim, 100, 0)
            .unwrap_err()
            .to_string()
            .contains("group limit")
    );
    let ex = GroupExecution {
        max_group_bytes: 1024,
        workers: 1,
        window: 1,
        max_group_count: 100,
        failure: None,
    };
    assert!(
        external::sort_outputs_configured(path, &cfg, path, &lim, ex)
            .unwrap_err()
            .to_string()
            .contains("read limit")
    );
    assert!(bounded_workflow::run_presorted_configured(
        path,
        &cfg,
        path,
        "x",
        path,
        (&lim, 1024, 1, 1, 100)
    )
    .unwrap_err()
    .to_string()
    .contains("read limit"));
    for n in [0, 51, 99, 101, 249, 251, 499, 501] {
        assert!(Request::select_experimental(u64::MAX, 1, 256, n).is_err());
    }
    for n in [100, 250, 500] {
        let (r, p) = Request::select_experimental(u64::MAX, 1, 256, n).unwrap();
        assert_eq!(r.max_group_records, n);
        assert!(r.plan().is_err());
        assert_eq!(p, r.plan_experimental().unwrap());
        let mut below = r;
        below.budget = p.total - 1;
        assert!(below.plan_experimental().is_err());
    }
}
