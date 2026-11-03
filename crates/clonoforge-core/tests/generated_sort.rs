use clonoforge_core::{
    consensus::Config,
    external::{self, Limits},
};
use std::fs;
#[test]
fn generated_sort_matches_checked_path_through_multiple_merges_and_duplicate_index() {
    let root = std::env::temp_dir().join(format!("cf-generated-sort-{}", std::process::id()));
    fs::create_dir(&root).unwrap();
    let lim = Limits {
        chunk_bytes: 10000,
        max_record_bytes: 100,
        fan_in: 2,
    };
    let cfg = Config::default();
    let text: String = (0..41)
        .map(|i| {
            format!(
                "@r{i}|BARCODE={}\nAC\n+\nIJ\n",
                if i % 3 == 0 { "T" } else { "A" }
            )
        })
        .collect();
    let input = root.join("input");
    fs::write(&input, &text).unwrap();
    for name in ["checked", "direct", "unique", "duplicate"] {
        fs::create_dir(root.join(name)).unwrap();
    }
    let checked = external::sort_input(&input, &cfg, &root.join("checked"), &lim).unwrap();
    let direct =
        external::sort_generated_by_barcode(&input, &cfg, &root.join("direct"), &lim).unwrap();
    assert_eq!(fs::read(checked).unwrap(), fs::read(direct).unwrap());
    external::check_unique_headers(&input, &cfg, &root.join("unique"), &lim).unwrap();
    fs::write(&input, format!("{text}@r0|BARCODE=T\nAC\n+\nIJ\n")).unwrap();
    assert!(
        external::check_unique_headers(&input, &cfg, &root.join("duplicate"), &lim)
            .unwrap_err()
            .to_string()
            .contains("duplicate complete header")
    );
    fs::remove_dir_all(root).unwrap();
}

#[test]
fn presorted_filter_preserves_order_and_excludes_whole_groups() {
    let root = std::env::temp_dir().join(format!("cf-generated-filter-{}", std::process::id()));
    fs::create_dir(&root).unwrap();
    let lim = Limits {
        chunk_bytes: 10000,
        max_record_bytes: 200,
        fan_in: 2,
    };
    let cfg = Config::default();
    let input = root.join("input");
    fs::write(&input,"@r0|BARCODE=B\nAC\n+\nII\n@r1|BARCODE=A\nAC\n+\nII\n@r2|BARCODE=B\nAC\n+\nII\n@r3|BARCODE=C\nA\n+\nI\n@r4|BARCODE=C\nAC\n+\nII\n@r5|BARCODE=B\nAC\n+\nII\n@r6|BARCODE=D\nAC\n+\nII\n@r7|BARCODE=D\nAC\n+\nII\n").unwrap();
    let sorted = external::sort_generated_by_barcode(&input, &cfg, &root, &lim).unwrap();
    let output =
        external::filter_paired_groups(&sorted, &cfg, &root.join("filtered"), &lim, 2, 3).unwrap();
    let mut observed = Vec::new();
    external::for_each_record(&output, &cfg, &lim, |r| {
        observed.push((r.id, r.ordinal));
        Ok(())
    })
    .unwrap();
    assert_eq!(
        observed,
        vec![
            ("r1".to_owned(), 1),
            ("r6".to_owned(), 6),
            ("r7".to_owned(), 7)
        ]
    );
    assert!(
        external::filter_paired_groups(&sorted, &cfg, &root.join("wrong-count"), &lim, 2, 4)
            .is_err()
    );
    assert!(
        external::filter_paired_groups(&sorted, &cfg, &root.join("bad-limit"), &lim, 51, 3)
            .is_err()
    );
    fs::remove_dir_all(root).unwrap();
}
