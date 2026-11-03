use clonoforge_core::{
    admission::{Plan, Request},
    bounded_workflow,
    consensus::Config,
    external::Limits,
    filter::Filter,
    pair::HeaderMode,
    paired_prepare::{self, Extensions, Options},
    publication,
};
use std::{
    collections::BTreeMap,
    fs::{self, File, OpenOptions},
    io::Write,
    path::{Path, PathBuf},
};

/// The caller has already reserved the additional paired raw-parser allowance.
/// Grow only sort chunks within the unchanged complete admission model; retain
/// the selected worker/window and every other limit from M3's selector.
fn admit_sort_chunk(request: Request, plan: Plan) -> (Request, Plan) {
    for chunk_mib in [64u64, 32, 16] {
        let candidate = Request {
            sort_chunk_bytes: chunk_mib * 1024 * 1024,
            ..request
        };
        if candidate.sort_chunk_bytes <= request.sort_chunk_bytes {
            continue;
        }
        let plan = if candidate.max_group_records > 50 {
            #[cfg(feature = "experimental-large")]
            {
                candidate.plan_experimental()
            }
            #[cfg(not(feature = "experimental-large"))]
            {
                candidate.plan()
            }
        } else {
            candidate.plan()
        };
        if let Ok(candidate_plan) = plan {
            return (candidate, candidate_plan);
        }
    }
    (request, plan)
}

fn inject(point: &str) -> Result<(), Box<dyn std::error::Error>> {
    if std::env::var("CLONOFORGE_M4_FAIL").ok().as_deref() == Some(point) {
        return Err(format!("paired: injected {point}").into());
    }
    Ok(())
}
pub fn run(args: &[String]) -> Result<(), Box<dyn std::error::Error>> {
    if args == ["--help"] {
        if cfg!(feature = "experimental-large") {
            println!("clonoforge preprocess-paired FIRST_INPUT SECOND_INPUT --output NEW_DIRECTORY --header identical|sra|slash --key-start N --key-len N --target first|second --trim N [--primer ACGT] [--min-length N] [--max-n N] [--min-mean-q N] [--max-group N] [--experimental-large-group 100|250|500] [--threads N] [--memory-mib N] [--max-record N] [--primer-iupac IUPAC --primer-anchor 5|3] [--assembly-min-overlap N --assembly-max-mismatches N --assembly-max-read-length N]\nExtract an exact barcode from the first input and group target reads for consensus. Put the barcode-bearing mate first; it need not be biological read 1. Header mode must match input order. Optional ungapped overlap assembly requires target first. Use a new output directory; completion.json records a successful run. The standard limit is 50 reads per group. Limits 100, 250 and 500 require this experimental build and explicit selection.");
        } else {
            println!("clonoforge preprocess-paired FIRST_INPUT SECOND_INPUT --output NEW_DIRECTORY --header identical|sra|slash --key-start N --key-len N --target first|second --trim N [--primer ACGT] [--min-length N] [--max-n N] [--min-mean-q N] [--max-group N] [--threads N] [--memory-mib N] [--max-record N] [--primer-iupac IUPAC --primer-anchor 5|3] [--assembly-min-overlap N --assembly-max-mismatches N --assembly-max-read-length N]\nExtract an exact barcode from the first input and group target reads for consensus. Put the barcode-bearing mate first; it need not be biological read 1. Header mode must match input order. Optional ungapped overlap assembly requires target first. Use a new output directory; completion.json records a successful run. The standard limit is 50 reads per group.");
        }
        return Ok(());
    }
    if args.len() < 2 {
        return Err("paired: two input paths required".into());
    }
    let first = PathBuf::from(&args[0]);
    let second = PathBuf::from(&args[1]);
    let mut kv = BTreeMap::new();
    let mut i = 2;
    while i < args.len() {
        let k = &args[i];
        if !k.starts_with("--")
            || i + 1 >= args.len()
            || kv.insert(k.as_str(), args[i + 1].as_str()).is_some()
        {
            return Err("paired: invalid or duplicate option".into());
        }
        i += 2;
    }
    let get = |k| {
        kv.get(k)
            .copied()
            .ok_or_else(|| format!("paired: missing {k}"))
    };
    let number = |k, default: &str| -> Result<usize, Box<dyn std::error::Error>> {
        Ok(kv.get(k).copied().unwrap_or(default).parse()?)
    };
    for k in kv.keys() {
        if ![
            "--output",
            "--header",
            "--key-start",
            "--key-len",
            "--target",
            "--trim",
            "--primer",
            "--min-length",
            "--max-n",
            "--min-mean-q",
            "--max-group",
            "--experimental-large-group",
            "--threads",
            "--memory-mib",
            "--max-record",
            "--primer-iupac",
            "--primer-anchor",
            "--assembly-min-overlap",
            "--assembly-max-mismatches",
            "--assembly-max-read-length",
        ]
        .contains(k)
        {
            return Err(format!("paired: unknown option {k}").into());
        }
    }
    let output = PathBuf::from(get("--output")?);
    let parent = output
        .parent()
        .filter(|p| !p.as_os_str().is_empty())
        .unwrap_or(Path::new("."));
    let name = output
        .file_name()
        .and_then(|x| x.to_str())
        .ok_or("paired: invalid output name")?;
    publication::validate_prefix(name)?;
    let mode = match get("--header")? {
        "identical" => HeaderMode::Identical,
        "sra" => HeaderMode::SraNumbered,
        "slash" => HeaderMode::SlashSuffix,
        _ => return Err("paired: invalid header mode".into()),
    };
    let target = match get("--target")? {
        "first" => true,
        "second" => false,
        _ => return Err("paired: invalid target".into()),
    };
    let max_record = number("--max-record", "16384")?;
    let compatibility_max = number("--max-group", "50")?;
    let max_group = if let Some(value) = kv.get("--experimental-large-group") {
        if !cfg!(feature = "experimental-large") {
            return Err("paired: experimental large-group build required".into());
        }
        if kv.contains_key("--max-group") {
            return Err(
                "paired: experimental and compatibility group options are mutually exclusive"
                    .into(),
            );
        }
        let n: usize = value.parse()?;
        if ![100, 250, 500].contains(&n) {
            return Err("paired: experimental large-group must be 100, 250 or 500".into());
        }
        n
    } else {
        compatibility_max
    };
    if !kv.contains_key("--experimental-large-group") && max_group > 50 {
        return Err("paired: compatibility maximum is 50".into());
    }
    if max_group > 50 && std::env::var_os("CLONOFORGE_FINAL_UNFUSED").is_some_and(|v| v == "1") {
        return Err("paired: experimental mode requires presorted path".into());
    }
    let threads = number("--threads", "1")?;
    let memory = number("--memory-mib", "1024")?;
    let options = Options {
        header_mode: mode,
        key_start: get("--key-start")?.parse()?,
        key_len: get("--key-len")?.parse()?,
        target_is_first: target,
        remove_target_prefix: get("--trim")?.parse()?,
        primer: kv.get("--primer").map(|s| s.as_bytes().to_vec()),
        filter: Filter {
            min_length: number("--min-length", "1")?,
            max_n: number("--max-n", "18446744073709551615")?,
            min_mean_phred: number("--min-mean-q", "0")?.try_into()?,
        },
        max_raw_record: max_record,
        max_group_count: max_group,
    };
    if !kv.contains_key("--assembly-min-overlap")
        && (kv.contains_key("--assembly-max-mismatches")
            || kv.contains_key("--assembly-max-read-length"))
    {
        return Err("paired: assembly options require --assembly-min-overlap".into());
    }
    if kv.contains_key("--primer-anchor") && !kv.contains_key("--primer-iupac") {
        return Err("paired: --primer-anchor requires --primer-iupac".into());
    }
    let extensions = Extensions {
        assembly: if kv.contains_key("--assembly-min-overlap") {
            Some(clonoforge_core::assembly::AssemblyConfig {
                min_overlap: number("--assembly-min-overlap", "20")?,
                max_mismatches: number("--assembly-max-mismatches", "0")?,
                max_read_length: number("--assembly-max-read-length", "4096")?,
                quality_cap: 60,
            })
        } else {
            None
        },
        primer_iupac: kv.get("--primer-iupac").map(|s| s.as_bytes().to_vec()),
        primer_anchor: match kv.get("--primer-anchor").copied().unwrap_or("5") {
            "5" => Some(clonoforge_core::primer_ext::Anchor::FivePrime),
            "3" => Some(clonoforge_core::primer_ext::Anchor::ThreePrime),
            _ => return Err("paired: primer anchor must be 5 or 3".into()),
        },
    };
    // Validate scientific choices before creating any staging directory.
    if let Some(cfg) = &extensions.assembly {
        cfg.validate()?;
        if !target {
            return Err("paired: assembly requires --target first".into());
        }
    }
    if let Some(primer) = &extensions.primer_iupac {
        clonoforge_core::primer_ext::validate_iupac(primer)?;
        if options.primer.is_some() {
            return Err("paired: primer modes are mutually exclusive".into());
        }
    }
    // M3 admission plus a conservative additional raw parsing allowance; modeled,
    // not an RSS guarantee.
    let budget = u64::try_from(memory)?
        .checked_mul(1024 * 1024)
        .ok_or("paired: budget overflow")?;
    let extra = 32u64
        .checked_mul(u64::try_from(max_record)?)
        .and_then(|n| n.checked_add(32 * 1024 * 1024))
        .ok_or("paired: memory overflow")?;
    let select = if max_group > 50 {
        #[cfg(feature = "experimental-large")]
        {
            Request::select_experimental
        }
        #[cfg(not(feature = "experimental-large"))]
        {
            Request::select
        }
    } else {
        Request::select
    };
    let (req, plan) = select(
        budget
            .checked_sub(extra)
            .ok_or("paired: insufficient memory")?,
        threads.try_into()?,
        max_record.try_into()?,
        max_group.try_into()?,
    )?;
    let (req, plan) = admit_sort_chunk(req, plan);
    let _ = writeln!(
        std::io::stderr().lock(),
        "paired: requested_workers={threads} actual_workers={} modeled_memory_bytes={}",
        req.workers,
        plan.total + extra
    );
    let lim = Limits {
        chunk_bytes: req.sort_chunk_bytes.try_into()?,
        max_record_bytes: max_record,
        fan_in: req.merge_fan_in.try_into()?,
    };
    let stage = (0..1000)
        .find_map(|attempt| {
            let path = parent.join(format!(
                ".{name}.m4-incomplete-{}-{attempt}",
                std::process::id()
            ));
            match fs::create_dir(&path) {
                Ok(()) => Some(Ok(path)),
                Err(e) if e.kind() == std::io::ErrorKind::AlreadyExists => None,
                Err(e) => Some(Err(e)),
            }
        })
        .ok_or("paired: no available stage")??;
    inject("stage")?;
    let work = stage.join("work");
    fs::create_dir(&work)?;
    let work_metadata = fs::symlink_metadata(&work)?;
    let profile_start = std::time::Instant::now();
    let before1 = publication::digest(&first)?;
    let before2 = publication::digest(&second)?;
    let hash_before_ns = profile_start.elapsed().as_nanos();
    let prepare_start = std::time::Instant::now();
    let prepare = if max_group > 50 {
        #[cfg(feature = "experimental-large")]
        {
            paired_prepare::prepare_experimental
        }
        #[cfg(not(feature = "experimental-large"))]
        {
            paired_prepare::prepare_extended
        }
    } else {
        paired_prepare::prepare_extended
    };
    let counts = prepare(&first, &second, &work, &options, &lim, &extensions)?;
    let prepare_ns = prepare_start.elapsed().as_nanos();
    inject("prepare")?;
    let input1 = publication::digest(&first)?;
    let input2 = publication::digest(&second)?;
    if before1 != input1 || before2 != input2 {
        return Err("paired: input changed during processing".into());
    }
    let qualified = stage.join("qualified.fastq");
    let source_qualified = publication::digest(&work.join("qualified.fastq"))?;
    fs::copy(work.join("qualified.fastq"), &qualified)?;
    if publication::digest(&qualified)? != source_qualified {
        return Err("paired: qualified copy mismatch".into());
    }
    File::open(&qualified)?.sync_all()?;
    let consensus_start = std::time::Instant::now();
    let cfg = Config::default();
    let resource = (
        &lim,
        plan.group_weight.try_into()?,
        req.workers.try_into()?,
        req.window.try_into()?,
        max_group,
    );
    let manifest = if max_group > 50 {
        #[cfg(feature = "experimental-large")]
        {
            bounded_workflow::run_presorted_experimental(
                &work.join("group-qualified.run"),
                &cfg,
                &stage,
                "m3",
                &work,
                resource,
            )?
        }
        #[cfg(not(feature = "experimental-large"))]
        {
            return Err("paired: experimental support unavailable".into());
        }
    } else if std::env::var_os("CLONOFORGE_FINAL_UNFUSED").is_some_and(|v| v == "1") {
        bounded_workflow::run_configured(&qualified, &cfg, &stage, "m3", &work, resource, None)?
    } else {
        bounded_workflow::run_presorted_configured(
            &work.join("group-qualified.run"),
            &cfg,
            &stage,
            "m3",
            &work,
            resource,
        )?
    };
    inject("consensus")?;
    let consensus_ns = consensus_start.elapsed().as_nanos();
    let prepared = publication::digest(&qualified)?;
    let m3 = publication::digest(&manifest)?;
    let m3_receipt: serde_json::Value = serde_json::from_reader(File::open(&manifest)?)?;
    if m3_receipt["success"] != true {
        return Err("paired: M3 completion invalid".into());
    }
    let entries = m3_receipt["files"]
        .as_array()
        .ok_or("paired: M3 files invalid")?;
    let mut output_groups = 0u64;
    for entry in entries {
        let name = entry["name"].as_str().ok_or("paired: M3 name invalid")?;
        if name != "m3_consensus-pass.fastq" {
            return Err("paired: unexpected M3 category".into());
        }
        let (size, hash) = publication::digest(&stage.join(name))?;
        if entry["bytes"].as_u64() != Some(size) || entry["sha256"].as_str() != Some(hash.as_str())
        {
            return Err("paired: M3 output hash mismatch".into());
        }
        output_groups = output_groups
            .checked_add(
                entry["records"]
                    .as_u64()
                    .ok_or("paired: invalid M3 count")?,
            )
            .ok_or("paired: count overflow")?;
    }
    if output_groups != counts.qualified_groups {
        return Err("paired: consensus group count mismatch".into());
    }
    let mut replay = clonoforge_core::io::Reader::fastq_limited(
        clonoforge_core::io::input(&qualified)?,
        max_record,
    )?;
    let mut replay_count = 0u64;
    while replay.next_record()?.is_some() {
        replay_count = replay_count
            .checked_add(1)
            .ok_or("paired: count overflow")?;
    }
    if replay_count != counts.qualified {
        return Err("paired: qualified count mismatch".into());
    }
    let pass = stage.join("m3_consensus-pass.fastq");
    let pass_hash = if pass.exists() {
        Some(publication::digest(&pass)?)
    } else {
        None
    };
    #[cfg(unix)]
    {
        use std::os::unix::fs::MetadataExt;
        let current = fs::symlink_metadata(&work)?;
        if !current.is_dir()
            || current.dev() != work_metadata.dev()
            || current.ino() != work_metadata.ino()
        {
            return Err("paired: scratch ownership changed".into());
        }
    }
    #[cfg(not(unix))]
    return Err("paired: scratch ownership verification unsupported".into());
    fs::remove_dir_all(&work)?;
    let mut receipt = serde_json::json!({"schema":"clonoforge-m4-preliminary-v1", "complete":true,"inputs":[{"path":first,"bytes":input1.0,"sha256":input1.1},{"path":second,"bytes":input2.0,"sha256":input2.1}],"configuration":{"header":get("--header")?,"key_start":options.key_start,"key_len":options.key_len,"target":get("--target")?,"trim":options.remove_target_prefix,"primer":kv.get("--primer"),"min_length":options.filter.min_length,"max_n":options.filter.max_n,"min_mean_q":options.filter.min_mean_phred,"max_group":max_group,"max_record":max_record},"counts":{"pairs":counts.pairs,"invalid_key":counts.invalid_key,"short_target":counts.short_target,"primer_miss":counts.primer_miss,"too_short":counts.too_short,"too_many_n":counts.too_many_n,"low_quality":counts.low_quality,"candidates":counts.retained,"group_excluded":counts.group_excluded,"group_excluded_groups":counts.group_excluded_groups,"oversize_records":counts.oversize_records,"oversize_groups":counts.oversize_groups,"unequal_records":counts.unequal_records,"unequal_groups":counts.unequal_groups,"qualified":counts.qualified,"qualified_groups":counts.qualified_groups},"files":{"qualified.fastq":{"bytes":prepared.0,"sha256":prepared.1},"m3_completion.json":{"bytes":m3.0,"sha256":m3.1},"m3_consensus-pass.fastq":pass_hash.map(|(bytes,sha256)| serde_json::json!({"bytes":bytes,"sha256":sha256}))},"target_only":true});
    if extensions.assembly.is_some() || extensions.primer_iupac.is_some() {
        receipt["schema"] = serde_json::json!("clonoforge-paired-v2");
        receipt["target_only"] = serde_json::json!(extensions.assembly.is_none());
        receipt["configuration"]["primer_iupac"] = serde_json::json!(kv.get("--primer-iupac"));
        receipt["configuration"]["primer_anchor"] =
            serde_json::json!(kv.get("--primer-anchor").copied().unwrap_or("5"));
        receipt["configuration"]["assembly"] = extensions
            .assembly
            .as_ref()
            .map(|cfg| {
                serde_json::json!({
                    "method":"unique-ungapped-suffix-prefix-v1", "min_overlap":cfg.min_overlap,
                    "max_mismatches":cfg.max_mismatches,"max_read_length":cfg.max_read_length,
                    "quality_cap":cfg.quality_cap,"orientation":"R1-plus-reverse-complement-R2",
                    "quality_model":"independent-symmetric-substitution-uniform-prior-v1"
                })
            })
            .unwrap_or(serde_json::Value::Null);
        receipt["counts"]["assembly_no_overlap"] = serde_json::json!(counts.assembly_no_overlap);
        receipt["counts"]["assembly_ambiguous"] = serde_json::json!(counts.assembly_ambiguous);
        receipt["counts"]["assembled"] = serde_json::json!(counts.assembled);
    }
    if kv.contains_key("--experimental-large-group") {
        receipt["schema"] = serde_json::json!("clonoforge-experimental-large-group-v1");
        receipt["experimental"] = serde_json::json!(true);
    }
    let mut file = OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(stage.join("completion.json"))?;
    writeln!(file, "{receipt}")?;
    file.sync_all()?;
    drop(file);
    let replay_receipt: serde_json::Value =
        serde_json::from_reader(File::open(stage.join("completion.json"))?)?;
    if replay_receipt != receipt {
        return Err("paired: completion receipt readback mismatch".into());
    }
    if publication::digest(&qualified)? != prepared || publication::digest(&manifest)? != m3 {
        return Err("paired: staged output changed".into());
    }
    File::open(&stage)?.sync_all()?;
    if std::env::var_os("CLONOFORGE_PROFILE").is_some_and(|v| v == "1") {
        let _ = writeln!(
            std::io::stderr().lock(),
            "{}",
            serde_json::json!({"requested_workers":threads,"actual_workers":req.workers,"sort_chunk_bytes":req.sort_chunk_bytes,"merge_fan_in":req.merge_fan_in,"modeled_bytes":plan.total+extra,"hash_before_ns":hash_before_ns,"prepare_ns":prepare_ns,"consensus_ns":consensus_ns,"prepublish_total_ns":profile_start.elapsed().as_nanos()})
        );
    }
    inject("precommit")?;
    publication::rename_no_replace(&stage, &output)?;
    let _ = writeln!(std::io::stdout().lock(), "{}", output.display());
    Ok(())
}

#[cfg(test)]
mod admission_tests {
    use super::*;
    const MIB: u64 = 1024 * 1024;

    #[test]
    fn growing_chunks_preserves_worker_selection_and_the_full_budget() {
        let raw_allowance = 32 * MIB + 32 * 16384;
        for total_mib in [600, 1024, 2048] {
            for requested_workers in [1, 2, 8] {
                let (base, plan) = Request::select(
                    total_mib * MIB - raw_allowance,
                    requested_workers,
                    16384,
                    50,
                )
                .unwrap();
                let (selected, admitted) = admit_sort_chunk(base, plan);
                let mut expected = base;
                expected.sort_chunk_bytes = selected.sort_chunk_bytes;
                assert_eq!(selected, expected, "only the sort chunk may change");
                assert!(selected.sort_chunk_bytes >= base.sort_chunk_bytes);
                assert!(selected.sort_chunk_bytes <= 64 * MIB);
                assert_eq!(selected.plan().unwrap(), admitted);
                assert!(admitted.total + raw_allowance <= total_mib * MIB);
                if total_mib == 2048 && requested_workers == 8 {
                    assert_eq!(selected.workers, 8);
                    assert_eq!(selected.sort_chunk_bytes, 64 * MIB);
                }
            }
        }
    }

    #[test]
    fn exact_budget_and_large_records_keep_the_original_admission() {
        let (mut base, plan) = Request::select(8 * 1024 * MIB, 1, 16384, 50).unwrap();
        base.budget = plan.total;
        let (selected, admitted) = admit_sort_chunk(base, plan.clone());
        assert_eq!(selected, base);
        assert_eq!(admitted, plan);
        let (large, plan) = Request::select(8 * 1024 * MIB, 2, 2 * MIB, 2).unwrap();
        assert!(large.sort_chunk_bytes > 64 * MIB);
        let (selected, admitted) = admit_sort_chunk(large, plan.clone());
        assert_eq!(selected, large);
        assert_eq!(admitted, plan);
    }
}
