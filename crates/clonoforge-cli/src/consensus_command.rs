use clonoforge_core::{
    admission::Request, bounded_workflow, consensus::Config, external::Limits, publication,
};
use std::{
    collections::BTreeSet,
    io::{BufRead, Write},
    path::PathBuf,
};
pub fn run(args: &[String]) -> Result<(), Box<dyn std::error::Error>> {
    if args == ["--help"] {
        println!("clonoforge build-consensus INPUT --outdir DIR [--outname PREFIX] [--bf KEY] [-n COUNT] [--freq FRACTION] [-q QUALITY] [--failed] [--force-spill] [--threads N] [--memory-mib MIB] [--max-record-bytes N] [--max-group-count N] [--workdir DIR]\nBuild consensus from a prepared FASTQ with barcode annotations. Plain or gzip input is detected by content. Reads in each group must have equal lengths and uppercase ACGTN bases. Standard groups contain 1 to 50 reads. Defaults: 1 worker, 512 MiB modeled memory, 16384 record bytes, and 50 reads per group. Temporary files default to the output directory. Modeled memory is not an exact RSS limit. The record limit counts header, sequence and quality bytes. See docs/USER_MANUAL.md.");
        return Ok(());
    }
    let input = args
        .first()
        .filter(|s| !s.starts_with('-'))
        .ok_or("configuration: INPUT path required")?;
    let mut config = Config::default();
    let mut output = None;
    let mut prefix = None;
    let mut threads = 1u64;
    let mut memory_mib = 512u64;
    let mut max_record = 16384u64;
    let mut max_group = 50u64;
    let mut workdir = None;
    let mut force_spill = false;
    let mut seen = BTreeSet::new();
    let mut i = 1;
    while i < args.len() {
        let flag = &args[i];
        if !seen.insert(flag) {
            return Err("configuration: repeated option".into());
        }
        if flag == "--failed" {
            config.failed = true;
            i += 1;
            continue;
        }
        if flag == "--force-spill" {
            force_spill = true;
            i += 1;
            continue;
        }
        if ![
            "--bf",
            "-n",
            "--freq",
            "-q",
            "--outdir",
            "--outname",
            "--threads",
            "--memory-mib",
            "--max-record-bytes",
            "--max-group-count",
            "--workdir",
        ]
        .contains(&flag.as_str())
        {
            return Err(format!("configuration: unsupported option {flag}").into());
        }
        let value = args
            .get(i + 1)
            .ok_or("configuration: option value missing")?;
        match flag.as_str() {
            "--bf" => config.barcode_field = value.clone(),
            "-n" => {
                config.min_count = value
                    .parse()
                    .map_err(|_| "configuration: invalid minimum count")?
            }
            "--freq" => {
                config.frequency = value
                    .parse()
                    .map_err(|_| "configuration: invalid frequency")?
            }
            "-q" => {
                config.min_quality = value
                    .parse()
                    .map_err(|_| "configuration: invalid minimum quality")?
            }
            "--outdir" => output = Some(PathBuf::from(value)),
            "--outname" => prefix = Some(value.clone()),
            "--threads" => {
                threads = value
                    .parse()
                    .map_err(|_| "configuration: invalid threads")?
            }
            "--memory-mib" => {
                memory_mib = value
                    .parse()
                    .map_err(|_| "configuration: invalid memory budget")?
            }
            "--max-record-bytes" => {
                max_record = value
                    .parse()
                    .map_err(|_| "configuration: invalid record limit")?
            }
            "--max-group-count" => {
                max_group = value
                    .parse()
                    .map_err(|_| "configuration: invalid group limit")?;
            }
            "--workdir" => workdir = Some(PathBuf::from(value)),
            _ => unreachable!(),
        }
        i += 2;
    }
    config.validate()?;
    let output = output.ok_or("configuration: --outdir required")?;
    let input = PathBuf::from(input);
    let prefix = prefix
        .or_else(|| {
            input
                .file_stem()
                .and_then(|s| s.to_str())
                .map(str::to_owned)
        })
        .ok_or("configuration: UTF-8 input stem or --outname required")?;
    publication::validate_prefix(&prefix)?;
    let budget = memory_mib
        .checked_mul(1024 * 1024)
        .ok_or("configuration: memory budget overflow")?;
    let profile_enabled = std::env::var_os("CLONOFORGE_PROFILE").is_some_and(|v| v == "1");
    let (mut request, mut plan) = Request::select(budget, threads, max_record, max_group)?;
    {
        let mut expanded = request;
        expanded.sort_chunk_bytes = expanded.sort_chunk_bytes.max(64 * 1024 * 1024);
        if let Ok(p) = expanded.plan() {
            request = expanded;
            plan = p;
        }
    }
    let mut probe = std::io::BufReader::new(std::fs::File::open(&input)?);
    let input_mode = if probe.fill_buf()?.starts_with(&[0x1f, 0x8b]) {
        "gzip"
    } else {
        "plain"
    };
    let limits = Limits {
        chunk_bytes: usize::try_from(request.sort_chunk_bytes)?,
        max_record_bytes: usize::try_from(request.max_record_bytes)?,
        fan_in: usize::try_from(request.merge_fan_in)?,
    };
    std::fs::create_dir_all(&output)?;
    let workdir = workdir.unwrap_or_else(|| output.clone());
    eprintln!(
        "input: {input_mode}; modeled memory {} bytes, limiting phase {}, workers {}",
        plan.total, plan.limiting_phase, request.workers
    );
    let mut timings = bounded_workflow::StageTimings::default();
    timings.input.detail_enabled = profile_enabled;
    timings.output.detail_enabled = profile_enabled;
    timings.publication.detail_enabled = profile_enabled;
    let (manifest, spilled) = if force_spill {
        (
            bounded_workflow::run_configured_profiled(
                &input,
                &config,
                &output,
                &prefix,
                &workdir,
                (
                    &limits,
                    usize::try_from(plan.group_weight)?,
                    usize::try_from(request.workers)?,
                    usize::try_from(request.window)?,
                    usize::try_from(max_group)?,
                ),
                None,
                Some(&mut timings),
            )?,
            true,
        )
    } else {
        bounded_workflow::run_adaptive_profiled(
            &input,
            &config,
            &output,
            &prefix,
            &workdir,
            (
                &limits,
                usize::try_from(plan.group_weight)?,
                usize::try_from(request.workers)?,
                usize::try_from(request.window)?,
                usize::try_from(max_group)?,
            ),
            None,
            Some(&mut timings),
        )?
    };
    let profile = if profile_enabled {
        Some(
            serde_json::json!({"execution_path":if force_spill {"forced-spill"} else if spilled {"seeded-spill"} else {"no-spill"},"consensus_measurement":if spilled {"one kernel execution per group"} else {"kernel elapsed sum includes write and readback replays"},"input_sort_ns":timings.input_sort_ns,"group_and_output_order_ns":timings.group_and_output_order_ns,"publication_ns":timings.publication_ns,"total_ns":timings.total_ns,"read_parse_header_spill_ns":timings.input.read_parse_header_spill_ns,"header_finalize_verify_ns":timings.input.header_finalize_verify_ns,"duplicate_scan_barcode_spill_ns":timings.input.duplicate_scan_barcode_spill_ns,"barcode_finalize_verify_ns":timings.input.barcode_finalize_verify_ns,"group_pipeline_ns":if spilled {Some(timings.output.group_pipeline_ns)} else {None},"output_finalize_merge_ns":timings.output.output_finalize_merge_ns,"staged_write_sync_validate_hash_ns":timings.publication.staged_write_sync_validate_hash_ns,"manifest_and_publish_ns":timings.publication.manifest_and_publish_ns,
"record_read_parse_inclusive_ns":timings.input.record_read_parse_inclusive_ns,"raw_file_read_ns":timings.input.raw_file_read_ns,"decoded_read_inclusive_ns":timings.input.decoded_read_inclusive_ns,"raw_bytes":timings.input.raw_bytes,"decoded_bytes":timings.input.decoded_bytes,"compressed":timings.input.compressed,
"input_initial_run_build_verify_ns":timings.input.spill.initial_run_build_verify_ns,"input_merge_pass_inclusive_ns":timings.input.spill.merge_pass_inclusive_ns,"input_spill_flush_ns":timings.input.spill.flush_ns,"input_spill_sync_ns":timings.input.spill.sync_ns,"input_initial_runs":timings.input.spill.initial_runs,"input_merge_passes":timings.input.spill.merge_passes,
"consensus_group_elapsed_sum_ns":timings.output.consensus_group_elapsed_sum_ns,"consensus_groups":timings.output.consensus_groups,"output_initial_run_build_verify_ns":timings.output.spill.initial_run_build_verify_ns,"output_merge_pass_inclusive_ns":timings.output.spill.merge_pass_inclusive_ns,"output_spill_flush_ns":timings.output.spill.flush_ns,"output_spill_sync_ns":timings.output.spill.sync_ns,"output_initial_runs":timings.output.spill.initial_runs,"output_merge_passes":timings.output.spill.merge_passes,
"serialization_elapsed_sum_ns":timings.publication.serialization_elapsed_sum_ns,"data_flush_ns":timings.publication.data_flush_ns,"data_fsync_ns":timings.publication.data_fsync_ns,"staged_readback_inclusive_ns":timings.publication.staged_readback_inclusive_ns,"staged_hash_ns":timings.publication.staged_hash_ns,"manifest_directory_sync_ns":timings.publication.manifest_directory_sync_ns,"precommit_cleanup_ns":timings.publication.precommit_cleanup_ns}),
        )
    } else {
        None
    };
    let report = serde_json::json!({"profile":profile,"spilled":spilled,"forced_spill":force_spill,"success":true,"manifest":manifest,"input_mode":input_mode,"workdir":workdir,"requested_workers":threads,"workers":request.workers,"memory_budget_bytes":budget,"modeled_memory_bytes":plan.total,"limiting_phase":plan.limiting_phase,"max_record_bytes":max_record,"max_group_count":max_group,"merge_fan_in":request.merge_fan_in,"sort_chunk_bytes":request.sort_chunk_bytes,"window":request.window,"spill_cleanup_error":timings.cleanup_error,"retained_spill":timings.retained_spill,"stage_cleanup_error":timings.publication.cleanup_error,"retained_stage":timings.publication.retained_stage});
    let _ = writeln!(std::io::stdout().lock(), "{report}");
    Ok(())
}
