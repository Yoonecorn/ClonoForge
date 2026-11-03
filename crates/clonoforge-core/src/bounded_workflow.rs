//! Experimental external-sort workflow; caller must apply an admission model first.
use crate::{
    consensus::Config,
    external::{self, Limits},
    model::invalid,
    publication::{self, Failure},
};
use std::{
    fs,
    path::{Path, PathBuf},
    time::Instant,
};

fn cleanup_owned(
    dir: &Path,
    metadata: &fs::Metadata,
    failure: Option<Failure>,
) -> std::io::Result<()> {
    publication::inject(failure, Failure::Cleanup)?;
    #[cfg(unix)]
    let still_owned = {
        use std::os::unix::fs::MetadataExt;
        fs::symlink_metadata(dir)
            .is_ok_and(|m| m.is_dir() && m.dev() == metadata.dev() && m.ino() == metadata.ino())
    };
    #[cfg(not(unix))]
    let still_owned = false;
    if !still_owned {
        return Err(invalid("spill directory ownership changed"));
    }
    fs::remove_dir_all(dir)
}
/// The selected work root must already exist. A unique child is exclusively owned;
/// it is retained on any failure and removed after successful publication.
pub fn run(
    input: &Path,
    cfg: &Config,
    out: &Path,
    prefix: &str,
    root: &Path,
    resource: (&Limits, usize, usize, usize),
    failure: Option<Failure>,
) -> std::io::Result<PathBuf> {
    let (limits, max_group_bytes, workers, window) = resource;
    run_configured(
        input,
        cfg,
        out,
        prefix,
        root,
        (limits, max_group_bytes, workers, window, 50),
        failure,
    )
}
/// Times coarse complete stages; values are not stored in the completion manifest.
#[derive(Clone, Debug, Default)]
pub struct StageTimings {
    pub input_sort_ns: u64,
    pub group_and_output_order_ns: u64,
    pub publication_ns: u64,
    pub total_ns: u64,
    pub input: external::InputTimings,
    pub output: external::OutputTimings,
    pub publication: publication::PublicationTimings,
    pub cleanup_error: Option<String>,
    pub retained_spill: Option<PathBuf>,
}
fn elapsed(t: Instant) -> u64 {
    u64::try_from(t.elapsed().as_nanos()).unwrap_or(u64::MAX)
}
pub fn run_configured(
    input: &Path,
    cfg: &Config,
    out: &Path,
    prefix: &str,
    root: &Path,
    resource: (&Limits, usize, usize, usize, usize),
    failure: Option<Failure>,
) -> std::io::Result<PathBuf> {
    run_configured_profiled(input, cfg, out, prefix, root, resource, failure, None)
}
/// Instrumentation is opt-in and intentionally execution metadata, not part of
/// deterministic scientific output or completion manifest.
#[allow(clippy::too_many_arguments)]
pub fn run_configured_profiled(
    input: &Path,
    cfg: &Config,
    out: &Path,
    prefix: &str,
    root: &Path,
    resource: (&Limits, usize, usize, usize, usize),
    failure: Option<Failure>,
    profile: Option<&mut StageTimings>,
) -> std::io::Result<PathBuf> {
    run_configured_impl(
        input, cfg, out, prefix, root, resource, failure, profile, None, false,
    )
}

/// Consume a verified barcode run from paired preprocessing after global identity
/// validation and whole-group eligibility. Raw and qualified input ordinals have
/// the same strict order, so output ordering and M3 publication bytes agree.
/// This is an internal trusted-input boundary; arbitrary raw FASTQ must use
/// `run_configured` and its complete duplicate checks instead.
#[allow(clippy::too_many_arguments)]
pub fn run_presorted_configured(
    sorted: &Path,
    cfg: &Config,
    out: &Path,
    prefix: &str,
    root: &Path,
    resource: (&Limits, usize, usize, usize, usize),
) -> std::io::Result<PathBuf> {
    run_configured_impl(
        sorted,
        cfg,
        out,
        prefix,
        root,
        resource,
        None,
        None,
        Some(sorted),
        false,
    )
}
#[cfg(feature = "experimental-large")]
pub fn run_presorted_experimental(
    sorted: &Path,
    cfg: &Config,
    out: &Path,
    prefix: &str,
    root: &Path,
    resource: (&Limits, usize, usize, usize, usize),
) -> std::io::Result<PathBuf> {
    if ![100, 250, 500].contains(&resource.4) {
        return Err(invalid("grouping: unsupported experimental limit"));
    }
    run_configured_impl(
        sorted,
        cfg,
        out,
        prefix,
        root,
        resource,
        None,
        None,
        Some(sorted),
        true,
    )
}

#[allow(clippy::too_many_arguments)]
fn run_configured_impl(
    input: &Path,
    cfg: &Config,
    out: &Path,
    prefix: &str,
    root: &Path,
    resource: (&Limits, usize, usize, usize, usize),
    failure: Option<Failure>,
    mut profile: Option<&mut StageTimings>,
    presorted: Option<&Path>,
    experimental: bool,
) -> std::io::Result<PathBuf> {
    let total_start = Instant::now();
    let (limits, max_group_bytes, workers, window, max_group_count) = resource;
    publication::inject(failure, Failure::Admission)?;
    if max_group_count == 0 || max_group_count > if experimental { 500 } else { 50 } {
        return Err(invalid("grouping: invalid configured read limit"));
    }
    cfg.validate()?;
    publication::validate_prefix(prefix)?;
    limits.validate()?;
    let mut owned = None;
    for attempt in 0u64..1000 {
        let path = root.join(format!(
            ".clonoforge-{prefix}-spill-incomplete-{}-{attempt}",
            std::process::id()
        ));
        match fs::create_dir(&path) {
            Ok(()) => {
                owned = Some(path);
                break;
            }
            Err(e) if e.kind() == std::io::ErrorKind::AlreadyExists => (),
            Err(e) => return Err(e),
        }
    }
    let dir = owned.ok_or_else(|| invalid("external: no available owned scratch directory"))?;
    let owned_metadata = fs::symlink_metadata(&dir)?;
    let result = (|| {
        let stage_start = Instant::now();
        let groups = if let Some(run) = presorted {
            run.to_path_buf()
        } else {
            external::sort_input_profiled(
                input,
                cfg,
                &dir,
                limits,
                failure,
                profile.as_deref_mut().map(|p| &mut p.input),
            )?
        };
        if let Some(p) = profile.as_deref_mut() {
            p.input_sort_ns = elapsed(stage_start);
        }
        let stage_start = Instant::now();
        let sort_outputs = if experimental {
            #[cfg(feature = "experimental-large")]
            {
                external::sort_outputs_experimental
            }
            #[cfg(not(feature = "experimental-large"))]
            {
                external::sort_outputs_profiled
            }
        } else {
            external::sort_outputs_profiled
        };
        let (pass, fail) = sort_outputs(
            &groups,
            cfg,
            &dir,
            limits,
            external::GroupExecution {
                max_group_bytes,
                workers,
                window,
                max_group_count,
                failure,
            },
            profile.as_deref_mut().map(|p| &mut p.output),
        )?;
        if let Some(p) = profile.as_deref_mut() {
            p.group_and_output_order_ns = elapsed(stage_start);
        }
        let stage_start = Instant::now();
        let result = publication::publish_stream_precommit(
            cfg,
            out,
            prefix,
            failure,
            |category, sink| {
                let run = if category == 0 { &pass } else { &fail };
                external::for_each_output(run, category == 0, cfg, limits, |record| sink(&record))
            },
            profile.as_deref_mut().map(|p| &mut p.publication),
            || cleanup_owned(&dir, &owned_metadata, failure),
        );
        if let Some(p) = profile.as_deref_mut() {
            p.publication_ns = elapsed(stage_start);
        }
        result
    })();
    match result {
        Ok(p) => {
            if let Some(t) = profile {
                t.total_ns = elapsed(total_start);
            }
            Ok(p)
        }
        Err(e) => Err(invalid(&format!(
            "bounded workflow: {e}; spill path {} (may have been cleaned before publication)",
            dir.display()
        ))),
    }
}

/// Experimental single-pass adaptive path.
pub fn run_adaptive(
    input: &Path,
    cfg: &Config,
    out: &Path,
    prefix: &str,
    root: &Path,
    resource: (&Limits, usize, usize, usize, usize),
) -> std::io::Result<(PathBuf, bool)> {
    run_adaptive_profiled(input, cfg, out, prefix, root, resource, None, None)
}
#[allow(clippy::too_many_arguments)]
pub fn run_adaptive_profiled(
    input: &Path,
    cfg: &Config,
    out: &Path,
    prefix: &str,
    root: &Path,
    resource: (&Limits, usize, usize, usize, usize),
    failure: Option<Failure>,
    mut profile: Option<&mut StageTimings>,
) -> std::io::Result<(PathBuf, bool)> {
    let total = Instant::now();
    let (limits, group_bytes, workers, window, max_count) = resource;
    publication::inject(failure, Failure::Admission)?;
    cfg.validate()?;
    publication::validate_prefix(prefix)?;
    limits.validate()?;
    publication::inject(failure, Failure::Parse)?;
    let trial_start = Instant::now();
    let detailed = profile.as_ref().is_some_and(|t| t.input.detail_enabled);
    let (source, metrics, compressed) = if detailed {
        let (s, m, c) = crate::io::input_timed(input)?;
        (s, Some(m), c)
    } else {
        (crate::io::input(input)?, None, false)
    };
    let source = crate::io::Reader::fastq_limited(source, limits.max_record_bytes)?;
    let mut trial_parse_ns = 0u64;
    let trial = crate::group::consume_capped(
        source,
        cfg,
        limits.chunk_bytes,
        max_count,
        group_bytes,
        failure,
        detailed.then_some(&mut trial_parse_ns),
    )?;
    if let Some(t) = profile.as_deref_mut() {
        t.input_sort_ns = elapsed(trial_start);
        t.input.record_read_parse_inclusive_ns = trial_parse_ns;
        t.input.read_parse_header_spill_ns = t.input_sort_ns;
    }
    match trial {
        crate::group::GroupAttempt::Complete(groups) => {
            if let Some(t) = profile.as_deref_mut() {
                t.output.consensus_groups = u64::try_from(
                    groups
                        .iter()
                        .filter(|g| g.records.len() >= cfg.min_count)
                        .count(),
                )
                .map_err(|_| invalid("group count overflow"))?;
                if let Some(m) = &metrics {
                    use std::sync::atomic::Ordering;
                    t.input.raw_bytes = m.raw_bytes.load(Ordering::Relaxed);
                    t.input.decoded_bytes = m.decoded_bytes.load(Ordering::Relaxed);
                    t.input.raw_file_read_ns = m.raw_ns.load(Ordering::Relaxed);
                    t.input.decoded_read_inclusive_ns = m.decoded_ns.load(Ordering::Relaxed);
                    t.input.compressed = compressed;
                }
            }
            publication::inject(failure, Failure::OutputOrder)?;
            let detailed = profile.as_ref().is_some_and(|t| t.output.detail_enabled);
            let kernel_sum = std::sync::atomic::AtomicU64::new(0);
            let stage = Instant::now();
            let manifest = publication::publish_stream_profiled(
                cfg,
                out,
                prefix,
                failure,
                |category, sink| {
                    let jobs = groups
                        .iter()
                        .map(|g| Ok((g.barcode.clone(), g.records.clone())));
                    crate::ordered_executor::ordered_map_inner(
                        jobs,
                        workers,
                        window,
                        |(barcode, records)| {
                            publication::inject(failure, Failure::Worker)?;
                            if records.len() >= cfg.min_count {
                                if category == 0 {
                                    let start = detailed.then(Instant::now);
                                    let record =
                                        crate::consensus::consensus(&records, &barcode, cfg)?;
                                    if let Some(start) = start {
                                        let ns = elapsed(start);
                                        let _ = kernel_sum.fetch_update(
                                            std::sync::atomic::Ordering::Relaxed,
                                            std::sync::atomic::Ordering::Relaxed,
                                            |v| Some(v.saturating_add(ns)),
                                        );
                                    }
                                    return Ok(vec![record]);
                                }
                            } else if category == 1 && cfg.failed {
                                return Ok(records);
                            }
                            Ok(Vec::new())
                        },
                        |records| {
                            for r in records {
                                sink(&r)?;
                            }
                            Ok(())
                        },
                        if failure == Some(Failure::WorkerEarlyExit) {
                            2
                        } else {
                            0
                        },
                    )
                },
                profile.as_deref_mut().map(|t| &mut t.publication),
            )?;
            if let Some(t) = profile {
                t.output.consensus_group_elapsed_sum_ns =
                    kernel_sum.load(std::sync::atomic::Ordering::Relaxed);
                t.publication_ns = elapsed(stage);
                t.total_ns = elapsed(total);
            }
            Ok((manifest, false))
        }
        crate::group::GroupAttempt::Spill {
            reader,
            groups,
            pending,
        } => {
            let mut owned = None;
            for attempt in 0..1000 {
                let p = root.join(format!(
                    ".clonoforge-{prefix}-spill-incomplete-{}-{attempt}",
                    std::process::id()
                ));
                match fs::create_dir(&p) {
                    Ok(()) => {
                        owned = Some(p);
                        break;
                    }
                    Err(e) if e.kind() == std::io::ErrorKind::AlreadyExists => (),
                    Err(e) => return Err(e),
                }
            }
            let dir =
                owned.ok_or_else(|| invalid("external: no available owned scratch directory"))?;
            let metadata = fs::symlink_metadata(&dir)?;
            let result = (|| {
                let stage = Instant::now();
                let sorted = external::sort_reader_with_seed(
                    reader,
                    groups,
                    pending,
                    cfg,
                    &dir,
                    limits,
                    failure,
                    profile.as_deref_mut().map(|t| &mut t.input),
                )?;
                if let Some(t) = profile.as_deref_mut() {
                    t.input_sort_ns = t.input_sort_ns.saturating_add(elapsed(stage));
                    if let Some(m) = &metrics {
                        use std::sync::atomic::Ordering;
                        t.input.raw_bytes = m.raw_bytes.load(Ordering::Relaxed);
                        t.input.decoded_bytes = m.decoded_bytes.load(Ordering::Relaxed);
                        t.input.raw_file_read_ns = m.raw_ns.load(Ordering::Relaxed);
                        t.input.decoded_read_inclusive_ns = m.decoded_ns.load(Ordering::Relaxed);
                        t.input.compressed = compressed;
                    }
                }
                let stage = Instant::now();
                let (pass, fail) = external::sort_outputs_profiled(
                    &sorted,
                    cfg,
                    &dir,
                    limits,
                    external::GroupExecution {
                        max_group_bytes: group_bytes,
                        workers,
                        window,
                        max_group_count: max_count,
                        failure,
                    },
                    profile.as_deref_mut().map(|t| &mut t.output),
                )?;
                if let Some(t) = profile.as_deref_mut() {
                    t.group_and_output_order_ns = elapsed(stage);
                }
                let stage = Instant::now();
                let result = publication::publish_stream_precommit(
                    cfg,
                    out,
                    prefix,
                    failure,
                    |category, sink| {
                        let run = if category == 0 { &pass } else { &fail };
                        external::for_each_output(run, category == 0, cfg, limits, |r| sink(&r))
                    },
                    profile.as_deref_mut().map(|t| &mut t.publication),
                    || cleanup_owned(&dir, &metadata, failure),
                );
                if let Some(t) = profile.as_deref_mut() {
                    t.publication_ns = elapsed(stage);
                }
                result
            })();
            match result {
                Ok(path)=>{if let Some(t)=profile{t.total_ns=elapsed(total);} Ok((path,true))},
                Err(e)=>Err(invalid(&format!("bounded workflow: {e}; spill path {} (may have been cleaned before publication)",dir.display()))),
            }
        }
    }
}
