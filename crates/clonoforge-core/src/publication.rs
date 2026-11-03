//! Manifest-last single-writer transaction. Incomplete staging is retained for recovery.
use crate::{
    consensus::{consensus, Config},
    group::Group,
    io::{write_record, Format, Reader},
    model::{invalid, ReadRecord},
};
use serde_json::{json, Value};
use sha2::{Digest, Sha256};
use std::{
    fs::{self, File, OpenOptions},
    io::{BufReader, BufWriter, Read, Write},
    path::{Path, PathBuf},
    time::Instant,
};
#[derive(Clone, Debug, Default)]
pub struct PublicationTimings {
    pub detail_enabled: bool,
    pub serialization_elapsed_sum_ns: u64,
    pub data_flush_ns: u64,
    pub data_fsync_ns: u64,
    pub staged_readback_inclusive_ns: u64,
    pub staged_hash_ns: u64,
    pub manifest_directory_sync_ns: u64,
    pub precommit_cleanup_ns: u64,
    pub staged_write_sync_validate_hash_ns: u64,
    pub manifest_and_publish_ns: u64,
    pub cleanup_error: Option<String>,
    pub retained_stage: Option<PathBuf>,
}
fn timed_ns(t: Instant) -> u64 {
    u64::try_from(t.elapsed().as_nanos()).unwrap_or(u64::MAX)
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Failure {
    Admission,
    Parse,
    SpillCreate,
    SpillWrite,
    SpillFlush,
    SpillSync,
    SpillValidate,
    Merge,
    MergePartial,
    DuplicateDetection,
    Group,
    Worker,
    WorkerEarlyExit,
    OutputOrder,
    Stage,
    WritePass,
    WriteFail,
    FlushPass,
    FlushFail,
    ValidatePass,
    ValidateFail,
    ManifestWrite,
    ManifestFlush,
    ManifestValidate,
    CollisionPass,
    CollisionFail,
    CollisionManifest,
    PublishPass,
    PublishFail,
    DirectorySync,
    ManifestPublish,
    Cleanup,
    StageCleanup,
}
pub(crate) fn inject(at: Option<Failure>, point: Failure) -> std::io::Result<()> {
    if at == Some(point) {
        Err(invalid(&format!("publication: injected {point:?}")))
    } else {
        Ok(())
    }
}
pub fn digest(path: &Path) -> std::io::Result<(u64, String)> {
    let mut f = File::open(path)?;
    let mut hash = Sha256::new();
    let mut buf = [0u8; 65536];
    let mut bytes = 0u64;
    loop {
        let n = f.read(&mut buf)?;
        if n == 0 {
            break;
        }
        bytes = bytes
            .checked_add(n as u64)
            .ok_or_else(|| invalid("publication: size overflow"))?;
        hash.update(&buf[..n]);
    }
    Ok((bytes, format!("{:x}", hash.finalize())))
}
#[cfg(test)]
fn same_file(a: &Path, b: &Path) -> bool {
    #[cfg(unix)]
    {
        use std::os::unix::fs::MetadataExt;
        match (fs::symlink_metadata(a), fs::symlink_metadata(b)) {
            (Ok(a), Ok(b)) => a.dev() == b.dev() && a.ino() == b.ino(),
            _ => false,
        }
    }
    #[cfg(not(unix))]
    {
        let _ = (a, b);
        false
    }
}
fn same_identity(expected: &fs::Metadata, path: &Path) -> bool {
    #[cfg(unix)]
    {
        use std::os::unix::fs::MetadataExt;
        fs::symlink_metadata(path)
            .is_ok_and(|m| m.is_file() && m.dev() == expected.dev() && m.ino() == expected.ino())
    }
    #[cfg(not(unix))]
    {
        let _ = (expected, path);
        false
    }
}
fn same_owned_directory(expected: &fs::Metadata, path: &Path) -> bool {
    #[cfg(unix)]
    {
        use std::os::unix::fs::MetadataExt;
        fs::symlink_metadata(path)
            .is_ok_and(|m| m.is_dir() && m.dev() == expected.dev() && m.ino() == expected.ino())
    }
    #[cfg(not(unix))]
    {
        let _ = (expected, path);
        false
    }
}
#[cfg(any(target_os = "macos", target_os = "linux"))]
pub fn rename_no_replace(from: &Path, to: &Path) -> std::io::Result<()> {
    use rustix::fs::{renameat_with, RenameFlags, CWD};
    renameat_with(CWD, from, CWD, to, RenameFlags::NOREPLACE).map_err(std::io::Error::from)
}
#[cfg(not(any(target_os = "macos", target_os = "linux")))]
pub fn rename_no_replace(_from: &Path, _to: &Path) -> std::io::Result<()> {
    Err(std::io::Error::new(
        std::io::ErrorKind::Unsupported,
        "exclusive marker rename unsupported on this platform",
    ))
}
pub fn validate_prefix(prefix: &str) -> std::io::Result<()> {
    if prefix.is_empty()
        || !prefix
            .bytes()
            .all(|b| b.is_ascii_alphanumeric() || b"_-.".contains(&b))
        || prefix == "."
        || prefix == ".."
    {
        return Err(invalid("configuration: unsafe output prefix"));
    }
    Ok(())
}
pub fn publish(
    groups: &[Group],
    config: &Config,
    out: &Path,
    prefix: &str,
    failure: Option<Failure>,
) -> std::io::Result<PathBuf> {
    publish_stream(config, out, prefix, failure, |index, sink| {
        for group in groups {
            if group.records.len() >= config.min_count {
                if index == 0 {
                    sink(&consensus(&group.records, &group.barcode, config)?)?;
                }
            } else if index == 1 && config.failed {
                for record in &group.records {
                    sink(record)?;
                }
            }
        }
        Ok(())
    })
}

/// `source` must replay the same category on every call; staged readback checks it.
pub fn publish_stream(
    config: &Config,
    out: &Path,
    prefix: &str,
    failure: Option<Failure>,
    source: impl FnMut(usize, &mut dyn FnMut(&ReadRecord) -> std::io::Result<()>) -> std::io::Result<()>,
) -> std::io::Result<PathBuf> {
    publish_stream_profiled(config, out, prefix, failure, source, None)
}
pub fn publish_stream_profiled(
    config: &Config,
    out: &Path,
    prefix: &str,
    failure: Option<Failure>,
    source: impl FnMut(usize, &mut dyn FnMut(&ReadRecord) -> std::io::Result<()>) -> std::io::Result<()>,
    timing: Option<&mut PublicationTimings>,
) -> std::io::Result<PathBuf> {
    publish_stream_precommit(config, out, prefix, failure, source, timing, || Ok(()))
}
/// Runs ownership-aware cleanup after staged categories are validated and hashed,
/// before writing or publishing any completion marker.
#[allow(clippy::too_many_arguments)]
pub fn publish_stream_precommit(
    config: &Config,
    out: &Path,
    prefix: &str,
    failure: Option<Failure>,
    mut source: impl FnMut(
        usize,
        &mut dyn FnMut(&ReadRecord) -> std::io::Result<()>,
    ) -> std::io::Result<()>,
    mut timing: Option<&mut PublicationTimings>,
    precommit: impl FnOnce() -> std::io::Result<()>,
) -> std::io::Result<PathBuf> {
    config.validate()?;
    validate_prefix(prefix)?;
    fs::create_dir_all(out)?;
    // Stable advisory coordination inode; never truncate or unlink it. Locks
    // release on process exit, including a crash, so abandoned staging is retryable.
    let lock = OpenOptions::new()
        .read(true)
        .write(true)
        .create(true)
        .truncate(false)
        .open(out.join(format!(".{prefix}.transaction.lock")))?;
    lock.try_lock()
        .map_err(|e| invalid(&format!("publication: transaction lock unavailable: {e}")))?;
    let names = [
        format!("{prefix}_consensus-pass.fastq"),
        format!("{prefix}_consensus-fail.fastq"),
        format!("{prefix}_completion.json"),
    ];
    for name in &names {
        match fs::symlink_metadata(out.join(name)) {
            Ok(_) => {
                return Err(invalid(
                    "publication: destination exists (including lazy category)",
                ))
            }
            Err(e) if e.kind() == std::io::ErrorKind::NotFound => {}
            Err(e) => return Err(e),
        }
    }
    let parent = out
        .parent()
        .filter(|p| !p.as_os_str().is_empty())
        .unwrap_or(Path::new("."));
    let mut stage = None;
    for attempt in 0..1000 {
        let p = parent.join(format!(
            ".clonoforge-{prefix}.incomplete.{}.{attempt}",
            std::process::id()
        ));
        match fs::create_dir(&p) {
            Ok(()) => {
                stage = Some(p);
                break;
            }
            Err(e) if e.kind() == std::io::ErrorKind::AlreadyExists => {}
            Err(e) => return Err(e),
        }
    }
    let stage = stage.ok_or_else(|| invalid("publication: no free staging name"))?;
    let stage_metadata = fs::symlink_metadata(&stage)?;
    let mut published: Vec<(String, fs::Metadata)> = Vec::new();
    let mut pending: Option<(PathBuf, fs::Metadata)> = None;
    let mut stage_removed = false;
    let result = (|| {
        inject(failure, Failure::Stage)?;
        let mut files: Vec<Value> = Vec::new();
        let mut absent = Vec::new();
        let stage_timer = Instant::now();
        let detailed = timing.as_ref().is_some_and(|t| t.detail_enabled);
        let mut serialization_ns = 0u64;
        let mut data_flush_ns = 0u64;
        let mut data_fsync_ns = 0u64;
        let mut readback_ns = 0u64;
        let mut hash_ns = 0u64;
        let mut commit_sync_ns = 0u64;
        for (index, name) in names[..2].iter().enumerate() {
            let path = stage.join(name);
            let mut file: Option<BufWriter<File>> = None;
            let mut count = 0u64;
            source(index, &mut |record| {
                if file.is_none() {
                    let created = OpenOptions::new()
                        .write(true)
                        .create_new(true)
                        .open(&path)?;
                    inject(
                        failure,
                        if index == 0 {
                            Failure::WritePass
                        } else {
                            Failure::WriteFail
                        },
                    )?;
                    file = Some(BufWriter::with_capacity(64 * 1024, created));
                }
                let start = detailed.then(Instant::now);
                write_record(file.as_mut().unwrap(), record, Format::Fastq)?;
                if let Some(start) = start {
                    serialization_ns = serialization_ns.saturating_add(timed_ns(start));
                }
                count = count
                    .checked_add(1)
                    .ok_or_else(|| invalid("publication: record count overflow"))?;
                Ok(())
            })?;
            let Some(mut file) = file else {
                absent.push(name.clone());
                continue;
            };
            inject(
                failure,
                if index == 0 {
                    Failure::FlushPass
                } else {
                    Failure::FlushFail
                },
            )?;
            let start = detailed.then(Instant::now);
            file.flush()?;
            if let Some(start) = start {
                data_flush_ns = data_flush_ns.saturating_add(timed_ns(start));
            }
            let start = detailed.then(Instant::now);
            file.get_ref().sync_all()?;
            if let Some(start) = start {
                data_fsync_ns = data_fsync_ns.saturating_add(timed_ns(start));
            }
            drop(file);
            inject(
                failure,
                if index == 0 {
                    Failure::ValidatePass
                } else {
                    Failure::ValidateFail
                },
            )?;
            let readback_start = detailed.then(Instant::now);
            let mut reader = Reader::staged_fastq(BufReader::new(File::open(&path)?));
            let mut checked = 0u64;
            source(index, &mut |record| {
                let parsed = reader
                    .next_record()?
                    .ok_or_else(|| invalid("publication: missing staged record"))?;
                if parsed.header() != record.header()
                    || parsed.bases != record.bases
                    || parsed.phred != record.phred
                {
                    return Err(invalid("publication: staged record differs"));
                }
                checked = checked
                    .checked_add(1)
                    .ok_or_else(|| invalid("publication: record count overflow"))?;
                Ok(())
            })?;
            if checked != count || reader.next_record()?.is_some() {
                return Err(invalid("publication: staged record count differs"));
            }
            if let Some(start) = readback_start {
                readback_ns = readback_ns.saturating_add(timed_ns(start));
            }
            let hash_start = detailed.then(Instant::now);
            let (bytes, sha256) = digest(&path)?;
            if let Some(start) = hash_start {
                hash_ns = hash_ns.saturating_add(timed_ns(start));
            }
            files.push(json!({"name":name,"bytes":bytes,"sha256":sha256,"records":count}));
        }
        if let Some(t) = timing.as_deref_mut() {
            t.staged_write_sync_validate_hash_ns = timed_ns(stage_timer);
            t.serialization_elapsed_sum_ns = serialization_ns;
            t.data_flush_ns = data_flush_ns;
            t.data_fsync_ns = data_fsync_ns;
            t.staged_readback_inclusive_ns = readback_ns;
            t.staged_hash_ns = hash_ns;
        }
        let cleanup_start = Instant::now();
        precommit()?;
        if let Some(t) = timing.as_deref_mut() {
            t.precommit_cleanup_ns = timed_ns(cleanup_start);
        }
        let stage_timer = Instant::now();
        let manifest = json!({"schema_version":1,"command":"clonoforge build-consensus","success":true,"files":files,"absent_categories":absent,"configuration":{"barcode_field":config.barcode_field,"min_count":config.min_count,"frequency":config.frequency,"min_quality":config.min_quality,"failed":config.failed,"prefix":prefix,"logical_processes":1}});
        inject(failure, Failure::ManifestWrite)?;
        let mut file = OpenOptions::new()
            .create_new(true)
            .write(true)
            .open(stage.join(&names[2]))?;
        serde_json::to_writer_pretty(&mut file, &manifest)?;
        file.write_all(b"\n")?;
        inject(failure, Failure::ManifestFlush)?;
        file.flush()?;
        let start = detailed.then(Instant::now);
        file.sync_all()?;
        if let Some(start) = start {
            commit_sync_ns = commit_sync_ns.saturating_add(timed_ns(start));
        }
        drop(file);
        inject(failure, Failure::ManifestValidate)?;
        let reread: Value = serde_json::from_reader(File::open(stage.join(&names[2]))?)?;
        if reread != manifest {
            return Err(invalid("publication: staged manifest differs"));
        }
        let start = detailed.then(Instant::now);
        File::open(&stage)?.sync_all()?;
        if let Some(start) = start {
            commit_sync_ns = commit_sync_ns.saturating_add(timed_ns(start));
        }
        for (index, name) in names[..2].iter().enumerate() {
            if stage.join(name).exists() {
                inject(
                    failure,
                    if index == 0 {
                        Failure::PublishPass
                    } else {
                        Failure::PublishFail
                    },
                )?;
                if failure
                    == Some(if index == 0 {
                        Failure::CollisionPass
                    } else {
                        Failure::CollisionFail
                    })
                {
                    OpenOptions::new()
                        .write(true)
                        .create_new(true)
                        .open(out.join(name))?
                        .write_all(b"external-user-file")?;
                }
                let metadata = fs::symlink_metadata(stage.join(name))?;
                fs::hard_link(stage.join(name), out.join(name))?;
                published.push((name.clone(), metadata));
            }
        }
        inject(failure, Failure::DirectorySync)?;
        let start = detailed.then(Instant::now);
        File::open(out)?.sync_all()?;
        if let Some(start) = start {
            commit_sync_ns = commit_sync_ns.saturating_add(timed_ns(start));
        }
        // Move the validated marker into a unique task-owned pending name, then
        // remove staging before atomically consuming that name into the final marker.
        let marker_meta = fs::symlink_metadata(stage.join(&names[2]))?;
        for attempt in 0..1000 {
            let path = out.join(format!(
                ".clonoforge-{prefix}.pending-manifest.{}.{attempt}",
                std::process::id()
            ));
            match fs::hard_link(stage.join(&names[2]), &path) {
                Ok(()) => {
                    pending = Some((path, marker_meta));
                    break;
                }
                Err(e) if e.kind() == std::io::ErrorKind::AlreadyExists => (),
                Err(e) => return Err(e),
            }
        }
        let sync_start = detailed.then(Instant::now);
        File::open(out)?.sync_all()?;
        if let Some(start) = sync_start {
            commit_sync_ns = commit_sync_ns.saturating_add(timed_ns(start));
        }
        let (pending_path, pending_identity) = pending
            .as_ref()
            .ok_or_else(|| invalid("publication: no free pending-marker name"))?;
        if !same_identity(pending_identity, pending_path) {
            return Err(invalid("publication: pending marker ownership changed"));
        }
        if failure == Some(Failure::StageCleanup) {
            return Err(invalid("publication: injected stage cleanup failure"));
        }
        if !same_owned_directory(&stage_metadata, &stage) {
            return Err(invalid("publication: stage ownership changed"));
        }
        fs::remove_dir_all(&stage)?;
        stage_removed = true;
        inject(failure, Failure::ManifestPublish)?;
        if failure == Some(Failure::CollisionManifest) {
            OpenOptions::new()
                .write(true)
                .create_new(true)
                .open(out.join(&names[2]))?
                .write_all(b"external-user-file")?;
        }
        rename_no_replace(pending_path, &out.join(&names[2]))?;
        if let Some(t) = timing {
            t.manifest_and_publish_ns = timed_ns(stage_timer);
            t.manifest_directory_sync_ns = commit_sync_ns;
        }
        // The exclusive rename is the final fallible operation; no cleanup follows.
        Ok(out.join(&names[2]))
    })();
    if result.is_err() {
        for (name, identity) in &published {
            if same_identity(identity, &out.join(name)) {
                let _ = fs::remove_file(out.join(name));
            }
        }
        if !stage_removed {
            if let Some((path, identity)) = &pending {
                if same_identity(identity, path) {
                    let _ = fs::remove_file(path);
                }
            }
        }
    }
    result.map_err(|e| {
        if stage_removed {
            if let Some((path, identity)) = &pending {
                if same_identity(identity, path) {
                    return invalid(&format!(
                        "publication: {e}; incomplete pending marker retained {}",
                        path.display()
                    ));
                }
            }
        }
        invalid(&format!(
            "publication: {e}; stage path {} (may already be removed)",
            stage.display()
        ))
    })
}

#[cfg(all(test, unix))]
mod ownership_tests {
    use super::*;
    #[test]
    fn replacement_symlink_is_not_owned_link() {
        let root = std::env::temp_dir().join(format!("clonoforge-owner-{}", std::process::id()));
        fs::create_dir(&root).unwrap();
        let a = root.join("stage");
        let b = root.join("public");
        fs::write(&a, b"data").unwrap();
        std::os::unix::fs::symlink(&a, &b).unwrap();
        assert!(!same_file(&a, &b));
        let metadata = fs::symlink_metadata(&a).unwrap();
        assert!(!same_identity(&metadata, &b));
        fs::remove_file(&b).unwrap();
        fs::hard_link(&a, &b).unwrap();
        assert!(same_identity(&metadata, &b));
        fs::remove_dir_all(root).unwrap();
    }
}
