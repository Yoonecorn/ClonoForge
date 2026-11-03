mod consensus_command;
mod paired_command;
use clonoforge_core::io::{input, write_record, Format, Reader, StagedOutput};
use std::path::Path;
fn run() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = std::env::args().skip(1).collect();
    if args.as_slice() == ["--version"] {
        println!("clonoforge {}", env!("CARGO_PKG_VERSION"));
        return Ok(());
    }
    if args.is_empty() || args.as_slice() == ["--help"] {
        println!("clonoforge: prepare paired reads and build barcode consensus\nPaired reads: clonoforge preprocess-paired --help\nPrepared FASTQ: clonoforge build-consensus --help\nSequence files: clonoforge <validate|inspect|normalize|compare> <fasta|fastq> INPUT [OUTPUT]\nnormalize requires a new OUTPUT path. These commands are not pRESTO command-compatible.");
        return Ok(());
    }
    if args.first().is_some_and(|s| s == "preprocess-paired") {
        return paired_command::run(&args[1..]);
    }
    if args.first().is_some_and(|s| s == "build-consensus") {
        return consensus_command::run(&args[1..]);
    }
    if args.len() < 3 {
        return Err("expected command, format, and input".into());
    }
    let format = match args[1].as_str() {
        "fasta" => Format::Fasta,
        "fastq" => Format::Fastq,
        _ => return Err("format must be fasta or fastq".into()),
    };
    if args[0] == "compare" {
        if args.len() != 4 {
            return Err("compare requires format and two input paths".into());
        }
        let differences = clonoforge_core::compare::compare(
            input(Path::new(&args[2]))?,
            input(Path::new(&args[3]))?,
            format,
        );
        for difference in &differences {
            println!("{}\t{}", difference.category, difference.id);
        }
        if !differences.is_empty() {
            return Err("structural differences found".into());
        }
        return Ok(());
    }
    let normalize = args[0] == "normalize";
    if !["normalize", "validate", "inspect"].contains(&args[0].as_str())
        || args.len() != if normalize { 4 } else { 3 }
    {
        return Err("invalid command or argument count; use --help".into());
    }
    let mut reader = Reader::new(input(Path::new(&args[2]))?, format);
    let mut output = if normalize {
        Some(StagedOutput::new(Path::new(&args[3]))?)
    } else {
        None
    };
    let mut count = 0u64;
    let mut bases = 0u64;
    while let Some(record) = reader.next_record()? {
        count = count.checked_add(1).ok_or("record count overflow")?;
        bases = bases
            .checked_add(record.bases.len() as u64)
            .ok_or("base count overflow")?;
        if let Some(out) = output.as_mut() {
            write_record(out, &record, format)?;
        }
    }
    if let Some(out) = output {
        out.commit()?;
    }
    println!("records={count}\tbases={bases}");
    Ok(())
}
fn main() {
    if let Err(error) = run() {
        eprintln!("clonoforge: {error}");
        std::process::exit(2);
    }
}
