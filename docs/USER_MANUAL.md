# ClonoForge 1.0.0 user manual

This guide describes ClonoForge 1.0.0. It has three parts: start with a small example, choose settings for your assay, and consult the technical reference when needed.

## 1. Start here

### What ClonoForge does

ClonoForge prepares paired FASTQ reads and builds consensus sequences from reads with the same exact barcode. A consensus summarizes the observed bases and qualities in a group. It can also build consensus from a prepared FASTQ. A barcode is a short sequence used as an identifier. A unique molecular identifier, or UMI, is a barcode intended to label an original molecule.

ClonoForge does not correct barcode errors, discover your assay layout, align reads generally, identify clones, or run a complete pRESTO or repertoire analysis. Identical barcodes do not prove that reads came from one molecule: barcode collisions can merge molecules, and barcode sequencing errors can split them. Single-read groups can pass the default settings; they do not provide agreement between independent reads.

### Install and try it

Use the [Quick Start](QUICK_START.md) to check Rust, Cargo and Python 3, build the source package, and run a synthetic example. The package has no prebuilt executable. The checked system was macOS on Apple ARM64; validate other systems for your work. A first build may need internet access to obtain locked dependencies.

The tutorial creates a fresh working directory. The paired verifier should report PASS, three qualified reads and two groups. That checks the files and their counts; it does not prove biological accuracy. If a step fails, stop and use [Troubleshooting](TROUBLESHOOTING.md).

### Before using your samples

Find your assay protocol and confirm where the barcode and target sequence occur, their orientation, the pair-header format, and any primer or quality filters. Put the barcode-bearing file first when you run the paired workflow, even if it is biological read 2. The two inputs must be synchronized. ClonoForge does not find matching mates by searching ahead.

For a first check, use the standard limit of 50 reads per group, one worker, and no assembly. Choose the barcode, target, trimming and filters from your assay, not from the tutorial. Test a small representative subset, inspect exclusions and verify the output before analyzing the full experiment. These steps are checks, not proof of biological validity.

### Recognize a successful paired run

The requested output directory contains `completion.json`, `qualified.fastq`, and `m3_completion.json`. A passing consensus file, `m3_consensus-pass.fastq`, is present when produced. `qualified.fastq` contains retained targets with barcode annotations. Review the outer receipt's counts for invalid keys, filtering and whole-group exclusions. Check the files using the supplied paired verifier, as shown in the Quick Start. The verifier checks standard and extended paired receipts, but does not support the experimental large-group receipt. A successful run and a verifier PASS do not establish biological accuracy.

## 2. Choose settings

### Match inputs to your assay

Provide two paired FASTQ files, plain or gzip. Each FASTQ record must have four physical lines and matching sequence and quality lengths. Put the barcode-bearing file first. The barcode is a substring of that file; its first base is position zero. The target can be either input. The program checks mate identity and rejects duplicates. It does not repair input order.

Choose the header mode that matches the files. `identical` requires identical full headers. In `slash` mode the first and second tokens end in `/1` and `/2`, with matching trailing descriptions. In `sra` mode the files need the supported common spot identifier and numbered mate descriptions. Simply swapping numbered biological mates may violate the selected mode; verify identities and prepare valid synchronized headers using your laboratory's established procedure.

The paired workflow excludes non-ACGT barcodes. Target bases for consensus must be uppercase A, C, G, T or N, with Phred+33 qualities. In the standard workflow a group contains at most 50 reads, and all target reads in a group must have the same length. Paired processing excludes whole groups that are oversized or have unequal target lengths; oversize takes precedence if both apply. Review the counts before drawing conclusions.

### Required choices

| Choice | What to decide |
|---|---|
| Output directory | Use a new directory name. Do not overwrite earlier results. |
| Header mode | Select the mode that describes the synchronized input headers. |
| Barcode start and length | Locate the exact barcode in the first input. Start counting at zero. |
| Target | Select the first or second input based on the assay. |
| Target trimming | Specify how many bases to remove from the beginning of the target. Zero is allowed. |

The technical reference lists the corresponding options and defaults. The synthetic command in the Quick Start is a tested example; replace its choices only after checking your assay. Use a new output path for every run.

### Optional choices

The default paired settings use one requested worker, a 1024 MiB modeled memory budget and a maximum of 50 reads per group. The actual worker count can be lower. Memory is a planning model, not an operating-system limit. Provide enough temporary disk space. The tutorial explicitly requests 2048 MiB for its example.

Exact anchored primer trimming and simple length, N-count and mean-quality filters are optional. Set thresholds for your study only after checking their consequences. Optional ungapped overlap assembly requires the first target and a single eligible overlap with the reverse-complemented second input. Use it only when both mates contain overlapping target sequence and you have validated orientation and thresholds. Its quality model assumes independent substitution errors and is not empirically calibrated; it is not reference-guided assembly.

Separately compiled experimental paired profiles allow explicit limits of 100, 250 or 500. These have finite qualification, are not the default and do not establish broad biological validity. The standard verifier does not validate their experimental receipt. See the reference before considering them.

### Record your work

Keep original inputs and checksums, source commit, version, build features, full command, assay layout, settings, completion receipts, exclusion counts and software environment. Protect temporary and output sequence files according to your institution's data rules. Store backups appropriate to your data.

For citation, identify ClonoForge 1.0.0 and the exact commit used. The package's `CITATION.cff` identifies Yoon Sing Yap and provides citation metadata. Cite your data sources and, where relevant, Vander Heiden et al., *Bioinformatics* (2014), doi:10.1093/bioinformatics/btu138, for pRESTO.

Adapt this Methods statement to your actual work: “FASTQ reads were processed with ClonoForge 1.0.0 (commit [commit]; build features [features]). Barcodes were extracted from [mate, zero-based start, length], using [header mode]. The [target] was trimmed by [bases]. [Primer and filtering settings, or none] were applied. Reads were grouped by exact barcode and consensus was generated using [settings], with a maximum group size of [limit]. [Assembly settings, or no assembly]. We reviewed the completion receipt and exclusions.” Replace every bracket and describe upstream and downstream processing separately.

## 3. Technical reference

### Command help

After completing the Quick Start in the same terminal, view the help for the actual binary:

```sh
"$CF" --help
"$CF" preprocess-paired --help
"$CF" build-consensus --help
```

The following are **syntax descriptions**, not commands to copy literally. Replace capitalized placeholders with your own paths or values. Options are case-sensitive. Unknown and repeated options fail. Quote paths that contain spaces.

### Paired preprocessing

```text
clonoforge preprocess-paired FIRST_INPUT SECOND_INPUT --output NEW_DIRECTORY --header MODE --key-start N --key-len N --target first|second --trim N [OPTIONS]
```

| Option | Meaning and default |
|---|---|
| `--output` | Required new output directory. |
| `--header` | Required: `identical`, `slash`, or `sra`. |
| `--key-start`, `--key-len` | Required zero-based barcode start and barcode length in the first input. |
| `--target` | Required target input: `first` or `second`. |
| `--trim` | Required number of bases removed from the start of the target. |
| `--primer` | Optional exact anchored 5-prime primer. Cannot be combined with the IUPAC option. |
| `--primer-iupac`, `--primer-anchor` | Optional IUPAC primer and anchor `5` or `3`; anchor defaults to `5` and requires the IUPAC option. |
| `--min-length`, `--max-n`, `--min-mean-q` | Filters after trimming; defaults 1, unlimited, and 0. |
| `--max-group` | Standard group limit; default and maximum 50. |
| `--threads`, `--memory-mib`, `--max-record` | Requested workers, modeled memory and record limit; defaults 1, 1024 MiB and 16384 bytes. |
| `--assembly-min-overlap` | Enables ungapped assembly; requires target `first`. |
| `--assembly-max-mismatches`, `--assembly-max-read-length` | Assembly bounds; defaults 0 and 4096. Require the minimum-overlap option. Maximum supported read length is 4096. |
| `--experimental-large-group` | Experimental build only: 100, 250 or 500. Cannot be combined with the standard group limit. |

Neither primer mode allows mismatches, indels or an internal search. Assembly compares a suffix of processed first input with a prefix of reverse-complemented second input. Pairs with no eligible overlap or multiple eligible overlaps are excluded. Output is in the processed first-input orientation.

### Prepared FASTQ consensus

```text
clonoforge build-consensus INPUT --outdir DIR [OPTIONS]
```

Prepared headers must contain the chosen annotation, for example `@read1|BARCODE=ACGT`. Plain and gzip input are detected by content. Groups must contain equal-length uppercase ACGTN reads and Phred+33 qualities. The standard contract supports 1–50 reads per group. Unlike paired preprocessing, unsupported groups in this workflow cause failure. Check the completion manifest and the files it lists.

| Option | Meaning and default |
|---|---|
| `--outdir` | Required output directory. |
| `--outname` | Output prefix; defaults to the input file stem. |
| `--bf` | Barcode header field; default `BARCODE`. |
| `-n` | Minimum group count; default 1. |
| `--freq` | Minimum base frequency; default 0.6, range 0–1. |
| `-q` | Minimum consensus quality; default 0, range 0–93. |
| `--failed` | Also write failed groups. |
| `--force-spill` | Force the external temporary-file path. |
| `--threads`, `--memory-mib` | Requested workers and modeled memory; defaults 1 and 512 MiB. |
| `--max-record-bytes`, `--max-group-count` | Record and group bounds; defaults 16384 bytes and 50. |
| `--workdir` | Temporary storage directory; defaults to the output directory. |

At each position, consensus uses observed bases and qualities. Frequency and minimum-quality filters can produce N. Singleton behavior is separate and retains the original quality. The source package's consensus implementation and tests give the algorithm and its tested bounds.

### Sequence file commands

```text
clonoforge validate FORMAT INPUT
clonoforge inspect FORMAT INPUT
clonoforge normalize FORMAT INPUT OUTPUT
clonoforge compare FORMAT LEFT RIGHT
```

FORMAT is `fasta` or `fastq`. Validate and inspect parse records and report counts. Normalize writes plain output to a new path; it does not overwrite. Compare reports structural differences. The general sequence reader supports a broader alphabet than consensus. These checks do not validate assay design or biological meaning. Success exits with status 0; argument, parse, I/O and comparison-difference errors exit with status 2.

### Record and output details

FASTQ uses four physical lines: an `@` header, sequence, a `+` line and qualities. The plus line may repeat the exact header. Qualities use ASCII 33–126, corresponding to Phred+33 scores 0–93. Wrapped FASTQ is unsupported. Input may be plain or gzip, detected by content. The general reader accepts multiline FASTA and a broader sequence alphabet; normalization writes plain output.

Prepared headers have an identifier and optional ordered `|KEY=VALUE` annotations. Identifiers and values cannot contain spaces, `|` or `=`. Keys use uppercase ASCII letters, digits or underscores. Empty or duplicate fields are errors. This is a strict subset, not a general pRESTO header compatibility claim.

Paired completion receipts record input hashes, counts and output file hashes. The standard schema is `clonoforge-m4-preliminary-v1`; assembly or IUPAC mode uses `clonoforge-paired-v2`. An experimental large-group run uses `clonoforge-experimental-large-group-v1`. Schemas describe structure, not biological validity. Inputs are hashed before and after processing. Interrupted runs may leave hidden incomplete directories containing sequence data. Successful publication does not guarantee survival of power loss or disk failure.

### Experimental build and evidence limits

The Cargo feature is `clonoforge-cli/experimental-large`; a separately compiled binary and explicit experimental option are required for paired limits of 100, 250 or 500. These profiles retain exact-key and equal-length rules, have finite test coverage and are not a general biological validation. The paired verifier supplied with the source does not validate the experimental receipt.

In the accepted comparison, prepared reads from the pRESTO example immunoglobulin dataset (SRR4026043) and a second mixed heavy/light immunoglobulin dataset (SRR4431771) yielded byte-identical consensus outputs to the tested legacy pRESTO reference. This was a prepared-FASTQ-to-consensus comparison on an Apple ARM64 host, not a full-workflow comparison or proof of biological equivalence. Preparation and exclusions differed between the two datasets. The [manuscript](ClonoForge-Manuscript.pdf) describes methods and limitations and identifies the evaluated binary and its separate benchmark campaign.
