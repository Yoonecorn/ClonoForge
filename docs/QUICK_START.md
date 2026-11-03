# ClonoForge Quick Start

This tutorial uses three artificial read pairs. It teaches the workflow; its barcode and read layout are not recommendations for your experiment. You need the ClonoForge 1.0.0 source package, Rust, Cargo, Python 3, and a macOS or Unix-like terminal. The qualified system was macOS on Apple ARM64. Validate other systems locally. A prebuilt executable is not distributed because a downloaded, quarantined macOS installation has not been validated.

Copy each block in order into the same Terminal window. Press Return after pasting. **Stop if a command reports an error.** Do not continue until you have fixed it. See [Troubleshooting](TROUBLESHOOTING.md).

## 1. Check the tools

```sh
rustc --version
cargo --version
python3 --version
```

Each line should print a version. If a tool is missing, install it before continuing. Cargo is normally installed with Rust. Python 3 runs the independent output verifier.

## 2. Open the source folder

In Terminal, change into your unpacked ClonoForge source folder. Ask your local support person if you do not know its location. It must contain `Cargo.toml` and `examples/verify_paired.py`. Check that both files are present:

```sh
test -f Cargo.toml && test -f examples/verify_paired.py && printf 'Source files found\n'
```

Continue only if you see `Source files found`. Stay in this Terminal window for the remaining steps.

## 3. Make a fresh workspace and build

```sh
SOURCE="$PWD"
WORK="$(mktemp -d "${TMPDIR:-/tmp}/clonoforge-tutorial.XXXXXX")"
cd "$WORK"
printf 'Tutorial directory: %s\n' "$WORK"
CARGO_TARGET_DIR="$WORK/build" cargo build --release --locked --manifest-path "$SOURCE/Cargo.toml"
```

The build may need internet access to obtain the locked dependencies. It can take several minutes. The tutorial directory contains the build and the example files; the source folder is not changed. Check the program:

```sh
CF="$WORK/build/release/clonoforge"
"$CF" --version
```

The version should be `clonoforge 1.0.0`.

## 4. Create the example reads

```sh
cat > R1.fastq <<'FASTQ'
@x
TT
+
II
@y
AA
+
II
@z
TT
+
II
FASTQ
cat > R2.fastq <<'FASTQ'
@x
AC
+
II
@y
TG
+
II
@z
AC
+
II
FASTQ
```

The first file contains two-base barcodes. The second contains target sequences. Two reads have barcode `TT`; one has `AA`.

## 5. Process the reads

```sh
"$CF" preprocess-paired R1.fastq R2.fastq --output result --header identical --key-start 0 --key-len 2 --target second --trim 0 --threads 1 --memory-mib 2048
```

This extracts the first two bases of the first input as the barcode and uses the second input as the target. It creates a new `result` directory.

## 6. Verify and find the results

```sh
python3 "$SOURCE/examples/verify_paired.py" result
printf 'Tutorial directory: %s\n' "$WORK"
```

The verifier should report `"status": "PASS"`, `"qualified": 3` and `"qualified_groups": 2`. It checks file hashes, sizes, FASTQ framing and count relationships. It does not establish biological accuracy. In `result`, find `completion.json`, `qualified.fastq`, `m3_completion.json`, and `m3_consensus-pass.fastq`. The last line gives the tutorial directory. A fresh run makes a different directory.

For your own samples, first read [Choose settings in the user manual](USER_MANUAL.md). Check your assay layout and keep the input files unchanged.
