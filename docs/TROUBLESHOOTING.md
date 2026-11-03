# Troubleshooting

Stop when a command fails. Save the full error, the program version, the settings and any completion receipts. Do not share patient information, sample identifiers or sequence data in a public issue. See [Quick Start](QUICK_START.md) and [User manual](USER_MANUAL.md).

| Symptom | Safe next action |
|---|---|
| A tool version check fails | Install the missing Rust toolchain or Python 3. Open a new terminal and repeat all three checks in Quick Start step 1. If Cargo cannot fetch dependencies, check your network and local policies. |
| Source files are not found | Find the unpacked source folder and change into it in Terminal. Repeat the source-file check in Quick Start step 2. Do not run the tutorial from inside a random data folder. |
| The executable cannot be found | Return to the same terminal used for the build. If it was closed, start the tutorial again from step 2. Check the build completed without an error. |
| Pair identity or duplicate error | Stop. Compare mate order, pair IDs and the selected header mode against the original files. Obtain correctly synchronized inputs; do not guess at header edits. Retry with a new output directory. |
| FASTQ parse error | Check that each record has four physical lines, a valid header and plus line, and matching sequence and quality lengths. Check that the file is complete and any gzip input is intact. Repair or reacquire the input using your laboratory's validated procedure. |
| Output already exists | Pick a new output directory. Keep the earlier output for review. The paired workflow will not overwrite it. |
| Few or no expected consensus reads | Open `completion.json`. Review invalid-key, filter, oversize and unequal-length counts. Check barcode location, target choice, primer and filtering settings against the assay protocol. Do not interpret excluded counts as biological recovery. |
| Prepared consensus rejects a group | Check that each header has the selected barcode field, target bases are uppercase ACGTN, and each group has equal-length reads and at most 50 reads for the standard build. Correct preparation only after checking the assay. |
| Memory or disk error | Check free memory and temporary disk. Request fewer workers and retry with a new output path. The memory setting is a model, not a hard memory limit. |
| Verification fails or a receipt is missing | Preserve the result directory and the error. Do not use that output for downstream analysis until the discrepancy is understood. The supplied paired verifier does not support experimental large-group receipts. |
| A run is interrupted | A hidden incomplete directory may contain sequences. Preserve it if needed for diagnosis. After addressing the cause, retry with a new output path. |

## When to report a problem

Stop and contact the maintainer if the synthetic Quick Start fails after the checks above, a verifier reports a hash or count mismatch, or a reproducible result contradicts the documented contract. Open an issue at the [ClonoForge issue tracker](https://github.com/Yoonecorn/ClonoForge/issues).

Provide the version, operating system, exact command with private paths removed, complete error text, and steps to reproduce using non-sensitive synthetic data. Include receipt details only after removing confidential information. Never upload raw reads or identifiable sample data to a public issue.
