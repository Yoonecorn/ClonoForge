# ClonoForge

ClonoForge is a Rust tool for bounded, deterministic FASTQ preprocessing and exact-key consensus. It supports paired target preprocessing and optional ungapped overlap assembly, as well as consensus on prepared FASTQ. It implements a defined subset of pRESTO-style consensus; it is not a complete pRESTO workflow or a repertoire analysis pipeline.

**ClonoForge 1.0.0** · [Quick Start](docs/QUICK_START.md) · [User manual](docs/USER_MANUAL.md) · [Manuscript](docs/ClonoForge-Manuscript.pdf) · [Troubleshooting](docs/TROUBLESHOOTING.md)

Repository: <https://github.com/Yoonecorn/ClonoForge>

## Install and try

Build from source with Rust and Cargo. A prebuilt macOS executable is not distributed because a downloaded, quarantined installation has not been validated:

```sh
cargo build --release --locked
cargo test --workspace --locked
./target/release/clonoforge --help
```

See the [Quick Start](docs/QUICK_START.md) and [user manual](docs/USER_MANUAL.md) for a reproducible synthetic example, input requirements, options, output verification and limitations. The standard contract permits at most **50 reads per group**. Separately compiled 100, 250 and 500 profiles are opt-in experimental methods with finite qualification; they are not defaults or broad biological validation.

## Matched consensus benchmark

In an earlier benchmark campaign at matched prepared-FASTQ-to-consensus boundaries on two complete public accessions, ClonoForge produced byte-identical outputs and was 1.58× faster than legacy pRESTO 0.7.9 on SRR4026043 and 2.11× faster on SRR4431771 on the tested Apple ARM64 host.

Downloading, SRA conversion, primer and quality filtering, paired-read preparation, output validation, assembly, and downstream biological analysis were outside the timed boundary. These measurements are not whole-workflow speedups or biological equivalence claims. The [manuscript](docs/ClonoForge-Manuscript.pdf) reports a later, separate ClonoForge campaign and describes data provenance and limitations.

## Documentation and citation

- [Quick Start](docs/QUICK_START.md)
- [User manual](docs/USER_MANUAL.md)
- [Manuscript](docs/ClonoForge-Manuscript.pdf)
- [Troubleshooting](docs/TROUBLESHOOTING.md)
- [Paired receipt schemas](schemas/)
- [Citation information](CITATION.cff)

Original code is licensed under MIT OR Apache-2.0. Reference software and datasets retain their own terms. Cite the exact version and the relevant pRESTO and dataset publications.
