# Third-party provenance and licenses

Original ClonoForge code is MIT OR Apache-2.0. Cargo.lock pins Rust dependency versions; exact resolved package license metadata is retained in the final environment evidence. Dependency copyright and license terms remain with their respective authors. No pRESTO implementation source is incorporated in the Rust product.

The measured Apple Silicon acceleration enables the `sha2` crate's optional `asm` feature. Its additional locked dependencies are sha2-asm 0.6.4 (MIT), cc 1.5.1, find-msvc-tools 0.1.14 and shlex 2.0.1 (each MIT OR Apache-2.0). Runtime CPU feature detection retains the SHA software fallback on the tested aarch64 path. Other platforms are not qualified by these measurements.

The separate pRESTO 0.7.9 reference runner is AGPL-3.0 and is invoked as a black box. It is not a runtime requirement. Its separately retained distribution archives, dependencies and generated reference evidence retain their own terms. Public sequencing data retain their original provenance and use conditions; local qualification is not a grant of redistribution rights. Local Git source bundles exclude downloaded sequence inputs and reference distributions.
