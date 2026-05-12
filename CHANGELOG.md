# Changelog

## Unreleased

### Added

- Rust workspace with `kanalyze` and `kestrel` crates.
- KAnalyze subset used by Kestrel: bases, k-mer utilities, counters, sequence readers, IKC I/O, count module, conditions, batching, and queue wrappers.
- Kestrel constants, logging levels, digest utilities, interval parsing, streamable output, reference readers, active-region data structures, variant data structures, filters, writers, haplotype writers, counters, alignment primitives, active-region detection, and a limited variant-emitting runner.
- Clap-based Rust CLI that accepts the VNtyper-style Kestrel argv shape and supports `-h`, `-hreader`, `--help=topics`, `reader`, `writer`, and `format` help topics.
- Criterion benchmark suite for Rust k-mer counting, the current Rust runner path, and optional Java Kestrel CLI comparison via `KESTREL_JAR`.
- Rust CI coverage enforcement with `cargo-llvm-cov` at 90% minimum line coverage, above the measured Java combined coverage baseline.

### Preserved Java Behavior

- `AlignNode.compareTo` equal-alignments quirk is tested and preserved.
- `VariantDeletion` currently preserves the Java type-constructor bug for compatibility while the final decision remains open.
- `VariantFilterRunner(0)` panics on first add, matching Java behavior.
- `CoverageVariantFilter` preserves the Java `coverage=X` / `depth=Y` parse bug; bare positional coverage/depth remains the primary working form.

### Incomplete

- Java-equivalent graph traversal for complex haplotype assembly.
- Java CLI parity matrix.
- VNtyper end-to-end parity verification.
- Publishing `kanalyze` to crates.io is deferred until the public API stabilizes.

### Licensing

- Added `COPYING`, `COPYING.LESSER`, and `COPYING.DOC` from the Java source tree.
