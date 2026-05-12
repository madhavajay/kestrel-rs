# kanalyze Publishing Decision

Decision: do not publish the `kanalyze` crate to crates.io yet.

Reasoning:

- The crate currently exposes only the subset of KAnalyze needed by the Rust Kestrel port, not the full Java KAnalyze API.
- Several public APIs are still shaped by porting needs and may change while full Kestrel runner and parity work is completed.
- Publishing before parity would make semver promises around APIs that are still being validated against Java behavior.

Revisit publishing after:

- The Rust Kestrel runner emits variants end to end.
- Java CLI parity and VNtyper parity have been verified or explicitly scoped.
- Public rustdoc coverage is completed for the exposed `kanalyze` API.
- Any intentionally preserved Java quirks are documented in the crate docs.

Until then, `kanalyze` remains a workspace crate consumed by `kestrel`.
