# kestrel-rs

Rust port of Kestrel, a k-mer-based variant caller, plus the subset of KAnalyze that Kestrel needs. The Java source is kept as the `kestrel/` submodule on branch `madhava/bioscript` and is the reference implementation for behavior and tests.

## Status

This port is under active development. The lower-level data structures, readers, counters, filters, writers, active-region detection, variant data types, and alignment primitives are implemented with Rust tests mirroring the Java behavior.

The Rust CLI parses Kestrel-style arguments and runs the VNtyper-facing pipeline: references and samples are read, k-mers are counted, active regions are detected, read-backed haplotypes are derived for spanning reads, variants are called and filtered, and VCF/SAM output is written. Java graph traversal is still the main remaining area where deeper parity work may be needed.

## Build

```sh
cargo build --workspace
cargo build --release -p kestrel --bin kestrel
```

The workspace pins the Rust stable toolchain via `rust-toolchain.toml`.

## Test

Rust unit, integration, lint, and docs checks:

```sh
cargo fmt --check
cargo clippy --workspace --all-targets -- -D warnings
cargo test --workspace --all-features
RUSTDOCFLAGS="-D missing-docs" cargo doc -p kestrel --no-deps
RUSTDOCFLAGS="-D missing-docs" cargo doc -p kanalyze --no-deps
```

Java reference tests from the fork:

```sh
git submodule update --init --recursive kestrel
cd kestrel
mkdir -p tools
curl -sL -o tools/cfr.jar https://repo1.maven.org/maven2/org/benf/cfr/0.152/cfr-0.152.jar
curl -sL -o tools/jacoco-agent.jar https://repo1.maven.org/maven2/org/jacoco/org.jacoco.agent/0.8.12/org.jacoco.agent-0.8.12-runtime.jar
curl -sL -o tools/jacoco-cli.jar https://repo1.maven.org/maven2/org/jacoco/org.jacoco.cli/0.8.12/org.jacoco.cli-0.8.12-nodeps.jar
scripts/test.sh
scripts/cli-smoke.sh
scripts/coverage-all.sh
```

Rust-vs-Java CLI parity:

```sh
kestrel/scripts/build-kestrel.sh
KESTREL_RUN_JAVA_PARITY=1 KESTREL_JAR=kestrel/lib/kestrel.jar \
  cargo test -p kestrel --test cli_parity -- --nocapture
```

VNtyper shim verification:

```sh
cargo build --release -p kestrel --bin kestrel
scripts/vntyper-verify.sh target/release/kestrel
```

Benchmarks use Criterion:

```sh
cargo bench -p kestrel --bench pipeline
```

Set `KESTREL_JAR=/path/to/kestrel.jar` to include the optional Java CLI comparison benchmark on the same generated fixture.

## CLI

The Rust binary is `kestrel`:

```sh
cargo run -p kestrel -- -h
cargo run -p kestrel -- -hreader
```

Accepted VNtyper-style argument shape:

```sh
cargo run -p kestrel -- \
  -k 25 --maxalignstates 40 --maxhapstates 40 \
  -r ref.fasta -o out.vcf -s sample \
  reads_1.fastq reads_2.fastq \
  --hapfmt sam -p haplotypes.sam \
  --logstderr --loglevel ERROR --temploc /tmp
```

The command line is implemented with `clap`; Java's custom `argparse` package is intentionally not ported.

## VNtyper Shim

VNtyper currently shells out through the Java Kestrel command shape:

```sh
java -jar dependencies/kestrel/kestrel.jar ...
```

`scripts/vntyper-verify.sh` installs the Rust release binary beside VNtyper's expected Kestrel dependency and builds a small `kestrel.jar` shim there. The shim preserves the Java invocation form but delegates execution to the sibling Rust `kestrel` binary, so VNtyper's existing tests exercise the Rust package through the same call path.

The verifier creates `.venv-vntyper`, installs VNtyper's Python test dependencies, runs the VNtyper unit test suite, and checks that `java -jar .../kestrel.jar -h` reaches the Rust CLI.

Large VNtyper integration data is not stored in this repository. The current CI job runs the VNtyper unit suite and the Rust shim probe.

## Java Kestrel Compatibility

The target is observable compatibility with Java Kestrel for VNtyper's usage. Known Java quirks are preserved where tests cover them, including:

- `AlignNode.compareTo` equal-alignment ordering behavior.
- `VariantDeletion` retaining the Java constructor/type bug.
- `VariantFilterRunner(0)` panicking on first add.
- `CoverageVariantFilter` preserving Java's `attribute=value` parse bug.

CI is split into four lanes:

- `java-kestrel`: clones the `madhavajay/kestrel` fork via the submodule and runs the Java JUnit, CLI smoke, and coverage scripts.
- `rust-unit`: runs Rust formatting, linting, unit/integration tests, rustdoc checks, and coverage.
- `rust-java-parity`: builds the Java reference jar and runs the opt-in Rust CLI parity tests.
- `vntyper-shim`: builds the Rust release binary and runs VNtyper's unit tests through the Java shim.

Coverage is enforced in Rust CI with `cargo-llvm-cov --fail-under-lines 90`.

See `TODO.md` for the current checklist.

## Publishing

The `kanalyze` crate is not ready for crates.io publication. It remains a workspace crate until runner parity, public rustdoc coverage, and compatibility decisions are complete. See `docs/kanalyze-publishing.md`.
