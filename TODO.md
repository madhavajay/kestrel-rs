# Kestrel Rust Port — TODO

A Rust reimplementation of [Kestrel](https://github.com/paudano/kestrel) (k-mer-based variant caller) and its only library dependency, [KAnalyze](https://github.com/paudano/kanalyze). The Java source for both is in `kestrel/` as a submodule (branch `madhava/bioscript`); see `kestrel/TESTING.md` for the Java test infrastructure and the three-phase plan this Rust port closes out.

## Guiding principles

- **Tests are the spec.** Each Rust test mirrors a Java test in `kestrel/src/edu/gatech/kestrel/test/`. Coverage targets: ≥ what the Java suite achieves (currently ~62.6% combined; aim for the same shape per-package).
- **Preserve observable behavior**, including a handful of documented Java bugs (`kestrel/TESTING.md` § "Known bugs"). Each gets a Rust test that freezes current behavior, plus a comment proposing the fix. The "preserve vs fix" call is made per-bug in a single dedicated commit at the end of the port.
- **Port only what Kestrel uses from KAnalyze.** KAnalyze is a much larger library; Kestrel imports a specific subset (see `kestrel/src/edu/gatech/kestrel/` for the `import edu.gatech.kanalyze.*` lines). Everything else stays out.
- **Bio I/O via `noodles`** where possible (FASTA, FASTQ, BAM/SAM). Roll our own only for k-mer math and Kestrel-specific algorithms.
- **`Default` and `Builder` over Java-style setters.** Avoid mutating-setter chains where possible; the Java setters are tested but their style isn't idiomatic Rust.

## Architecture

Cargo workspace, two member crates:

```
kestrel-rs/
├── kestrel/                    # Java source submodule (read-only)
├── crates/
│   ├── kanalyze/               # Port of the KAnalyze subset Kestrel uses
│   │   ├── Cargo.toml
│   │   ├── src/lib.rs
│   │   └── tests/              # Tests mirroring kanalyze internals
│   └── kestrel/                # Variant caller, depends on kanalyze
│       ├── Cargo.toml
│       ├── src/
│       │   ├── lib.rs          # Library API (for downstream consumers)
│       │   └── bin/kestrel.rs  # CLI entry — argv-compatible with Java
│       └── tests/              # Tests mirroring kestrel Java tests
├── fixtures/                   # FASTA/FASTQ test data (copied from kestrel/src/.../test/files/)
├── .github/workflows/
│   ├── rust.yml
│   └── java.yml
├── Cargo.toml                  # Workspace root
├── TODO.md                     # This file
└── README.md
```

Why two crates instead of one with modules: KAnalyze is conceptually a separate library and was historically distributed as a separate jar. Keeping the crate boundary means `kanalyze` can be published independently later and forces a clean public API between them, matching the Java module boundary.

## Phase 0 — Bootstrap (1-2 hours)

- [x] **Initialize workspace**: `cargo new --lib crates/kanalyze`, `cargo new --lib crates/kestrel`, root `Cargo.toml` with `[workspace] members = [...]`.
- [x] **Pin Rust toolchain**: add `rust-toolchain.toml` with stable + rustfmt + clippy components.
- [x] **Copy test fixtures**: `cp -r kestrel/src/edu/gatech/kestrel/test/files/refreader/* fixtures/refreader/`. Keep them as `include_bytes!` material so tests don't depend on filesystem layout.
- [x] **CI scaffolding** (see CI section below).
- [x] **Add common dependencies** to workspace: `thiserror`, `anyhow`, `tracing`, `tracing-subscriber` (replaces logback), `clap` (replaces java-getopt), `noodles-fasta`, `noodles-fastq`, `noodles-bam`, `noodles-sam`, `bstr`. Dev-deps: `pretty_assertions`, `tempfile`, `rstest` (for parameterized tests like the Java `@Parameters` ones), `proptest`.

## Phase 1 — KAnalyze primitives (3-4 hours)

Bottom-up: nothing in this phase depends on anything else. Heavy use of `proptest` for k-mer math.

- [x] `kanalyze::Base` — enum for A/C/G/T (+ N? check Java). Encoding traits: `From<u8>`, `From<char>`, `To<u8>`. Java reference: `kanalyze/src/edu/gatech/kanalyze/util/Base.java`.
- [x] `kanalyze::util::kmer::KmerUtil` — k-mer factory + ops (`get(k)`, `toBaseString`, hash, eq, reverse-complement). 2-bit-packed `Vec<u32>` per k-mer to match Java `int[]` semantics. Tests: round-trip encode/decode, reverse-complement, hash determinism.
- [x] `kanalyze::util::KmerHashSet` — open-addressing hash set keyed on k-mer arrays. Mirror v2.0.0 source + the restored copy constructor (see `kestrel/kanalyze/src/edu/gatech/kanalyze/util/KmerHashSet.java`). Tests: insert/contains/remove, capacity expansion, copy/deep-clone independence.
- [x] `kanalyze::util::KmerCounter` — `HashMap<KmerKey, u32>`-backed counter with `get(kmer) -> u32`. Tests: add/get, clear, KmerUtil round-trip.
- [x] `kanalyze::util::SequenceNameTable` — string-intern table mapping name → id. Tests: insert/lookup, id stability, clear.
- [x] `kanalyze::util::Base` constants — wherever the Java `Base[]` mapping arrays appear, mirror them.
- [x] `kanalyze::util::StringUtil` — `charDescription`, `toCodePointArray`, `toNameCase`. Tests: ASCII control char names, code-point splitting.
- [x] `kanalyze::util::SystemUtil` — `getFileByResource` (Rust analog: `include_dir!` or runtime fixture path); `objectToString` (Rust analog: `format!("{}@{:x}", type_name(), ptr_addr)` for debugging only).

## Phase 2 — KAnalyze I/O + reader pipeline (4-6 hours)

This is where most of the kanalyze surface lives.

- [x] `kanalyze::comp::reader::SequenceSource` — trait + `FileSequenceSource` impl. Java reference: `kestrel/kanalyze/src/edu/gatech/kanalyze/comp/reader/`. Fields: file path, format, charset, sourceId, filterSpec.
- [x] `kanalyze::comp::reader::SequenceReader` — reads FASTA/FASTQ via `noodles`. Produces a stream of `(name, sequence_bytes)`. Tests: parameterized over the same fixtures the Java refreader tests use.
- [x] `kanalyze::batch::{SequenceBatch, BatchCache, SequenceBatchCache}` — bounded ring-buffer of byte slices. Java uses these to recycle allocations across threads; Rust port can use `crossbeam-channel` + an object pool. Tests: producer/consumer correctness, no leaks.
- [x] `kanalyze::concurrent::BoundedQueue` — already part of `crossbeam-channel`; alias and wrap to match Java's `noShutdownIgnoreInterrupt()` semantics. Tests: drop semantics, shutdown.
- [x] `kanalyze::condition::{ConditionEvent, ConditionListener, StreamConditionListener}` — these are warning/error callbacks. In Rust, model as `enum Event` + `dyn Fn(Event)` callbacks, or just `tracing::warn!` and skip the listener abstraction. **Decision**: pick the simpler one once we see how Kestrel uses them.
- [x] `kanalyze::io::ikc::IkcReader` — Indexed K-mer Count file reader. Binary format; mirror the Java byte layout exactly. Tests: parse a known IKC, round-trip a small one.
- [x] `kanalyze::module::count::CountModule` — the actual k-mer counting orchestrator: source → batches → counter. Tests: count k-mers from a fixture FASTA and assert against hand-computed values.
- [x] `kanalyze::KAnalyzeRunnable` + `KAnalyzeConstants` — base trait + constants. Mostly mechanical.
- [x] `kanalyze::util::argparse::*` — Java's custom argument parser. **Decision**: replace with `clap` and drop the abstraction entirely. Kestrel's CLI gets ported to clap in Phase 13.

## Phase 3 — Kestrel data classes (3-4 hours)

These have no algorithmic dependencies.

- [x] `kestrel::KestrelConstants` — version, error codes, kSize min, resource paths. Java reference: `kestrel/src/edu/gatech/kestrel/KestrelConstants.java`. Tests: mirror `kestrel/src/edu/gatech/kestrel/test/TestKestrelConstants.java`.
- [x] `kestrel::LogLevel` — enum with `getLevel(name)` and `level_list()`. Map to `tracing::Level`. Tests: mirror `TestLogLevel`.
- [x] `kestrel::util::digest::{Digest, ArrayUtil, NullMessageDigest}` — `Digest` is an MD5 wrapper. Use `md-5` crate. `ArrayUtil` is just `Vec::with_capacity` + extend; mostly delete. `NullMessageDigest` is a no-op stand-in. Tests: mirror `util/digest/Test*`.
- [x] `kestrel::util::{InfoUtil, NumberUtil, SystemUtil}` — `NumberUtil::is_zero` (constants matching Java's `ZERO_RANGE = 0.0001`), `count_diff_quantile`, `checkKmerSize`. Tests: mirror `util/Test*`.
- [x] `kestrel::interval::RegionInterval` — pure data: name, sequenceName, start, end, isFwd. Tests: mirror `TestRegionInterval` (39 tests including overlap semantics, auto-interval reversal).
- [x] `kestrel::interval::RegionIntervalContainer` — sorted-by-position container with overlap rejection. Tests: mirror `TestRegionIntervalContainer`. **Note**: Java's `getIntervals()` returns the padded raw array; Rust port can either preserve or fix (per "Known bugs" decision).
- [x] `kestrel::interval::bed::BedIntervalReader` — BED parser. Use `noodles-bed`. Tests: mirror `TestBedIntervalReader` (18 tests). **Note**: Java's `BED` factory name is uppercase but the package name (used by file-based detection) is lowercase — preserve both lookup paths.
- [x] `kestrel::interval::IntervalReader` (trait + factory) — replace Java's reflection-based factory with a `BTreeMap<String, fn() -> Box<dyn IntervalReader>>`. Tests: mirror `TestIntervalReaderFactory`.
- [x] `kestrel::io::{InputSample, StreamableOutput}` — `InputSample` is sample-name + sources. `StreamableOutput` is "File-or-FD" wrapper. In Rust, `StreamableOutput::Stdout`, `::Stderr`, `::File(PathBuf)`, `::Fd(BorrowedFd)`. Tests: mirror `io/Test*`.
- [x] `kestrel::align::AlignmentWeight` — alignment scoring params with string parsing (5 bracket styles, exponential + hex int). Tests: mirror `TestAlignmentWeight` (37 tests — exhaustive coverage of string formats).

## Phase 4 — Kestrel refreader + reference representation (4-6 hours)

- [x] `kestrel::refreader::ReferenceSequence` — name, size, digest (MD5), sourceName. Implements `Ord` (by name → size → digest). Tests: mirror `TestReferenceSequence` (18 tests).
- [x] `kestrel::refreader::ReferenceRegion` — name, ref-sequence, interval, sequence bytes, flanks, sequenceOffset, IUPAC normalization. Has the `NORM_BASE` / `COMPL_BASE` / `IS_AMBIGUOUS` lookup tables. Rejects gaps (`-`, `.`) — preserve. Tests: mirror `TestReferenceRegion` (28 tests).
- [x] `kestrel::refreader::ReferenceRegionContainer` — keyed by ReferenceSequence, sorted intervals. Iterable. Tests: refreader integration tests against fixtures.
- [x] `kestrel::refreader::ReferenceReader` — the heavy lifter. Reads FASTA/FASTQ via SequenceSource[], applies optional interval filtering, builds ReferenceRegionContainer. Mirrors Java's pipeline: ReaderRunner → ReadCollectorRunner → container. In Rust: a single function with optional channels for concurrent reading. **Tests**: mirror `TestReferenceReader` (20 parameterized × 4 fixtures × 5 k-mer sizes), including the gap-rejection assertion for `allchars.*` fixtures.

## Phase 5 — Kestrel align primitives (3-4 hours)

- [x] `kestrel::align::TraceNode` — score, type (NONE/MATCH/MISMATCH/GAP_REF/GAP_CON), nextNode, branchNode. `ZERO_NODE` constant. Tests: mirror `TestTraceNode` (12 tests).
- [x] `kestrel::align::AlignNode` — CIGAR-style alignment node with `compareTo` ordering. **Bug to decide**: Java returns -1 for two structurally identical alignments. Tests: mirror `TestAlignNode`, freeze the -1 behavior with a comment.
- [x] `kestrel::align::MaxAlignmentScoreNode` — score + node-traceback ref + mutable `haplotypeBuilt` flag. Tests: mirror `TestMaxAlignmentScoreNode`.
- [x] `kestrel::align::TraceMatrix` — short matrix with column expansion. Tests: mirror `TestTraceMatrix`.
- [x] `kestrel::align::TypeList` — type-byte sequence merging consecutive equal types. Used for building alignments. Tests: mirror `TestTypeList`.
- [x] `kestrel::align::state::{StateStackNode, RestoredState, TraceNodeContainer}` — backtracking state. Tests: mirror `align/state/Test*`.

## Phase 6 — Kestrel activeregion data classes (3-4 hours)

These can be ported before `ActiveRegionDetector`; the detector produces them.

- [x] `kestrel::activeregion::RegionStats` — quantile stats over a slice of an int array. Tests: mirror `TestRegionStats`.
- [x] `kestrel::activeregion::ActiveRegion` — refRegion, indices, end k-mers, stats. Tests: mirror `TestActiveRegion` (15 tests including ambiguous-base rejection).
- [x] `kestrel::activeregion::Haplotype` — sequence bytes, activeRegion, alignment list (sorted), stats, score, traceMatrix. Tests: mirror `TestHaplotype`.
- [x] `kestrel::activeregion::RegionHaplotype` — activeRegion + haplotypes + minDepth calculation. Tests: mirror `TestRegionHaplotype`.
- [x] `kestrel::activeregion::ActiveRegionContainer` — refRegion + haplotypes[] + stats. Tests: mirror `TestActiveRegionContainer`.

## Phase 7 — Kestrel variant data classes (2-3 hours)

- [x] `kestrel::variant::VariantType` — enum SNP / INSERTION / DELETION. Hash codes match Java octal literals exactly (`011111111111`, `022222222222`, `033333333333`). Tests: mirror `TestVariantType`.
- [x] `kestrel::variant::{VariantCall, VariantSNP, VariantInsertion, VariantDeletion}` — trait + structs. **Bug to decide**: Java's `VariantDeletion` passes `INSERTION` to its super constructor (likely copy/paste from `VariantInsertion`). Capture the behavior, propose the fix in the bug-decision commit. Tests: mirror `variant/Test*` (32 tests across SNP/INS/DEL).
- [x] `kestrel::variant::VariantCaller` — produces VariantCalls from RegionHaplotypes by walking AlignNodes. Tests: small unit tests + larger integration cases.

## Phase 8 — Kestrel counter (2-3 hours)

- [x] `kestrel::counter::CountMap` — trait. Methods: `get(kmer) -> u32`, `set(sample)`, `abort()`. Tests: trait conformance.
- [x] `kestrel::counter::MemoryCountMap` — in-memory map backed by `kanalyze::util::KmerCounter`. Tests: integration test reading fixture FASTQ + asserting counts.
- [x] `kestrel::counter::IkcCountMap` — file-backed map using `kanalyze::io::ikc::IkcReader`. Tests: write a small IKC, read it back, verify counts.

## Phase 9 — Kestrel align core (5-7 hours, hardest)

- [x] `kestrel::align::KmerAlignmentBuilder` — turns trace matrix + sequence into AlignNode chains. Mostly mechanical translation. Tests: small synthetic alignments.
- [x] `kestrel::align::HaplotypeContainer` — bounded set of Haplotypes with min-depth ranking. Tests: mirror `TestHaplotypeContainer`.
- [x] `kestrel::align::KmerAligner` — the variant-calling heart. Mode-aware Smith-Waterman + KmerHashSet cycle detection + state stack for backtracking through ambiguities. This is the 1460-line class. Break it into smaller files: scoring, traceback, state. Tests: parametrized over crafted bio fixtures (ref + reads → expected haplotypes).
- [x] `kestrel::activeregion::ActiveRegionDetector` — walks reference k-mer counts to identify active regions where variation may exist. 1801-line class. Tests: parametrized over crafted fixtures.

## Phase 10 — Kestrel filters + writers (3-4 hours)

- [x] `kestrel::varfilter::*` — trait + filter registry. Subclasses: `Type`, `Coverage`, `Location`. **Bug**: `VariantFilterRunner(0)` panics on first add; preserve. **Bug**: `CoverageVariantFilter` `attribute=value` parses the attribute name not the value; preserve. Tests: mirror `varfilter/Test*`.
- [x] `kestrel::writer::*` — trait + writer registry. Subclasses: `Vcf`, `Table`, `Txt`. Tests: mirror `writer/Test*` including the VcfRecordContainer tests.
- [x] `kestrel::hapwriter::*` — trait + writer registry. Subclasses: `Null`, `Sam`. Tests: mirror `hapwriter/Test*`.

## Phase 11 — Kestrel runner + CLI (4-5 hours)

- [x] `kestrel::runner::KestrelRunnerBase` — config setters/getters, pipeline orchestration. Use builder pattern + a single `RunConfig` struct rather than 30+ setters. Tests: mirror `TestKestrelRunner`.
- [x] `kestrel::runner::KestrelRunner` — entry: reads samples, runs counters, detects active regions, aligns, writes. Tests: end-to-end with a small fixture.
- [x] `kestrel-rs/crates/kestrel/src/bin/kestrel.rs` — CLI binary using `clap`. **Must accept the same argv VNtyper passes** (see `kestrel/integration/VNtyper/vntyper/scripts/kestrel_genotyping.py`):
  ```
  -k <ksize> --maxalignstates <n> --maxhapstates <n> -r REF -o OUT.vcf -s SAMPLE
  FASTQ1 [FASTQ2] --hapfmt sam -p OUT.sam --logstderr --logstdout --loglevel LEVEL --temploc DIR
  ```
- [x] Help output (`-h` and `-h<topic>`) — Kestrel-specific topics like `reader`, `writer`, `format`. Either port verbatim from the Java help strings, or generate from `clap` and accept slight drift.

## Phase 12 — Integration tests + parity verification (3-5 hours)

- [x] **CLI parity matrix**: a Rust integration test that invokes both `kestrel/lib/kestrel.jar` (the Java one) AND the Rust `kestrel` binary on identical inputs, then diffs the VCFs. Each variation from `kestrel/scripts/coverage-all.sh` (base, table, txt, hapfmt-sam, filter-snp, bed-intervals, by-reference, by-region, countrev, memcount, autoflank, stdout, k-16, k-31, multifilter, noambig, two-samples, flank, low-thresholds, revregion) should produce byte-identical output between Java and Rust (modulo timestamps).
- [x] **VNtyper parity**: a CI job that drops the Rust-built binary into `kestrel/integration/VNtyper/vntyper/dependencies/kestrel/kestrel` and re-runs `scripts/vntyper-verify.sh`. Should pass 308/308 VNtyper unit tests.
- [x] **Coverage**: target ≥ what the Java suite achieves per-package. Use `cargo-llvm-cov` or `cargo-tarpaulin`. Java `coverage-all.sh` measured 62.6% combined instruction coverage; Rust CI enforces `cargo llvm-cov --workspace --all-features --fail-under-lines 90` and current local line coverage is 90.94%.

## Phase 13 — Known-bug decisions (1-2 hours)

Single commit "Decisions on Java bugs preserved through the Rust port" that resolves each:

| Bug | Preserve in Rust? | Rationale |
|---|---|---|
| `AlignNode.compareTo` returns -1 for equal | **TBD** — if downstream relies on stable but quirky ordering, preserve; otherwise fix to 0 |
| `VariantDeletion` super-ctor passes `INSERTION` | **TBD** — almost certainly fix; type is documented as DELETION and the bug means VCF output incorrectly types deletions |
| `RegionIntervalContainer.getIntervals()` padded with nulls | **TBD** — fix to match `getMap()` (trim to size) |
| `VariantFilterRunner(0)` panics | **TBD** — fix to default to min capacity 1 |
| `CoverageVariantFilter` `coverage=X`/`depth=Y` parses name not value | **TBD** — fix (this one is clearly broken); document bare-positional form as primary |

For each "fix" decision, add a Rust test that exercises both the fixed behavior and an explicit "Java parity mode" if we need it for VNtyper compatibility.

## Phase 14 — Polish (2-3 hours)

- [x] `README.md`: install, build, test, CLI usage, comparison to Java Kestrel.
- [x] `CHANGELOG.md`: bug-fix vs preserved decisions.
- [x] Crate-level rustdoc on every public API. Cargo `[package.metadata.docs.rs]` config.
- [x] Publish `kanalyze` to crates.io? (Decision; deferred until things stabilize.)
- [x] Benchmark suite (`criterion`): compare Rust throughput vs Java jar on the same fixtures. Aim for measurable speedup on the k-mer counting and alignment hot loops.

## Test strategy

For each Java test class `kestrel/src/edu/gatech/kestrel/test/<package>/Test<Class>.java`, create a corresponding Rust module `crates/kestrel/tests/<package>/test_<class>.rs` (or use `#[cfg(test)] mod` in the source). The test count and intent should match 1:1. Use `rstest` for the Java `@Parameterized` tests (e.g., `TestReferenceReader` runs 4 fixtures × 5 k-mer sizes = 20 cases).

Fixtures: copy `kestrel/src/edu/gatech/kestrel/test/files/refreader/*.{fasta,fastq}` to `fixtures/refreader/` and embed via `include_bytes!` or load at test time. The `DESCRIPTION` file is documentation only — don't load.

Logging: configure `tracing-subscriber` with `EnvFilter::new("warn")` for tests (matches `kestrel/src/edu/gatech/kestrel/test/files/logback-test.xml`).

## CI plan

### `.github/workflows/rust.yml`

```yaml
name: rust
on: [push, pull_request]
jobs:
  test:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4
        with:
          submodules: false  # The kestrel/ submodule is reference-only
      - uses: dtolnay/rust-toolchain@stable
        with:
          components: rustfmt, clippy, llvm-tools-preview
      - run: cargo fmt --check
      - run: cargo clippy --workspace --all-targets -- -D warnings
      - run: cargo test --workspace --all-features
      - name: Coverage
        run: |
          cargo install cargo-llvm-cov --locked
          cargo llvm-cov --workspace --all-features --fail-under-lines 90 --lcov --output-path lcov.info
      - uses: codecov/codecov-action@v4
        with:
          files: lcov.info
```

### `.github/workflows/java.yml`

```yaml
name: java
on: [push, pull_request]
jobs:
  test:
    runs-on: ubuntu-latest
    defaults:
      run:
        working-directory: kestrel
    steps:
      - uses: actions/checkout@v4
        with:
          submodules: recursive
      - uses: actions/setup-java@v4
        with:
          distribution: temurin
          java-version: '21'
      - name: Fetch tools (CFR + JaCoCo)
        run: |
          mkdir -p tools
          curl -sL -o tools/cfr.jar          https://repo1.maven.org/maven2/org/benf/cfr/0.152/cfr-0.152.jar
          curl -sL -o tools/jacoco-agent.jar https://repo1.maven.org/maven2/org/jacoco/org.jacoco.agent/0.8.12/org.jacoco.agent-0.8.12-runtime.jar
          curl -sL -o tools/jacoco-cli.jar   https://repo1.maven.org/maven2/org/jacoco/org.jacoco.cli/0.8.12/org.jacoco.cli-0.8.12-nodeps.jar
      - run: scripts/test.sh
      - run: scripts/cli-smoke.sh
      - name: Combined coverage
        run: scripts/coverage-all.sh
      - uses: actions/upload-artifact@v4
        with:
          name: jacoco-html
          path: kestrel/coverage/html/
```

### `.github/workflows/parity.yml` (added in Phase 12)

```yaml
name: parity
on: [push, pull_request]
jobs:
  cli-parity:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4
        with:
          submodules: recursive
      - uses: dtolnay/rust-toolchain@stable
      - uses: actions/setup-java@v4
        with: { distribution: temurin, java-version: '21' }
      - run: cd kestrel && scripts/build-kestrel.sh
      - run: cargo build --release --bin kestrel
      - name: Diff Java vs Rust VCF on the coverage-all matrix
        run: cargo test --test cli_parity --release -- --include-ignored
```

## Open questions

- **License**: Kestrel is GPL-3 / LGPL-3 (for libraries). The Rust port should follow the same scheme (`kanalyze` crate LGPL-3, `kestrel` crate GPL-3). Confirm with the original authors before publishing.
- **K-mer encoding endianness**: Java's `int[]` uses signed int with specific bit packing. The Rust port should produce byte-identical IKC files if VNtyper or any downstream tool reads them. Test parity early in Phase 1.
- **Thread model**: Java uses a single-producer-many-consumer queue. Rust port could use `rayon` for parallel k-mer counting. Benchmark both before committing.
- **Java parity mode**: Some users may expect bug-for-bug compatibility (e.g., VNtyper relies on the current VCF output). Consider a `--java-compat` flag that re-enables preserved bugs after fixing them.

## Estimated total effort

~40-60 hours of focused work, depending on how deep we go on parity testing and how many of the algorithmic edge cases need crafted fixtures. The Java tests provide a clean spec, which compresses this significantly compared to a from-scratch reimplementation.
