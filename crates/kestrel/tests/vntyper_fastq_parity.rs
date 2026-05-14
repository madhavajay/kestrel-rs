use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

use kanalyze::comp::reader::FileSequenceSource;
use kestrel::io::{InputSample, StreamableOutput};
use kestrel::runner::KestrelRunner;
use tempfile::TempDir;

#[test]
fn vntyper_positive_fastq_matches_java_expected_vcf() {
    run_case_if_enabled(VntyperCase::new("positive", "example_6449_hg19_subset"));
}

#[test]
fn vntyper_negative_fastq_matches_java_expected_vcf() {
    run_case_if_enabled(VntyperCase::new("negative", "example_66bf_hg19_subset"));
}

fn run_case_if_enabled(case: VntyperCase) {
    let Some(root) = enabled_bioscript_root() else {
        eprintln!("set KESTREL_RUN_VNTYPER_FASTQ_PARITY=1 to run VNtyper FASTQ parity");
        return;
    };
    let reference = root
        .join("ports/vntyper/vntyper/reference")
        .join("All_Pairwise_and_Self_Merged_MUC1_motifs_filtered.fa");
    let temp = TempDir::new().unwrap();
    let work_dir = parity_output_dir(case.label).unwrap_or_else(|| temp.path().to_path_buf());
    fs::create_dir_all(&work_dir)
        .unwrap_or_else(|err| panic!("failed to create {}: {err}", work_dir.display()));
    let actual_vcf = work_dir.join(format!("{}.rust.vcf", case.label));
    let fastq_1 = decompress_gzip(
        &root
            .join("ports/vntyper/test-data")
            .join(format!("{}_R1.fastq.gz", case.sample)),
        &work_dir,
    );
    let fastq_2 = decompress_gzip(
        &root
            .join("ports/vntyper/test-data")
            .join(format!("{}_R2.fastq.gz", case.sample)),
        &work_dir,
    );

    run_kestrel(
        &reference,
        &[fastq_1, fastq_2],
        case.sample,
        &actual_vcf,
        &work_dir,
    );

    let expected_vcf = root
        .join("ports/vntyper/test-data/expected")
        .join(case.label)
        .join("kestrel/output.vcf");
    if let Some(output_dir) = std::env::var_os("KESTREL_VNTYPER_PARITY_OUT") {
        let copied_expected = PathBuf::from(output_dir)
            .join(case.label)
            .join(format!("{}.java-expected.vcf", case.label));
        fs::copy(&expected_vcf, &copied_expected).unwrap_or_else(|err| {
            panic!(
                "failed to copy expected VCF to {}: {err}",
                copied_expected.display()
            )
        });
        eprintln!(
            "kept VNtyper Kestrel parity files in {}",
            work_dir.display()
        );
    }

    assert_vcf_records_match(&expected_vcf, &actual_vcf, case.label);
}

struct VntyperCase {
    label: &'static str,
    sample: &'static str,
}

impl VntyperCase {
    const fn new(label: &'static str, sample: &'static str) -> Self {
        Self { label, sample }
    }
}

fn run_kestrel(reference: &Path, fastqs: &[PathBuf], sample: &str, output: &Path, temp: &Path) {
    let mut runner = KestrelRunner::new();
    runner.set_k_size(20).unwrap();
    runner.set_output_file(Some(StreamableOutput::from_path(output, None)));
    runner.set_output_format("vcf").unwrap();
    runner.set_log_file(Some(StreamableOutput::stderr()));
    runner.set_temp_dir_name(Some(&temp.display().to_string()));
    runner.set_kmer_count_in_memory(true);
    runner.set_count_reverse_kmers(true);
    runner.set_minimum_difference(5).unwrap();
    runner.set_difference_quantile(0.90).unwrap();
    runner.set_anchor_both_ends(true);
    runner.set_decay_minimum(0.55).unwrap();
    runner.set_decay_alpha(0.80).unwrap();
    runner.set_peak_scan_length(7).unwrap();
    runner.set_scan_limit_factor(7.0).unwrap();
    runner.set_call_ambiguous_regions(true);
    // Use Java's default `DEFAULT_MIN_KMER_COUNT = 5`, which translates to
    // `kmercount:5` post-count filter in kanalyze. The Java CLI used to
    // generate the expected fixture passes no --mincount, so 5 is the
    // effective threshold.
    runner.set_min_kmer_count(5).unwrap();
    runner
        .set_max_haplotypes(parity_usize_env("KESTREL_VNTYPER_MAX_HAPLOTYPES", 2) as i32)
        .unwrap();
    runner.set_max_repeat_count(0).unwrap();
    runner
        .set_max_aligner_state(parity_usize_env("KESTREL_VNTYPER_MAX_ALIGNER_STATES", 2) as i32)
        .unwrap();
    runner.add_reference(FileSequenceSource::from_path(reference, 1).unwrap());
    let sources = fastqs
        .iter()
        .enumerate()
        .map(|(index, path)| FileSequenceSource::from_path(path, i32::try_from(index + 1).unwrap()))
        .collect::<Result<Vec<_>, _>>()
        .unwrap();
    runner.add_sample(InputSample::new(Some(sample), sources).unwrap());
    runner.run().unwrap();
}

fn assert_vcf_records_match(expected: &Path, actual: &Path, label: &str) {
    let expected_records = vcf_records(expected);
    let actual_records = vcf_records(actual);
    let context = vcf_record_diff_context(&expected_records, &actual_records);
    assert_eq!(
        actual_records.len(),
        expected_records.len(),
        "{label} VCF record count differs: {context}"
    );
    for (index, (actual, expected)) in actual_records
        .iter()
        .zip(expected_records.iter())
        .enumerate()
    {
        assert_eq!(
            actual, expected,
            "{label} VCF record {index} differs: {context}"
        );
    }
}

fn vcf_record_diff_context(expected_records: &[String], actual_records: &[String]) -> String {
    let expected = expected_records
        .iter()
        .map(|record| ParsedVcfRecord::parse(record))
        .collect::<Vec<_>>();
    let actual = actual_records
        .iter()
        .map(|record| ParsedVcfRecord::parse(record))
        .collect::<Vec<_>>();
    let expected_keys = expected
        .iter()
        .map(|record| &record.key)
        .collect::<BTreeSet<_>>();
    let actual_keys = actual
        .iter()
        .map(|record| &record.key)
        .collect::<BTreeSet<_>>();
    let expected_by_key = expected
        .iter()
        .map(|record| (&record.key, record))
        .collect::<BTreeMap<_, _>>();
    let actual_by_key = actual
        .iter()
        .map(|record| (&record.key, record))
        .collect::<BTreeMap<_, _>>();
    let missing = expected
        .iter()
        .filter(|record| !actual_keys.contains(&record.key))
        .take(5)
        .map(|record| record.line.clone())
        .collect::<Vec<_>>();
    let extra = actual
        .iter()
        .filter(|record| !expected_keys.contains(&record.key))
        .take(5)
        .map(|record| record.line.clone())
        .collect::<Vec<_>>();
    let shared_keys = expected_keys.intersection(&actual_keys).collect::<Vec<_>>();
    let shared_count = shared_keys.len();
    let same_gdp_count = shared_keys
        .iter()
        .filter(|key| expected_by_key[***key].gdp == actual_by_key[***key].gdp)
        .count();
    let same_dp_count = shared_keys
        .iter()
        .filter(|key| expected_by_key[***key].dp == actual_by_key[***key].dp)
        .count();
    format!(
        "shared={shared_count}; missing={}; extra={}; same_gdp={same_gdp_count}; same_dp={same_dp_count}; expected_types={:?}; actual_types={:?}; expected_gdp_buckets={:?}; actual_gdp_buckets={:?}; missing_examples={missing:?}; extra_examples={extra:?}",
        expected_keys.len() - shared_count,
        actual_keys.len() - shared_count,
        type_counts(&expected),
        type_counts(&actual),
        gdp_buckets(&expected),
        gdp_buckets(&actual),
    )
}

#[derive(Debug)]
struct ParsedVcfRecord {
    line: String,
    key: (String, String, String, String, String),
    gdp: Option<i64>,
    dp: Option<i64>,
}

impl ParsedVcfRecord {
    fn parse(line: &str) -> Self {
        let fields = line.split('\t').collect::<Vec<_>>();
        let sample = fields.get(9).copied().unwrap_or_default();
        let sample_fields = sample.split(':').collect::<Vec<_>>();
        Self {
            line: line.to_owned(),
            key: (
                fields.first().copied().unwrap_or_default().to_owned(),
                fields.get(1).copied().unwrap_or_default().to_owned(),
                fields.get(3).copied().unwrap_or_default().to_owned(),
                fields.get(4).copied().unwrap_or_default().to_owned(),
                fields.get(8).copied().unwrap_or_default().to_owned(),
            ),
            gdp: sample_fields.get(1).and_then(|value| value.parse().ok()),
            dp: sample_fields.get(2).and_then(|value| value.parse().ok()),
        }
    }

    fn variant_kind(&self) -> &'static str {
        match self.key.3.len().cmp(&self.key.2.len()) {
            std::cmp::Ordering::Greater => "ins",
            std::cmp::Ordering::Less => "del",
            std::cmp::Ordering::Equal => "snp",
        }
    }
}

fn type_counts(records: &[ParsedVcfRecord]) -> BTreeMap<&'static str, usize> {
    let mut counts = BTreeMap::new();
    for record in records {
        *counts.entry(record.variant_kind()).or_insert(0) += 1;
    }
    counts
}

fn gdp_buckets(records: &[ParsedVcfRecord]) -> BTreeMap<&'static str, usize> {
    let mut counts = BTreeMap::new();
    for record in records {
        let bucket = match record.gdp {
            Some(0) => "0",
            Some(1) => "1",
            Some(2..=5) => "2-5",
            Some(6..=20) => "6-20",
            Some(21..=100) => "21-100",
            Some(101..) => ">100",
            Some(..=-1) | None => "missing",
        };
        *counts.entry(bucket).or_insert(0) += 1;
    }
    counts
}

fn vcf_records(path: &Path) -> Vec<String> {
    fs::read_to_string(path)
        .unwrap_or_else(|err| panic!("failed to read {}: {err}", path.display()))
        .lines()
        .filter(|line| !line.starts_with('#') && !line.trim().is_empty())
        .map(str::to_owned)
        .collect()
}

fn decompress_gzip(input: &Path, temp: &Path) -> PathBuf {
    let output = temp.join(
        input
            .file_name()
            .and_then(|name| name.to_str())
            .unwrap()
            .trim_end_matches(".gz"),
    );
    let output_file = fs::File::create(&output).unwrap();
    let status = Command::new("gzip")
        .arg("-dc")
        .arg(input)
        .stdout(Stdio::from(output_file))
        .status()
        .unwrap_or_else(|err| panic!("failed to run gzip for {}: {err}", input.display()));
    assert!(status.success(), "gzip failed for {}", input.display());
    output
}

fn default_bioscript_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .ancestors()
        .nth(5)
        .unwrap()
        .to_path_buf()
}

fn enabled_bioscript_root() -> Option<PathBuf> {
    std::env::var_os("KESTREL_RUN_VNTYPER_FASTQ_PARITY")?;
    Some(
        std::env::var_os("BIOSCRIPT_REPO_ROOT")
            .map(PathBuf::from)
            .unwrap_or_else(default_bioscript_root),
    )
}

fn parity_output_dir(label: &str) -> Option<PathBuf> {
    std::env::var_os("KESTREL_VNTYPER_PARITY_OUT").map(|root| PathBuf::from(root).join(label))
}

fn parity_usize_env(name: &str, default: usize) -> usize {
    std::env::var(name)
        .ok()
        .and_then(|value| value.parse().ok())
        .unwrap_or(default)
}
