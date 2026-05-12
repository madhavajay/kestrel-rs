use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;

use tempfile::TempDir;

const REF_SEQ: &str = "GCTAAAGACAATTACATAACATACACGTCAGCACGAAACTTGTTGGCCCAGTGTGAATCGCTTAAGGGTTAAGTAAGTGTGATGCATACGCCTTTACTTGCTGTGTCCACCCCATCGGACTGGCATTTTTATTACACTCAGAAACAGAACTCGGGTAATTTTGACAGGTCACGCAGAGGCGCGCCCTCCTGAA";

#[test]
fn java_and_rust_cli_outputs_match_for_coverage_matrix() {
    if std::env::var_os("KESTREL_RUN_JAVA_PARITY").is_none() {
        eprintln!("set KESTREL_RUN_JAVA_PARITY=1 to run Java/Rust CLI parity");
        return;
    }

    let java_jar = workspace_path("kestrel/lib/kestrel.jar");
    if !java_jar.exists() {
        panic!(
            "missing {}; run `git submodule update --init --recursive` and `cd kestrel && scripts/build-kestrel.sh`",
            java_jar.display()
        );
    }

    let rust_bin = PathBuf::from(env!("CARGO_BIN_EXE_kestrel"));
    let temp = TempDir::new().unwrap();
    let fixture = Fixture::new(temp.path());

    for case in parity_cases(&fixture) {
        run_case(&case, &java_jar, &rust_bin, temp.path());
    }
}

fn run_case(case: &ParityCase, java_jar: &Path, rust_bin: &Path, temp: &Path) {
    let java_out = temp.join(format!("java-{}", case.output_name));
    let rust_out = temp.join(format!("rust-{}", case.output_name));
    let java_stdout = temp.join(format!("java-{}.stdout", case.name));
    let rust_stdout = temp.join(format!("rust-{}.stdout", case.name));

    let java_hap = case
        .haplotype_output_name
        .map(|name| temp.join(format!("java-{name}")));
    let rust_hap = case
        .haplotype_output_name
        .map(|name| temp.join(format!("rust-{name}")));
    let java_args = case.args(&java_out, java_hap.as_deref());
    let rust_args = case.args(&rust_out, rust_hap.as_deref());

    let java_status = Command::new("java")
        .arg("-cp")
        .arg(java_jar)
        .arg("edu.gatech.kestrel.clui.Main")
        .args(&java_args)
        .stdout(fs::File::create(&java_stdout).unwrap())
        .stderr(fs::File::create(temp.join(format!("java-{}.stderr", case.name))).unwrap())
        .status()
        .unwrap_or_else(|err| panic!("failed to run Java Kestrel for {}: {err}", case.name));
    assert!(
        java_status.success(),
        "Java Kestrel failed for {}",
        case.name
    );

    let rust_status = Command::new(rust_bin)
        .args(&rust_args)
        .stdout(fs::File::create(&rust_stdout).unwrap())
        .stderr(fs::File::create(temp.join(format!("rust-{}.stderr", case.name))).unwrap())
        .status()
        .unwrap_or_else(|err| panic!("failed to run Rust Kestrel for {}: {err}", case.name));
    assert!(
        rust_status.success(),
        "Rust Kestrel failed for {}",
        case.name
    );

    let (java_bytes, rust_bytes) = if case.stdout {
        (
            read_output(&java_stdout, case.name, "Java stdout"),
            read_output(&rust_stdout, case.name, "Rust stdout"),
        )
    } else {
        (
            read_output(&java_out, case.name, "Java output"),
            read_output(&rust_out, case.name, "Rust output"),
        )
    };
    assert_eq!(
        java_bytes, rust_bytes,
        "Java/Rust CLI output differs for {}",
        case.name
    );

    if let (Some(java_hap), Some(rust_hap)) = (java_hap, rust_hap) {
        let java_hap = fs::read(java_hap).unwrap();
        let rust_hap = fs::read(rust_hap).unwrap();
        assert_eq!(
            java_hap, rust_hap,
            "Java/Rust haplotype output differs for {}",
            case.name
        );
    }
}

fn read_output(path: &Path, case_name: &str, label: &str) -> Vec<u8> {
    fs::read(path).unwrap_or_else(|err| {
        panic!(
            "failed to read {label} for {case_name} at {}: {err}",
            path.display()
        )
    })
}

struct Fixture {
    ref_path: PathBuf,
    reads_path: PathBuf,
    regions_path: PathBuf,
}

impl Fixture {
    fn new(temp: &Path) -> Self {
        let ref_path = temp.join("ref.fasta");
        let reads_path = temp.join("reads.fastq");
        let regions_path = temp.join("regions.bed");
        fs::write(&ref_path, format!(">smoke_ref\n{REF_SEQ}\n")).unwrap();
        let quality = "I".repeat(REF_SEQ.len());
        let mut reads = String::new();
        for index in 1..=50 {
            reads.push_str(&format!("@read_{index}\n{REF_SEQ}\n+\n{quality}\n"));
        }
        fs::write(&reads_path, reads).unwrap();
        fs::write(&regions_path, "smoke_ref\t0\t90\nsmoke_ref\t100\t190\n").unwrap();
        Self {
            ref_path,
            reads_path,
            regions_path,
        }
    }
}

struct ParityCase {
    name: &'static str,
    output_name: &'static str,
    extra_args: Vec<String>,
    stdout: bool,
    haplotype_output_name: Option<&'static str>,
}

impl ParityCase {
    fn args(&self, output_path: &Path, haplotype_path: Option<&Path>) -> Vec<String> {
        let mut args = self
            .extra_args
            .iter()
            .map(|arg| {
                if arg == "{HAP}" {
                    haplotype_path.unwrap().display().to_string()
                } else {
                    arg.clone()
                }
            })
            .collect::<Vec<_>>();
        if self.stdout {
            args.push("--stdout".to_owned());
        } else {
            args.push("-o".to_owned());
            args.push(output_path.display().to_string());
        }
        args
    }
}

fn parity_cases(fixture: &Fixture) -> Vec<ParityCase> {
    let ref_path = fixture.ref_path.display().to_string();
    let reads_path = fixture.reads_path.display().to_string();
    let regions_path = fixture.regions_path.display().to_string();
    let common = |k_size: &str| {
        vec![
            "-k".to_owned(),
            k_size.to_owned(),
            "-r".to_owned(),
            ref_path.clone(),
            "-ssmoke".to_owned(),
        ]
    };
    let finish = |mut args: Vec<String>| {
        args.push("--logstderr".to_owned());
        args.push("--loglevel".to_owned());
        args.push("ERROR".to_owned());
        args.push(reads_path.clone());
        args
    };

    let mut cases = Vec::new();

    add_case(
        &mut cases,
        "base",
        "base.vcf",
        finish(with_args(
            common("25"),
            ["--maxalignstates", "40", "--maxhapstates", "40"],
        )),
    );
    add_case(
        &mut cases,
        "table",
        "table.txt",
        finish(with_args(common("25"), ["-m", "TABLE"])),
    );
    add_case(
        &mut cases,
        "txt",
        "txt.txt",
        finish(with_args(common("25"), ["-m", "TXT"])),
    );
    cases.push(ParityCase {
        name: "hapfmt-sam",
        output_name: "hapfmt.vcf",
        extra_args: finish(with_args(common("25"), ["--hapfmt", "sam", "-p", "{HAP}"])),
        stdout: false,
        haplotype_output_name: Some("hapfmt.sam"),
    });
    add_case(
        &mut cases,
        "filter-snp",
        "filter-snp.vcf",
        finish(with_args(common("25"), ["--varfilter", "TYPE:snp"])),
    );
    add_case(
        &mut cases,
        "bed-intervals",
        "bed.vcf",
        finish(with_owned_args(
            common("25"),
            ["-i".to_owned(), regions_path.clone()],
        )),
    );
    add_case(
        &mut cases,
        "by-reference",
        "by-reference.vcf",
        finish(with_args(common("25"), ["--byreference"])),
    );
    add_case(
        &mut cases,
        "by-region",
        "by-region.vcf",
        finish(with_owned_args(
            common("25"),
            [
                "-i".to_owned(),
                regions_path.clone(),
                "--byregion".to_owned(),
            ],
        )),
    );
    add_case(
        &mut cases,
        "countrev",
        "countrev.vcf",
        finish(with_args(common("25"), ["--countrev"])),
    );
    add_case(
        &mut cases,
        "memcount",
        "memcount.vcf",
        finish(with_args(common("25"), ["--memcount"])),
    );
    add_case(
        &mut cases,
        "autoflank",
        "autoflank.vcf",
        finish(with_owned_args(
            common("25"),
            [
                "-i".to_owned(),
                regions_path.clone(),
                "--autoflank".to_owned(),
            ],
        )),
    );
    cases.push(ParityCase {
        name: "stdout",
        output_name: "stdout.vcf",
        extra_args: finish(common("25")),
        stdout: true,
        haplotype_output_name: None,
    });
    add_case(&mut cases, "k-16", "k16.vcf", finish(common("16")));
    add_case(&mut cases, "k-31", "k31.vcf", finish(common("31")));
    add_case(
        &mut cases,
        "multifilter",
        "multifilter.vcf",
        finish(with_args(
            common("25"),
            [
                "--varfilter",
                "TYPE:snp",
                "--varfilter",
                "LOCATION:length=10",
            ],
        )),
    );
    add_case(
        &mut cases,
        "noambig",
        "noambig.vcf",
        finish(with_args(common("25"), ["--noambigregions", "--noambivar"])),
    );
    cases.push(ParityCase {
        name: "two-samples",
        output_name: "two-samples.vcf",
        extra_args: vec![
            "-k".to_owned(),
            "25".to_owned(),
            "-r".to_owned(),
            ref_path.clone(),
            "-ssample_a".to_owned(),
            reads_path.clone(),
            "-ssample_b".to_owned(),
            reads_path.clone(),
            "--logstderr".to_owned(),
            "--loglevel".to_owned(),
            "ERROR".to_owned(),
        ],
        stdout: false,
        haplotype_output_name: None,
    });
    add_case(
        &mut cases,
        "flank",
        "flank.vcf",
        finish(with_owned_args(
            common("25"),
            [
                "-i".to_owned(),
                regions_path.clone(),
                "--flank=10".to_owned(),
            ],
        )),
    );
    add_case(
        &mut cases,
        "low-thresholds",
        "low-thresholds.vcf",
        finish(with_args(common("25"), ["--mincount=2", "--mindiff=2"])),
    );
    add_case(
        &mut cases,
        "revregion",
        "revregion.vcf",
        finish(with_args(common("25"), ["--revregion"])),
    );
    cases
}

fn add_case(
    cases: &mut Vec<ParityCase>,
    name: &'static str,
    output_name: &'static str,
    extra_args: Vec<String>,
) {
    cases.push(ParityCase {
        name,
        output_name,
        extra_args,
        stdout: false,
        haplotype_output_name: None,
    });
}

fn with_args<const N: usize>(mut args: Vec<String>, extra: [&str; N]) -> Vec<String> {
    args.extend(extra.into_iter().map(str::to_owned));
    args
}

fn with_owned_args<const N: usize>(mut args: Vec<String>, extra: [String; N]) -> Vec<String> {
    args.extend(extra);
    args
}

fn workspace_path(relative: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../..")
        .join(relative)
}
