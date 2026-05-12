use std::path::{Path, PathBuf};
use std::process::Command;

use criterion::{Criterion, criterion_group, criterion_main};
use kanalyze::comp::reader::FileSequenceSource;
use kanalyze::module::count::CountModule;
use kanalyze::util::KmerUtil;
use kestrel::io::{InputSample, StreamableOutput};
use kestrel::runner::KestrelRunner;
use tempfile::TempDir;

const K_SIZE: usize = 25;

fn bench_pipeline(c: &mut Criterion) {
    let fixture = Fixture::new();

    c.bench_function("rust_count_module", |b| {
        let source = FileSequenceSource::from_path(&fixture.reads_path, 1).unwrap();
        let kmer_util = KmerUtil::new(K_SIZE).unwrap();
        b.iter(|| {
            CountModule::new(source.clone(), kmer_util.clone())
                .count()
                .unwrap()
        });
    });

    c.bench_function("rust_runner_no_variant", |b| {
        b.iter(|| {
            let runner = configured_runner(&fixture, fixture.rust_output_path());
            runner.run().unwrap();
        });
    });

    if let Some(jar_path) = java_jar_path() {
        c.bench_function("java_kestrel_cli_no_variant", |b| {
            b.iter(|| {
                run_java_kestrel(&jar_path, &fixture).unwrap();
            });
        });
    }
}

fn configured_runner(fixture: &Fixture, output_path: PathBuf) -> KestrelRunner {
    let mut runner = KestrelRunner::new();
    runner.set_k_size(K_SIZE).unwrap();
    runner.set_minimizer_size(0).unwrap();
    runner.set_kmer_count_in_memory(true);
    runner.set_output_file(Some(StreamableOutput::from_path(output_path, None)));
    runner.add_reference(FileSequenceSource::from_path(&fixture.ref_path, 1).unwrap());
    runner.add_sample(
        InputSample::new(
            Some("bench"),
            vec![FileSequenceSource::from_path(&fixture.reads_path, 1).unwrap()],
        )
        .unwrap(),
    );
    runner
}

fn run_java_kestrel(jar_path: &Path, fixture: &Fixture) -> std::io::Result<()> {
    let status = Command::new("java")
        .arg("-jar")
        .arg(jar_path)
        .arg("-k")
        .arg(K_SIZE.to_string())
        .arg("--maxalignstates")
        .arg("40")
        .arg("--maxhapstates")
        .arg("40")
        .arg("-r")
        .arg(&fixture.ref_path)
        .arg("-o")
        .arg(fixture.java_output_path())
        .arg("-sbench")
        .arg("--temploc")
        .arg(fixture.temp.path())
        .arg("--logstderr")
        .arg("--loglevel")
        .arg("ERROR")
        .arg(&fixture.reads_path)
        .status()?;

    assert!(status.success(), "Java Kestrel benchmark command failed");
    Ok(())
}

fn java_jar_path() -> Option<PathBuf> {
    std::env::var_os("KESTREL_JAR")
        .map(PathBuf::from)
        .filter(|path| path.is_file())
        .or_else(|| {
            let path = PathBuf::from("kestrel/lib/kestrel.jar");
            path.is_file().then_some(path)
        })
}

struct Fixture {
    temp: TempDir,
    ref_path: PathBuf,
    reads_path: PathBuf,
}

impl Fixture {
    fn new() -> Self {
        let temp = TempDir::new().unwrap();
        let ref_path = temp.path().join("bench_ref.fasta");
        let reads_path = temp.path().join("bench_reads.fastq");
        let reference = repeated_reference();
        let reads = repeated_reads(&reference, 48);

        std::fs::write(&ref_path, format!(">bench_ref\n{reference}\n")).unwrap();
        std::fs::write(&reads_path, reads).unwrap();

        Self {
            temp,
            ref_path,
            reads_path,
        }
    }

    fn rust_output_path(&self) -> PathBuf {
        self.temp.path().join("rust-out.vcf")
    }

    fn java_output_path(&self) -> PathBuf {
        self.temp.path().join("java-out.vcf")
    }
}

fn repeated_reference() -> String {
    "ACGT".repeat(128)
}

fn repeated_reads(reference: &str, copies: usize) -> String {
    let quality = "F".repeat(reference.len());
    (0..copies)
        .map(|index| format!("@read{index}\n{reference}\n+\n{quality}\n"))
        .collect()
}

criterion_group!(benches, bench_pipeline);
criterion_main!(benches);
