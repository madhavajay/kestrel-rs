use std::path::PathBuf;

use clap::{ArgAction, CommandFactory, Parser};
use kanalyze::comp::reader::{FileSequenceSource, ReaderError};
use thiserror::Error;

use crate::io::{InputSample, InputSampleError, StreamableOutput};
use crate::log_level::LogLevelError;
use crate::runner::{KestrelRunner, RunnerError};

/// Errors returned while parsing CLI arguments or creating a runner.
#[derive(Debug, Error)]
pub enum CliError {
    /// Error reported by the argument parser.
    #[error(transparent)]
    Clap(#[from] clap::Error),
    /// Error reported while configuring or running Kestrel.
    #[error(transparent)]
    Runner(#[from] RunnerError),
    /// Error reported while creating an input reader.
    #[error(transparent)]
    Reader(#[from] ReaderError),
    /// Error reported while creating an input sample.
    #[error(transparent)]
    InputSample(#[from] InputSampleError),
    /// Error reported while parsing the log level.
    #[error(transparent)]
    LogLevel(#[from] LogLevelError),
    /// Multiple sample names were provided without matching input files.
    #[error("multiple sample names require an equal number of input files")]
    AmbiguousSampleInputs,
}

/// Parsed command-line arguments for the `kestrel` binary.
#[derive(Clone, Debug, Parser, PartialEq)]
#[command(
    name = "kestrel",
    disable_help_flag = true,
    version,
    about = "Rust port of the Kestrel variant caller"
)]
pub struct CliArgs {
    /// Optional help topic requested with `-h` or `--help`.
    #[arg(short = 'h', long = "help", num_args = 0..=1, default_missing_value = "main")]
    pub help: Option<Option<String>>,

    /// K-mer size.
    #[arg(short = 'k')]
    pub k_size: Option<usize>,

    /// Reference FASTA files.
    #[arg(short = 'r', action = ArgAction::Append)]
    pub references: Vec<PathBuf>,

    /// Variant output path.
    #[arg(short = 'o')]
    pub output: Option<PathBuf>,

    /// Variant output format.
    #[arg(short = 'm', long = "format")]
    pub output_format: Option<String>,

    /// Write variant output to stdout.
    #[arg(long = "stdout", action = ArgAction::SetTrue)]
    pub stdout: bool,

    /// Sample names.
    #[arg(short = 's', long = "sample", action = ArgAction::Append)]
    pub sample_names: Vec<String>,

    /// Haplotype output format.
    #[arg(long = "hapfmt")]
    pub haplotype_format: Option<String>,

    /// Haplotype output path.
    #[arg(short = 'p')]
    pub haplotype_output: Option<PathBuf>,

    /// Write logs to stderr.
    #[arg(long = "logstderr", action = ArgAction::SetTrue)]
    pub log_stderr: bool,

    /// Write logs to stdout.
    #[arg(long = "logstdout", action = ArgAction::SetTrue)]
    pub log_stdout: bool,

    /// Log level name.
    #[arg(long = "loglevel")]
    pub log_level: Option<String>,

    /// Temporary directory path.
    #[arg(long = "temploc")]
    pub temp_dir: Option<String>,

    /// Maximum number of saved aligner states.
    #[arg(long = "maxalignstates")]
    pub max_align_states: Option<i32>,

    /// Maximum number of haplotypes retained per region.
    #[arg(long = "maxhapstates")]
    pub max_haplotypes: Option<i32>,

    /// Minimum k-mer count.
    #[arg(long = "mincount")]
    pub min_count: Option<i32>,

    /// Minimum active-region count difference.
    #[arg(long = "mindiff")]
    pub min_difference: Option<i32>,

    /// Variant filter specifications.
    #[arg(long = "varfilter", action = ArgAction::Append)]
    pub variant_filters: Vec<String>,

    /// Interval BED files.
    #[arg(short = 'i', action = ArgAction::Append)]
    pub interval_files: Vec<PathBuf>,

    /// Emit calls grouped by reference.
    #[arg(long = "byreference", action = ArgAction::SetTrue)]
    pub by_reference: bool,

    /// Emit calls grouped by region.
    #[arg(long = "byregion", action = ArgAction::SetTrue)]
    pub by_region: bool,

    /// Count reverse-complement k-mers.
    #[arg(long = "countrev", action = ArgAction::SetTrue)]
    pub count_reverse: bool,

    /// Keep k-mer counts in memory.
    #[arg(long = "memcount", action = ArgAction::SetTrue)]
    pub memory_count: bool,

    /// Derive flank length automatically.
    #[arg(long = "autoflank", action = ArgAction::SetTrue)]
    pub auto_flank: bool,

    /// Explicit variant flank length.
    #[arg(long = "flank")]
    pub flank: Option<i32>,

    /// Suppress active regions containing ambiguous bases.
    #[arg(long = "noambigregions", action = ArgAction::SetTrue)]
    pub no_ambiguous_regions: bool,

    /// Suppress variants containing ambiguous bases.
    #[arg(long = "noambivar", action = ArgAction::SetTrue)]
    pub no_ambiguous_variant: bool,

    /// Reverse-complement references on negative-strand intervals.
    #[arg(long = "revregion", action = ArgAction::SetTrue)]
    pub reverse_complement_negative_region: bool,

    /// Input read files.
    #[arg(value_name = "FASTQ")]
    pub inputs: Vec<PathBuf>,
}

/// Parsed CLI action.
#[derive(Clone, Debug, PartialEq)]
pub enum CliCommand {
    /// Print a help topic.
    Help(String),
    /// Run a configured Kestrel runner.
    Run(Box<KestrelRunner>),
}

/// Parses command-line arguments into a help or run command.
pub fn parse_command<I, T>(args: I) -> Result<CliCommand, CliError>
where
    I: IntoIterator<Item = T>,
    T: Into<std::ffi::OsString> + Clone,
{
    let args = CliArgs::try_parse_from(args)?;
    if let Some(topic) = args.help {
        return Ok(CliCommand::Help(topic.unwrap_or_else(|| "main".to_owned())));
    }

    Ok(CliCommand::Run(Box::new(runner_from_args(args)?)))
}

/// Returns help text for a named topic.
pub fn help_topic(topic: &str) -> String {
    match topic.trim().to_ascii_lowercase().as_str() {
        "" | "main" => CliArgs::command().render_long_help().to_string(),
        "topics" => "Kestrel help topics: reader, writer, format\n".to_owned(),
        "reader" => "reader: reference and sample sequence input options\n".to_owned(),
        "writer" => "writer: variant and haplotype output options\n".to_owned(),
        "format" => "format: supported output formats include vcf, table, and txt\n".to_owned(),
        other => format!("unknown Kestrel help topic: {other}\n"),
    }
}

fn runner_from_args(args: CliArgs) -> Result<KestrelRunner, CliError> {
    let mut runner = KestrelRunner::new();

    if let Some(k_size) = args.k_size {
        runner.set_k_size(k_size)?;
    }
    for reference in args.references {
        let source_id = runner.config().references.len() as i32 + 1;
        runner.add_reference(FileSequenceSource::from_path(reference, source_id)?);
    }
    if args.stdout {
        runner.set_output_file(Some(StreamableOutput::stdout()));
    }
    if let Some(output) = args.output {
        runner.set_output_path(output);
    }
    if let Some(output_format) = args.output_format {
        runner.set_output_format(&output_format)?;
    }
    if args.log_stderr {
        runner.set_log_file(Some(StreamableOutput::stderr()));
    }
    if args.log_stdout {
        runner.set_log_file(Some(StreamableOutput::stdout()));
    }
    if let Some(log_level) = args.log_level {
        runner.set_log_level_name(&log_level)?;
    }
    if let Some(temp_dir) = args.temp_dir {
        runner.set_temp_dir_name(Some(&temp_dir));
    }
    if let Some(max_align_states) = args.max_align_states {
        runner.set_max_aligner_state(max_align_states)?;
    }
    if let Some(max_haplotypes) = args.max_haplotypes {
        runner.set_max_haplotypes(max_haplotypes)?;
    }
    if let Some(min_count) = args.min_count {
        runner.set_min_kmer_count(min_count)?;
    }
    if let Some(min_difference) = args.min_difference {
        runner.set_minimum_difference(min_difference)?;
    }
    if let Some(haplotype_format) = args.haplotype_format {
        runner.set_haplotype_output_format(&haplotype_format)?;
    }
    if let Some(haplotype_output) = args.haplotype_output {
        runner.set_haplotype_output_file(Some(StreamableOutput::from_path(haplotype_output, None)));
    }
    if args.by_reference {
        runner.set_variant_call_by_reference();
    }
    if args.by_region {
        runner.set_variant_call_by_region();
    }
    if args.count_reverse {
        runner.set_count_reverse_kmers(true);
    }
    if args.memory_count {
        runner.set_kmer_count_in_memory(true);
    }
    if args.auto_flank {
        runner.set_default_flank_length();
    }
    if let Some(flank) = args.flank {
        runner.set_flank_length(flank)?;
    }
    if args.no_ambiguous_regions {
        runner.set_call_ambiguous_regions(false);
    }
    if args.no_ambiguous_variant {
        runner.set_call_ambiguous_variant(false);
    }
    if args.reverse_complement_negative_region {
        runner.set_rev_complement_neg_reference_strand(true);
    }
    for filter in args.variant_filters {
        runner.add_variant_filter_spec(filter);
    }
    for interval_file in args.interval_files {
        runner.add_interval_file(interval_file);
    }

    add_samples(&mut runner, &args.sample_names, args.inputs)?;
    Ok(runner)
}

fn add_samples(
    runner: &mut KestrelRunner,
    sample_names: &[String],
    inputs: Vec<PathBuf>,
) -> Result<(), CliError> {
    if inputs.is_empty() {
        return Ok(());
    }

    if sample_names.len() <= 1 {
        let sources = input_sources(inputs, 1)?;
        let name = sample_names.first().map(String::as_str);
        runner.add_sample(InputSample::new(name, sources)?);
        return Ok(());
    }

    if inputs.len() != sample_names.len() {
        return Err(CliError::AmbiguousSampleInputs);
    }

    for (index, (name, input)) in sample_names.iter().zip(inputs).enumerate() {
        let source = FileSequenceSource::from_path(input, index as i32 + 1)?;
        runner.add_sample(InputSample::new(Some(name), vec![source])?);
    }

    Ok(())
}

fn input_sources(
    inputs: Vec<PathBuf>,
    first_source_id: i32,
) -> Result<Vec<FileSequenceSource>, ReaderError> {
    inputs
        .into_iter()
        .enumerate()
        .map(|(index, input)| FileSequenceSource::from_path(input, first_source_id + index as i32))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_vntyper_shaped_command_line() {
        let command = parse_command([
            "kestrel",
            "-k",
            "25",
            "--maxalignstates",
            "40",
            "--maxhapstates",
            "30",
            "-r",
            "ref.fasta",
            "-o",
            "out.vcf",
            "-s",
            "sample",
            "reads_1.fastq",
            "reads_2.fastq",
            "--hapfmt",
            "sam",
            "-p",
            "out.sam",
            "--logstderr",
            "--loglevel",
            "ERROR",
            "--temploc",
            "/tmp",
        ])
        .unwrap();

        let CliCommand::Run(runner) = command else {
            panic!("expected run command");
        };
        let config = runner.config();
        assert_eq!(config.k_size, 25);
        assert_eq!(config.max_aligner_state, 40);
        assert_eq!(config.max_haplotypes, 30);
        assert_eq!(config.references.len(), 1);
        assert_eq!(config.samples.len(), 1);
        assert_eq!(config.samples[0].name, "sample");
        assert_eq!(config.samples[0].sources.len(), 2);
        assert_eq!(config.haplotype_output_format, "sam");
        assert!(matches!(
            config.haplotype_output_file,
            Some(StreamableOutput::File { ref path, .. }) if path == &PathBuf::from("out.sam")
        ));
        assert_eq!(config.log_level, crate::LogLevel::Error);
        assert_eq!(config.temp_dir_name.as_deref(), Some("/tmp"));
    }

    #[test]
    fn accepts_coverage_script_options_without_running_pipeline() {
        let command = parse_command([
            "kestrel",
            "-k",
            "25",
            "-r",
            "ref.fasta",
            "-m",
            "TABLE",
            "--stdout",
            "-ssmoke",
            "--varfilter",
            "TYPE:snp",
            "-i",
            "regions.bed",
            "--byregion",
            "--memcount",
            "--noambigregions",
            "--noambivar",
            "--flank=10",
            "reads.fastq",
        ])
        .unwrap();

        let CliCommand::Run(runner) = command else {
            panic!("expected run command");
        };
        let config = runner.config();
        assert_eq!(config.output_format, "TABLE");
        assert!(config.variant_call_by_region);
        assert!(config.kmer_count_in_memory);
        assert!(!config.call_ambiguous_regions);
        assert!(!config.call_ambiguous_variant);
        assert_eq!(config.flank_length, 10);
        assert_eq!(config.variant_filter_specs, vec!["TYPE:snp"]);
        assert_eq!(config.interval_files, vec![PathBuf::from("regions.bed")]);
    }

    #[test]
    fn parses_help_topics() {
        assert_eq!(
            parse_command(["kestrel", "-hreader"]).unwrap(),
            CliCommand::Help("reader".to_owned())
        );
        assert_eq!(
            parse_command(["kestrel", "--help=topics"]).unwrap(),
            CliCommand::Help("topics".to_owned())
        );
        assert!(help_topic("writer").contains("writer"));
    }

    #[test]
    fn multiple_sample_names_split_equal_inputs() {
        let command =
            parse_command(["kestrel", "-ssample_a", "a.fastq", "-ssample_b", "b.fastq"]).unwrap();

        let CliCommand::Run(runner) = command else {
            panic!("expected run command");
        };
        assert_eq!(runner.config().samples.len(), 2);
        assert_eq!(runner.config().samples[0].name, "sample_a");
        assert_eq!(runner.config().samples[1].name, "sample_b");
    }
}
