use std::path::{Path, PathBuf};
use std::str::FromStr;

use kanalyze::comp::reader::{FileSequenceSource, SequenceReader};
use kanalyze::util::{Base, KmerError, KmerUtil};
use thiserror::Error;

use crate::activeregion::{
    ActiveRegion, ActiveRegionDetector, ActiveRegionDetectorError, Haplotype,
};
use crate::align::{
    AlignmentWeight, AlignmentWeightError, KmerAligner, KmerAlignerError, KmerAlignmentBuilder,
};
use crate::constants::MIN_KMER_SIZE;
use crate::counter::{CountMap, CountMapError, IkcCountMap, MemoryCountMap};
use crate::hapwriter::{self, HaplotypeWriter, HaplotypeWriterError};
use crate::interval::{self, RegionIntervalContainer, RegionIntervalReadError};
use crate::io::{InputSample, StreamableOutput, StreamableOutputError};
use crate::log_level::{LogLevel, LogLevelError};
use crate::refreader::{ReferenceIntervalMap, ReferenceReader, ReferenceSequenceError};
use crate::varfilter::{VariantFilterError, VariantFilterKind, VariantFilterRunner};
use crate::variant::{VariantCaller, VariantError};
use crate::writer::{self, VariantWriter, VariantWriterError};

/// Errors returned while configuring or running Kestrel.
#[derive(Debug, Error, PartialEq)]
pub enum RunnerError {
    /// K-mer size was below the supported minimum.
    #[error("k-mer size must not be less than {minimum}: {actual}")]
    InvalidKSize {
        /// Minimum allowed k-mer size.
        minimum: usize,
        /// Actual requested k-mer size.
        actual: usize,
    },
    /// Output format string was empty.
    #[error("output format may not be empty")]
    EmptyOutputFormat,
    /// Haplotype output format string was empty.
    #[error("haplotype output format may not be empty")]
    EmptyHaplotypeOutputFormat,
    /// Minimizer size was negative.
    #[error("minimizer size may not be negative: {0}")]
    NegativeMinimizerSize(i32),
    /// Minimum k-mer count was negative.
    #[error("minimum k-mer count may not be negative: {0}")]
    NegativeMinKmerCount(i32),
    /// Minimum active-region difference was invalid.
    #[error("minimum difference must be at least 1: {0}")]
    InvalidMinimumDifference(i32),
    /// Active-region difference quantile was invalid.
    #[error("difference quantile must be in [0, 1): {0}")]
    InvalidDifferenceQuantile(f64),
    /// Peak scan length was negative.
    #[error("peak scan length may not be negative: {0}")]
    NegativePeakScanLength(i32),
    /// Scan limit factor was negative.
    #[error("scan limit factor may not be negative: {0}")]
    NegativeScanLimitFactor(f64),
    /// Decay minimum was outside the valid range.
    #[error("decay minimum must be in [0, 1]: {0}")]
    InvalidDecayMinimum(f64),
    /// Decay alpha was outside the valid range.
    #[error("decay alpha must be in (0, 1): {0}")]
    InvalidDecayAlpha(f64),
    /// Maximum repeat count was negative.
    #[error("maximum repeat count may not be negative: {0}")]
    NegativeMaxRepeatCount(i32),
    /// Flank length was invalid.
    #[error("flank length must be -1 or greater: {0}")]
    InvalidFlankLength(i32),
    /// Maximum aligner-state count was invalid.
    #[error("maximum aligner states must be at least 1: {0}")]
    InvalidMaxAlignerState(i32),
    /// Maximum haplotype count was invalid.
    #[error("maximum haplotypes must be at least 1: {0}")]
    InvalidMaxHaplotypes(i32),
    /// Required inputs were missing for the runner pipeline.
    #[error("the Rust Kestrel runner pipeline is not implemented yet")]
    PipelineNotImplemented,
    /// Count-map error.
    #[error("count map error: {0}")]
    CountMap(String),
    /// Reference-reader error.
    #[error("reference reader error: {0}")]
    Reference(String),
    /// Variant-writer error.
    #[error("variant writer error: {0}")]
    Writer(String),
    /// K-mer utility error.
    #[error("k-mer error: {0}")]
    Kmer(String),
    /// Active-region detector error.
    #[error("active region detector error: {0}")]
    ActiveRegion(String),
    /// Interval-reader error.
    #[error("interval reader error: {0}")]
    Interval(String),
    /// Variant-filter error.
    #[error("variant filter error: {0}")]
    VariantFilter(String),
    /// Haplotype-writer error.
    #[error("haplotype writer error: {0}")]
    HaplotypeWriter(String),
    /// Sequence-reader error.
    #[error("sequence reader error: {0}")]
    Reader(String),
    /// Aligner error.
    #[error("aligner error: {0}")]
    Aligner(String),
    /// Variant-caller error.
    #[error("variant caller error: {0}")]
    Variant(String),
    /// Output target error.
    #[error(transparent)]
    Output(#[from] StreamableOutputError),
    /// Log-level parsing error.
    #[error(transparent)]
    LogLevel(#[from] LogLevelError),
    /// Alignment-weight error.
    #[error(transparent)]
    AlignmentWeight(#[from] AlignmentWeightError),
}

impl From<CountMapError> for RunnerError {
    fn from(value: CountMapError) -> Self {
        Self::CountMap(value.to_string())
    }
}

impl From<ReferenceSequenceError> for RunnerError {
    fn from(value: ReferenceSequenceError) -> Self {
        Self::Reference(value.to_string())
    }
}

impl From<VariantWriterError> for RunnerError {
    fn from(value: VariantWriterError) -> Self {
        Self::Writer(value.to_string())
    }
}

impl From<KmerError> for RunnerError {
    fn from(value: KmerError) -> Self {
        Self::Kmer(value.to_string())
    }
}

impl From<ActiveRegionDetectorError> for RunnerError {
    fn from(value: ActiveRegionDetectorError) -> Self {
        Self::ActiveRegion(value.to_string())
    }
}

impl From<RegionIntervalReadError> for RunnerError {
    fn from(value: RegionIntervalReadError) -> Self {
        Self::Interval(value.to_string())
    }
}

impl From<VariantFilterError> for RunnerError {
    fn from(value: VariantFilterError) -> Self {
        Self::VariantFilter(value.to_string())
    }
}

impl From<HaplotypeWriterError> for RunnerError {
    fn from(value: HaplotypeWriterError) -> Self {
        Self::HaplotypeWriter(value.to_string())
    }
}

impl From<KmerAlignerError> for RunnerError {
    fn from(value: KmerAlignerError) -> Self {
        Self::Aligner(value.to_string())
    }
}

impl From<VariantError> for RunnerError {
    fn from(value: VariantError) -> Self {
        Self::Variant(value.to_string())
    }
}

/// Runtime configuration for a Kestrel run.
#[derive(Clone, Debug, PartialEq)]
pub struct RunConfig {
    /// K-mer size.
    pub k_size: usize,
    /// Minimizer size.
    pub minimizer_size: i32,
    /// Minimizer mask.
    pub minimizer_mask: u32,
    /// Variant output target.
    pub output_file: StreamableOutput,
    /// Variant output format.
    pub output_format: String,
    /// Log output target.
    pub log_file: StreamableOutput,
    /// Log level.
    pub log_level: LogLevel,
    /// Alignment scoring weights.
    pub alignment_weight: AlignmentWeight,
    /// Temporary directory name.
    pub temp_dir_name: Option<String>,
    /// Minimum k-mer count.
    pub min_kmer_count: i32,
    /// Whether to keep k-mer counts in memory.
    pub kmer_count_in_memory: bool,
    /// Whether to free intermediate resources.
    pub free_resources: bool,
    /// Whether to remove temporary IKC files.
    pub remove_ikc: bool,
    /// Whether active regions require both anchors.
    pub anchor_both_ends: bool,
    /// Whether reverse-complement k-mers are counted.
    pub count_reverse_kmers: bool,
    /// Minimum active-region count difference.
    pub minimum_difference: i32,
    /// Active-region count-difference quantile.
    pub difference_quantile: f64,
    /// Active-region peak scan length.
    pub peak_scan_length: i32,
    /// Active-region scan limit factor.
    pub scan_limit_factor: f64,
    /// Active-region decay alpha.
    pub decay_alpha: f64,
    /// Active-region decay minimum.
    pub decay_minimum: f64,
    /// Maximum repeat count.
    pub max_repeat_count: i32,
    /// Whether ambiguous active regions are called.
    pub call_ambiguous_regions: bool,
    /// Whether ambiguous variants are called.
    pub call_ambiguous_variant: bool,
    /// Flank length, or `-1` to use the default formula.
    pub flank_length: i32,
    /// Whether variants are emitted by region instead of by reference.
    pub variant_call_by_region: bool,
    /// Maximum saved aligner states.
    pub max_aligner_state: i32,
    /// Maximum haplotypes per active region.
    pub max_haplotypes: i32,
    /// Whether reference descriptions are stripped from names.
    pub remove_reference_sequence_description: bool,
    /// Whether negative-strand intervals are reverse-complemented.
    pub reverse_complement_negative_strand: bool,
    /// Optional haplotype output target.
    pub haplotype_output_file: Option<StreamableOutput>,
    /// Haplotype output format.
    pub haplotype_output_format: String,
    /// Input samples.
    pub samples: Vec<InputSample>,
    /// Reference sequence sources.
    pub references: Vec<FileSequenceSource>,
    /// Extra library paths retained for Java compatibility.
    pub libraries: Vec<PathBuf>,
    /// Variant filter specifications.
    pub variant_filter_specs: Vec<String>,
    /// Interval files.
    pub interval_files: Vec<PathBuf>,
}

impl Default for RunConfig {
    fn default() -> Self {
        Self {
            k_size: DEFAULT_KSIZE,
            minimizer_size: DEFAULT_MINIMIZER_SIZE,
            minimizer_mask: DEFAULT_MINIMIZER_MASK,
            output_file: DEFAULT_OUTPUT_FILE,
            output_format: DEFAULT_OUTPUT_FORMAT.to_owned(),
            log_file: DEFAULT_LOG_FILE,
            log_level: DEFAULT_LOG_LEVEL,
            alignment_weight: AlignmentWeight::defaults(),
            temp_dir_name: None,
            min_kmer_count: DEFAULT_MIN_KMER_COUNT,
            kmer_count_in_memory: DEFAULT_KMER_COUNT_IN_MEMORY,
            free_resources: DEFAULT_FREE_RESOURCES,
            remove_ikc: DEFAULT_REMOVE_IKC,
            anchor_both_ends: ActiveRegionDetector::DEFAULT_ANCHOR_BOTH_ENDS,
            count_reverse_kmers: DEFAULT_COUNT_REV_KMER,
            minimum_difference: ActiveRegionDetector::DEFAULT_MINIMUM_DIFFERENCE,
            difference_quantile: ActiveRegionDetector::DEFAULT_DIFFERENCE_QUANTILE,
            peak_scan_length: ActiveRegionDetector::DEFAULT_PEAK_SCAN_LENGTH,
            scan_limit_factor: ActiveRegionDetector::DEFAULT_SCAN_LIMIT_FACTOR,
            decay_alpha: ActiveRegionDetector::DEFAULT_DECAY_ALPHA,
            decay_minimum: ActiveRegionDetector::DEFAULT_DECAY_MINIMUM,
            max_repeat_count: ActiveRegionDetector::DEFAULT_MAX_REPEAT_COUNT,
            call_ambiguous_regions: ActiveRegionDetector::DEFAULT_CALL_AMBIGUOUS_REGIONS,
            call_ambiguous_variant: VariantCaller::DEFAULT_CALL_AMBIGUOUS_VARIANT,
            flank_length: -1,
            variant_call_by_region: false,
            max_aligner_state: KmerAligner::DEFAULT_MAX_STATE,
            max_haplotypes: KmerAlignmentBuilder::DEFAULT_MAX_HAPLOTYPES,
            remove_reference_sequence_description:
                ReferenceReader::DEFAULT_REMOVE_SEQUENCE_DESCRIPTION,
            reverse_complement_negative_strand:
                ReferenceReader::DEFAULT_REVERSE_COMPLEMENT_NEGATIVE_STRAND,
            haplotype_output_file: None,
            haplotype_output_format: DEFAULT_HAPLOTYPE_OUTPUT_FORMAT.to_owned(),
            samples: Vec::new(),
            references: Vec::new(),
            libraries: Vec::new(),
            variant_filter_specs: Vec::new(),
            interval_files: Vec::new(),
        }
    }
}

/// Configurable Kestrel runner.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct KestrelRunner {
    config: RunConfig,
}

impl KestrelRunner {
    /// Creates a runner with default configuration.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Returns the current run configuration.
    #[must_use]
    pub fn config(&self) -> &RunConfig {
        &self.config
    }

    /// Runs Kestrel with the current configuration.
    pub fn run(&self) -> Result<(), RunnerError> {
        run_pipeline(&self.config)
    }

    /// Adds an input sample.
    pub fn add_sample(&mut self, sample: InputSample) {
        self.config.samples.push(sample);
    }

    /// Removes all input samples.
    pub fn clear_samples(&mut self) {
        self.config.samples.clear();
    }

    /// Adds a reference source.
    pub fn add_reference(&mut self, reference: FileSequenceSource) {
        self.config.references.push(reference);
    }

    /// Removes all reference sources.
    pub fn clear_reference(&mut self) {
        self.config.references.clear();
    }

    /// Adds a variant filter specification.
    pub fn add_variant_filter_spec(&mut self, filter_spec: impl Into<String>) {
        self.config.variant_filter_specs.push(filter_spec.into());
    }

    /// Removes all variant filter specifications.
    pub fn clear_variant_filter_specs(&mut self) {
        self.config.variant_filter_specs.clear();
    }

    /// Adds an interval file.
    pub fn add_interval_file(&mut self, interval_file: impl Into<PathBuf>) {
        self.config.interval_files.push(interval_file.into());
    }

    /// Removes all interval files.
    pub fn clear_interval_files(&mut self) {
        self.config.interval_files.clear();
    }

    /// Adds a library path if it is not already present.
    pub fn add_library_file(&mut self, path: impl Into<PathBuf>) {
        let path = path.into();
        if !self.config.libraries.contains(&path) {
            self.config.libraries.push(path);
        }
    }

    /// Removes all library paths.
    pub fn clear_libraries(&mut self) {
        self.config.libraries.clear();
    }

    /// Returns configured library paths.
    #[must_use]
    pub fn libraries(&self) -> &[PathBuf] {
        &self.config.libraries
    }

    /// Sets the k-mer size.
    pub fn set_k_size(&mut self, k_size: usize) -> Result<(), RunnerError> {
        if k_size < MIN_KMER_SIZE {
            return Err(RunnerError::InvalidKSize {
                minimum: MIN_KMER_SIZE,
                actual: k_size,
            });
        }
        self.config.k_size = k_size;
        Ok(())
    }

    /// Returns the k-mer size.
    #[must_use]
    pub fn k_size(&self) -> usize {
        self.config.k_size
    }

    /// Sets the variant output target.
    pub fn set_output_file(&mut self, output_file: Option<StreamableOutput>) {
        self.config.output_file = output_file.unwrap_or(DEFAULT_OUTPUT_FILE);
    }

    /// Sets the variant output path.
    pub fn set_output_path(&mut self, path: impl AsRef<Path>) {
        self.config.output_file = StreamableOutput::from_path(path.as_ref(), None);
    }

    /// Returns the variant output target.
    #[must_use]
    pub fn output_file(&self) -> &StreamableOutput {
        &self.config.output_file
    }

    /// Sets the variant output format.
    pub fn set_output_format(&mut self, output_format: &str) -> Result<(), RunnerError> {
        let output_format = output_format.trim();
        if output_format.is_empty() {
            return Err(RunnerError::EmptyOutputFormat);
        }
        self.config.output_format = output_format.to_owned();
        Ok(())
    }

    /// Returns the variant output format.
    #[must_use]
    pub fn output_format(&self) -> &str {
        &self.config.output_format
    }

    /// Sets the log output target.
    pub fn set_log_file(&mut self, log_file: Option<StreamableOutput>) {
        self.config.log_file = log_file.unwrap_or(DEFAULT_LOG_FILE);
    }

    /// Returns the log output target.
    #[must_use]
    pub fn log_file(&self) -> &StreamableOutput {
        &self.config.log_file
    }

    /// Sets the log level.
    pub fn set_log_level(&mut self, log_level: Option<LogLevel>) {
        self.config.log_level = log_level.unwrap_or(DEFAULT_LOG_LEVEL);
    }

    /// Parses and sets the log level.
    pub fn set_log_level_name(&mut self, log_level: &str) -> Result<(), RunnerError> {
        self.config.log_level = LogLevel::from_str(log_level)?;
        Ok(())
    }

    /// Returns the log level.
    #[must_use]
    pub fn log_level(&self) -> LogLevel {
        self.config.log_level
    }

    /// Sets the minimizer size.
    pub fn set_minimizer_size(&mut self, minimizer_size: i32) -> Result<(), RunnerError> {
        if minimizer_size < 0 {
            return Err(RunnerError::NegativeMinimizerSize(minimizer_size));
        }
        self.config.minimizer_size = minimizer_size;
        Ok(())
    }

    /// Returns the minimizer size.
    #[must_use]
    pub fn minimizer_size(&self) -> i32 {
        self.config.minimizer_size
    }

    /// Sets the minimizer mask.
    pub fn set_minimizer_mask(&mut self, minimizer_mask: u32) {
        self.config.minimizer_mask = minimizer_mask;
    }

    /// Returns the minimizer mask.
    #[must_use]
    pub fn minimizer_mask(&self) -> u32 {
        self.config.minimizer_mask
    }

    /// Parses and sets all alignment weights.
    pub fn set_alignment_weight(&mut self, value: Option<&str>) -> Result<(), RunnerError> {
        self.config.alignment_weight = AlignmentWeight::parse(value)?;
        Ok(())
    }

    /// Sets the alignment match score.
    pub fn set_alignment_weight_match(&mut self, value: f32) -> Result<(), RunnerError> {
        self.config.alignment_weight = self.config.alignment_weight.with_match(value)?;
        Ok(())
    }

    /// Sets the alignment mismatch score.
    pub fn set_alignment_weight_mismatch(&mut self, value: f32) -> Result<(), RunnerError> {
        self.config.alignment_weight = self.config.alignment_weight.with_mismatch(value)?;
        Ok(())
    }

    /// Sets the alignment gap-open score.
    pub fn set_alignment_weight_gap_open(&mut self, value: f32) -> Result<(), RunnerError> {
        self.config.alignment_weight = self.config.alignment_weight.with_gap_open(value)?;
        Ok(())
    }

    /// Sets the alignment gap-extend score.
    pub fn set_alignment_weight_gap_extend(&mut self, value: f32) -> Result<(), RunnerError> {
        self.config.alignment_weight = self.config.alignment_weight.with_gap_extend(value)?;
        Ok(())
    }

    /// Returns alignment scoring weights.
    #[must_use]
    pub fn alignment_weight(&self) -> &AlignmentWeight {
        &self.config.alignment_weight
    }

    /// Sets the temporary directory name.
    pub fn set_temp_dir_name(&mut self, temp_dir_name: Option<&str>) {
        self.config.temp_dir_name = temp_dir_name
            .map(str::trim)
            .filter(|value| !value.is_empty())
            .map(str::to_owned);
    }

    /// Returns the temporary directory name.
    #[must_use]
    pub fn temp_dir_name(&self) -> Option<&str> {
        self.config.temp_dir_name.as_deref()
    }

    /// Sets the minimum k-mer count.
    pub fn set_min_kmer_count(&mut self, min_kmer_count: i32) -> Result<(), RunnerError> {
        if min_kmer_count < 0 {
            return Err(RunnerError::NegativeMinKmerCount(min_kmer_count));
        }
        self.config.min_kmer_count = min_kmer_count;
        Ok(())
    }

    /// Returns the minimum k-mer count.
    #[must_use]
    pub fn min_kmer_count(&self) -> i32 {
        self.config.min_kmer_count
    }

    /// Sets the minimum active-region count difference.
    pub fn set_minimum_difference(&mut self, minimum_difference: i32) -> Result<(), RunnerError> {
        if minimum_difference < 1 {
            return Err(RunnerError::InvalidMinimumDifference(minimum_difference));
        }
        self.config.minimum_difference = minimum_difference;
        Ok(())
    }

    /// Returns the minimum active-region count difference.
    #[must_use]
    pub fn minimum_difference(&self) -> i32 {
        self.config.minimum_difference
    }

    /// Sets the active-region difference quantile.
    pub fn set_difference_quantile(&mut self, difference_quantile: f64) -> Result<(), RunnerError> {
        if !(0.0..1.0).contains(&difference_quantile) {
            return Err(RunnerError::InvalidDifferenceQuantile(difference_quantile));
        }
        self.config.difference_quantile = difference_quantile;
        Ok(())
    }

    /// Returns the active-region difference quantile.
    #[must_use]
    pub fn difference_quantile(&self) -> f64 {
        self.config.difference_quantile
    }

    /// Sets whether active regions must have both anchors.
    pub fn set_anchor_both_ends(&mut self, anchor_both_ends: bool) {
        self.config.anchor_both_ends = anchor_both_ends;
    }

    /// Returns whether active regions must have both anchors.
    #[must_use]
    pub fn anchor_both_ends(&self) -> bool {
        self.config.anchor_both_ends
    }

    /// Sets whether k-mer counts are kept in memory.
    pub fn set_kmer_count_in_memory(&mut self, kmer_count_in_memory: bool) {
        self.config.kmer_count_in_memory = kmer_count_in_memory;
    }

    /// Returns whether k-mer counts are kept in memory.
    #[must_use]
    pub fn kmer_count_in_memory(&self) -> bool {
        self.config.kmer_count_in_memory
    }

    /// Sets whether intermediate resources are freed.
    pub fn set_free_resources(&mut self, free_resources: bool) {
        self.config.free_resources = free_resources;
    }

    /// Returns whether intermediate resources are freed.
    #[must_use]
    pub fn free_resources(&self) -> bool {
        self.config.free_resources
    }

    /// Sets whether temporary IKC files are removed.
    pub fn set_remove_ikc(&mut self, remove_ikc: bool) {
        self.config.remove_ikc = remove_ikc;
    }

    /// Returns whether temporary IKC files are removed.
    #[must_use]
    pub fn remove_ikc(&self) -> bool {
        self.config.remove_ikc
    }

    /// Sets whether reverse-complement k-mers are counted.
    pub fn set_count_reverse_kmers(&mut self, count_reverse_kmers: bool) {
        self.config.count_reverse_kmers = count_reverse_kmers;
    }

    /// Returns whether reverse-complement k-mers are counted.
    #[must_use]
    pub fn count_reverse_kmers(&self) -> bool {
        self.config.count_reverse_kmers
    }

    /// Sets whether ambiguous active regions are called.
    pub fn set_call_ambiguous_regions(&mut self, call_ambiguous_regions: bool) {
        self.config.call_ambiguous_regions = call_ambiguous_regions;
    }

    /// Returns whether ambiguous active regions are called.
    #[must_use]
    pub fn call_ambiguous_regions(&self) -> bool {
        self.config.call_ambiguous_regions
    }

    /// Sets whether ambiguous variants are called.
    pub fn set_call_ambiguous_variant(&mut self, call_ambiguous_variant: bool) {
        self.config.call_ambiguous_variant = call_ambiguous_variant;
    }

    /// Returns whether ambiguous variants are called.
    #[must_use]
    pub fn call_ambiguous_variant(&self) -> bool {
        self.config.call_ambiguous_variant
    }

    /// Sets the active-region peak scan length.
    pub fn set_peak_scan_length(&mut self, peak_scan_length: i32) -> Result<(), RunnerError> {
        if peak_scan_length < 0 {
            return Err(RunnerError::NegativePeakScanLength(peak_scan_length));
        }
        self.config.peak_scan_length = peak_scan_length;
        Ok(())
    }

    /// Returns the active-region peak scan length.
    #[must_use]
    pub fn peak_scan_length(&self) -> i32 {
        self.config.peak_scan_length
    }

    /// Sets the active-region scan limit factor.
    pub fn set_scan_limit_factor(&mut self, scan_limit_factor: f64) -> Result<(), RunnerError> {
        if scan_limit_factor < 0.0 {
            return Err(RunnerError::NegativeScanLimitFactor(scan_limit_factor));
        }
        self.config.scan_limit_factor = scan_limit_factor;
        Ok(())
    }

    /// Returns the active-region scan limit factor.
    #[must_use]
    pub fn scan_limit_factor(&self) -> f64 {
        self.config.scan_limit_factor
    }

    /// Sets the active-region decay minimum.
    pub fn set_decay_minimum(&mut self, decay_minimum: f64) -> Result<(), RunnerError> {
        if !(0.0..=1.0).contains(&decay_minimum) {
            return Err(RunnerError::InvalidDecayMinimum(decay_minimum));
        }
        self.config.decay_minimum = decay_minimum;
        Ok(())
    }

    /// Returns the active-region decay minimum.
    #[must_use]
    pub fn decay_minimum(&self) -> f64 {
        self.config.decay_minimum
    }

    /// Sets the active-region decay alpha.
    pub fn set_decay_alpha(&mut self, decay_alpha: f64) -> Result<(), RunnerError> {
        if decay_alpha <= 0.0 || decay_alpha >= 1.0 {
            return Err(RunnerError::InvalidDecayAlpha(decay_alpha));
        }
        self.config.decay_alpha = decay_alpha;
        Ok(())
    }

    /// Returns the active-region decay alpha.
    #[must_use]
    pub fn decay_alpha(&self) -> f64 {
        self.config.decay_alpha
    }

    /// Sets the maximum repeat count.
    pub fn set_max_repeat_count(&mut self, max_repeat_count: i32) -> Result<(), RunnerError> {
        if max_repeat_count < 0 {
            return Err(RunnerError::NegativeMaxRepeatCount(max_repeat_count));
        }
        self.config.max_repeat_count = max_repeat_count;
        Ok(())
    }

    /// Returns the maximum repeat count.
    #[must_use]
    pub fn max_repeat_count(&self) -> i32 {
        self.config.max_repeat_count
    }

    /// Sets the flank length.
    pub fn set_flank_length(&mut self, flank_length: i32) -> Result<(), RunnerError> {
        if flank_length < -1 {
            return Err(RunnerError::InvalidFlankLength(flank_length));
        }
        self.config.flank_length = flank_length;
        Ok(())
    }

    /// Resets flank length to the default formula.
    pub fn set_default_flank_length(&mut self) {
        self.config.flank_length = -1;
    }

    /// Returns the flank length, or `-1` when defaulted.
    #[must_use]
    pub fn flank_length(&self) -> i32 {
        self.config.flank_length
    }

    /// Emits variants grouped by reference.
    pub fn set_variant_call_by_reference(&mut self) {
        self.config.variant_call_by_region = false;
    }

    /// Returns true when variants are grouped by reference.
    #[must_use]
    pub fn variant_call_by_reference(&self) -> bool {
        !self.config.variant_call_by_region
    }

    /// Emits variants grouped by region.
    pub fn set_variant_call_by_region(&mut self) {
        self.config.variant_call_by_region = true;
    }

    /// Returns true when variants are grouped by region.
    #[must_use]
    pub fn variant_call_by_region(&self) -> bool {
        self.config.variant_call_by_region
    }

    /// Sets the maximum number of saved aligner states.
    pub fn set_max_aligner_state(&mut self, max_aligner_state: i32) -> Result<(), RunnerError> {
        if max_aligner_state < 1 {
            return Err(RunnerError::InvalidMaxAlignerState(max_aligner_state));
        }
        self.config.max_aligner_state = max_aligner_state;
        Ok(())
    }

    /// Returns the maximum number of saved aligner states.
    #[must_use]
    pub fn max_aligner_state(&self) -> i32 {
        self.config.max_aligner_state
    }

    /// Sets the maximum number of haplotypes.
    pub fn set_max_haplotypes(&mut self, max_haplotypes: i32) -> Result<(), RunnerError> {
        if max_haplotypes < 1 {
            return Err(RunnerError::InvalidMaxHaplotypes(max_haplotypes));
        }
        self.config.max_haplotypes = max_haplotypes;
        Ok(())
    }

    /// Returns the maximum number of haplotypes.
    #[must_use]
    pub fn max_haplotypes(&self) -> i32 {
        self.config.max_haplotypes
    }

    /// Sets whether reference descriptions are removed from names.
    pub fn set_remove_reference_description(&mut self, remove: bool) {
        self.config.remove_reference_sequence_description = remove;
    }

    /// Returns whether reference descriptions are removed from names.
    #[must_use]
    pub fn remove_reference_description(&self) -> bool {
        self.config.remove_reference_sequence_description
    }

    /// Sets whether negative-strand reference intervals are reverse-complemented.
    pub fn set_rev_complement_neg_reference_strand(&mut self, reverse_complement: bool) {
        self.config.reverse_complement_negative_strand = reverse_complement;
    }

    /// Returns whether negative-strand reference intervals are reverse-complemented.
    #[must_use]
    pub fn rev_complement_neg_reference_strand(&self) -> bool {
        self.config.reverse_complement_negative_strand
    }

    /// Sets the haplotype output target.
    pub fn set_haplotype_output_file(&mut self, output_file: Option<StreamableOutput>) {
        self.config.haplotype_output_file = output_file;
    }

    /// Returns the haplotype output target.
    #[must_use]
    pub fn haplotype_output_file(&self) -> Option<&StreamableOutput> {
        self.config.haplotype_output_file.as_ref()
    }

    /// Sets the haplotype output format.
    pub fn set_haplotype_output_format(
        &mut self,
        haplotype_output_format: &str,
    ) -> Result<(), RunnerError> {
        let haplotype_output_format = haplotype_output_format.trim();
        if haplotype_output_format.is_empty() {
            return Err(RunnerError::EmptyHaplotypeOutputFormat);
        }
        self.config.haplotype_output_format = haplotype_output_format.to_owned();
        Ok(())
    }

    /// Returns the haplotype output format.
    #[must_use]
    pub fn haplotype_output_format(&self) -> &str {
        &self.config.haplotype_output_format
    }
}

/// Default generic format name.
pub const DEFAULT_FORMAT: &str = "auto";
/// Default k-mer size.
pub const DEFAULT_KSIZE: usize = 31;
/// Default variant output target.
pub const DEFAULT_OUTPUT_FILE: StreamableOutput = StreamableOutput::Stdout;
/// Default variant output format.
pub const DEFAULT_OUTPUT_FORMAT: &str = "vcf";
/// Default log output target.
pub const DEFAULT_LOG_FILE: StreamableOutput = StreamableOutput::Stderr;
/// Default log level.
pub const DEFAULT_LOG_LEVEL: LogLevel = LogLevel::Warn;
/// Default character set name.
pub const DEFAULT_CHARSET: &str = "UTF-8";
/// Default minimizer size.
pub const DEFAULT_MINIMIZER_SIZE: i32 = 15;
/// Default minimizer mask.
pub const DEFAULT_MINIMIZER_MASK: u32 = 0x0000_0000;
/// Default minimum k-mer count.
pub const DEFAULT_MIN_KMER_COUNT: i32 = 5;
/// Default sequence reader buffer size.
pub const DEFAULT_READER_SEQUENCE_BUFFER_SIZE: usize = 1024;
/// Default cache size.
pub const DEFAULT_CACHE_SIZE: usize = 100;
/// Default in-memory counting setting.
pub const DEFAULT_KMER_COUNT_IN_MEMORY: bool = false;
/// Default free-resources setting.
pub const DEFAULT_FREE_RESOURCES: bool = false;
/// Default temporary IKC removal setting.
pub const DEFAULT_REMOVE_IKC: bool = true;
/// Default reverse-complement k-mer counting setting.
pub const DEFAULT_COUNT_REV_KMER: bool = true;
/// Default multiplier for deriving flank length.
pub const DEFAULT_FLANK_LENGTH_MULTIPLIER: f64 = 3.5;
/// Default haplotype output format.
pub const DEFAULT_HAPLOTYPE_OUTPUT_FORMAT: &str = "sam";

fn run_pipeline(config: &RunConfig) -> Result<(), RunnerError> {
    if config.samples.is_empty() || config.references.is_empty() {
        return Err(RunnerError::PipelineNotImplemented);
    }
    let kmer_util = KmerUtil::new(config.k_size)?;
    let filter_runner = build_filter_runner(config)?;
    let mut reference_reader = ReferenceReader::new(kmer_util.clone());
    reference_reader.set_remove_description(config.remove_reference_sequence_description);
    reference_reader.set_rev_complement_neg_strand(config.reverse_complement_negative_strand);
    reference_reader.set_flank_length(resolved_flank_length(config))?;
    let interval_map = read_interval_map(config)?;
    let references = reference_reader.read(&config.references, interval_map.as_ref())?;

    let mut writer = writer::get_writer(
        Some(&config.output_format),
        Some(config.output_file.clone()),
        references.reference_sequence_array(),
        config.variant_call_by_region,
    )?;
    let mut haplotype_writer = if let Some(output) = &config.haplotype_output_file {
        Some(hapwriter::get_writer(
            Some(&config.haplotype_output_format),
            Some(output.clone()),
            references.reference_sequence_array(),
        )?)
    } else {
        None
    };

    for sample in &config.samples {
        let mut counter = count_map(config, kmer_util.clone())?;
        counter.set(sample.clone())?;
        let read_sequences = sample_read_sequences(sample)?;

        let mut detector = ActiveRegionDetector::new(kmer_util.clone())?;
        detector.set_anchor_both_ends(config.anchor_both_ends);
        detector.set_count_reverse_kmers(config.count_reverse_kmers);
        detector.set_minimum_difference(config.minimum_difference)?;
        detector.set_difference_quantile(config.difference_quantile)?;
        detector.set_peak_scan_length(config.peak_scan_length)?;
        detector.set_scan_limit_factor(config.scan_limit_factor)?;
        detector.set_decay_alpha(config.decay_alpha)?;
        detector.set_decay_minimum(config.decay_minimum)?;
        detector.set_max_repeat_count(config.max_repeat_count)?;
        detector.set_call_ambiguous_regions(config.call_ambiguous_regions);

        writer.set_sample_name(Some(&sample.name))?;
        if let Some(haplotype_writer) = &mut haplotype_writer {
            haplotype_writer.set_sample_name(Some(&sample.name))?;
        }
        for ref_region in references.iter() {
            writer.set_reference_region(ref_region.clone())?;
            if let Some(haplotype_writer) = &mut haplotype_writer {
                haplotype_writer.set_reference_region(ref_region.clone())?;
            }
            let regions = detector.detect_active_regions(ref_region, counter.as_ref())?;
            let emit_context = EmitContext {
                config,
                kmer_util: &kmer_util,
                counter: counter.as_ref(),
                read_sequences: &read_sequences,
                filter_runner: &filter_runner,
            };
            emit_region_variants(
                &emit_context,
                &regions,
                writer.as_mut(),
                &mut haplotype_writer,
            )?;
        }
    }

    writer.flush()?;
    if let Some(haplotype_writer) = &mut haplotype_writer {
        haplotype_writer.flush()?;
    }
    Ok(())
}

struct EmitContext<'a> {
    config: &'a RunConfig,
    kmer_util: &'a KmerUtil,
    counter: &'a dyn CountMap,
    read_sequences: &'a [Vec<u8>],
    filter_runner: &'a VariantFilterRunner,
}

fn emit_region_variants(
    context: &EmitContext<'_>,
    regions: &[ActiveRegion],
    writer: &mut dyn VariantWriter,
    haplotype_writer: &mut Option<Box<dyn HaplotypeWriter>>,
) -> Result<(), RunnerError> {
    for region in regions {
        let haplotypes = read_backed_haplotypes(
            context.config,
            context.kmer_util,
            context.counter,
            context.read_sequences,
            region,
        )?;
        if haplotypes.is_empty() {
            continue;
        }

        let mut caller = VariantCaller::new();
        caller.init(region.clone());
        if context.config.variant_call_by_region {
            caller.set_variant_call_by_region();
        } else {
            caller.set_variant_call_by_reference();
        }
        caller.set_call_ambiguous_variant(context.config.call_ambiguous_variant);

        for haplotype in &haplotypes {
            if let Some(haplotype_writer) = haplotype_writer.as_mut() {
                haplotype_writer.add(Some(haplotype))?;
            }
            caller.add(haplotype.clone())?;
        }

        for variant in caller.variants() {
            if context.filter_runner.filter(Some(&variant)).is_some() {
                writer.write_variant(Some(&variant))?;
            }
        }
    }
    Ok(())
}

fn sample_read_sequences(sample: &InputSample) -> Result<Vec<Vec<u8>>, RunnerError> {
    let mut sequences = Vec::new();
    for source in &sample.sources {
        let records = SequenceReader::new(source.clone())
            .read_all()
            .map_err(|err| RunnerError::Reader(err.to_string()))?;
        sequences.extend(
            records
                .into_iter()
                .filter_map(|record| normalize_read_sequence(&record.sequence)),
        );
    }
    Ok(sequences)
}

fn normalize_read_sequence(sequence: &[u8]) -> Option<Vec<u8>> {
    sequence
        .iter()
        .map(|base| match base {
            b'A' | b'a' => Some(b'A'),
            b'C' | b'c' => Some(b'C'),
            b'G' | b'g' => Some(b'G'),
            b'T' | b't' | b'U' | b'u' => Some(b'T'),
            _ => None,
        })
        .collect()
}

fn read_backed_haplotypes(
    config: &RunConfig,
    kmer_util: &KmerUtil,
    counter: &dyn CountMap,
    read_sequences: &[Vec<u8>],
    region: &ActiveRegion,
) -> Result<Vec<Haplotype>, RunnerError> {
    let mut haplotypes = Vec::new();
    let mut seen = std::collections::BTreeSet::new();

    for read_sequence in read_sequences {
        let Some(consensus) = candidate_consensus(read_sequence, region, kmer_util.k_size()) else {
            continue;
        };
        if !seen.insert(consensus.clone()) {
            continue;
        }

        let mut aligner = KmerAligner::new(
            kmer_util.clone(),
            config.alignment_weight,
            config.haplotype_output_file.is_some(),
        )?;
        aligner.set_max_state(config.max_aligner_state)?;
        aligner.init(region.clone())?;

        let mut valid = true;
        for base in consensus.iter().skip(kmer_util.k_size()) {
            let Some(base) = Base::from_char(char::from(*base)) else {
                valid = false;
                break;
            };
            aligner.add_base(base)?;
        }
        if !valid {
            continue;
        }

        haplotypes.extend(aligner.get_haplotypes(counter, config.count_reverse_kmers)?);
        if haplotypes.len() >= config.max_haplotypes as usize {
            haplotypes.truncate(config.max_haplotypes as usize);
            break;
        }
    }

    Ok(haplotypes)
}

fn candidate_consensus(
    read_sequence: &[u8],
    region: &ActiveRegion,
    k_size: usize,
) -> Option<Vec<u8>> {
    let start = usize::try_from(region.start_index).ok()?;
    let end = usize::try_from(region.end_index).ok()?.checked_add(1)?;
    let ref_sequence = &region.ref_region.sequence;
    let ref_slice = ref_sequence.get(start..end)?;

    if read_sequence.len() == region.ref_region.sequence.len() {
        return read_sequence.get(start..end).map(<[u8]>::to_vec);
    }

    if ref_slice.len() < k_size {
        return None;
    }

    let left_anchor = ref_slice.get(..k_size)?;
    let right_anchor = ref_slice.get(ref_slice.len().checked_sub(k_size)?..)?;
    let left_pos = find_subslice(read_sequence, left_anchor)?;
    let right_search_start = left_pos.checked_add(k_size)?;
    let right_pos =
        right_search_start + find_subslice(read_sequence.get(right_search_start..)?, right_anchor)?;
    read_sequence
        .get(left_pos..right_pos.checked_add(k_size)?)
        .map(<[u8]>::to_vec)
}

fn find_subslice(haystack: &[u8], needle: &[u8]) -> Option<usize> {
    if needle.is_empty() || haystack.len() < needle.len() {
        return None;
    }
    haystack
        .windows(needle.len())
        .position(|window| window == needle)
}

fn count_map(config: &RunConfig, kmer_util: KmerUtil) -> Result<Box<dyn CountMap>, CountMapError> {
    if config.kmer_count_in_memory {
        return Ok(Box::new(MemoryCountMap::new(kmer_util)));
    }

    let minimizer_size = if config.minimizer_size <= 0 {
        kmer_util.k_size()
    } else {
        config.minimizer_size as usize
    };
    IkcCountMap::new(kmer_util, minimizer_size, config.minimizer_mask)
        .map(|map| Box::new(map) as Box<dyn CountMap>)
}

fn resolved_flank_length(config: &RunConfig) -> i32 {
    if config.flank_length >= 0 {
        config.flank_length
    } else {
        (config.k_size as f64 * DEFAULT_FLANK_LENGTH_MULTIPLIER) as i32
    }
}

fn read_interval_map(config: &RunConfig) -> Result<Option<ReferenceIntervalMap>, RunnerError> {
    if config.interval_files.is_empty() {
        return Ok(None);
    }

    let mut container = RegionIntervalContainer::new();
    for path in &config.interval_files {
        for interval in interval::read_path(path, None)? {
            container
                .add(interval)
                .map_err(|error| RunnerError::Interval(error.to_string()))?;
        }
    }
    Ok(Some(container.get_map()))
}

fn build_filter_runner(config: &RunConfig) -> Result<VariantFilterRunner, RunnerError> {
    let mut runner = VariantFilterRunner::default();
    for spec in &config.variant_filter_specs {
        runner.add_filter(Some(VariantFilterKind::get_filter(spec)?))?;
    }
    Ok(runner)
}

#[cfg(test)]
mod tests {
    use super::*;
    use kanalyze::comp::reader::FileSequenceSource;
    use tempfile::tempdir;

    #[test]
    fn defaults_match_java_runner_base() {
        let runner = KestrelRunner::new();
        assert_eq!(runner.k_size(), DEFAULT_KSIZE);
        assert_eq!(runner.output_file(), &StreamableOutput::Stdout);
        assert_eq!(runner.output_format(), "vcf");
        assert_eq!(runner.log_file(), &StreamableOutput::Stderr);
        assert_eq!(runner.log_level(), LogLevel::Warn);
        assert_eq!(runner.minimizer_size(), DEFAULT_MINIMIZER_SIZE);
        assert_eq!(runner.minimizer_mask(), DEFAULT_MINIMIZER_MASK);
        assert!(runner.alignment_weight().is_default());
        assert_eq!(runner.temp_dir_name(), None);
        assert_eq!(runner.min_kmer_count(), DEFAULT_MIN_KMER_COUNT);
        assert!(!runner.kmer_count_in_memory());
        assert!(!runner.free_resources());
        assert!(runner.remove_ikc());
        assert!(runner.anchor_both_ends());
        assert!(runner.count_reverse_kmers());
        assert!(runner.call_ambiguous_regions());
        assert!(runner.call_ambiguous_variant());
        assert_eq!(runner.flank_length(), -1);
        assert!(runner.variant_call_by_reference());
        assert_eq!(runner.max_aligner_state(), KmerAligner::DEFAULT_MAX_STATE);
        assert_eq!(
            runner.max_haplotypes(),
            KmerAlignmentBuilder::DEFAULT_MAX_HAPLOTYPES
        );
        assert!(runner.remove_reference_description());
        assert!(!runner.rev_complement_neg_reference_strand());
        assert_eq!(runner.haplotype_output_file(), None);
        assert_eq!(runner.haplotype_output_format(), "sam");
        assert!(runner.libraries().is_empty());
    }

    #[test]
    fn setter_surface_matches_java_runner_tests() {
        let mut runner = KestrelRunner::new();
        runner.set_k_size(20).unwrap();
        assert_eq!(runner.k_size(), 20);
        assert_eq!(
            runner.set_k_size(1),
            Err(RunnerError::InvalidKSize {
                minimum: MIN_KMER_SIZE,
                actual: 1
            })
        );

        runner.set_output_format("table").unwrap();
        assert_eq!(runner.output_format(), "table");
        assert_eq!(
            runner.set_output_format(" "),
            Err(RunnerError::EmptyOutputFormat)
        );

        runner.set_output_path("/tmp/out.vcf");
        assert!(matches!(
            runner.output_file(),
            StreamableOutput::File { path, .. } if path == &PathBuf::from("/tmp/out.vcf")
        ));
        runner.set_output_file(Some(StreamableOutput::from_fd(1, None)));
        assert_eq!(runner.output_file(), &StreamableOutput::Stdout);

        runner.set_log_level(Some(LogLevel::Debug));
        assert_eq!(runner.log_level(), LogLevel::Debug);
        runner.set_log_level_name("error").unwrap();
        assert_eq!(runner.log_level(), LogLevel::Error);
        assert!(matches!(
            runner.set_log_level_name("not-a-level"),
            Err(RunnerError::LogLevel(LogLevelError::Unknown(_)))
        ));

        runner.set_minimizer_size(10).unwrap();
        assert_eq!(runner.minimizer_size(), 10);
        runner.set_minimizer_mask(0x0f);
        assert_eq!(runner.minimizer_mask(), 0x0f);

        runner.set_alignment_weight(Some("5,5,40,4")).unwrap();
        assert_eq!(runner.alignment_weight().match_score, 5.0);
        runner.set_alignment_weight_match(20.0).unwrap();
        assert_eq!(runner.alignment_weight().match_score, 20.0);
        runner.set_alignment_weight_mismatch(20.0).unwrap();
        assert_eq!(runner.alignment_weight().mismatch, -20.0);
        runner.set_alignment_weight_gap_open(50.0).unwrap();
        assert_eq!(runner.alignment_weight().gap_open, -50.0);
        runner.set_alignment_weight_gap_extend(5.0).unwrap();
        assert_eq!(runner.alignment_weight().gap_extend, -5.0);

        runner.set_temp_dir_name(Some("/tmp/foo"));
        assert_eq!(runner.temp_dir_name(), Some("/tmp/foo"));
        runner.clear_samples();
        runner.clear_reference();
        runner.clear_libraries();
        assert!(runner.libraries().is_empty());
    }

    #[test]
    fn validates_extended_config_values() {
        let mut runner = KestrelRunner::new();
        assert_eq!(
            runner.set_minimizer_size(-1),
            Err(RunnerError::NegativeMinimizerSize(-1))
        );
        assert_eq!(
            runner.set_min_kmer_count(-1),
            Err(RunnerError::NegativeMinKmerCount(-1))
        );
        assert_eq!(
            runner.set_minimum_difference(0),
            Err(RunnerError::InvalidMinimumDifference(0))
        );
        assert_eq!(
            runner.set_difference_quantile(1.0),
            Err(RunnerError::InvalidDifferenceQuantile(1.0))
        );
        assert_eq!(
            runner.set_peak_scan_length(-1),
            Err(RunnerError::NegativePeakScanLength(-1))
        );
        assert_eq!(
            runner.set_scan_limit_factor(-1.0),
            Err(RunnerError::NegativeScanLimitFactor(-1.0))
        );
        assert_eq!(
            runner.set_decay_minimum(1.1),
            Err(RunnerError::InvalidDecayMinimum(1.1))
        );
        assert_eq!(
            runner.set_decay_alpha(0.0),
            Err(RunnerError::InvalidDecayAlpha(0.0))
        );
        assert_eq!(
            runner.set_max_repeat_count(-1),
            Err(RunnerError::NegativeMaxRepeatCount(-1))
        );
        assert_eq!(
            runner.set_flank_length(-2),
            Err(RunnerError::InvalidFlankLength(-2))
        );
        assert_eq!(
            runner.set_max_aligner_state(0),
            Err(RunnerError::InvalidMaxAlignerState(0))
        );
        assert_eq!(
            runner.set_max_haplotypes(0),
            Err(RunnerError::InvalidMaxHaplotypes(0))
        );
        assert_eq!(
            runner.set_haplotype_output_format(" "),
            Err(RunnerError::EmptyHaplotypeOutputFormat)
        );
    }

    #[test]
    fn run_requires_configured_inputs() {
        assert_eq!(
            KestrelRunner::new().run(),
            Err(RunnerError::PipelineNotImplemented)
        );
    }

    #[test]
    fn run_reads_counts_detects_and_writes_header_for_no_variant_fixture() {
        let temp = tempdir().unwrap();
        let ref_path = temp.path().join("ref.fasta");
        let reads_path = temp.path().join("reads.fastq");
        let out_path = temp.path().join("out.vcf");

        std::fs::write(&ref_path, b">chr1\nACGTACGTACGT\n").unwrap();
        std::fs::write(&reads_path, b"@r1\nACGTACGTACGT\n+\nFFFFFFFFFFFF\n").unwrap();

        let mut runner = KestrelRunner::new();
        runner.set_k_size(4).unwrap();
        runner.set_minimizer_size(0).unwrap();
        runner.set_kmer_count_in_memory(true);
        runner.set_output_path(&out_path);
        runner.add_reference(FileSequenceSource::from_path(&ref_path, 1).unwrap());
        runner.add_sample(
            InputSample::new(
                Some("sample"),
                vec![FileSequenceSource::from_path(&reads_path, 1).unwrap()],
            )
            .unwrap(),
        );

        runner.run().unwrap();

        let output = std::fs::read_to_string(out_path).unwrap();
        assert!(output.contains("##fileformat=VCFv4.2"));
        assert!(output.contains("##contig=<ID=chr1,length=12"));
        assert!(output.contains("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample"));
    }

    #[test]
    fn run_detects_read_backed_snp_and_writes_variant() {
        let temp = tempdir().unwrap();
        let ref_path = temp.path().join("ref.fasta");
        let reads_path = temp.path().join("reads.fastq");
        let out_path = temp.path().join("out.vcf");

        std::fs::write(&ref_path, b">chr1\nAAAACCCCGGGGTTTT\n").unwrap();
        std::fs::write(&reads_path, b"@r1\nAAAATCCCGGGGTTTT\n+\nFFFFFFFFFFFFFFFF\n").unwrap();

        let mut runner = KestrelRunner::new();
        runner.set_k_size(4).unwrap();
        runner.set_minimizer_size(0).unwrap();
        runner.set_minimum_difference(1).unwrap();
        runner.set_count_reverse_kmers(false);
        runner.set_kmer_count_in_memory(true);
        runner.set_output_path(&out_path);
        runner.add_reference(FileSequenceSource::from_path(&ref_path, 1).unwrap());
        runner.add_sample(
            InputSample::new(
                Some("sample"),
                vec![FileSequenceSource::from_path(&reads_path, 1).unwrap()],
            )
            .unwrap(),
        );

        runner.run().unwrap();

        let output = std::fs::read_to_string(out_path).unwrap();
        assert!(output.contains("chr1\t5\t.\tC\tT"));
    }

    #[test]
    fn run_applies_bed_intervals_to_reference_regions() {
        let temp = tempdir().unwrap();
        let ref_path = temp.path().join("ref.fasta");
        let reads_path = temp.path().join("reads.fastq");
        let bed_path = temp.path().join("regions.bed");
        let out_path = temp.path().join("out.vcf");

        std::fs::write(&ref_path, b">chr1\nACGTACGTACGT\n").unwrap();
        std::fs::write(&reads_path, b"@r1\nACGTACGTACGT\n+\nFFFFFFFFFFFF\n").unwrap();
        std::fs::write(&bed_path, b"chr1\t2\t10\tinner\n").unwrap();

        let mut runner = KestrelRunner::new();
        runner.set_k_size(4).unwrap();
        runner.set_minimizer_size(0).unwrap();
        runner.set_flank_length(0).unwrap();
        runner.set_kmer_count_in_memory(true);
        runner.set_output_path(&out_path);
        runner.add_reference(FileSequenceSource::from_path(&ref_path, 1).unwrap());
        runner.add_interval_file(&bed_path);
        runner.add_sample(
            InputSample::new(
                Some("sample"),
                vec![FileSequenceSource::from_path(&reads_path, 1).unwrap()],
            )
            .unwrap(),
        );

        runner.run().unwrap();

        let output = std::fs::read_to_string(out_path).unwrap();
        assert!(output.contains("##contig=<ID=chr1,length=12"));
        assert!(output.contains("\tsample"));
    }

    #[test]
    fn run_initializes_variant_filters_from_specs() {
        let temp = tempdir().unwrap();
        let ref_path = temp.path().join("ref.fasta");
        let reads_path = temp.path().join("reads.fastq");
        let out_path = temp.path().join("out.vcf");

        std::fs::write(&ref_path, b">chr1\nACGTACGTACGT\n").unwrap();
        std::fs::write(&reads_path, b"@r1\nACGTACGTACGT\n+\nFFFFFFFFFFFF\n").unwrap();

        let mut runner = KestrelRunner::new();
        runner.set_k_size(4).unwrap();
        runner.set_minimizer_size(0).unwrap();
        runner.set_kmer_count_in_memory(true);
        runner.set_output_path(&out_path);
        runner.add_variant_filter_spec("TYPE:snp");
        runner.add_reference(FileSequenceSource::from_path(&ref_path, 1).unwrap());
        runner.add_sample(
            InputSample::new(
                Some("sample"),
                vec![FileSequenceSource::from_path(&reads_path, 1).unwrap()],
            )
            .unwrap(),
        );

        runner.run().unwrap();
    }

    #[test]
    fn run_initializes_and_flushes_haplotype_writer() {
        let temp = tempdir().unwrap();
        let ref_path = temp.path().join("ref.fasta");
        let reads_path = temp.path().join("reads.fastq");
        let out_path = temp.path().join("out.vcf");
        let hap_path = temp.path().join("hap.sam");

        std::fs::write(&ref_path, b">chr1\nACGTACGTACGT\n").unwrap();
        std::fs::write(&reads_path, b"@r1\nACGTACGTACGT\n+\nFFFFFFFFFFFF\n").unwrap();

        let mut runner = KestrelRunner::new();
        runner.set_k_size(4).unwrap();
        runner.set_minimizer_size(0).unwrap();
        runner.set_kmer_count_in_memory(true);
        runner.set_output_path(&out_path);
        runner.set_haplotype_output_format("sam").unwrap();
        runner.set_haplotype_output_file(Some(StreamableOutput::from_path(&hap_path, None)));
        runner.add_reference(FileSequenceSource::from_path(&ref_path, 1).unwrap());
        runner.add_sample(
            InputSample::new(
                Some("sample"),
                vec![FileSequenceSource::from_path(&reads_path, 1).unwrap()],
            )
            .unwrap(),
        );

        runner.run().unwrap();

        let hap_output = std::fs::read_to_string(hap_path).unwrap();
        assert!(hap_output.contains("@HD\tVN:1.5\tSO:coordinate"));
        assert!(hap_output.contains("@SQ\tSN:chr1\tLN:12"));
        assert!(hap_output.contains("@PG\tID:Kestrel"));
    }

    #[test]
    fn run_rejects_invalid_variant_filter_specs() {
        let temp = tempdir().unwrap();
        let ref_path = temp.path().join("ref.fasta");
        let reads_path = temp.path().join("reads.fastq");
        let out_path = temp.path().join("out.vcf");

        std::fs::write(&ref_path, b">chr1\nACGTACGTACGT\n").unwrap();
        std::fs::write(&reads_path, b"@r1\nACGTACGTACGT\n+\nFFFFFFFFFFFF\n").unwrap();

        let mut runner = KestrelRunner::new();
        runner.set_k_size(4).unwrap();
        runner.set_minimizer_size(0).unwrap();
        runner.set_kmer_count_in_memory(true);
        runner.set_output_path(&out_path);
        runner.add_variant_filter_spec("NOPE:snp");
        runner.add_reference(FileSequenceSource::from_path(&ref_path, 1).unwrap());
        runner.add_sample(
            InputSample::new(
                Some("sample"),
                vec![FileSequenceSource::from_path(&reads_path, 1).unwrap()],
            )
            .unwrap(),
        );

        assert!(matches!(runner.run(), Err(RunnerError::VariantFilter(_))));
    }
}
