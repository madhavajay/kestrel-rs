use std::collections::{BTreeMap, BTreeSet};
use std::io::{self, BufWriter, Write};

use thiserror::Error;

use crate::constants::{FORMAT_TYPE_PATTERN, VERSION};
use crate::io::StreamableOutput;
use crate::refreader::{ReferenceRegion, ReferenceSequence, ReferenceSequenceError};
use crate::variant::{Variant, VariantCall};

/// Default sample name used when no sample is provided.
pub const DEFAULT_SAMPLE_NAME: &str = "UnknownSample";

/// Errors returned by variant writer selection and output.
#[derive(Debug, Error)]
pub enum VariantWriterError {
    /// Writer specification was absent.
    #[error("cannot get variant writer with specification: null")]
    NullSpec,
    /// Writer specification was empty.
    #[error("cannot get variant writer with an empty specification")]
    EmptySpec,
    /// Writer name did not match the expected format pattern.
    #[error("writer name does not match regular expression \"{FORMAT_TYPE_PATTERN}\": {0}")]
    InvalidWriterName(String),
    /// No writer exists for the requested name.
    #[error("cannot find class for variant writer: {0}")]
    UnknownWriter(String),
    /// Writer initialization requires at least one reference sequence.
    #[error("reference sequence array is empty")]
    EmptyReferenceArray,
    /// Reference array contained an absent reference.
    #[error("reference sequence array contains a null reference at index {0}")]
    NullReference(usize),
    /// Variant argument was absent.
    #[error("cannot write variant: null")]
    NullVariant,
    /// Reference region was absent.
    #[error("cannot set reference region: null")]
    NullReferenceRegion,
    /// Sample name was absent.
    #[error("cannot set new sample by name: null")]
    NullSampleName,
    /// Sample name was empty.
    #[error("sample name is empty")]
    EmptySampleName,
    /// Sample name contained whitespace.
    #[error("sample name contains whitespace: \"{0}\"")]
    WhitespaceSampleName(String),
    /// Sample name was already present in a VCF container.
    #[error("sample name already used in this VCF record container: {0}")]
    DuplicateSampleName(String),
    /// A VCF variant was added before declaring the first sample.
    #[error(
        "cannot add variants to VCF container before the first sample (new_sample()) is declared"
    )]
    MissingSample,
    /// I/O error while writing variants.
    #[error(transparent)]
    Io(#[from] io::Error),
    /// Reference sequence error while rendering variants.
    #[error(transparent)]
    Reference(#[from] ReferenceSequenceError),
}

/// Writer interface for variant calls.
pub trait VariantWriter {
    /// Returns the writer name.
    fn name(&self) -> &str;
    /// Returns a short writer description.
    fn description(&self) -> &str;
    /// Initializes the writer.
    fn init(
        &mut self,
        writer_spec: Option<&str>,
        output: Option<StreamableOutput>,
        reference_sequences: Vec<ReferenceSequence>,
        by_region: bool,
    ) -> Result<(), VariantWriterError>;
    /// Sets the active sample name.
    fn set_sample_name(&mut self, sample_name: Option<&str>) -> Result<(), VariantWriterError>;
    /// Sets the active reference region.
    fn set_reference_region(
        &mut self,
        reference_region: ReferenceRegion,
    ) -> Result<(), VariantWriterError>;
    /// Writes one variant.
    fn write_variant(&mut self, variant: Option<&Variant>) -> Result<(), VariantWriterError>;
    /// Flushes pending output.
    fn flush(&mut self) -> Result<(), VariantWriterError>;
}

/// Creates and initializes a variant writer from a writer specification.
pub fn get_writer(
    writer_spec: Option<&str>,
    output: Option<StreamableOutput>,
    reference_sequences: Vec<ReferenceSequence>,
    by_region: bool,
) -> Result<Box<dyn VariantWriter>, VariantWriterError> {
    let writer_spec = writer_spec.ok_or(VariantWriterError::NullSpec)?.trim();
    if writer_spec.is_empty() {
        return Err(VariantWriterError::EmptySpec);
    }

    let (writer_name, writer_args) = writer_spec
        .split_once(':')
        .map_or((writer_spec, ""), |(name, args)| (name.trim(), args.trim()));
    let mut writer = get_writer_class(writer_name)?;
    writer.init(Some(writer_args), output, reference_sequences, by_region)?;
    Ok(writer)
}

/// Creates a variant writer implementation by name.
pub fn get_writer_class(writer_name: &str) -> Result<Box<dyn VariantWriter>, VariantWriterError> {
    let writer_name = writer_name.trim();
    if !matches_format_type(writer_name) {
        return Err(VariantWriterError::InvalidWriterName(
            writer_name.to_owned(),
        ));
    }

    match writer_name.to_ascii_lowercase().as_str() {
        "vcf" => Ok(Box::new(VcfVariantWriter::new())),
        "table" => Ok(Box::new(TableVariantWriter::new())),
        "txt" => Ok(Box::new(TxtVariantWriter::new())),
        _ => Err(VariantWriterError::UnknownWriter(writer_name.to_owned())),
    }
}

/// Returns the description for a variant writer name.
pub fn get_writer_description(writer_name: &str) -> Option<&'static str> {
    match writer_name.trim().to_ascii_lowercase().as_str() {
        "vcf" => Some("Writes variant call format (VCF) files."),
        "table" => Some("Tab-delimited table of variant information"),
        "txt" => Some("Plain-text list of variants"),
        _ => None,
    }
}

/// Lists supported variant writer names.
pub fn list_writers() -> Vec<&'static str> {
    vec!["table", "txt", "vcf"]
}

fn matches_format_type(name: &str) -> bool {
    let mut chars = name.chars();
    let Some(first) = chars.next() else {
        return false;
    };
    first.is_ascii_alphabetic()
        && chars.all(|ch| ch.is_ascii_alphanumeric() || matches!(ch, '_' | '-'))
}

struct WriterBase {
    name: String,
    output: StreamableOutput,
    sample_name: String,
    reference_region: Option<ReferenceRegion>,
    reference_sequences: Vec<ReferenceSequence>,
    by_region: bool,
}

impl WriterBase {
    fn new(name: &str) -> Self {
        let name = normalize_writer_name(name);
        Self {
            name,
            output: StreamableOutput::stdout(),
            sample_name: DEFAULT_SAMPLE_NAME.to_owned(),
            reference_region: None,
            reference_sequences: Vec::new(),
            by_region: false,
        }
    }

    fn init(
        &mut self,
        output: Option<StreamableOutput>,
        reference_sequences: Vec<ReferenceSequence>,
        by_region: bool,
    ) -> Result<(), VariantWriterError> {
        if reference_sequences.is_empty() {
            return Err(VariantWriterError::EmptyReferenceArray);
        }

        self.output = output.unwrap_or_else(StreamableOutput::stdout);
        self.reference_sequences = reference_sequences;
        self.by_region = by_region;
        self.reference_region = None;
        Ok(())
    }

    fn set_sample_name(&mut self, sample_name: Option<&str>) -> bool {
        let sample_name = sample_name
            .map(str::trim)
            .filter(|name| !name.is_empty())
            .unwrap_or(DEFAULT_SAMPLE_NAME);
        if sample_name == self.sample_name {
            return false;
        }
        self.sample_name = sample_name.to_owned();
        true
    }

    fn set_reference_region(&mut self, reference_region: ReferenceRegion) -> bool {
        if self.reference_region.as_ref() == Some(&reference_region) {
            return false;
        }
        self.reference_region = Some(reference_region);
        true
    }
}

fn normalize_writer_name(name: &str) -> String {
    let normalized = name
        .trim()
        .chars()
        .map(|ch| if ch.is_whitespace() { '_' } else { ch })
        .collect::<String>();
    if normalized.is_empty() {
        "UnknownVariantWriter".to_owned()
    } else {
        normalized
    }
}

pub(crate) fn open_output(output: &StreamableOutput) -> Result<Box<dyn Write>, io::Error> {
    match output {
        StreamableOutput::Stdout => Ok(Box::new(io::stdout())),
        StreamableOutput::Stderr => Ok(Box::new(io::stderr())),
        StreamableOutput::File { .. } => Ok(Box::new(BufWriter::new(
            output
                .create_file()?
                .expect("file outputs always create a file handle"),
        ))),
        StreamableOutput::Fd { .. } => Err(io::Error::new(
            io::ErrorKind::Unsupported,
            "raw file descriptor outputs are not supported yet",
        )),
    }
}

/// Variant writer that emits a tab-delimited table.
pub struct TableVariantWriter {
    base: WriterBase,
    out: Option<Box<dyn Write>>,
    reference_name: String,
    region_name: String,
}

impl TableVariantWriter {
    const DEFAULT_REFERENCE_NAME: &'static str = "UnknownRef";
    const DEFAULT_REGION_NAME: &'static str = "UnknownRegion";

    /// Creates a table variant writer.
    #[must_use]
    pub fn new() -> Self {
        Self {
            base: WriterBase::new("table"),
            out: None,
            reference_name: Self::DEFAULT_REFERENCE_NAME.to_owned(),
            region_name: Self::DEFAULT_REGION_NAME.to_owned(),
        }
    }

    fn out(&mut self) -> &mut dyn Write {
        self.out
            .as_deref_mut()
            .expect("variant writer has been initialized")
    }
}

impl Default for TableVariantWriter {
    fn default() -> Self {
        Self::new()
    }
}

impl VariantWriter for TableVariantWriter {
    fn name(&self) -> &str {
        &self.base.name
    }

    fn description(&self) -> &str {
        "Tab-delimited table of variant information"
    }

    fn init(
        &mut self,
        _writer_spec: Option<&str>,
        output: Option<StreamableOutput>,
        reference_sequences: Vec<ReferenceSequence>,
        by_region: bool,
    ) -> Result<(), VariantWriterError> {
        self.base.init(output, reference_sequences, by_region)?;
        self.reference_name = Self::DEFAULT_REFERENCE_NAME.to_owned();
        self.region_name = Self::DEFAULT_REGION_NAME.to_owned();
        self.out = Some(open_output(&self.base.output)?);
        writeln!(
            self.out(),
            "sample\treference\tregion\tlocus\tref\talt\tvar_depth\tregion_depth"
        )?;
        Ok(())
    }

    fn set_sample_name(&mut self, sample_name: Option<&str>) -> Result<(), VariantWriterError> {
        self.base.set_sample_name(sample_name);
        Ok(())
    }

    fn set_reference_region(
        &mut self,
        reference_region: ReferenceRegion,
    ) -> Result<(), VariantWriterError> {
        if self.base.set_reference_region(reference_region) {
            let reference_region = self.base.reference_region.as_ref().unwrap();
            self.reference_name = reference_region.reference_sequence.name.clone();
            self.region_name = reference_region.name.clone();
        }
        Ok(())
    }

    fn write_variant(&mut self, variant: Option<&Variant>) -> Result<(), VariantWriterError> {
        let variant = variant.ok_or(VariantWriterError::NullVariant)?;
        let data = variant.data();
        let sample_name = self.base.sample_name.clone();
        let reference_name = self.reference_name.clone();
        let region_name = self.region_name.clone();
        writeln!(
            self.out(),
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            sample_name,
            reference_name,
            region_name,
            data.start,
            data.ref_allele,
            data.alt_allele,
            data.variant_depth,
            data.locus_depth
        )?;
        Ok(())
    }

    fn flush(&mut self) -> Result<(), VariantWriterError> {
        if let Some(out) = &mut self.out {
            out.flush()?;
        }
        Ok(())
    }
}

/// Variant writer that emits a plain text report.
pub struct TxtVariantWriter {
    base: WriterBase,
    out: Option<Box<dyn Write>>,
    first_sample: bool,
    last_region: Option<ReferenceRegion>,
}

impl TxtVariantWriter {
    /// Creates a text variant writer.
    #[must_use]
    pub fn new() -> Self {
        Self {
            base: WriterBase::new("txt"),
            out: None,
            first_sample: true,
            last_region: None,
        }
    }

    fn out(&mut self) -> &mut dyn Write {
        self.out
            .as_deref_mut()
            .expect("variant writer has been initialized")
    }
}

impl Default for TxtVariantWriter {
    fn default() -> Self {
        Self::new()
    }
}

impl VariantWriter for TxtVariantWriter {
    fn name(&self) -> &str {
        &self.base.name
    }

    fn description(&self) -> &str {
        "Plain-text list of variants"
    }

    fn init(
        &mut self,
        _writer_spec: Option<&str>,
        output: Option<StreamableOutput>,
        reference_sequences: Vec<ReferenceSequence>,
        by_region: bool,
    ) -> Result<(), VariantWriterError> {
        self.base.init(output, reference_sequences, by_region)?;
        self.out = Some(open_output(&self.base.output)?);
        self.first_sample = true;
        self.last_region = None;
        Ok(())
    }

    fn set_sample_name(&mut self, sample_name: Option<&str>) -> Result<(), VariantWriterError> {
        if self.base.set_sample_name(sample_name) {
            if self.first_sample {
                self.first_sample = false;
            } else {
                writeln!(self.out())?;
                writeln!(self.out())?;
            }
            let sample_name = self.base.sample_name.clone();
            writeln!(self.out(), "** Sample: {sample_name}")?;
        }
        Ok(())
    }

    fn set_reference_region(
        &mut self,
        reference_region: ReferenceRegion,
    ) -> Result<(), VariantWriterError> {
        if self.base.set_reference_region(reference_region) {
            let reference_region = self.base.reference_region.as_ref().unwrap().clone();
            if self.base.by_region {
                writeln!(
                    self.out(),
                    "\n* Reference region: {}",
                    reference_region.name
                )?;
            } else if self
                .last_region
                .as_ref()
                .is_none_or(|last| last.reference_sequence != reference_region.reference_sequence)
            {
                writeln!(
                    self.out(),
                    "\n* Reference: {}",
                    reference_region.reference_sequence.name
                )?;
            }
            self.last_region = Some(reference_region);
        }
        Ok(())
    }

    fn write_variant(&mut self, variant: Option<&Variant>) -> Result<(), VariantWriterError> {
        let variant = variant.ok_or(VariantWriterError::NullVariant)?;
        let data = variant.data();
        writeln!(
            self.out(),
            "{} ({}/{})",
            variant.hgvs(),
            data.variant_depth,
            data.locus_depth
        )?;
        Ok(())
    }

    fn flush(&mut self) -> Result<(), VariantWriterError> {
        if let Some(out) = &mut self.out {
            out.flush()?;
        }
        self.last_region = None;
        Ok(())
    }
}

/// Variant writer that emits VCF records.
pub struct VcfVariantWriter {
    base: WriterBase,
    out: Option<Box<dyn Write>>,
    record_container: VcfRecordContainer,
}

impl VcfVariantWriter {
    /// Creates a VCF variant writer.
    #[must_use]
    pub fn new() -> Self {
        Self {
            base: WriterBase::new("VCF"),
            out: None,
            record_container: VcfRecordContainer::new(),
        }
    }

    fn out(&mut self) -> &mut dyn Write {
        self.out
            .as_deref_mut()
            .expect("variant writer has been initialized")
    }
}

impl Default for VcfVariantWriter {
    fn default() -> Self {
        Self::new()
    }
}

impl VariantWriter for VcfVariantWriter {
    fn name(&self) -> &str {
        &self.base.name
    }

    fn description(&self) -> &str {
        "Writes variant call format (VCF) files."
    }

    fn init(
        &mut self,
        _writer_spec: Option<&str>,
        output: Option<StreamableOutput>,
        reference_sequences: Vec<ReferenceSequence>,
        by_region: bool,
    ) -> Result<(), VariantWriterError> {
        self.base.init(output, reference_sequences, by_region)?;
        self.record_container.clear();
        self.out = Some(open_output(&self.base.output)?);
        writeln!(self.out(), "##fileformat=VCF4.2")?;
        writeln!(self.out(), "##source=Kestrel{VERSION}")?;
        Ok(())
    }

    fn set_sample_name(&mut self, sample_name: Option<&str>) -> Result<(), VariantWriterError> {
        if self.base.set_sample_name(sample_name) {
            self.record_container.new_sample(&self.base.sample_name)?;
        }
        Ok(())
    }

    fn set_reference_region(
        &mut self,
        reference_region: ReferenceRegion,
    ) -> Result<(), VariantWriterError> {
        self.base.set_reference_region(reference_region);
        Ok(())
    }

    fn write_variant(&mut self, variant: Option<&Variant>) -> Result<(), VariantWriterError> {
        let variant = variant.ok_or(VariantWriterError::NullVariant)?;
        self.record_container.add_variant(variant.clone())
    }

    fn flush(&mut self) -> Result<(), VariantWriterError> {
        let reference_sequences = self.base.reference_sequences.clone();
        for ref_sequence in &reference_sequences {
            let digest = ref_sequence
                .digest
                .as_ref()
                .map_or_else(|| ".".to_owned(), ToString::to_string);
            writeln!(
                self.out(),
                "##contig=<ID={},length={},md5={}>",
                ref_sequence.name,
                ref_sequence.size,
                digest
            )?;
        }

        for format_line in VcfRecordContainer::FORMAT_HEADERS {
            writeln!(self.out(), "{format_line}")?;
        }

        write!(
            self.out(),
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"
        )?;
        for sample_name in self.record_container.sample_names().to_vec() {
            write!(self.out(), "\t{sample_name}")?;
        }
        writeln!(self.out())?;

        for line in self.record_container.vcf_lines()? {
            writeln!(self.out(), "{line}")?;
        }

        self.out().flush()?;
        Ok(())
    }
}

/// Container that groups variants into VCF records across samples.
#[derive(Clone, Debug)]
pub struct VcfRecordContainer {
    records: BTreeMap<VcfRecordKey, VcfVariantRecord>,
    sample_name: Option<String>,
    sample_names: Vec<String>,
    sample_name_set: BTreeSet<String>,
}

impl VcfRecordContainer {
    const FORMAT_HEADERS: [&'static str; 3] = [
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
        "##FORMAT=<ID=GDP,Number=A,Type=Integer,Description=\"Estimated depth of all haplotypes supporting the alternate variant\">",
        "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Estimated depth of all haplotypes in the variant active region\">",
    ];

    /// Creates an empty VCF record container.
    #[must_use]
    pub fn new() -> Self {
        Self {
            records: BTreeMap::new(),
            sample_name: None,
            sample_names: Vec::new(),
            sample_name_set: BTreeSet::new(),
        }
    }

    /// Starts recording variants for a new sample.
    pub fn new_sample(&mut self, sample_name: &str) -> Result<(), VariantWriterError> {
        let sample_name = sample_name.trim();
        if sample_name.is_empty() {
            return Err(VariantWriterError::EmptySampleName);
        }
        if sample_name.chars().any(char::is_whitespace) {
            return Err(VariantWriterError::WhitespaceSampleName(
                sample_name.to_owned(),
            ));
        }
        if self.sample_name_set.contains(sample_name) {
            return Err(VariantWriterError::DuplicateSampleName(
                sample_name.to_owned(),
            ));
        }

        self.sample_name = Some(sample_name.to_owned());
        self.sample_names.push(sample_name.to_owned());
        self.sample_name_set.insert(sample_name.to_owned());
        Ok(())
    }

    /// Clears all samples and records.
    pub fn clear(&mut self) {
        self.records.clear();
        self.sample_name = None;
        self.sample_names.clear();
        self.sample_name_set.clear();
    }

    /// Adds one variant for the current sample.
    pub fn add_variant(&mut self, variant: Variant) -> Result<(), VariantWriterError> {
        let sample_name = self
            .sample_name
            .clone()
            .ok_or(VariantWriterError::MissingSample)?;
        let key = VcfRecordKey::from_variant(&variant)?;
        self.records
            .entry(key.clone())
            .and_modify(|record| record.add_variant(sample_name.clone(), variant.clone()))
            .or_insert_with(|| VcfVariantRecord::new(key, sample_name, variant));
        Ok(())
    }

    /// Returns sample names in output order.
    #[must_use]
    pub fn sample_names(&self) -> &[String] {
        &self.sample_names
    }

    /// Returns VCF FORMAT header lines.
    #[must_use]
    pub fn format_headers(&self) -> [&'static str; 3] {
        Self::FORMAT_HEADERS
    }

    /// Renders grouped VCF records as lines.
    pub fn vcf_lines(&self) -> Result<Vec<String>, VariantWriterError> {
        self.records
            .values()
            .map(|record| record.to_vcf_line(&self.sample_names))
            .collect()
    }
}

impl Default for VcfRecordContainer {
    fn default() -> Self {
        Self::new()
    }
}

#[derive(Clone, Debug, Eq, Ord, PartialEq, PartialOrd)]
struct VcfRecordKey {
    chrom: String,
    pos: i32,
    ref_allele: String,
    alt_allele: String,
}

impl VcfRecordKey {
    fn from_variant(variant: &Variant) -> Result<Self, VariantWriterError> {
        Ok(Self {
            chrom: variant
                .data()
                .active_region
                .ref_region
                .reference_sequence
                .name
                .clone(),
            pos: variant.vcf_pos(),
            ref_allele: variant.vcf_ref()?,
            alt_allele: variant.vcf_alt()?,
        })
    }
}

#[derive(Clone, Debug)]
struct VcfVariantRecord {
    key: VcfRecordKey,
    sample_infos: Vec<VcfVariantSampleInfo>,
}

impl VcfVariantRecord {
    fn new(key: VcfRecordKey, sample_name: String, variant: Variant) -> Self {
        Self {
            key,
            sample_infos: vec![VcfVariantSampleInfo {
                sample_name,
                variant,
            }],
        }
    }

    fn add_variant(&mut self, sample_name: String, variant: Variant) {
        self.sample_infos.push(VcfVariantSampleInfo {
            sample_name,
            variant,
        });
    }

    fn to_vcf_line(&self, sample_names: &[String]) -> Result<String, VariantWriterError> {
        let mut record = format!(
            "{}\t{}\t.\t{}\t{}\t.\t.\t.\tGT:GDP:DP",
            self.key.chrom, self.key.pos, self.key.ref_allele, self.key.alt_allele
        );

        let first_sample_info = self
            .sample_infos
            .first()
            .expect("VCF records are only created with a sample info");
        for sample_name in sample_names {
            if sample_name == &first_sample_info.sample_name {
                let data = first_sample_info.variant.data();
                record.push_str(&format!("\t1:{}:{}", data.variant_depth, data.locus_depth));
            } else {
                record.push_str("\t0:.:.");
            }
        }

        Ok(record)
    }
}

#[derive(Clone, Debug)]
struct VcfVariantSampleInfo {
    sample_name: String,
    variant: Variant,
}

#[cfg(test)]
mod tests {
    use std::fs;

    use kanalyze::util::KmerUtil;
    use tempfile::NamedTempFile;

    use crate::activeregion::{ActiveRegion, Haplotype, RegionStats};
    use crate::align::AlignNode;
    use crate::refreader::{ReferenceRegion, ReferenceSequence};
    use crate::util::digest::Digest;
    use crate::variant::{Variant, VariantSnp};

    use super::*;

    #[test]
    fn writer_factory_resolves_known_writers() {
        assert_eq!(get_writer_class("VCF").unwrap().name(), "VCF");
        assert_eq!(get_writer_class("TABLE").unwrap().name(), "table");
        assert_eq!(get_writer_class("TXT").unwrap().name(), "txt");
        assert!(matches!(
            get_writer_class("DOES_NOT_EXIST"),
            Err(VariantWriterError::UnknownWriter(_))
        ));
        assert!(get_writer_description("VCF").unwrap().contains("VCF"));
        assert_eq!(get_writer_description("NOPE"), None);
        assert_eq!(list_writers(), vec!["table", "txt", "vcf"]);
    }

    #[test]
    fn table_writer_writes_header_and_variant() {
        let out = NamedTempFile::new().unwrap();
        let mut writer = get_writer(
            Some("TABLE"),
            Some(StreamableOutput::from_path(out.path(), None)),
            vec![make_ref()],
            false,
        )
        .unwrap();
        writer
            .write_variant(Some(&snp_variant(5, "C", "A")))
            .unwrap();
        writer.flush().unwrap();

        let content = fs::read_to_string(out.path()).unwrap();
        let lines = content.lines().collect::<Vec<_>>();
        assert!(lines[0].contains("sample"));
        assert!(lines[0].contains("ref"));
        assert!(lines[1].contains("\tC\t"));
        assert!(lines[1].contains("\tA\t"));
        assert!(lines[1].contains("\t7\t"));
    }

    #[test]
    fn txt_writer_writes_variant_with_depth_info() {
        let out = NamedTempFile::new().unwrap();
        let mut writer = get_writer(
            Some("TXT"),
            Some(StreamableOutput::from_path(out.path(), None)),
            vec![make_ref()],
            false,
        )
        .unwrap();
        writer
            .write_variant(Some(&snp_variant(5, "C", "A")))
            .unwrap();
        writer.flush().unwrap();

        let content = fs::read_to_string(out.path()).unwrap();
        assert!(content.contains("7/12"));
        assert!(content.contains("5C>A"));
    }

    #[test]
    fn vcf_record_container_validates_samples_and_formats_line() {
        let mut container = VcfRecordContainer::new();
        assert_eq!(container.format_headers().len(), 3);
        assert!(container.sample_names().is_empty());
        assert!(matches!(
            container.add_variant(snp_variant(5, "C", "A")),
            Err(VariantWriterError::MissingSample)
        ));
        assert!(matches!(
            container.new_sample("space here"),
            Err(VariantWriterError::WhitespaceSampleName(_))
        ));

        container.new_sample("s1").unwrap();
        assert!(matches!(
            container.new_sample("s1"),
            Err(VariantWriterError::DuplicateSampleName(_))
        ));
        container.add_variant(snp_variant(5, "C", "A")).unwrap();
        let lines = container.vcf_lines().unwrap();
        assert_eq!(lines.len(), 1);
        assert!(lines[0].contains("chr1"));
        assert!(lines[0].contains("\t5\t"));
        assert!(lines[0].contains("\tC\tA\t"));

        container.clear();
        assert!(container.sample_names().is_empty());
        assert!(container.vcf_lines().unwrap().is_empty());
    }

    #[test]
    fn vcf_writer_flushes_header_samples_and_records() {
        let out = NamedTempFile::new().unwrap();
        let mut writer = get_writer(
            Some("VCF"),
            Some(StreamableOutput::from_path(out.path(), None)),
            vec![make_ref()],
            false,
        )
        .unwrap();
        writer.set_sample_name(Some("s1")).unwrap();
        writer
            .write_variant(Some(&snp_variant(5, "C", "A")))
            .unwrap();
        writer.flush().unwrap();

        let content = fs::read_to_string(out.path()).unwrap();
        assert!(content.contains("##fileformat=VCF4.2"));
        assert!(content.contains("##source=Kestrel1.0.2"));
        assert!(content.contains("##contig=<ID=chr1,length=16,md5="));
        assert!(content.contains("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1"));
        assert!(content.contains("chr1\t5\t.\tC\tA\t.\t.\t.\tGT:GDP:DP\t1:7:12"));
    }

    #[test]
    fn base_methods_reject_invalid_initialization_and_null_variant() {
        let mut writer = TableVariantWriter::new();
        assert!(matches!(
            writer.init(None, Some(StreamableOutput::stdout()), Vec::new(), false),
            Err(VariantWriterError::EmptyReferenceArray)
        ));
        let out = NamedTempFile::new().unwrap();
        writer
            .init(
                None,
                Some(StreamableOutput::from_path(out.path(), None)),
                vec![make_ref()],
                false,
            )
            .unwrap();
        assert!(matches!(
            writer.write_variant(None),
            Err(VariantWriterError::NullVariant)
        ));
        writer.set_sample_name(None).unwrap();
        writer.set_sample_name(Some("   ")).unwrap();
        writer
            .set_reference_region(reference_region("AAAACCCCGGGGTTTT"))
            .unwrap();
    }

    fn snp_variant(start: i32, ref_allele: &str, alt_allele: &str) -> Variant {
        let haplotype = haplotype("AAAACCCCGGGGTTTT");
        Variant::Snp(
            VariantSnp::new(start, 7, 12, ref_allele, alt_allele, &[haplotype], 1, true).unwrap(),
        )
    }

    fn haplotype(sequence: &str) -> Haplotype {
        let ref_region = reference_region(sequence);
        let count = vec![10; sequence.len() - 4 + 1];
        let kmer_util = KmerUtil::new(4).unwrap();
        let active_region = ActiveRegion::new(ref_region, 0, 12, &count, &kmer_util).unwrap();
        let stats = RegionStats::from_counts(
            &count,
            active_region.start_kmer_index,
            active_region.end_kmer_index,
        )
        .unwrap();
        Haplotype::new(
            sequence.as_bytes().to_vec(),
            active_region,
            vec![AlignNode::new(AlignNode::MATCH, sequence.len() as i32, None).unwrap()],
            100.0,
            None,
            stats,
        )
        .unwrap()
    }

    fn make_ref() -> ReferenceSequence {
        ReferenceSequence::new("chr1", 16, Some(digest()), Some("test")).unwrap()
    }

    fn reference_region(sequence: &str) -> ReferenceRegion {
        ReferenceRegion::whole(make_ref(), sequence.as_bytes(), 0).unwrap()
    }

    fn digest() -> Digest {
        Digest::new([0_u8; 16], "MD5").unwrap()
    }
}
