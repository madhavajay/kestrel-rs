use std::io::{self, Write};

use thiserror::Error;

use crate::activeregion::Haplotype;
use crate::constants::{FORMAT_TYPE_PATTERN, VERSION};
use crate::io::StreamableOutput;
use crate::refreader::{ReferenceRegion, ReferenceSequence};
use crate::writer::open_output;

/// Default sample name used when no sample name is provided.
pub const DEFAULT_SAMPLE_NAME: &str = "UnknownSample";

/// Errors returned by haplotype writer selection and output.
#[derive(Debug, Error)]
pub enum HaplotypeWriterError {
    /// Writer specification was absent.
    #[error("cannot get haplotype writer with specification: null")]
    NullSpec,
    /// Writer specification was empty.
    #[error("cannot get haplotype writer with an empty specification")]
    EmptySpec,
    /// Writer name did not match the expected format pattern.
    #[error(
        "haplotype writer name does not match regular expression \"{FORMAT_TYPE_PATTERN}\": {0}"
    )]
    InvalidWriterName(String),
    /// No writer exists for the requested name.
    #[error("cannot find class for haplotype writer: {0}")]
    UnknownWriter(String),
    /// Writer initialization requires at least one reference sequence.
    #[error("reference sequence array is empty")]
    EmptyReferenceArray,
    /// Reference region was absent.
    #[error("cannot set reference region: null")]
    NullReferenceRegion,
    /// Haplotype argument was absent.
    #[error("cannot add haplotype: null")]
    NullHaplotype,
    /// I/O error while writing haplotypes.
    #[error(transparent)]
    Io(#[from] io::Error),
}

/// Writer interface for resolved haplotypes.
pub trait HaplotypeWriter {
    /// Returns the writer name.
    fn name(&self) -> &str;
    /// Returns a short writer description.
    fn description(&self) -> &str;
    /// Initializes the writer with a spec, output target, and references.
    fn init(
        &mut self,
        writer_spec: Option<&str>,
        output: Option<StreamableOutput>,
        reference_sequences: Vec<ReferenceSequence>,
    ) -> Result<(), HaplotypeWriterError>;
    /// Sets the active sample name.
    fn set_sample_name(&mut self, sample_name: Option<&str>) -> Result<(), HaplotypeWriterError>;
    /// Sets the active reference region.
    fn set_reference_region(
        &mut self,
        reference_region: ReferenceRegion,
    ) -> Result<(), HaplotypeWriterError>;
    /// Adds a haplotype to the writer.
    fn add(&mut self, haplotype: Option<&Haplotype>) -> Result<(), HaplotypeWriterError>;
    /// Flushes pending output.
    fn flush(&mut self) -> Result<(), HaplotypeWriterError>;
}

/// Creates and initializes a haplotype writer from a writer specification.
pub fn get_writer(
    writer_spec: Option<&str>,
    output: Option<StreamableOutput>,
    reference_sequences: Vec<ReferenceSequence>,
) -> Result<Box<dyn HaplotypeWriter>, HaplotypeWriterError> {
    let writer_spec = writer_spec.ok_or(HaplotypeWriterError::NullSpec)?.trim();
    if writer_spec.is_empty() {
        return Err(HaplotypeWriterError::EmptySpec);
    }

    let (writer_name, writer_args) = writer_spec
        .split_once(':')
        .map_or((writer_spec, ""), |(name, args)| (name.trim(), args.trim()));
    let mut writer = get_writer_class(writer_name)?;
    writer.init(Some(writer_args), output, reference_sequences)?;
    Ok(writer)
}

/// Creates a haplotype writer implementation by name.
pub fn get_writer_class(
    writer_name: &str,
) -> Result<Box<dyn HaplotypeWriter>, HaplotypeWriterError> {
    let writer_name = writer_name.trim();
    if !matches_format_type(writer_name) {
        return Err(HaplotypeWriterError::InvalidWriterName(
            writer_name.to_owned(),
        ));
    }

    match writer_name.to_ascii_lowercase().as_str() {
        "sam" => Ok(Box::new(SamHaplotypeWriter::new())),
        "null" => Ok(Box::new(NullHaplotypeWriter::new())),
        _ => Err(HaplotypeWriterError::UnknownWriter(writer_name.to_owned())),
    }
}

/// Returns the description for a haplotype writer name.
pub fn get_writer_description(writer_name: &str) -> Option<&'static str> {
    match writer_name.trim().to_ascii_lowercase().as_str() {
        "sam" => Some("Write resolved haplotype sequences in SAM format."),
        "null" => Some("Discards haplotypes without writing them."),
        _ => None,
    }
}

/// Lists supported haplotype writer names.
pub fn list_writers() -> Vec<&'static str> {
    vec!["null", "sam"]
}

fn matches_format_type(name: &str) -> bool {
    let mut chars = name.chars();
    let Some(first) = chars.next() else {
        return false;
    };
    first.is_ascii_alphabetic()
        && chars.all(|ch| ch.is_ascii_alphanumeric() || matches!(ch, '_' | '-'))
}

#[derive(Clone, Debug)]
struct HaplotypeWriterBase {
    name: String,
    output: StreamableOutput,
    sample_name: String,
    reference_region: Option<ReferenceRegion>,
    reference_sequences: Vec<ReferenceSequence>,
}

impl HaplotypeWriterBase {
    fn new(name: &str) -> Self {
        Self {
            name: normalize_writer_name(name),
            output: StreamableOutput::stdout(),
            sample_name: DEFAULT_SAMPLE_NAME.to_owned(),
            reference_region: None,
            reference_sequences: Vec::new(),
        }
    }

    fn init(
        &mut self,
        output: Option<StreamableOutput>,
        mut reference_sequences: Vec<ReferenceSequence>,
    ) -> Result<(), HaplotypeWriterError> {
        if reference_sequences.is_empty() {
            return Err(HaplotypeWriterError::EmptyReferenceArray);
        }
        reference_sequences.sort();
        self.output = output.unwrap_or_else(StreamableOutput::stdout);
        self.reference_sequences = reference_sequences;
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
        "UnknownHaplotypeWriter".to_owned()
    } else {
        normalized
    }
}

/// Haplotype writer that discards all records.
pub struct NullHaplotypeWriter {
    base: HaplotypeWriterBase,
}

impl NullHaplotypeWriter {
    /// Creates a null haplotype writer.
    #[must_use]
    pub fn new() -> Self {
        Self {
            base: HaplotypeWriterBase::new("EmptyHaplotypeWriter"),
        }
    }
}

impl Default for NullHaplotypeWriter {
    fn default() -> Self {
        Self::new()
    }
}

impl HaplotypeWriter for NullHaplotypeWriter {
    fn name(&self) -> &str {
        &self.base.name
    }

    fn description(&self) -> &str {
        "Discards haplotypes without writing them."
    }

    fn init(
        &mut self,
        _writer_spec: Option<&str>,
        output: Option<StreamableOutput>,
        reference_sequences: Vec<ReferenceSequence>,
    ) -> Result<(), HaplotypeWriterError> {
        self.base.init(output, reference_sequences)
    }

    fn set_sample_name(&mut self, sample_name: Option<&str>) -> Result<(), HaplotypeWriterError> {
        self.base.set_sample_name(sample_name);
        Ok(())
    }

    fn set_reference_region(
        &mut self,
        reference_region: ReferenceRegion,
    ) -> Result<(), HaplotypeWriterError> {
        self.base.set_reference_region(reference_region);
        Ok(())
    }

    fn add(&mut self, _haplotype: Option<&Haplotype>) -> Result<(), HaplotypeWriterError> {
        Ok(())
    }

    fn flush(&mut self) -> Result<(), HaplotypeWriterError> {
        Ok(())
    }
}

/// Haplotype writer that emits SAM records.
pub struct SamHaplotypeWriter {
    base: HaplotypeWriterBase,
    out: Option<Box<dyn Write>>,
    records: Vec<SamRecord>,
}

impl SamHaplotypeWriter {
    /// Creates a SAM haplotype writer.
    #[must_use]
    pub fn new() -> Self {
        Self {
            base: HaplotypeWriterBase::new("SamHaplotypeWriter"),
            out: None,
            records: Vec::new(),
        }
    }

    fn out(&mut self) -> &mut dyn Write {
        self.out
            .as_deref_mut()
            .expect("haplotype writer has been initialized")
    }
}

impl Default for SamHaplotypeWriter {
    fn default() -> Self {
        Self::new()
    }
}

impl HaplotypeWriter for SamHaplotypeWriter {
    fn name(&self) -> &str {
        &self.base.name
    }

    fn description(&self) -> &str {
        "Write resolved haplotype sequences in SAM format."
    }

    fn init(
        &mut self,
        _writer_spec: Option<&str>,
        output: Option<StreamableOutput>,
        reference_sequences: Vec<ReferenceSequence>,
    ) -> Result<(), HaplotypeWriterError> {
        self.base.init(output, reference_sequences)?;
        self.records.clear();
        self.out = Some(open_output(&self.base.output)?);
        Ok(())
    }

    fn set_sample_name(&mut self, sample_name: Option<&str>) -> Result<(), HaplotypeWriterError> {
        self.base.set_sample_name(sample_name);
        Ok(())
    }

    fn set_reference_region(
        &mut self,
        reference_region: ReferenceRegion,
    ) -> Result<(), HaplotypeWriterError> {
        self.base.set_reference_region(reference_region);
        Ok(())
    }

    fn add(&mut self, haplotype: Option<&Haplotype>) -> Result<(), HaplotypeWriterError> {
        let haplotype = haplotype.ok_or(HaplotypeWriterError::NullHaplotype)?;
        self.records.push(SamRecord::new(
            self.base.sample_name.clone(),
            haplotype.clone(),
        ));
        Ok(())
    }

    fn flush(&mut self) -> Result<(), HaplotypeWriterError> {
        let reference_sequences = self.base.reference_sequences.clone();
        writeln!(self.out(), "@HD\tVN:1.5\tSO:coordinate")?;
        for reference in &reference_sequences {
            let digest = reference
                .digest
                .as_ref()
                .map_or_else(|| ".".to_owned(), ToString::to_string);
            writeln!(
                self.out(),
                "@SQ\tSN:{}\tLN:{}\tM5:{}",
                reference.name,
                reference.size,
                digest
            )?;
        }
        writeln!(self.out(), "@PG\tID:Kestrel\tVN:{VERSION}")?;

        let mut records = self.records.clone();
        records.sort();
        for record in &records {
            let sequence = String::from_utf8_lossy(&record.haplotype.sequence);
            writeln!(
                self.out(),
                "{}-{}-{}-{}\t0\t{}\t{}\t255\t{}\t*\t0\t0\t{}\t*\tXD:i:{}\tXN:Z:{}\tXL:i:{}\tXR:i:{}",
                record.sample_name,
                record
                    .haplotype
                    .active_region
                    .ref_region
                    .reference_sequence
                    .name,
                record.haplotype.active_region.start_index + 1,
                record.haplotype.length,
                record.ref_name,
                record.pos,
                record.haplotype.alignment.cigar_string(),
                sequence,
                record.haplotype.stats.min,
                record.haplotype.active_region.ref_region.name,
                i32::from(!record.haplotype.active_region.left_end),
                i32::from(!record.haplotype.active_region.right_end),
            )?;
        }
        self.out().flush()?;
        Ok(())
    }
}

#[derive(Clone, Debug)]
struct SamRecord {
    sample_name: String,
    haplotype: Haplotype,
    ref_name: String,
    pos: i32,
}

impl SamRecord {
    fn new(sample_name: String, haplotype: Haplotype) -> Self {
        let ref_region = &haplotype.active_region.ref_region;
        let pos = haplotype.active_region.start_index - ref_region.left_flank_length
            + ref_region.interval.start;
        Self {
            sample_name,
            ref_name: ref_region.reference_sequence.name.clone(),
            haplotype,
            pos,
        }
    }
}

impl Eq for SamRecord {}

impl PartialEq for SamRecord {
    fn eq(&self, other: &Self) -> bool {
        self.ref_name == other.ref_name
            && self.pos == other.pos
            && self.haplotype.length == other.haplotype.length
    }
}

impl Ord for SamRecord {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.ref_name
            .cmp(&other.ref_name)
            .then_with(|| self.pos.cmp(&other.pos))
            .then_with(|| self.haplotype.length.cmp(&other.haplotype.length))
    }
}

impl PartialOrd for SamRecord {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
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

    use super::*;

    #[test]
    fn factory_resolves_sam_and_null_writers() {
        assert_eq!(
            get_writer_class("SAM").unwrap().name(),
            "SamHaplotypeWriter"
        );
        assert_eq!(
            get_writer_class("NULL").unwrap().name(),
            "EmptyHaplotypeWriter"
        );
        assert!(matches!(
            get_writer_class("DOES_NOT_EXIST"),
            Err(HaplotypeWriterError::UnknownWriter(_))
        ));
        assert!(get_writer_description("SAM").unwrap().contains("SAM"));
        assert_eq!(get_writer_description("NOPE"), None);
        assert_eq!(list_writers(), vec!["null", "sam"]);
    }

    #[test]
    fn null_writer_ignores_null_haplotypes() {
        let mut writer = NullHaplotypeWriter::new();
        assert_eq!(writer.name(), "EmptyHaplotypeWriter");
        writer
            .init(None, Some(StreamableOutput::stdout()), vec![make_ref()])
            .unwrap();
        writer.add(None).unwrap();
        writer.add(Some(&haplotype())).unwrap();
        writer.flush().unwrap();
    }

    #[test]
    fn base_methods_validate_and_accept_state_changes() {
        let mut writer = NullHaplotypeWriter::new();
        assert!(matches!(
            writer.init(None, Some(StreamableOutput::stdout()), Vec::new()),
            Err(HaplotypeWriterError::EmptyReferenceArray)
        ));
        writer
            .init(None, Some(StreamableOutput::stdout()), vec![make_ref()])
            .unwrap();
        writer.set_sample_name(None).unwrap();
        writer.set_sample_name(Some("  good  ")).unwrap();
        writer.set_reference_region(reference_region()).unwrap();
        writer.set_reference_region(reference_region()).unwrap();
    }

    #[test]
    fn sam_writer_rejects_null_haplotype() {
        let mut writer = SamHaplotypeWriter::new();
        assert!(matches!(
            writer.add(None),
            Err(HaplotypeWriterError::NullHaplotype)
        ));
    }

    #[test]
    fn sam_writer_flushes_headers_and_records() {
        let out = NamedTempFile::new().unwrap();
        let mut writer = SamHaplotypeWriter::new();
        writer
            .init(
                None,
                Some(StreamableOutput::from_path(out.path(), None)),
                vec![make_ref()],
            )
            .unwrap();
        writer.set_sample_name(Some("s1")).unwrap();
        writer.add(Some(&haplotype())).unwrap();
        writer.flush().unwrap();

        let content = fs::read_to_string(out.path()).unwrap();
        assert!(content.contains("@HD\tVN:1.5"));
        assert!(content.contains("@SQ\tSN:chr1"));
        assert!(content.contains("@PG\tID:Kestrel"));
        assert!(content.contains("s1-chr1-1-16\t0\tchr1\t1\t255\t16="));
        assert!(content.contains("\tXD:i:10\tXN:Z:chr1"));
    }

    #[test]
    fn factory_gets_initialized_sam_writer() {
        let out = NamedTempFile::new().unwrap();
        let mut writer = get_writer(
            Some("SAM"),
            Some(StreamableOutput::from_path(out.path(), None)),
            vec![make_ref()],
        )
        .unwrap();
        writer.add(Some(&haplotype())).unwrap();
        writer.flush().unwrap();

        let content = fs::read_to_string(out.path()).unwrap();
        assert!(content.contains("@HD"));
    }

    fn haplotype() -> Haplotype {
        let ref_region = reference_region();
        let count = vec![10; 13];
        let kmer_util = KmerUtil::new(4).unwrap();
        let active_region = ActiveRegion::new(ref_region, 0, 12, &count, &kmer_util).unwrap();
        let stats = RegionStats::from_counts(
            &count,
            active_region.start_kmer_index,
            active_region.end_kmer_index,
        )
        .unwrap();
        Haplotype::new(
            b"AAAACCCCGGGGTTTT".to_vec(),
            active_region,
            vec![AlignNode::new(AlignNode::MATCH, 16, None).unwrap()],
            100.0,
            None,
            stats,
        )
        .unwrap()
    }

    fn reference_region() -> ReferenceRegion {
        ReferenceRegion::whole(make_ref(), b"AAAACCCCGGGGTTTT", 0).unwrap()
    }

    fn make_ref() -> ReferenceSequence {
        ReferenceSequence::new("chr1", 16, Some(digest()), Some("test")).unwrap()
    }

    fn digest() -> Digest {
        Digest::new([0_u8; 16], "MD5").unwrap()
    }
}
