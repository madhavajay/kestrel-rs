use std::cmp::Ordering;
use std::collections::BTreeMap;
use std::fmt;

use kanalyze::comp::reader::{FileSequenceSource, ReaderError, SequenceReader, SequenceSource};
use kanalyze::util::KmerUtil;
use md5::{Digest as Md5Digest, Md5};
use thiserror::Error;

use crate::interval::RegionInterval;
use crate::util::digest::Digest;

/// Errors returned by reference sequence and region handling.
#[derive(Clone, Debug, Error, Eq, PartialEq)]
pub enum ReferenceSequenceError {
    /// Reference name was empty.
    #[error("sequence name is empty")]
    EmptyName,
    /// Reference size was invalid.
    #[error("sequence size is less than 1: {0}")]
    InvalidSize(i32),
    /// Buffer offset was negative.
    #[error("buffer offset is negative: {0}")]
    NegativeBufferOffset(i32),
    /// Left flank length was negative.
    #[error("left flank length is negative: {0}")]
    NegativeLeftFlank(i32),
    /// Right flank length was negative.
    #[error("right flank length is negative: {0}")]
    NegativeRightFlank(i32),
    /// Input buffer was too small for the requested region.
    #[error("buffer is not large enough to contain this reference region")]
    BufferTooSmall,
    /// Gap character was present in a reference region.
    #[error("found a gap in reference region {region} at location {location}: {base}")]
    Gap {
        /// Region name.
        region: String,
        /// One-based location.
        location: i32,
        /// Offending base.
        base: char,
    },
    /// Non-IUPAC character was present in a reference region.
    #[error(
        "found non-IUPAC character in reference region {region} at location {location}: {base}"
    )]
    NonIupac {
        /// Region name.
        region: String,
        /// One-based location.
        location: i32,
        /// Offending base.
        base: char,
    },
    /// Range start and end were invalid.
    #[error("invalid range: start={start}, end={end}")]
    InvalidRange {
        /// Range start.
        start: i32,
        /// Range end.
        end: i32,
    },
    /// Requested base location was outside the region.
    #[error("base location {0} is out of bounds for interval with flanks")]
    BaseOutOfBounds(i32),
    /// Reference sequence did not have regions in the container.
    #[error("reference sequence does not have regions in this container: {0}")]
    MissingReference(String),
    /// Reference metadata conflicted with an existing entry.
    #[error("reference sequence with the same name has already been added: {0}")]
    ReferenceMismatch(String),
    /// New region overlapped an existing region.
    #[error("new region overlaps with an existing region: new={new}, existing={existing}")]
    RegionOverlap {
        /// New region description.
        new: String,
        /// Existing region description.
        existing: String,
    },
    /// Flank length was negative.
    #[error("flank length is negative: {0}")]
    NegativeFlankLength(i32),
    /// Sequence reader error.
    #[error("error reading sequence source: {0}")]
    Reader(String),
    /// Requested interval was outside the reference.
    #[error(
        "missing interval for reference sequence {reference}: reference was too short for interval {interval} (reference size = {size})"
    )]
    MissingInterval {
        /// Reference name.
        reference: String,
        /// Interval description.
        interval: String,
        /// Reference size.
        size: i32,
    },
}

/// Reference sequence metadata.
#[derive(Clone, Debug, Eq)]
pub struct ReferenceSequence {
    /// Reference name.
    pub name: String,
    /// Reference length in bases.
    pub size: i32,
    /// Optional sequence digest.
    pub digest: Option<Digest>,
    /// Source name used to read this sequence.
    pub source_name: String,
}

impl PartialEq for ReferenceSequence {
    fn eq(&self, other: &Self) -> bool {
        self.name == other.name && self.size == other.size && self.digest == other.digest
    }
}

impl ReferenceSequence {
    /// Creates validated reference sequence metadata.
    pub fn new(
        name: &str,
        size: i32,
        digest: Option<Digest>,
        source_name: Option<&str>,
    ) -> Result<Self, ReferenceSequenceError> {
        let name = name.trim();
        if name.is_empty() {
            return Err(ReferenceSequenceError::EmptyName);
        }

        if size < 1 {
            return Err(ReferenceSequenceError::InvalidSize(size));
        }

        let source_name = source_name
            .map(str::trim)
            .filter(|source_name| !source_name.is_empty())
            .unwrap_or("<UNKNOWN_SOURCE>")
            .to_owned();

        Ok(Self {
            name: name.to_owned(),
            size,
            digest,
            source_name,
        })
    }

    /// Returns an interval covering the whole reference sequence.
    pub fn reference_interval(&self) -> RegionInterval {
        RegionInterval::new(None, &self.name, 1, self.size, true)
            .expect("validated reference sequence names and sizes create valid intervals")
    }
}

impl Ord for ReferenceSequence {
    fn cmp(&self, other: &Self) -> Ordering {
        self.name
            .cmp(&other.name)
            .then_with(|| self.size.cmp(&other.size))
            .then_with(|| self.digest.cmp(&other.digest))
    }
}

impl PartialOrd for ReferenceSequence {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl fmt::Display for ReferenceSequence {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "ReferenceSequence[name={}, size={}",
            self.name, self.size
        )?;
        if let Some(digest) = &self.digest {
            write!(f, ", digest={digest}")?;
        }
        f.write_str("]")
    }
}

/// Concrete reference region sequence and interval metadata.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ReferenceRegion {
    /// Region name.
    pub name: String,
    /// Parent reference sequence.
    pub reference_sequence: ReferenceSequence,
    /// Interval represented by this region.
    pub interval: RegionInterval,
    /// Region sequence including flanks.
    pub sequence: Vec<u8>,
    /// Region length including flanks.
    pub size: i32,
    /// Number of left-flank bases.
    pub left_flank_length: i32,
    /// First index after the non-right-flank portion.
    pub right_flank_index: i32,
    /// Offset from reference coordinate to sequence index.
    pub sequence_offset: i32,
    ambiguous_regions: Vec<(i32, i32)>,
}

impl ReferenceRegion {
    /// Creates a reference region from a buffer and interval.
    pub fn new(
        reference_sequence: ReferenceSequence,
        interval: Option<RegionInterval>,
        buffer: &[u8],
        buffer_offset: i32,
        left_flank_length: i32,
        right_flank_length: i32,
    ) -> Result<Self, ReferenceSequenceError> {
        if buffer_offset < 0 {
            return Err(ReferenceSequenceError::NegativeBufferOffset(buffer_offset));
        }
        if left_flank_length < 0 {
            return Err(ReferenceSequenceError::NegativeLeftFlank(left_flank_length));
        }
        if right_flank_length < 0 {
            return Err(ReferenceSequenceError::NegativeRightFlank(
                right_flank_length,
            ));
        }

        let interval = interval.unwrap_or_else(|| reference_sequence.reference_interval());
        let size = interval.end - interval.start + 1 + left_flank_length + right_flank_length;
        if size < 1 || buffer.len() < buffer_offset as usize + size as usize {
            return Err(ReferenceSequenceError::BufferTooSmall);
        }

        let sequence = copy_sequence_from_buffer(&interval, buffer, buffer_offset as usize, size)?;
        let ambiguous_regions = ambiguous_regions(&sequence);
        let right_flank_index = size - right_flank_length;
        let sequence_offset = -(interval.start - left_flank_length);

        Ok(Self {
            name: interval.name.clone(),
            reference_sequence,
            interval,
            sequence,
            size,
            left_flank_length,
            right_flank_index,
            sequence_offset,
            ambiguous_regions,
        })
    }

    /// Creates a region covering an entire reference sequence.
    pub fn whole(
        reference_sequence: ReferenceSequence,
        buffer: &[u8],
        buffer_offset: i32,
    ) -> Result<Self, ReferenceSequenceError> {
        Self::new(reference_sequence, None, buffer, buffer_offset, 0, 0)
    }

    /// Reverse-complements this region in place.
    pub fn reverse_complement(&mut self) {
        self.sequence.reverse();
        for base in &mut self.sequence {
            *base = complement_base(*base);
        }
        self.ambiguous_regions = ambiguous_regions(&self.sequence);
    }

    /// Returns true when an index range overlaps ambiguous bases.
    pub fn contains_ambiguous_by_index(
        &self,
        start_index: i32,
        end_index: i32,
    ) -> Result<bool, ReferenceSequenceError> {
        if start_index < 0 || end_index < start_index {
            return Err(ReferenceSequenceError::InvalidRange {
                start: start_index,
                end: end_index,
            });
        }

        Ok(self
            .ambiguous_regions
            .iter()
            .any(|(start, end)| *start <= end_index && *end >= start_index))
    }

    /// Returns true when a one-based coordinate range overlaps ambiguous bases.
    pub fn contains_ambiguous_by_base_coordinate(
        &self,
        start_coordinate: i32,
        end_coordinate: i32,
    ) -> Result<bool, ReferenceSequenceError> {
        if start_coordinate < 1 || end_coordinate < start_coordinate {
            return Err(ReferenceSequenceError::InvalidRange {
                start: start_coordinate,
                end: end_coordinate,
            });
        }
        self.contains_ambiguous_by_index(start_coordinate - 1, end_coordinate - 1)
    }

    /// Returns true when an index range touches flank bases.
    pub fn is_flank_by_index(&self, start: i32, end: i32) -> Result<bool, ReferenceSequenceError> {
        if start < 0 || end < start {
            return Err(ReferenceSequenceError::InvalidRange { start, end });
        }

        Ok(start < self.left_flank_length || end >= self.right_flank_index)
    }

    /// Returns true when a one-based coordinate range touches flank bases.
    pub fn is_flank_by_coordinate(
        &self,
        start: i32,
        end: i32,
    ) -> Result<bool, ReferenceSequenceError> {
        if start < 1 || end < start {
            return Err(ReferenceSequenceError::InvalidRange { start, end });
        }

        Ok(start <= self.left_flank_length || (end - 1) >= self.right_flank_index)
    }

    /// Returns a base by one-based reference coordinate.
    pub fn get_base(&self, location: i32) -> Result<u8, ReferenceSequenceError> {
        let seq_index = location + self.sequence_offset;
        if seq_index < 0 || seq_index >= self.size {
            return Err(ReferenceSequenceError::BaseOutOfBounds(location));
        }
        Ok(self.sequence[seq_index as usize])
    }
}

impl Ord for ReferenceRegion {
    fn cmp(&self, other: &Self) -> Ordering {
        self.interval.cmp(&other.interval)
    }
}

impl PartialOrd for ReferenceRegion {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl fmt::Display for ReferenceRegion {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "ReferenceRegion[name={}, size={}]", self.name, self.size)
    }
}

fn copy_sequence_from_buffer(
    interval: &RegionInterval,
    buffer: &[u8],
    offset: usize,
    size: i32,
) -> Result<Vec<u8>, ReferenceSequenceError> {
    let mut sequence = Vec::with_capacity(size as usize);
    for index in 0..size as usize {
        let original = buffer[offset + index];
        let Some(base) = normalize_base(original) else {
            let location = index as i32 + 1;
            let base = char::from(original);
            if matches!(original, b'.' | b'-') {
                return Err(ReferenceSequenceError::Gap {
                    region: interval.name.clone(),
                    location,
                    base,
                });
            }
            return Err(ReferenceSequenceError::NonIupac {
                region: interval.name.clone(),
                location,
                base,
            });
        };
        sequence.push(base);
    }
    Ok(sequence)
}

fn normalize_base(base: u8) -> Option<u8> {
    match base {
        b'A' | b'a' => Some(b'A'),
        b'C' | b'c' => Some(b'C'),
        b'G' | b'g' => Some(b'G'),
        b'T' | b't' => Some(b'T'),
        b'U' | b'u' => Some(b'U'),
        b'R' | b'r' => Some(b'R'),
        b'Y' | b'y' => Some(b'Y'),
        b'S' | b's' => Some(b'S'),
        b'W' | b'w' => Some(b'W'),
        b'K' | b'k' => Some(b'K'),
        b'M' | b'm' => Some(b'M'),
        b'B' | b'b' => Some(b'B'),
        b'D' | b'd' => Some(b'D'),
        b'H' | b'h' => Some(b'H'),
        b'V' | b'v' => Some(b'V'),
        b'N' | b'n' => Some(b'N'),
        _ => None,
    }
}

fn complement_base(base: u8) -> u8 {
    match base {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' | b'U' => b'A',
        b'R' => b'Y',
        b'Y' => b'R',
        b'S' => b'S',
        b'W' => b'W',
        b'K' => b'M',
        b'M' => b'K',
        b'B' => b'V',
        b'D' => b'H',
        b'H' => b'D',
        b'V' => b'B',
        _ => b'N',
    }
}

fn is_ambiguous(base: u8) -> bool {
    !matches!(base, b'A' | b'C' | b'G' | b'T' | b'U')
}

fn ambiguous_regions(sequence: &[u8]) -> Vec<(i32, i32)> {
    let mut regions = Vec::new();
    let mut index = 0;
    while index < sequence.len() {
        if is_ambiguous(sequence[index]) {
            let start = index;
            while index < sequence.len() && is_ambiguous(sequence[index]) {
                index += 1;
            }
            regions.push((start as i32, index as i32 - 1));
        } else {
            index += 1;
        }
    }
    regions
}

/// Container for reference regions grouped by reference sequence.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct ReferenceRegionContainer {
    references: Vec<ReferenceSequence>,
    regions: BTreeMap<String, Vec<ReferenceRegion>>,
    region_count: usize,
}

impl ReferenceRegionContainer {
    /// Creates an empty region container.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Adds a reference region, rejecting metadata mismatches and overlaps.
    pub fn add(&mut self, ref_region: ReferenceRegion) -> Result<(), ReferenceSequenceError> {
        let reference_name = ref_region.reference_sequence.name.clone();
        if let Some(existing) = self
            .references
            .iter()
            .find(|reference| reference.name == reference_name)
        {
            if *existing != ref_region.reference_sequence {
                return Err(ReferenceSequenceError::ReferenceMismatch(reference_name));
            }
        } else {
            self.references.push(ref_region.reference_sequence.clone());
        }

        let regions = self.regions.entry(reference_name).or_default();
        let add_index = regions
            .binary_search_by(|existing| existing.cmp(&ref_region))
            .unwrap_or_else(|index| index);

        if add_index > 0
            && regions[add_index - 1]
                .interval
                .overlaps(&ref_region.interval)
        {
            return Err(ReferenceSequenceError::RegionOverlap {
                new: ref_region.to_string(),
                existing: regions[add_index - 1].to_string(),
            });
        }
        if add_index < regions.len() && regions[add_index].interval.overlaps(&ref_region.interval) {
            return Err(ReferenceSequenceError::RegionOverlap {
                new: ref_region.to_string(),
                existing: regions[add_index].to_string(),
            });
        }

        regions.insert(add_index, ref_region);
        self.region_count += 1;
        Ok(())
    }

    /// Adds all present reference regions from an iterator.
    pub fn add_all<I>(&mut self, ref_regions: I) -> Result<(), ReferenceSequenceError>
    where
        I: IntoIterator<Item = Option<ReferenceRegion>>,
    {
        for ref_region in ref_regions.into_iter().flatten() {
            self.add(ref_region)?;
        }
        Ok(())
    }

    /// Returns regions for a reference sequence.
    pub fn get(
        &self,
        ref_sequence: &ReferenceSequence,
    ) -> Result<Vec<ReferenceRegion>, ReferenceSequenceError> {
        self.regions
            .get(&ref_sequence.name)
            .cloned()
            .ok_or_else(|| ReferenceSequenceError::MissingReference(ref_sequence.to_string()))
    }

    /// Iterates over regions in reference order.
    pub fn iter(&self) -> impl Iterator<Item = &ReferenceRegion> {
        self.references.iter().flat_map(|reference| {
            self.regions
                .get(&reference.name)
                .into_iter()
                .flat_map(|regions| regions.iter())
        })
    }

    /// Returns reference sequences in container order.
    pub fn reference_sequences(&self) -> &[ReferenceSequence] {
        &self.references
    }

    /// Returns cloned reference sequences in container order.
    #[must_use]
    pub fn reference_sequence_array(&self) -> Vec<ReferenceSequence> {
        self.references.clone()
    }

    /// Sorts references by their natural ordering.
    pub fn sort_references(&mut self) {
        self.references.sort();
    }

    /// Sorts references using a caller-provided comparison.
    pub fn sort_references_by<F>(&mut self, compare: F)
    where
        F: FnMut(&ReferenceSequence, &ReferenceSequence) -> Ordering,
    {
        self.references.sort_by(compare);
    }

    /// Returns the number of regions.
    #[must_use]
    pub fn len(&self) -> usize {
        self.region_count
    }

    /// Returns true when no regions are stored.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.region_count == 0
    }

    /// Returns the number of reference sequences.
    #[must_use]
    pub fn reference_len(&self) -> usize {
        self.references.len()
    }

    /// Clears all reference sequences and regions.
    pub fn clear(&mut self) {
        self.references.clear();
        self.regions.clear();
        self.region_count = 0;
    }
}

impl fmt::Display for ReferenceRegionContainer {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "ReferenceRegionContainer[references={}, regions={}]",
            self.references.len(),
            self.region_count
        )
    }
}

/// Map from reference name to requested intervals.
pub type ReferenceIntervalMap = BTreeMap<String, Vec<RegionInterval>>;

/// Reads FASTA references into reference regions.
#[derive(Clone, Debug)]
pub struct ReferenceReader {
    /// K-mer utility used by downstream region processing.
    pub kmer_util: KmerUtil,
    flank_length: i32,
    remove_sequence_description: bool,
    reverse_complement_negative_strand: bool,
}

impl ReferenceReader {
    /// Default setting for removing sequence descriptions from names.
    pub const DEFAULT_REMOVE_SEQUENCE_DESCRIPTION: bool = true;
    /// Default setting for reverse-complementing negative-strand intervals.
    pub const DEFAULT_REVERSE_COMPLEMENT_NEGATIVE_STRAND: bool = false;
    /// Digest algorithm used for reference sequence metadata.
    pub const DIGEST_ALGORITHM: &'static str = "MD5";

    /// Creates a reference reader.
    #[must_use]
    pub fn new(kmer_util: KmerUtil) -> Self {
        Self {
            kmer_util,
            flank_length: 0,
            remove_sequence_description: Self::DEFAULT_REMOVE_SEQUENCE_DESCRIPTION,
            reverse_complement_negative_strand: Self::DEFAULT_REVERSE_COMPLEMENT_NEGATIVE_STRAND,
        }
    }

    /// Sets the flank length added around requested intervals.
    pub fn set_flank_length(&mut self, flank_length: i32) -> Result<(), ReferenceSequenceError> {
        if flank_length < 0 {
            return Err(ReferenceSequenceError::NegativeFlankLength(flank_length));
        }
        self.flank_length = flank_length;
        Ok(())
    }

    /// Returns the configured flank length.
    #[must_use]
    pub fn flank_length(&self) -> i32 {
        self.flank_length
    }

    /// Sets whether sequence descriptions are removed from reference names.
    pub fn set_remove_description(&mut self, remove_sequence_description: bool) {
        self.remove_sequence_description = remove_sequence_description;
    }

    /// Returns whether sequence descriptions are removed from reference names.
    #[must_use]
    pub fn remove_description(&self) -> bool {
        self.remove_sequence_description
    }

    /// Sets whether negative-strand intervals are reverse-complemented.
    pub fn set_rev_complement_neg_strand(&mut self, reverse_complement_negative_strand: bool) {
        self.reverse_complement_negative_strand = reverse_complement_negative_strand;
    }

    /// Returns whether negative-strand intervals are reverse-complemented.
    #[must_use]
    pub fn rev_complement_neg_strand(&self) -> bool {
        self.reverse_complement_negative_strand
    }

    /// Reads reference regions from sequence sources and optional interval maps.
    pub fn read(
        &self,
        sources: &[FileSequenceSource],
        interval_map: Option<&ReferenceIntervalMap>,
    ) -> Result<ReferenceRegionContainer, ReferenceSequenceError> {
        let mut container = ReferenceRegionContainer::new();
        for source in sources {
            let records = SequenceReader::new(source.clone())
                .read_all()
                .map_err(reader_error)?;
            for record in records {
                let reference_name = self.reference_name(&record.name);
                let reference_size = i32::try_from(record.sequence.len())
                    .map_err(|_| ReferenceSequenceError::InvalidSize(i32::MAX))?;
                let digest = Digest::new(md5_bytes(&record.sequence), Self::DIGEST_ALGORITHM)
                    .expect("MD5 digest bytes and algorithm are non-empty");
                let reference_sequence = ReferenceSequence::new(
                    &reference_name,
                    reference_size,
                    Some(digest),
                    Some(&source.name()),
                )?;
                let canonical_sequence =
                    canonicalize_reader_sequence(&reference_name, &record.sequence)?;

                if let Some(interval_map) = interval_map {
                    let intervals = interval_map
                        .get(&reference_name)
                        .map_or(&[][..], Vec::as_slice);
                    for interval in intervals {
                        let mut region = self.region_for_interval(
                            reference_sequence.clone(),
                            interval.clone(),
                            &canonical_sequence,
                        )?;
                        if self.reverse_complement_negative_strand && !region.interval.is_fwd {
                            region.reverse_complement();
                        }
                        container.add(region)?;
                    }
                } else {
                    container.add(ReferenceRegion::whole(
                        reference_sequence,
                        &canonical_sequence,
                        0,
                    )?)?;
                }
            }
        }
        Ok(container)
    }

    fn reference_name(&self, name: &str) -> String {
        let name = name.trim();
        if self.remove_sequence_description {
            name.split_whitespace().next().unwrap_or(name).to_owned()
        } else {
            name.to_owned()
        }
    }

    fn region_for_interval(
        &self,
        reference_sequence: ReferenceSequence,
        interval: RegionInterval,
        sequence: &[u8],
    ) -> Result<ReferenceRegion, ReferenceSequenceError> {
        if interval.end > reference_sequence.size {
            return Err(ReferenceSequenceError::MissingInterval {
                reference: reference_sequence.name.clone(),
                interval: interval.to_string(),
                size: reference_sequence.size,
            });
        }

        let left_start = (interval.start - self.flank_length).max(1);
        let right_end = (interval.end + self.flank_length).min(reference_sequence.size);
        let left_flank_length = interval.start - left_start;
        let right_flank_length = right_end - interval.end;

        ReferenceRegion::new(
            reference_sequence,
            Some(interval),
            sequence,
            left_start - 1,
            left_flank_length,
            right_flank_length,
        )
    }

    /// Requests that reading stop.
    pub fn stop(&self) {}
}

fn reader_error(error: ReaderError) -> ReferenceSequenceError {
    ReferenceSequenceError::Reader(error.to_string())
}

fn md5_bytes(sequence: &[u8]) -> Vec<u8> {
    let mut hasher = Md5::new();
    hasher.update(sequence);
    hasher.finalize().to_vec()
}

fn canonicalize_reader_sequence(
    reference_name: &str,
    sequence: &[u8],
) -> Result<Vec<u8>, ReferenceSequenceError> {
    let mut canonical = Vec::with_capacity(sequence.len());
    for (index, base) in sequence.iter().copied().enumerate() {
        let base = match base {
            b'A' | b'a' => b'A',
            b'C' | b'c' => b'C',
            b'G' | b'g' => b'G',
            b'T' | b't' | b'U' | b'u' => b'T',
            b'-' | b'.' => {
                return Err(ReferenceSequenceError::Gap {
                    region: reference_name.to_owned(),
                    location: index as i32 + 1,
                    base: char::from(base),
                });
            }
            _ => b'N',
        };
        canonical.push(base);
    }
    Ok(canonical)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn digest(seed: u8) -> Digest {
        let bytes = (0..16).map(|offset| seed + offset).collect::<Vec<_>>();
        Digest::new(bytes, "MD5").unwrap()
    }

    #[test]
    fn construct_stores_and_normalizes_fields() {
        let d = digest(1);
        let rs = ReferenceSequence::new("  chr1  ", 100, Some(d.clone()), Some("ref.fa")).unwrap();

        assert_eq!(rs.name, "chr1");
        assert_eq!(rs.size, 100);
        assert_eq!(rs.digest, Some(d));
        assert_eq!(rs.source_name, "ref.fa");
        assert_eq!(
            ReferenceSequence::new("chr1", 10, None, None)
                .unwrap()
                .source_name,
            "<UNKNOWN_SOURCE>"
        );
        assert_eq!(
            ReferenceSequence::new("chr1", 10, None, Some("   "))
                .unwrap()
                .source_name,
            "<UNKNOWN_SOURCE>"
        );
    }

    #[test]
    fn rejects_invalid_arguments() {
        assert_eq!(
            ReferenceSequence::new("   ", 100, Some(digest(1)), Some("ref.fa")),
            Err(ReferenceSequenceError::EmptyName)
        );
        assert_eq!(
            ReferenceSequence::new("chr1", 0, Some(digest(1)), Some("ref.fa")),
            Err(ReferenceSequenceError::InvalidSize(0))
        );
        assert_eq!(
            ReferenceSequence::new("chr1", -1, Some(digest(1)), Some("ref.fa")),
            Err(ReferenceSequenceError::InvalidSize(-1))
        );
    }

    #[test]
    fn reference_interval_spans_whole_sequence() {
        let rs = ReferenceSequence::new("chr1", 100, Some(digest(1)), Some("ref.fa")).unwrap();
        let interval = rs.reference_interval();

        assert_eq!(interval.sequence_name, "chr1");
        assert_eq!(interval.start, 1);
        assert_eq!(interval.end, 100);
        assert!(interval.is_fwd);
    }

    #[test]
    fn display_includes_name_size_and_digest() {
        let d = digest(1);
        let rendered = ReferenceSequence::new("chr1", 100, Some(d.clone()), Some("ref.fa"))
            .unwrap()
            .to_string();

        assert!(rendered.contains("chr1"));
        assert!(rendered.contains("100"));
        assert!(rendered.contains(&d.to_string()));
    }

    #[test]
    fn equality_ignores_source_name() {
        let a = ReferenceSequence::new("chr1", 100, Some(digest(1)), Some("ref-a.fa")).unwrap();
        let b = ReferenceSequence::new("chr1", 100, Some(digest(1)), Some("ref-b.fa")).unwrap();
        let c = ReferenceSequence::new("chr2", 100, Some(digest(1)), Some("ref-a.fa")).unwrap();
        let d = ReferenceSequence::new("chr1", 200, Some(digest(1)), Some("ref-a.fa")).unwrap();

        assert_eq!(a, b);
        assert_ne!(a, c);
        assert_ne!(a, d);
    }

    #[test]
    fn ordering_is_name_then_size_then_digest() {
        let a = ReferenceSequence::new("chr1", 100, Some(digest(1)), Some("ref.fa")).unwrap();
        let b = ReferenceSequence::new("chr2", 100, Some(digest(1)), Some("ref.fa")).unwrap();
        let c = ReferenceSequence::new("chr1", 200, Some(digest(1)), Some("ref.fa")).unwrap();
        let d = ReferenceSequence::new("chr1", 100, Some(digest(100)), Some("ref.fa")).unwrap();

        assert!(a < b);
        assert!(a < c);
        assert!(a < d);
    }

    #[test]
    fn reference_region_constructs_whole_interval_and_flanks() {
        let region = ReferenceRegion::whole(seq("chr1", 8), b"ACGTACGT", 0).unwrap();
        assert_eq!(region.sequence, b"ACGTACGT");
        assert_eq!(region.size, 8);
        assert_eq!(region.left_flank_length, 0);
        assert_eq!(region.right_flank_index, 8);

        let interval = RegionInterval::forward(Some("x"), "chr1", 3, 6).unwrap();
        let region =
            ReferenceRegion::new(seq("chr1", 8), Some(interval), b"ACGTACGT", 0, 0, 0).unwrap();
        assert_eq!(region.size, 4);
        assert_eq!(region.sequence, b"ACGT");

        let interval = RegionInterval::forward(Some("x"), "chr1", 4, 9).unwrap();
        let region =
            ReferenceRegion::new(seq("chr1", 14), Some(interval), b"AACCGGTTAACCGG", 0, 2, 2)
                .unwrap();
        assert_eq!(region.size, 10);
        assert_eq!(region.left_flank_length, 2);
        assert_eq!(region.right_flank_index, 8);
    }

    #[test]
    fn reference_region_normalizes_iupac_and_rejects_bad_bases() {
        assert_eq!(
            ReferenceRegion::whole(seq("chr1", 8), b"acgtACGT", 0)
                .unwrap()
                .sequence,
            b"ACGTACGT"
        );
        assert_eq!(
            ReferenceRegion::whole(seq("chr1", 6), b"ACGNRY", 0)
                .unwrap()
                .sequence,
            b"ACGNRY"
        );
        assert!(matches!(
            ReferenceRegion::whole(seq("chr1", 4), b"AC-G", 0),
            Err(ReferenceSequenceError::Gap { .. })
        ));
        assert!(matches!(
            ReferenceRegion::whole(seq("chr1", 4), b"AC.G", 0),
            Err(ReferenceSequenceError::Gap { .. })
        ));
        assert!(matches!(
            ReferenceRegion::whole(seq("chr1", 4), b"ACZG", 0),
            Err(ReferenceSequenceError::NonIupac { .. })
        ));
    }

    #[test]
    fn reference_region_validates_constructor_arguments() {
        let interval = RegionInterval::forward(Some("x"), "chr1", 1, 4).unwrap();
        assert_eq!(
            ReferenceRegion::new(seq("chr1", 4), Some(interval.clone()), b"ACGT", -1, 0, 0),
            Err(ReferenceSequenceError::NegativeBufferOffset(-1))
        );
        assert_eq!(
            ReferenceRegion::new(seq("chr1", 4), None, b"ACGT", 0, -1, 0),
            Err(ReferenceSequenceError::NegativeLeftFlank(-1))
        );
        assert_eq!(
            ReferenceRegion::new(seq("chr1", 4), None, b"ACGT", 0, 0, -1),
            Err(ReferenceSequenceError::NegativeRightFlank(-1))
        );
        assert_eq!(
            ReferenceRegion::new(seq("chr1", 4), Some(interval), b"AC", 0, 0, 0),
            Err(ReferenceSequenceError::BufferTooSmall)
        );
    }

    #[test]
    fn reference_region_reverse_complements() {
        let mut region = ReferenceRegion::whole(seq("chr1", 4), b"ACGT", 0).unwrap();
        region.reverse_complement();
        assert_eq!(region.sequence, b"ACGT");

        let mut region = ReferenceRegion::whole(seq("chr1", 5), b"AGGCC", 0).unwrap();
        region.reverse_complement();
        assert_eq!(region.sequence, b"GGCCT");

        let mut region = ReferenceRegion::whole(seq("chr1", 4), b"RYKM", 0).unwrap();
        region.reverse_complement();
        assert_eq!(region.sequence, b"KMRY");
    }

    #[test]
    fn reference_region_ambiguous_and_flank_queries_match_java_tests() {
        let region = ReferenceRegion::whole(seq("chr1", 5), b"ACNGT", 0).unwrap();
        assert!(region.contains_ambiguous_by_index(2, 2).unwrap());
        assert!(!region.contains_ambiguous_by_index(0, 1).unwrap());
        assert!(region.contains_ambiguous_by_base_coordinate(3, 3).unwrap());
        assert!(!region.contains_ambiguous_by_base_coordinate(1, 2).unwrap());
        assert!(matches!(
            region.contains_ambiguous_by_index(-1, 2),
            Err(ReferenceSequenceError::InvalidRange { .. })
        ));

        let interval = RegionInterval::forward(Some("x"), "chr1", 4, 9).unwrap();
        let flanked =
            ReferenceRegion::new(seq("chr1", 14), Some(interval), b"AACCGGTTAACCGG", 0, 2, 2)
                .unwrap();
        assert!(flanked.is_flank_by_index(0, 0).unwrap());
        assert!(flanked.is_flank_by_index(1, 1).unwrap());
        assert!(flanked.is_flank_by_index(8, 8).unwrap());
        assert!(flanked.is_flank_by_index(9, 9).unwrap());
        assert!(!flanked.is_flank_by_index(2, 7).unwrap());
        assert!(flanked.is_flank_by_coordinate(1, 1).unwrap());
        assert!(flanked.is_flank_by_coordinate(2, 2).unwrap());
        assert!(flanked.is_flank_by_coordinate(9, 9).unwrap());
        assert!(!flanked.is_flank_by_coordinate(3, 8).unwrap());
    }

    #[test]
    fn reference_region_get_base_display_and_ordering() {
        let region = ReferenceRegion::whole(seq("chr1", 4), b"ACGT", 0).unwrap();
        assert_eq!(region.get_base(1).unwrap(), b'A');
        assert_eq!(region.get_base(2).unwrap(), b'C');
        assert_eq!(region.get_base(3).unwrap(), b'G');
        assert_eq!(region.get_base(4).unwrap(), b'T');
        assert_eq!(
            region.get_base(100),
            Err(ReferenceSequenceError::BaseOutOfBounds(100))
        );
        assert!(region.to_string().contains("chr1"));

        let interval_1 = RegionInterval::forward(Some("x"), "chr1", 1, 4).unwrap();
        let interval_2 = RegionInterval::forward(Some("y"), "chr1", 5, 8).unwrap();
        let region_1 =
            ReferenceRegion::new(seq("chr1", 8), Some(interval_1), b"ACGTACGT", 0, 0, 0).unwrap();
        let region_2 =
            ReferenceRegion::new(seq("chr1", 8), Some(interval_2), b"ACGTACGT", 4, 0, 0).unwrap();
        assert!(region_1 < region_2);
        assert!(region_2 > region_1);
    }

    #[test]
    fn reference_region_container_adds_gets_and_iterates_regions() {
        let mut container = ReferenceRegionContainer::new();
        assert!(container.is_empty());

        let chr1 = seq("chr1", 12);
        let chr2 = seq("chr2", 8);
        let chr1_late = region(chr1.clone(), "late", 9, 12, b"TTTT", 0);
        let chr1_early = region(chr1.clone(), "early", 1, 4, b"AAAA", 0);
        let chr2_region = region(chr2.clone(), "other", 1, 4, b"CCCC", 0);

        container.add(chr1_late).unwrap();
        container.add(chr1_early.clone()).unwrap();
        container.add(chr2_region.clone()).unwrap();

        assert_eq!(container.len(), 3);
        assert_eq!(container.reference_len(), 2);
        assert_eq!(container.get(&chr1).unwrap()[0].name, "early");
        assert_eq!(container.get(&chr1).unwrap()[1].name, "late");
        assert_eq!(container.get(&chr2).unwrap(), vec![chr2_region]);

        let names = container
            .iter()
            .map(|region| region.name.as_str())
            .collect::<Vec<_>>();
        assert_eq!(names, vec!["early", "late", "other"]);
        assert!(container.to_string().contains("regions=3"));
    }

    #[test]
    fn reference_region_container_sorts_references_and_returns_copies() {
        let mut container = ReferenceRegionContainer::new();
        let chr_b = seq("chrB", 4);
        let chr_a = seq("chrA", 4);
        container
            .add(region(chr_b.clone(), "b", 1, 4, b"ACGT", 0))
            .unwrap();
        container
            .add(region(chr_a.clone(), "a", 1, 4, b"ACGT", 0))
            .unwrap();

        assert_eq!(container.reference_sequences()[0].name, "chrB");
        container.sort_references();
        assert_eq!(container.reference_sequences()[0].name, "chrA");

        let mut copied = container.reference_sequence_array();
        copied.clear();
        assert_eq!(container.reference_len(), 2);

        container.sort_references_by(|left, right| right.name.cmp(&left.name));
        assert_eq!(container.reference_sequences()[0].name, "chrB");
    }

    #[test]
    fn reference_region_container_rejects_bad_adds_and_clears() {
        let mut container = ReferenceRegionContainer::new();
        let chr1 = seq("chr1", 12);
        container
            .add(region(chr1.clone(), "a", 1, 6, b"AAAAAA", 0))
            .unwrap();
        assert!(matches!(
            container.add(region(chr1.clone(), "overlap", 5, 8, b"CCCC", 0)),
            Err(ReferenceSequenceError::RegionOverlap { .. })
        ));

        let mismatch = ReferenceSequence::new("chr1", 99, Some(digest(1)), Some("test")).unwrap();
        assert!(matches!(
            container.add(region(mismatch, "mismatch", 7, 9, b"GGG", 0)),
            Err(ReferenceSequenceError::ReferenceMismatch(_))
        ));

        let missing = seq("missing", 4);
        assert!(matches!(
            container.get(&missing),
            Err(ReferenceSequenceError::MissingReference(_))
        ));

        container
            .add_all(vec![
                None,
                Some(region(chr1.clone(), "b", 7, 9, b"GGG", 0)),
                Some(region(chr1, "c", 10, 12, b"TTT", 0)),
            ])
            .unwrap();
        assert_eq!(container.len(), 3);
        container.clear();
        assert!(container.is_empty());
        assert_eq!(container.reference_len(), 0);
    }

    #[test]
    fn reference_reader_reads_whole_fasta_and_fastq_fixtures() {
        for file_name in ["general.us-ascii.fasta", "general.us-ascii.fastq"] {
            let source = FileSequenceSource::from_path(fixture_path(file_name), 1).unwrap();
            let reader = ReferenceReader::new(KmerUtil::new(21).unwrap());
            let container = reader.read(&[source], None).unwrap();

            assert_eq!(container.reference_len(), 10);
            assert_eq!(container.len(), 10);
            let first_reference = &container.reference_sequences()[0];
            assert_eq!(first_reference.name, "Seq-1");
            assert_eq!(first_reference.size, 3000);
            assert_eq!(
                first_reference.digest.as_ref().unwrap().algorithm(),
                ReferenceReader::DIGEST_ALGORITHM
            );
            let regions = container.get(first_reference).unwrap();
            assert_eq!(regions.len(), 1);
            assert_eq!(regions[0].sequence.len(), 3000);
            assert!(
                regions[0]
                    .sequence
                    .iter()
                    .all(|base| matches!(base, b'A' | b'C' | b'G' | b'T' | b'N'))
            );
        }
    }

    #[test]
    fn reference_reader_rejects_gap_containing_fixtures() {
        for file_name in ["allchars.us-ascii.fasta", "allchars.us-ascii.fastq"] {
            let source = FileSequenceSource::from_path(fixture_path(file_name), 1).unwrap();
            let reader = ReferenceReader::new(KmerUtil::new(4).unwrap());
            assert!(matches!(
                reader.read(&[source], None),
                Err(ReferenceSequenceError::Gap { .. })
            ));
        }
    }

    #[test]
    fn reference_reader_matches_java_fixture_parameter_matrix() {
        let fixture_cases = [
            ("general.us-ascii.fasta", false),
            ("general.us-ascii.fastq", false),
            ("allchars.us-ascii.fasta", true),
            ("allchars.us-ascii.fastq", true),
        ];
        let k_sizes = [1, 2, 21, 32, 64];

        for (source_id, (file_name, contains_gap)) in fixture_cases.iter().enumerate() {
            for k_size in k_sizes {
                let source =
                    FileSequenceSource::from_path(fixture_path(file_name), source_id as i32 + 1)
                        .unwrap();
                let reader = ReferenceReader::new(KmerUtil::new(k_size).unwrap());
                let result = reader.read(&[source], None);

                if *contains_gap {
                    assert!(
                        matches!(result, Err(ReferenceSequenceError::Gap { .. })),
                        "expected gap rejection for {file_name} k={k_size}"
                    );
                } else {
                    let container = result
                        .unwrap_or_else(|error| panic!("failed {file_name} k={k_size}: {error}"));
                    assert_eq!(
                        container.reference_len(),
                        10,
                        "reference count mismatch for {file_name} k={k_size}"
                    );
                    assert_eq!(
                        container.len(),
                        10,
                        "region count mismatch for {file_name} k={k_size}"
                    );
                    for reference in container.reference_sequences() {
                        let regions = container.get(reference).unwrap();
                        assert_eq!(regions.len(), 1);
                        assert_eq!(regions[0].sequence.len() as i32, reference.size);
                    }
                }
            }
        }
    }

    #[test]
    fn reference_reader_extracts_intervals_with_flanks_and_orientation() {
        let mut reader = ReferenceReader::new(KmerUtil::new(4).unwrap());
        reader.set_flank_length(2).unwrap();
        reader.set_rev_complement_neg_strand(true);

        let mut intervals = ReferenceIntervalMap::new();
        intervals.insert(
            "Seq-1".to_owned(),
            vec![
                RegionInterval::forward(Some("left"), "Seq-1", 1, 4).unwrap(),
                RegionInterval::auto(Some("rev"), "Seq-1", 10, 6).unwrap(),
            ],
        );
        let source =
            FileSequenceSource::from_path(fixture_path("general.us-ascii.fasta"), 1).unwrap();
        let container = reader.read(&[source], Some(&intervals)).unwrap();

        assert_eq!(container.reference_len(), 1);
        assert_eq!(container.len(), 2);
        let reference = &container.reference_sequences()[0];
        let regions = container.get(reference).unwrap();
        assert_eq!(regions[0].name, "left");
        assert_eq!(regions[0].left_flank_length, 0);
        assert_eq!(regions[0].size, 6);
        assert!(regions[0].interval.is_fwd);
        assert_eq!(regions[1].name, "rev");
        assert_eq!(regions[1].left_flank_length, 2);
        assert_eq!(regions[1].size, 9);
        assert!(!regions[1].interval.is_fwd);
    }

    #[test]
    fn reference_reader_properties_match_java_defaults() {
        let mut reader = ReferenceReader::new(KmerUtil::new(4).unwrap());
        assert_eq!(reader.flank_length(), 0);
        assert!(reader.remove_description());
        assert!(!reader.rev_complement_neg_strand());
        assert_eq!(
            reader.set_flank_length(-1),
            Err(ReferenceSequenceError::NegativeFlankLength(-1))
        );

        reader.set_remove_description(false);
        reader.set_rev_complement_neg_strand(true);
        reader.set_flank_length(3).unwrap();
        assert_eq!(reader.flank_length(), 3);
        assert!(!reader.remove_description());
        assert!(reader.rev_complement_neg_strand());
    }

    fn region(
        reference_sequence: ReferenceSequence,
        name: &str,
        start: i32,
        end: i32,
        sequence: &[u8],
        buffer_offset: i32,
    ) -> ReferenceRegion {
        let interval =
            RegionInterval::forward(Some(name), &reference_sequence.name, start, end).unwrap();
        ReferenceRegion::new(
            reference_sequence,
            Some(interval),
            sequence,
            buffer_offset,
            0,
            0,
        )
        .unwrap()
    }

    fn seq(name: &str, size: i32) -> ReferenceSequence {
        ReferenceSequence::new(name, size, Some(digest(1)), Some("test")).unwrap()
    }

    fn fixture_path(name: &str) -> std::path::PathBuf {
        std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("../..")
            .join("fixtures/refreader")
            .join(name)
    }
}
