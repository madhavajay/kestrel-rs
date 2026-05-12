use std::cmp::Ordering;
use std::collections::BTreeMap;
use std::fmt;
use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::path::Path;

use thiserror::Error;

/// Errors returned while creating, storing, or resolving region intervals.
#[derive(Clone, Debug, Error, Eq, PartialEq)]
pub enum RegionIntervalError {
    /// Sequence name was empty.
    #[error("sequence name is empty")]
    EmptySequenceName,
    /// Sequence name contained a tab.
    #[error("sequence name contains tab characters: {0}")]
    TabInSequenceName(String),
    /// Interval name contained a tab.
    #[error("interval name contains tab characters: {0}")]
    TabInName(String),
    /// Start coordinate must be positive.
    #[error("start position is not a positive number: {0}")]
    InvalidStart(i32),
    /// End coordinate must be positive.
    #[error("end position is not a positive number: {0}")]
    InvalidEnd(i32),
    /// Start coordinate was greater than end coordinate.
    #[error("start position ({start}) is greater than the end position ({end})")]
    StartAfterEnd {
        /// Start coordinate.
        start: i32,
        /// End coordinate.
        end: i32,
    },
    /// Reference name was empty.
    #[error("reference name is empty")]
    EmptyReferenceName,
    /// New interval overlapped an existing interval.
    #[error("new interval overlaps with an existing interval: new={new}, existing={existing}")]
    Overlap {
        /// New interval description.
        new: String,
        /// Existing interval description.
        existing: String,
    },
    /// Interval reader name was empty.
    #[error("Cannot get interval reader with an empty name")]
    EmptyReaderName,
    /// No interval reader exists for the requested name.
    #[error("Cannot find class for interval reader: {0}")]
    UnknownReader(String),
}

/// One named reference interval.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct RegionInterval {
    /// Interval name.
    pub name: String,
    /// Reference sequence name.
    pub sequence_name: String,
    /// One-based inclusive start coordinate.
    pub start: i32,
    /// One-based inclusive end coordinate.
    pub end: i32,
    /// True when the interval is on the forward strand.
    pub is_fwd: bool,
}

impl RegionInterval {
    /// Creates a validated region interval.
    pub fn new(
        name: Option<&str>,
        sequence_name: &str,
        start: i32,
        end: i32,
        is_fwd: bool,
    ) -> Result<Self, RegionIntervalError> {
        let sequence_name = sequence_name.trim();

        if sequence_name.is_empty() {
            return Err(RegionIntervalError::EmptySequenceName);
        }

        if sequence_name.contains('\t') {
            return Err(RegionIntervalError::TabInSequenceName(
                sequence_name.to_owned(),
            ));
        }

        if start < 1 {
            return Err(RegionIntervalError::InvalidStart(start));
        }

        if end < 1 {
            return Err(RegionIntervalError::InvalidEnd(end));
        }

        if start > end {
            return Err(RegionIntervalError::StartAfterEnd { start, end });
        }

        let name = name.unwrap_or("").trim();
        let name = if name.is_empty() {
            format!(
                "{}_{start}-{end}",
                sequence_name
                    .split_whitespace()
                    .collect::<Vec<_>>()
                    .join("_")
            )
        } else {
            if name.contains('\t') {
                return Err(RegionIntervalError::TabInName(name.to_owned()));
            }
            name.to_owned()
        };

        Ok(Self {
            name,
            sequence_name: sequence_name.to_owned(),
            start,
            end,
            is_fwd,
        })
    }

    /// Creates a forward-strand region interval.
    pub fn forward(
        name: Option<&str>,
        sequence_name: &str,
        start: i32,
        end: i32,
    ) -> Result<Self, RegionIntervalError> {
        Self::new(name, sequence_name, start, end, true)
    }

    /// Creates an interval, swapping start and end to represent reverse-strand intervals.
    pub fn auto(
        name: Option<&str>,
        sequence_name: &str,
        start: i32,
        end: i32,
    ) -> Result<Self, RegionIntervalError> {
        if start > end {
            Self::new(name, sequence_name, end, start, false)
        } else {
            Self::new(name, sequence_name, start, end, true)
        }
    }

    /// Returns true if two intervals overlap on the same sequence.
    #[must_use]
    pub fn overlaps(&self, other: &Self) -> bool {
        self.sequence_name == other.sequence_name
            && self.start < other.end
            && other.start < self.end
    }
}

impl Ord for RegionInterval {
    fn cmp(&self, other: &Self) -> Ordering {
        self.sequence_name
            .cmp(&other.sequence_name)
            .then_with(|| self.start.cmp(&other.start))
            .then_with(|| self.end.cmp(&other.end))
    }
}

impl PartialOrd for RegionInterval {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl fmt::Display for RegionInterval {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "RegionInterval[name={}, start={}, end={}, fwd={}, seq={}]",
            self.name, self.start, self.end, self.is_fwd, self.sequence_name
        )
    }
}

/// Sorted non-overlapping interval collection grouped by reference sequence.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct RegionIntervalContainer {
    references: BTreeMap<String, Vec<RegionInterval>>,
    size: usize,
}

impl RegionIntervalContainer {
    /// Creates an empty interval container.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Adds an interval, rejecting overlaps on the same sequence.
    pub fn add(&mut self, interval: RegionInterval) -> Result<(), RegionIntervalError> {
        let intervals = self
            .references
            .entry(interval.sequence_name.clone())
            .or_default();
        let add_index = intervals
            .binary_search_by(|existing| existing.cmp(&interval))
            .unwrap_or_else(|index| index);

        if add_index > 0 && intervals[add_index - 1].overlaps(&interval) {
            return Err(RegionIntervalError::Overlap {
                new: interval.to_string(),
                existing: intervals[add_index - 1].to_string(),
            });
        }

        if add_index < intervals.len() && intervals[add_index].overlaps(&interval) {
            return Err(RegionIntervalError::Overlap {
                new: interval.to_string(),
                existing: intervals[add_index].to_string(),
            });
        }

        intervals.insert(add_index, interval);
        self.size += 1;
        Ok(())
    }

    /// Returns intervals for one reference sequence.
    pub fn get_intervals(
        &self,
        reference_name: &str,
    ) -> Result<Vec<RegionInterval>, RegionIntervalError> {
        let reference_name = reference_name.trim();

        if reference_name.is_empty() {
            return Err(RegionIntervalError::EmptyReferenceName);
        }

        Ok(self
            .references
            .get(reference_name)
            .cloned()
            .unwrap_or_default())
    }

    /// Returns a cloned map of all intervals.
    #[must_use]
    pub fn get_map(&self) -> BTreeMap<String, Vec<RegionInterval>> {
        self.references.clone()
    }

    /// Removes all intervals.
    pub fn clear(&mut self) {
        self.references.clear();
        self.size = 0;
    }

    /// Removes intervals for one reference sequence.
    pub fn clear_reference(&mut self, reference_name: Option<&str>) {
        if let Some(reference_name) = reference_name
            && let Some(removed) = self.references.remove(reference_name)
        {
            self.size = self.size.saturating_sub(removed.len());
        }
    }

    /// Returns the number of intervals.
    #[must_use]
    pub fn len(&self) -> usize {
        self.size
    }

    /// Returns true when no intervals are stored.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.size == 0
    }
}

/// Reader interface for interval files.
pub trait IntervalReader {
    /// Returns the reader name.
    fn name(&self) -> &str;
    /// Reads intervals from a path.
    fn read_path(&self, path: &Path) -> io::Result<Vec<RegionInterval>>;
    /// Returns a short reader description.
    fn description(&self) -> &'static str;
    /// Returns true if this reader recognizes a file name.
    fn matches_file_name(&self, file_name: &str) -> bool;
}

/// Creates an interval reader by name.
pub fn get_reader(reader_name: &str) -> Result<Box<dyn IntervalReader>, RegionIntervalError> {
    let reader_name = reader_name.trim();
    if reader_name.is_empty() {
        return Err(RegionIntervalError::EmptyReaderName);
    }

    if reader_name.eq_ignore_ascii_case(bed::BedIntervalReader::NAME) {
        return Ok(Box::new(bed::BedIntervalReader::new()));
    }

    Err(RegionIntervalError::UnknownReader(reader_name.to_owned()))
}

/// Returns the description for an interval reader name.
pub fn get_reader_description(reader_name: &str) -> Option<&'static str> {
    get_reader(reader_name)
        .ok()
        .map(|reader| reader.description())
}

/// Resolves an interval reader name from a file path.
#[must_use]
pub fn resolve_reader_name_by_file(path: &Path) -> Option<&'static str> {
    let file_name = path.file_name()?.to_str()?;
    let reader = bed::BedIntervalReader::new();
    reader.matches_file_name(file_name).then_some("bed")
}

/// Reads intervals from a path using an explicit or auto-detected reader.
pub fn read_path(
    path: &Path,
    reader_name: Option<&str>,
) -> Result<Vec<RegionInterval>, RegionIntervalReadError> {
    let reader_name = reader_name
        .map(str::trim)
        .filter(|name| !name.is_empty())
        .map(str::to_owned)
        .or_else(|| resolve_reader_name_by_file(path).map(str::to_owned))
        .ok_or_else(|| RegionIntervalReadError::AutoResolve(path.display().to_string()))?;
    let reader = get_reader(&reader_name).map_err(RegionIntervalReadError::Reader)?;
    reader.read_path(path).map_err(RegionIntervalReadError::Io)
}

/// Errors returned while reading an interval file.
#[derive(Debug, Error)]
pub enum RegionIntervalReadError {
    /// Interval reader error.
    #[error(transparent)]
    Reader(#[from] RegionIntervalError),
    /// I/O error while reading the interval file.
    #[error(transparent)]
    Io(#[from] io::Error),
    /// Reader could not be auto-detected from the file name.
    #[error("Reader for file could not be automatically determined by the file name: {0}")]
    AutoResolve(String),
}

/// BED interval reader implementation.
pub mod bed {
    use super::*;

    /// Interval reader for BED files.
    #[derive(Clone, Debug, Eq, PartialEq)]
    pub struct BedIntervalReader {
        tab_only: bool,
    }

    impl BedIntervalReader {
        /// Reader name.
        pub const NAME: &'static str = "BED";

        /// Creates a BED reader that splits fields on tabs.
        #[must_use]
        pub fn new() -> Self {
            Self { tab_only: true }
        }

        /// Creates a BED reader with configurable tab-only parsing.
        #[must_use]
        pub fn with_tab_only(tab_only: bool) -> Self {
            Self { tab_only }
        }

        fn split_fields<'a>(&self, line: &'a str) -> Vec<&'a str> {
            if self.tab_only {
                line.split('\t').collect()
            } else {
                line.split_whitespace().collect()
            }
        }
    }

    impl Default for BedIntervalReader {
        fn default() -> Self {
            Self::new()
        }
    }

    impl IntervalReader for BedIntervalReader {
        fn name(&self) -> &str {
            Self::NAME
        }

        fn read_path(&self, path: &Path) -> io::Result<Vec<RegionInterval>> {
            let file = File::open(path).map_err(|err| {
                if err.kind() == io::ErrorKind::NotFound {
                    io::Error::new(
                        io::ErrorKind::NotFound,
                        format!("Interval BED file was not found: {}", path.display()),
                    )
                } else {
                    err
                }
            })?;
            let reader = BufReader::new(file);
            let mut intervals = Vec::new();

            for (index, line) in reader.lines().enumerate() {
                let line_num = index + 1;
                let line = line?;
                let line = line.trim();
                if line.is_empty() || is_ignored_bed_line(line) {
                    continue;
                }

                let fields = self.split_fields(line);
                if fields.len() < 3 {
                    return Err(bed_error(
                        path,
                        line_num,
                        format!(
                            "Must contain at least 3 fields, but only found {}",
                            fields.len()
                        ),
                    ));
                }

                let chrom = fields[0].trim();
                if chrom.is_empty() {
                    return Err(bed_error(
                        path,
                        line_num,
                        "Chromosome name (first field) is empty",
                    ));
                }

                let (chrom_start, chrom_end) =
                    parse_bed3_with_noodles(path, line_num, chrom, fields[1], fields[2])?;
                if chrom_start < 0 {
                    return Err(bed_error(
                        path,
                        line_num,
                        format!("Start location (second field) is negative: {chrom_start}"),
                    ));
                }

                if chrom_end < 0 {
                    return Err(bed_error(
                        path,
                        line_num,
                        format!("End location (third field) is negative: {chrom_end}"),
                    ));
                }
                if chrom_end <= chrom_start {
                    return Err(bed_error(
                        path,
                        line_num,
                        format!(
                            "End location ({chrom_end}) must be greater than the start location ({chrom_start})"
                        ),
                    ));
                }

                let name = fields
                    .get(3)
                    .map(|name| name.trim())
                    .filter(|name| !name.is_empty());
                let is_fwd = match fields.get(5).map(|strand| strand.trim()) {
                    Some("+") | None => true,
                    Some("-") => false,
                    Some(strand) => {
                        return Err(bed_error(
                            path,
                            line_num,
                            format!("Strand (sixth field) must be \"+\" or \"-\": {strand}"),
                        ));
                    }
                };

                intervals.push(
                    RegionInterval::new(name, chrom, chrom_start + 1, chrom_end, is_fwd)
                        .map_err(|err| bed_error(path, line_num, err.to_string()))?,
                );
            }

            Ok(intervals)
        }

        fn description(&self) -> &'static str {
            "BED File."
        }

        fn matches_file_name(&self, file_name: &str) -> bool {
            file_name.to_ascii_lowercase().ends_with(".bed")
        }
    }

    fn is_ignored_bed_line(line: &str) -> bool {
        line.starts_with('#') || line.starts_with("browser ") || line.starts_with("track ")
    }

    fn parse_bed3_with_noodles(
        path: &Path,
        line_num: usize,
        chrom: &str,
        start: &str,
        end: &str,
    ) -> io::Result<(i32, i32)> {
        let line = format!("{chrom}\t{start}\t{end}\n");
        let mut reader = noodles_bed::io::Reader::<3, _>::new(line.as_bytes());
        let mut record = noodles_bed::Record::<3>::default();
        reader.read_record(&mut record).map_err(|err| {
            bed_error(
                path,
                line_num,
                format!("Unable to parse BED3 fields with noodles-bed: {err}"),
            )
        })?;

        let feature_start = record.feature_start().map_err(|err| {
            bed_error(
                path,
                line_num,
                format!("Start location (second field) is not an integer: {start}: {err}"),
            )
        })?;
        let feature_end = record.feature_end().transpose().map_err(|err| {
            bed_error(
                path,
                line_num,
                format!("End location (third field) is not an integer: {end}: {err}"),
            )
        })?;

        let chrom_start = i32::try_from(feature_start.get() - 1).map_err(|_| {
            bed_error(
                path,
                line_num,
                format!("Start location (second field) is out of range: {start}"),
            )
        })?;
        let chrom_end =
            i32::try_from(feature_end.map_or(0, |position| position.get())).map_err(|_| {
                bed_error(
                    path,
                    line_num,
                    format!("End location (third field) is out of range: {end}"),
                )
            })?;

        Ok((chrom_start, chrom_end))
    }

    fn bed_error(path: &Path, line_num: usize, message: impl fmt::Display) -> io::Error {
        io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "Bad BED record on line {line_num} in file {}: {message}",
                path.display()
            ),
        )
    }
}

#[cfg(test)]
mod tests {
    use std::io::Write;

    use tempfile::NamedTempFile;

    use super::*;

    #[test]
    fn construct_stores_fields_and_defaults_forward() {
        let r = RegionInterval::new(Some("exon1"), "chr1", 100, 200, true).unwrap();
        assert_eq!(r.name, "exon1");
        assert_eq!(r.sequence_name, "chr1");
        assert_eq!(r.start, 100);
        assert_eq!(r.end, 200);
        assert!(r.is_fwd);

        assert!(
            RegionInterval::forward(Some("x"), "chrA", 1, 10)
                .unwrap()
                .is_fwd
        );
    }

    #[test]
    fn derives_names_like_java() {
        assert_eq!(
            RegionInterval::forward(None, "chr X", 5, 9).unwrap().name,
            "chr_X_5-9"
        );
        assert_eq!(
            RegionInterval::forward(Some(""), "chrY", 100, 200)
                .unwrap()
                .name,
            "chrY_100-200"
        );
        assert_eq!(
            RegionInterval::forward(None, "chr  with   spaces", 1, 2)
                .unwrap()
                .name,
            "chr_with_spaces_1-2"
        );
    }

    #[test]
    fn rejects_invalid_arguments() {
        assert_eq!(
            RegionInterval::forward(Some("x"), "  ", 1, 2),
            Err(RegionIntervalError::EmptySequenceName)
        );
        assert_eq!(
            RegionInterval::forward(Some("x"), "chr\t1", 1, 2),
            Err(RegionIntervalError::TabInSequenceName("chr\t1".to_owned()))
        );
        assert_eq!(
            RegionInterval::forward(Some("bad\tname"), "chr1", 1, 2),
            Err(RegionIntervalError::TabInName("bad\tname".to_owned()))
        );
        assert_eq!(
            RegionInterval::forward(Some("x"), "chr1", 0, 10),
            Err(RegionIntervalError::InvalidStart(0))
        );
        assert_eq!(
            RegionInterval::forward(Some("x"), "chr1", 1, -1),
            Err(RegionIntervalError::InvalidEnd(-1))
        );
        assert_eq!(
            RegionInterval::forward(Some("x"), "chr1", 10, 5),
            Err(RegionIntervalError::StartAfterEnd { start: 10, end: 5 })
        );
    }

    #[test]
    fn ordering_ignores_interval_name() {
        let a = RegionInterval::forward(Some("a"), "chr1", 1, 10).unwrap();
        let b = RegionInterval::forward(Some("b"), "chr1", 1, 10).unwrap();
        let c = RegionInterval::forward(None, "chr1", 1, 11).unwrap();
        let d = RegionInterval::forward(None, "chr2", 1, 2).unwrap();

        assert_eq!(a.cmp(&b), Ordering::Equal);
        assert!(a < c);
        assert!(d > c);
    }

    #[test]
    fn overlap_uses_java_half_open_check() {
        let a = RegionInterval::forward(None, "chr1", 100, 200).unwrap();
        let b = RegionInterval::forward(None, "chr1", 150, 250).unwrap();
        let adjacent = RegionInterval::forward(None, "chr1", 200, 300).unwrap();
        let other_seq = RegionInterval::forward(None, "chr2", 150, 250).unwrap();

        assert!(a.overlaps(&b));
        assert!(b.overlaps(&a));
        assert!(!a.overlaps(&adjacent));
        assert!(!a.overlaps(&other_seq));
    }

    #[test]
    fn auto_interval_swaps_reverse() {
        let reverse = RegionInterval::auto(Some("x"), "chr1", 200, 100).unwrap();
        let forward = RegionInterval::auto(Some("x"), "chr1", 100, 200).unwrap();

        assert_eq!(reverse.start, 100);
        assert_eq!(reverse.end, 200);
        assert!(!reverse.is_fwd);
        assert!(forward.is_fwd);
    }

    #[test]
    fn display_contains_all_fields() {
        let rendered = RegionInterval::new(Some("name"), "chr1", 5, 9, false)
            .unwrap()
            .to_string();

        assert!(rendered.contains("name"));
        assert!(rendered.contains("chr1"));
        assert!(rendered.contains("5"));
        assert!(rendered.contains("9"));
        assert!(rendered.contains("fwd=false"));
    }

    #[test]
    fn container_adds_sorts_and_rejects_overlaps() {
        let mut container = RegionIntervalContainer::new();
        assert!(container.is_empty());

        container
            .add(RegionInterval::forward(None, "chr1", 100, 200).unwrap())
            .unwrap();
        container
            .add(RegionInterval::forward(None, "chr1", 1, 10).unwrap())
            .unwrap();
        container
            .add(RegionInterval::forward(None, "chr1", 50, 60).unwrap())
            .unwrap();

        assert_eq!(container.len(), 3);
        let intervals = container.get_intervals("chr1").unwrap();
        assert_eq!(intervals[0].start, 1);
        assert_eq!(intervals[1].start, 50);
        assert_eq!(intervals[2].start, 100);

        assert!(matches!(
            container.add(RegionInterval::forward(None, "chr1", 150, 250).unwrap()),
            Err(RegionIntervalError::Overlap { .. })
        ));
    }

    #[test]
    fn container_getters_and_clear_match_java_intent() {
        let mut container = RegionIntervalContainer::new();
        container
            .add(RegionInterval::forward(None, "chr1", 1, 10).unwrap())
            .unwrap();
        container
            .add(RegionInterval::forward(None, "chr2", 1, 10).unwrap())
            .unwrap();
        container
            .add(RegionInterval::forward(None, "chr1", 100, 200).unwrap())
            .unwrap();

        assert_eq!(container.get_intervals("chrUnknown").unwrap().len(), 0);
        assert_eq!(
            container.get_intervals("   "),
            Err(RegionIntervalError::EmptyReferenceName)
        );
        assert_eq!(container.get_intervals("  chr1  ").unwrap().len(), 2);

        let mut map = container.get_map();
        assert_eq!(map.len(), 2);
        map.get_mut("chr1").unwrap().clear();
        assert_eq!(container.get_map().get("chr1").unwrap().len(), 2);

        container.clear_reference(Some("chr1"));
        assert_eq!(container.get_intervals("chr1").unwrap().len(), 0);
        assert_eq!(container.get_intervals("chr2").unwrap().len(), 1);
        assert_eq!(container.len(), 1);
        container.clear_reference(None);
        assert_eq!(container.len(), 1);

        container.clear();
        assert!(container.is_empty());
    }

    #[test]
    fn container_grows_beyond_java_default_capacity() {
        let mut container = RegionIntervalContainer::new();

        for index in 0..20 {
            container
                .add(
                    RegionInterval::forward(None, "chr1", index * 100 + 1, index * 100 + 50)
                        .unwrap(),
                )
                .unwrap();
        }

        let intervals = container.get_intervals("chr1").unwrap();
        assert_eq!(intervals.len(), 20);
        assert!(intervals.windows(2).all(|window| window[0] < window[1]));
    }

    #[test]
    fn bed_reader_metadata_and_factory_match_java_tests() {
        let reader = bed::BedIntervalReader::new();
        assert_eq!(bed::BedIntervalReader::NAME, "BED");
        assert!(reader.matches_file_name("foo.bed"));
        assert!(reader.matches_file_name("foo.BED"));
        assert!(!reader.matches_file_name("foo.txt"));
        assert!(!reader.description().is_empty());

        assert_eq!(get_reader("BED").unwrap().name(), "BED");
        assert_eq!(get_reader("bed").unwrap().name(), "BED");
        assert_eq!(get_reader_description("BED"), Some("BED File."));
        assert_eq!(get_reader_description("NOPE"), None);
        assert_eq!(
            resolve_reader_name_by_file(Path::new("regions.bed")),
            Some("bed")
        );
        assert_eq!(resolve_reader_name_by_file(Path::new("regions.txt")), None);
        assert!(matches!(
            get_reader(""),
            Err(RegionIntervalError::EmptyReaderName)
        ));
        assert!(matches!(
            get_reader("DOES_NOT_EXIST"),
            Err(RegionIntervalError::UnknownReader(_))
        ));
    }

    #[test]
    fn bed_reader_reads_minimal_name_strand_and_ignores_lines() {
        let path = write_temp(
            "# comment line\n\
             browser hide all\n\
             track name=test\n\
             \n\
             chr1\t0\t10\n\
             chr1\t10\t20\texon1\n\
             chr1\t100\t200\tminus\t0\t-\n",
        );

        let intervals = bed::BedIntervalReader::new()
            .read_path(path.path())
            .unwrap();
        assert_eq!(intervals.len(), 3);
        assert_eq!(intervals[0].start, 1);
        assert_eq!(intervals[0].end, 10);
        assert_eq!(intervals[0].sequence_name, "chr1");
        assert!(intervals[0].is_fwd);
        assert_eq!(intervals[1].name, "exon1");
        assert_eq!(intervals[2].start, 101);
        assert_eq!(intervals[2].end, 200);
        assert!(!intervals[2].is_fwd);
    }

    #[test]
    fn bed_reader_supports_whitespace_mode_and_static_read() {
        let path = write_temp("chr1 0 10\nchr2 50 75\n");
        let intervals = bed::BedIntervalReader::with_tab_only(false)
            .read_path(path.path())
            .unwrap();
        assert_eq!(intervals.len(), 2);
        assert_eq!(intervals[1].sequence_name, "chr2");

        let path = write_temp("chr1\t0\t100\n");
        let intervals = read_path(path.path(), Some("BED")).unwrap();
        assert_eq!(intervals.len(), 1);
    }

    #[test]
    fn bed_reader_rejects_bad_records_and_handles_empty_or_missing_files() {
        for content in [
            "chr1\t10\n",
            "\t0\t10\n",
            "chr1\tNOT_INT\t10\n",
            "chr1\t10\tNOT_INT\n",
            "chr1\t-1\t10\n",
            "chr1\t10\t10\n",
            "chr1\t0\t10\tn\t0\t*\n",
        ] {
            let path = write_temp(content);
            assert!(
                bed::BedIntervalReader::new()
                    .read_path(path.path())
                    .is_err()
            );
        }

        assert_eq!(
            bed::BedIntervalReader::new()
                .read_path(Path::new("does-not-exist.bed"))
                .unwrap_err()
                .kind(),
            io::ErrorKind::NotFound
        );

        let empty = write_temp("");
        assert!(
            bed::BedIntervalReader::new()
                .read_path(empty.path())
                .unwrap()
                .is_empty()
        );
    }

    fn write_temp(content: &str) -> NamedTempFile {
        let mut file = NamedTempFile::new().unwrap();
        file.write_all(content.as_bytes()).unwrap();
        file
    }
}
