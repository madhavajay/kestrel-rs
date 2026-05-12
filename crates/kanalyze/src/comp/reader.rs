use std::fs::File;
use std::io::BufReader;
use std::path::{Path, PathBuf};

use bstr::ByteSlice;
use thiserror::Error;

/// Describes a sequence input source.
pub trait SequenceSource {
    /// Returns the filesystem path for the source.
    fn path(&self) -> &Path;
    /// Returns the source sequence format.
    fn format(&self) -> SequenceFormat;
    /// Returns the one-based source identifier.
    fn source_id(&self) -> i32;
    /// Returns an optional source filter specification.
    fn filter_spec(&self) -> Option<&str>;
    /// Returns the display name for the source.
    fn name(&self) -> String;
}

/// Supported sequence input formats.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum SequenceFormat {
    /// FASTA sequence records.
    Fasta,
    /// FASTQ sequence records.
    Fastq,
}

/// File-backed sequence source.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct FileSequenceSource {
    /// Path to the input file.
    pub path: PathBuf,
    /// Sequence format for the file.
    pub format: SequenceFormat,
    /// One-based source identifier.
    pub source_id: i32,
    /// Optional filter specification inherited from KAnalyze.
    pub filter_spec: Option<String>,
    /// Optional character set name.
    pub charset: Option<String>,
}

impl FileSequenceSource {
    /// Creates a file sequence source with an explicit format.
    pub fn new(
        path: impl Into<PathBuf>,
        format: SequenceFormat,
        source_id: i32,
    ) -> Result<Self, ReaderError> {
        if source_id < 1 {
            return Err(ReaderError::InvalidSourceId(source_id));
        }

        Ok(Self {
            path: path.into(),
            format,
            source_id,
            filter_spec: None,
            charset: None,
        })
    }

    /// Creates a file sequence source by inferring the format from its extension.
    pub fn from_path(path: impl Into<PathBuf>, source_id: i32) -> Result<Self, ReaderError> {
        let path = path.into();
        let format = SequenceFormat::from_path(&path)
            .ok_or_else(|| ReaderError::UnknownFormat(path.clone()))?;

        Self::new(path, format, source_id)
    }

    /// Sets the optional filter specification.
    #[must_use]
    pub fn with_filter_spec(mut self, filter_spec: Option<&str>) -> Self {
        self.filter_spec = filter_spec
            .map(str::trim)
            .filter(|value| !value.is_empty())
            .map(str::to_owned);
        self
    }

    /// Sets the optional character set name.
    #[must_use]
    pub fn with_charset(mut self, charset: Option<&str>) -> Self {
        self.charset = charset
            .map(str::trim)
            .filter(|value| !value.is_empty())
            .map(str::to_owned);
        self
    }
}

impl SequenceSource for FileSequenceSource {
    fn path(&self) -> &Path {
        &self.path
    }

    fn format(&self) -> SequenceFormat {
        self.format
    }

    fn source_id(&self) -> i32 {
        self.source_id
    }

    fn filter_spec(&self) -> Option<&str> {
        self.filter_spec.as_deref()
    }

    fn name(&self) -> String {
        self.path
            .file_name()
            .and_then(|name| name.to_str())
            .map_or_else(|| self.path.display().to_string(), str::to_owned)
    }
}

impl SequenceFormat {
    /// Infers a sequence format from a path extension.
    #[must_use]
    pub fn from_path(path: &Path) -> Option<Self> {
        let name = path.file_name()?.to_str()?.to_ascii_lowercase();

        if name.ends_with(".fa") || name.ends_with(".fasta") || name.ends_with(".fna") {
            Some(Self::Fasta)
        } else if name.ends_with(".fq") || name.ends_with(".fastq") {
            Some(Self::Fastq)
        } else {
            None
        }
    }
}

/// Sequence record read from a FASTA or FASTQ source.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct SequenceRecord {
    /// Record name without the format-specific marker.
    pub name: String,
    /// Raw sequence bytes.
    pub sequence: Vec<u8>,
}

/// Errors produced while constructing or reading sequence sources.
#[derive(Debug, Error)]
pub enum ReaderError {
    /// Source identifiers must be positive.
    #[error("source ID must be greater than 0: {0}")]
    InvalidSourceId(i32),
    /// The input path extension did not identify a supported format.
    #[error("could not infer sequence format from path: {0}")]
    UnknownFormat(PathBuf),
    /// An I/O or parser error occurred while reading the source.
    #[error("I/O error reading {path}: {source}")]
    Io {
        /// Path that failed to read.
        path: PathBuf,
        /// Underlying I/O error.
        #[source]
        source: std::io::Error,
    },
}

/// Reader for a sequence source.
pub struct SequenceReader<S = FileSequenceSource> {
    source: S,
}

impl<S> SequenceReader<S>
where
    S: SequenceSource,
{
    /// Creates a reader for `source`.
    #[must_use]
    pub fn new(source: S) -> Self {
        Self { source }
    }

    /// Reads every record from the source.
    pub fn read_all(&self) -> Result<Vec<SequenceRecord>, ReaderError> {
        match self.source.format() {
            SequenceFormat::Fasta => self.read_fasta(),
            SequenceFormat::Fastq => self.read_fastq(),
        }
    }

    fn read_fasta(&self) -> Result<Vec<SequenceRecord>, ReaderError> {
        let file = self.open()?;
        let mut reader = noodles_fasta::io::Reader::new(BufReader::new(file));
        let mut records = Vec::new();

        for result in reader.records() {
            let record = result.map_err(|source| ReaderError::Io {
                path: self.source.path().to_owned(),
                source,
            })?;
            records.push(SequenceRecord {
                name: record.name().as_bstr().to_string(),
                sequence: record.sequence().as_ref().to_vec(),
            });
        }

        Ok(records)
    }

    fn read_fastq(&self) -> Result<Vec<SequenceRecord>, ReaderError> {
        let file = self.open()?;
        let mut reader = noodles_fastq::io::Reader::new(BufReader::new(file));
        let mut records = Vec::new();

        for result in reader.records() {
            let record = result.map_err(|source| ReaderError::Io {
                path: self.source.path().to_owned(),
                source,
            })?;
            records.push(SequenceRecord {
                name: record.name().to_string(),
                sequence: record.sequence().to_vec(),
            });
        }

        Ok(records)
    }

    fn open(&self) -> Result<File, ReaderError> {
        File::open(self.source.path()).map_err(|source| ReaderError::Io {
            path: self.source.path().to_owned(),
            source,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reads_fasta_fixture() {
        let source =
            FileSequenceSource::from_path(fixture_path("general.us-ascii.fasta"), 1).unwrap();
        let records = SequenceReader::new(source).read_all().unwrap();

        assert_eq!(records.len(), 10);
        assert_eq!(records[0].name, "Seq-1");
        assert!(!records[0].sequence.is_empty());
        assert!(records[0].sequence.iter().all(u8::is_ascii_uppercase));
    }

    #[test]
    fn reads_fastq_fixture() {
        let source =
            FileSequenceSource::from_path(fixture_path("general.us-ascii.fastq"), 1).unwrap();
        let records = SequenceReader::new(source).read_all().unwrap();

        assert_eq!(records.len(), 10);
        assert_eq!(records[0].name, "Seq-1");
        assert!(!records[0].sequence.is_empty());
        assert!(records[0].sequence.iter().all(u8::is_ascii_uppercase));
    }

    #[test]
    fn validates_file_sources() {
        assert!(matches!(
            FileSequenceSource::new("reads.txt", SequenceFormat::Fasta, 0),
            Err(ReaderError::InvalidSourceId(0))
        ));
        assert!(matches!(
            FileSequenceSource::from_path("reads.txt", 1),
            Err(ReaderError::UnknownFormat(_))
        ));

        let source = FileSequenceSource::new("reads.fasta", SequenceFormat::Fasta, 1)
            .unwrap()
            .with_filter_spec(Some(" filter "))
            .with_charset(Some(" US-ASCII "));
        assert_eq!(source.filter_spec(), Some("filter"));
        assert_eq!(source.charset.as_deref(), Some("US-ASCII"));
        assert_eq!(source.source_id(), 1);
    }

    fn fixture_path(name: &str) -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("../..")
            .join("fixtures/refreader")
            .join(name)
    }
}
