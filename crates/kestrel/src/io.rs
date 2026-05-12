use std::fs::File;
use std::io;
use std::path::{Path, PathBuf};

use kanalyze::comp::reader::{FileSequenceSource, SequenceSource};
use thiserror::Error;

/// One logical input sample and its sequence sources.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct InputSample {
    /// Sample name.
    pub name: String,
    /// Sequence sources for this sample.
    pub sources: Vec<FileSequenceSource>,
}

/// Errors returned while creating input samples.
#[derive(Debug, Error, Eq, PartialEq)]
pub enum InputSampleError {
    /// At least one sequence source is required.
    #[error("cannot create sample with no input sources")]
    EmptySources,
}

impl InputSample {
    /// Creates an input sample, deriving the name from the first source when absent.
    pub fn new(
        name: Option<&str>,
        sources: Vec<FileSequenceSource>,
    ) -> Result<Self, InputSampleError> {
        if sources.is_empty() {
            return Err(InputSampleError::EmptySources);
        }

        let name = name
            .map(str::trim)
            .filter(|name| !name.is_empty())
            .map_or_else(|| sources[0].name(), str::to_owned);

        Ok(Self { name, sources })
    }
}

/// Output target that may refer to screen streams, files, or file descriptors.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum StreamableOutput {
    /// Standard output.
    Stdout,
    /// Standard error.
    Stderr,
    /// File-system output path.
    File {
        /// File path.
        path: PathBuf,
        /// Display name for diagnostics.
        name: String,
    },
    /// Existing file descriptor.
    Fd {
        /// File descriptor number.
        fd: i32,
        /// Display name for diagnostics.
        name: String,
    },
}

/// Errors returned while creating streamable output targets.
#[derive(Debug, Error, Eq, PartialEq)]
pub enum StreamableOutputError {
    /// File name was empty.
    #[error("cannot set output for an empty file name")]
    EmptyFileName,
    /// Output object is unsupported.
    #[error("unsupported output object")]
    Unsupported,
}

impl StreamableOutput {
    /// Returns a stdout output target.
    #[must_use]
    pub fn stdout() -> Self {
        Self::Stdout
    }

    /// Returns a stderr output target.
    #[must_use]
    pub fn stderr() -> Self {
        Self::Stderr
    }

    /// Creates a file output target from a file name.
    pub fn from_file_name(
        file_name: &str,
        name: Option<&str>,
    ) -> Result<Self, StreamableOutputError> {
        let file_name = file_name.trim();
        if file_name.is_empty() {
            return Err(StreamableOutputError::EmptyFileName);
        }

        Ok(Self::from_path(PathBuf::from(file_name), name))
    }

    /// Creates a file output target from a path.
    #[must_use]
    pub fn from_path(path: impl Into<PathBuf>, name: Option<&str>) -> Self {
        let path = path.into();
        let name = normalized_name(name).unwrap_or_else(|| file_name_or_display(&path));

        Self::File { path, name }
    }

    /// Creates an output target from a file descriptor.
    #[must_use]
    pub fn from_fd(fd: i32, name: Option<&str>) -> Self {
        match fd {
            1 => Self::Stdout,
            2 => Self::Stderr,
            _ => Self::Fd {
                fd,
                name: normalized_name(name)
                    .unwrap_or_else(|| "<UNKNOWN_FILE_DESCRIPTOR>".to_owned()),
            },
        }
    }

    /// Creates the file for file-backed outputs.
    pub fn create_file(&self) -> io::Result<Option<File>> {
        match self {
            Self::File { path, .. } => File::create(path).map(Some),
            Self::Stdout | Self::Stderr | Self::Fd { .. } => Ok(None),
        }
    }

    /// Returns true for stdout or stderr.
    #[must_use]
    pub fn is_screen(&self) -> bool {
        matches!(self, Self::Stdout | Self::Stderr)
    }

    /// Returns the diagnostic name for this output target.
    #[must_use]
    pub fn name(&self) -> &str {
        match self {
            Self::Stdout => "<STDOUT>",
            Self::Stderr => "<STDERR>",
            Self::File { name, .. } | Self::Fd { name, .. } => name,
        }
    }
}

fn normalized_name(name: Option<&str>) -> Option<String> {
    name.map(str::trim)
        .filter(|name| !name.is_empty())
        .map(str::to_owned)
}

fn file_name_or_display(path: &Path) -> String {
    path.file_name()
        .and_then(|name| name.to_str())
        .map_or_else(|| path.display().to_string(), str::to_owned)
}

#[cfg(test)]
mod tests {
    use tempfile::NamedTempFile;

    use super::*;
    use kanalyze::comp::reader::{FileSequenceSource, SequenceFormat};

    #[test]
    fn input_sample_stores_or_derives_name() {
        let source = FileSequenceSource::new("reads.fasta", SequenceFormat::Fasta, 1).unwrap();

        let explicit = InputSample::new(Some("  sample  "), vec![source.clone()]).unwrap();
        assert_eq!(explicit.name, "sample");

        let derived = InputSample::new(None, vec![source]).unwrap();
        assert_eq!(derived.name, "reads.fasta");
        assert_eq!(
            InputSample::new(Some("x"), Vec::new()),
            Err(InputSampleError::EmptySources)
        );
    }

    #[test]
    fn streamable_output_handles_screen_and_files() {
        assert!(StreamableOutput::stdout().is_screen());
        assert!(StreamableOutput::stderr().is_screen());
        assert_eq!(StreamableOutput::stdout().name(), "<STDOUT>");
        assert_eq!(StreamableOutput::stderr().name(), "<STDERR>");

        let out = StreamableOutput::from_file_name("/tmp/foo.txt", None).unwrap();
        assert_eq!(out.name(), "foo.txt");
        assert!(!out.is_screen());

        let out = StreamableOutput::from_file_name("/tmp/foo.txt", Some(" custom ")).unwrap();
        assert_eq!(out.name(), "custom");

        assert_eq!(
            StreamableOutput::from_file_name("   ", None),
            Err(StreamableOutputError::EmptyFileName)
        );
    }

    #[test]
    fn streamable_output_handles_fds_and_file_creation() {
        assert_eq!(StreamableOutput::from_fd(1, None), StreamableOutput::Stdout);
        assert_eq!(StreamableOutput::from_fd(2, None), StreamableOutput::Stderr);
        assert_eq!(
            StreamableOutput::from_fd(9, None).name(),
            "<UNKNOWN_FILE_DESCRIPTOR>"
        );

        let tmp = NamedTempFile::new().unwrap();
        let out = StreamableOutput::from_path(tmp.path(), None);
        let file = out.create_file().unwrap();
        assert!(file.is_some());
    }
}
