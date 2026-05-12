use std::collections::HashMap;

use thiserror::Error;

/// Default prefix for unnamed sequences.
pub const DEFAULT_SEQUENCE_NAME: &str = "Sequence";
/// Default prefix for unnamed sequence sources.
pub const DEFAULT_SOURCE_NAME: &str = "SequenceSource";

/// Errors from sequence and source name table operations.
#[derive(Debug, Error, Eq, PartialEq)]
pub enum SequenceNameError {
    /// Source identifiers must be positive.
    #[error("source ID must be greater than 0: {0}")]
    InvalidSourceId(i32),
    /// Sequence identifiers must be positive.
    #[error("sequence ID must be greater than 0: {0}")]
    InvalidSequenceId(i64),
    /// Character buffer range is invalid.
    #[error("invalid character range {start}..{end} for buffer length {len}")]
    InvalidRange {
        /// Start offset.
        start: usize,
        /// End offset.
        end: usize,
        /// Buffer length.
        len: usize,
    },
}

/// Stores optional display names for source and sequence identifiers.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct SequenceNameTable {
    names: HashMap<TableKey, String>,
}

impl SequenceNameTable {
    /// Creates an empty name table.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Adds or replaces a sequence name.
    pub fn add_sequence(
        &mut self,
        source_id: i32,
        sequence_id: i64,
        sequence_name: impl AsRef<str>,
    ) -> Result<(), SequenceNameError> {
        Self::validate_source_id(source_id)?;
        Self::validate_sequence_id(sequence_id)?;

        let sequence_name = normalize_name(
            sequence_name.as_ref(),
            DEFAULT_SEQUENCE_NAME,
            sequence_id.to_string(),
        );
        self.names.insert(
            TableKey {
                source_id,
                sequence_id,
            },
            sequence_name,
        );

        Ok(())
    }

    /// Adds a sequence name from a character buffer range.
    pub fn add_sequence_from_chars(
        &mut self,
        source_id: i32,
        sequence_id: i64,
        sequence_buffer: &[char],
        start: usize,
        end: usize,
    ) -> Result<(), SequenceNameError> {
        if end < start || end > sequence_buffer.len() {
            return Err(SequenceNameError::InvalidRange {
                start,
                end,
                len: sequence_buffer.len(),
            });
        }

        let name = sequence_buffer[start..end].iter().collect::<String>();
        self.add_sequence(source_id, sequence_id, name)
    }

    /// Adds or replaces a sequence source name.
    pub fn add_sequence_source(
        &mut self,
        source_id: i32,
        source_name: impl AsRef<str>,
    ) -> Result<(), SequenceNameError> {
        Self::validate_source_id(source_id)?;

        let source_name = normalize_name(source_name.as_ref(), DEFAULT_SOURCE_NAME, source_id);
        self.names.insert(
            TableKey {
                source_id,
                sequence_id: 0,
            },
            source_name,
        );

        Ok(())
    }

    /// Returns a sequence name if present.
    pub fn get_sequence_name(
        &self,
        source_id: i32,
        sequence_id: i64,
    ) -> Result<Option<&str>, SequenceNameError> {
        Self::validate_source_id(source_id)?;
        Self::validate_sequence_id(sequence_id)?;

        Ok(self
            .names
            .get(&TableKey {
                source_id,
                sequence_id,
            })
            .map(String::as_str))
    }

    /// Returns a sequence name or a Java-compatible default.
    pub fn get_sequence_name_with_default(
        &self,
        source_id: i32,
        sequence_id: i64,
    ) -> Result<String, SequenceNameError> {
        Ok(self.get_sequence_name(source_id, sequence_id)?.map_or_else(
            || format!("{DEFAULT_SEQUENCE_NAME}{sequence_id}"),
            str::to_owned,
        ))
    }

    /// Returns a source name if present.
    pub fn get_source_name(&self, source_id: i32) -> Result<Option<&str>, SequenceNameError> {
        Self::validate_source_id(source_id)?;

        Ok(self
            .names
            .get(&TableKey {
                source_id,
                sequence_id: 0,
            })
            .map(String::as_str))
    }

    /// Returns a source name or a Java-compatible default.
    pub fn get_source_name_with_default(
        &self,
        source_id: i32,
    ) -> Result<String, SequenceNameError> {
        Ok(self.get_source_name(source_id)?.map_or_else(
            || format!("{DEFAULT_SOURCE_NAME}{source_id}"),
            str::to_owned,
        ))
    }

    /// Removes a sequence name.
    pub fn remove_sequence(&mut self, source_id: i32, sequence_id: i64) -> Option<String> {
        if source_id < 1 || sequence_id < 1 {
            return None;
        }

        self.names.remove(&TableKey {
            source_id,
            sequence_id,
        })
    }

    /// Removes a source name.
    pub fn remove_source(&mut self, source_id: i32) -> Option<String> {
        if source_id < 1 {
            return None;
        }

        self.names.remove(&TableKey {
            source_id,
            sequence_id: 0,
        })
    }

    fn validate_source_id(source_id: i32) -> Result<(), SequenceNameError> {
        if source_id < 1 {
            return Err(SequenceNameError::InvalidSourceId(source_id));
        }

        Ok(())
    }

    fn validate_sequence_id(sequence_id: i64) -> Result<(), SequenceNameError> {
        if sequence_id < 1 {
            return Err(SequenceNameError::InvalidSequenceId(sequence_id));
        }

        Ok(())
    }
}

#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
struct TableKey {
    source_id: i32,
    sequence_id: i64,
}

fn normalize_name(name: &str, default_prefix: &str, id: impl std::fmt::Display) -> String {
    let name = name.trim();

    if name.is_empty() {
        format!("{default_prefix}{id}")
    } else {
        name.to_owned()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sequence_names_are_trimmed_and_defaulted() {
        let mut table = SequenceNameTable::new();

        table.add_sequence(1, 1, " chr1 ").unwrap();
        table.add_sequence(1, 2, "  ").unwrap();

        assert_eq!(table.get_sequence_name(1, 1).unwrap(), Some("chr1"));
        assert_eq!(table.get_sequence_name(1, 2).unwrap(), Some("Sequence2"));
        assert_eq!(
            table.get_sequence_name_with_default(1, 3).unwrap(),
            "Sequence3"
        );
    }

    #[test]
    fn source_names_are_trimmed_and_defaulted() {
        let mut table = SequenceNameTable::new();

        table.add_sequence_source(1, " reads.fastq ").unwrap();
        table.add_sequence_source(2, "").unwrap();

        assert_eq!(table.get_source_name(1).unwrap(), Some("reads.fastq"));
        assert_eq!(table.get_source_name(2).unwrap(), Some("SequenceSource2"));
        assert_eq!(
            table.get_source_name_with_default(3).unwrap(),
            "SequenceSource3"
        );
    }

    #[test]
    fn names_can_be_removed() {
        let mut table = SequenceNameTable::new();

        table.add_sequence(1, 1, "chr1").unwrap();
        table.add_sequence_source(1, "source").unwrap();

        assert_eq!(table.remove_sequence(1, 1), Some("chr1".to_owned()));
        assert_eq!(table.remove_sequence(1, 1), None);
        assert_eq!(table.remove_source(1), Some("source".to_owned()));
        assert_eq!(table.remove_source(1), None);
        assert_eq!(table.remove_sequence(0, 1), None);
        assert_eq!(table.remove_source(0), None);
    }

    #[test]
    fn invalid_ids_are_errors() {
        let mut table = SequenceNameTable::new();

        assert_eq!(
            table.add_sequence(0, 1, "chr1"),
            Err(SequenceNameError::InvalidSourceId(0))
        );
        assert_eq!(
            table.add_sequence(1, 0, "chr1"),
            Err(SequenceNameError::InvalidSequenceId(0))
        );
        assert_eq!(
            table.add_sequence_source(0, "source"),
            Err(SequenceNameError::InvalidSourceId(0))
        );
    }

    #[test]
    fn can_add_name_from_char_range() {
        let mut table = SequenceNameTable::new();
        let chars = ['x', 'c', 'h', 'r', '1', 'y'];

        table.add_sequence_from_chars(1, 1, &chars, 1, 5).unwrap();

        assert_eq!(table.get_sequence_name(1, 1).unwrap(), Some("chr1"));
        assert_eq!(
            table.add_sequence_from_chars(1, 1, &chars, 5, 7),
            Err(SequenceNameError::InvalidRange {
                start: 5,
                end: 7,
                len: 6
            })
        );
    }
}
