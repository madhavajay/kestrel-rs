use std::any::type_name_of_val;

use thiserror::Error;

use crate::constants::MIN_KMER_SIZE;

/// Errors returned by k-mer size validation.
#[derive(Debug, Error, Eq, PartialEq)]
pub enum KmerSizeError {
    /// K-mer size was below the minimum.
    #[error("K-mer size is less than the minimum ({minimum}): {actual}")]
    TooSmall {
        /// Minimum supported size.
        minimum: usize,
        /// Actual size.
        actual: usize,
    },
    /// K-mer size was below the minimum with contextual label.
    #[error("{context}: K-mer size is less than the minimum ({minimum}): {actual}")]
    TooSmallWithContext {
        /// Validation context.
        context: String,
        /// Minimum supported size.
        minimum: usize,
        /// Actual size.
        actual: usize,
    },
}

/// Formats an object identity string, or `null` for absent values.
#[must_use]
pub fn object_to_string<T: ?Sized>(value: Option<&T>) -> String {
    match value {
        Some(value) => format!(
            "{}@{:x}",
            type_name_of_val(value),
            value as *const T as *const () as usize
        ),
        None => "null".to_owned(),
    }
}

/// Checks whether a k-mer size meets Kestrel's minimum.
pub fn check_kmer_size(k_size: usize, context: Option<&str>) -> Result<(), KmerSizeError> {
    if k_size >= MIN_KMER_SIZE {
        return Ok(());
    }

    if let Some(context) = context {
        return Err(KmerSizeError::TooSmallWithContext {
            context: context.to_owned(),
            minimum: MIN_KMER_SIZE,
            actual: k_size,
        });
    }

    Err(KmerSizeError::TooSmall {
        minimum: MIN_KMER_SIZE,
        actual: k_size,
    })
}

/// Returns true when a k-mer size is valid for Kestrel.
#[must_use]
pub fn is_valid_kmer_size(k_size: usize) -> bool {
    k_size >= MIN_KMER_SIZE
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn object_to_string_formats_identity_and_null() {
        let value = 7_u32;
        let rendered = object_to_string(Some(&value));

        assert!(rendered.starts_with("u32@"));
        assert_eq!(object_to_string::<u32>(None), "null");
    }

    #[test]
    fn kmer_size_checks_match_java() {
        assert!(check_kmer_size(MIN_KMER_SIZE, None).is_ok());
        assert!(check_kmer_size(MIN_KMER_SIZE + 1, Some("extra")).is_ok());
        assert_eq!(
            check_kmer_size(MIN_KMER_SIZE - 1, None),
            Err(KmerSizeError::TooSmall {
                minimum: MIN_KMER_SIZE,
                actual: MIN_KMER_SIZE - 1
            })
        );
        assert!(matches!(
            check_kmer_size(MIN_KMER_SIZE - 1, Some("from reference")),
            Err(KmerSizeError::TooSmallWithContext { context, .. }) if context == "from reference"
        ));
        assert!(is_valid_kmer_size(MIN_KMER_SIZE));
        assert!(!is_valid_kmer_size(MIN_KMER_SIZE - 1));
    }
}
