use std::cmp::Ordering;
use std::fmt;

use thiserror::Error;

use crate::constants::{ARRAY_EXPAND_FACTOR, MAX_ARRAY_SIZE};

/// Errors returned while constructing digest values.
#[derive(Clone, Debug, Error, Eq, PartialEq)]
pub enum DigestError {
    /// Digest byte array was empty.
    #[error("digest bytes is an empty array")]
    EmptyBytes,
    /// Digest algorithm name was empty.
    #[error("digest algorithm is empty")]
    EmptyAlgorithm,
}

/// Sequence digest bytes and algorithm name.
#[derive(Clone, Debug, Eq, Hash, PartialEq)]
pub struct Digest {
    bytes: Vec<u8>,
    algorithm: String,
}

impl Digest {
    /// Creates a digest from nonempty bytes and a nonempty algorithm name.
    pub fn new(bytes: impl AsRef<[u8]>, algorithm: impl AsRef<str>) -> Result<Self, DigestError> {
        let bytes = bytes.as_ref();
        if bytes.is_empty() {
            return Err(DigestError::EmptyBytes);
        }

        let algorithm = algorithm.as_ref().trim();
        if algorithm.is_empty() {
            return Err(DigestError::EmptyAlgorithm);
        }

        Ok(Self {
            bytes: bytes.to_vec(),
            algorithm: algorithm.to_owned(),
        })
    }

    /// Returns the digest algorithm name.
    #[must_use]
    pub fn algorithm(&self) -> &str {
        &self.algorithm
    }

    /// Returns the digest bytes.
    #[must_use]
    pub fn bytes(&self) -> &[u8] {
        &self.bytes
    }
}

impl Ord for Digest {
    fn cmp(&self, other: &Self) -> Ordering {
        self.algorithm
            .cmp(&other.algorithm)
            .then_with(|| self.bytes.cmp(&other.bytes))
    }
}

impl PartialOrd for Digest {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl fmt::Display for Digest {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for byte in &self.bytes {
            write!(f, "{byte:02x}")?;
        }

        Ok(())
    }
}

/// Errors returned while expanding arrays.
#[derive(Clone, Debug, Error, Eq, PartialEq)]
pub enum ArrayExpandError {
    /// Array is already at the maximum capacity.
    #[error("cannot expand array: array is at its maximum size: {0}")]
    AtMaximum(usize),
}

/// Returns a clone of `values` with expanded capacity.
pub fn expand_array<T: Clone>(values: &[T]) -> Result<Vec<T>, ArrayExpandError> {
    let new_capacity = ((values.len() as f32) * ARRAY_EXPAND_FACTOR) as usize;
    let new_capacity = if new_capacity > MAX_ARRAY_SIZE {
        if values.len() == MAX_ARRAY_SIZE {
            return Err(ArrayExpandError::AtMaximum(MAX_ARRAY_SIZE));
        }
        MAX_ARRAY_SIZE
    } else {
        new_capacity
    };

    let mut expanded = Vec::with_capacity(new_capacity);
    expanded.extend_from_slice(values);
    Ok(expanded)
}

/// No-op message digest implementation.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub struct NullMessageDigest;

impl NullMessageDigest {
    /// Algorithm name for the null digest.
    pub const ALGORITHM: &'static str = "NO_DIGEST";

    /// Returns the null digest algorithm name.
    #[must_use]
    pub fn algorithm(&self) -> &'static str {
        Self::ALGORITHM
    }

    /// Returns the fixed null digest value.
    #[must_use]
    pub fn digest(&self) -> [u8; 1] {
        [0]
    }

    /// Accepts bytes without updating any state.
    pub fn update(&mut self, _bytes: impl AsRef<[u8]>) {}

    /// Resets the digest, which is a no-op.
    pub fn reset(&mut self) {}
}

#[cfg(test)]
mod tests {
    use super::*;

    const MD5_BYTES: [u8; 16] = [
        0xd4, 0x1d, 0x8c, 0xd9, 0x8f, 0x00, 0xb2, 0x04, 0xe9, 0x80, 0x09, 0x98, 0xec, 0xf8, 0x42,
        0x7e,
    ];
    const MD5_HEX: &str = "d41d8cd98f00b204e9800998ecf8427e";

    #[test]
    fn digest_constructs_and_formats_hex() {
        let digest = Digest::new(MD5_BYTES, "  MD5  ").unwrap();

        assert_eq!(digest.algorithm(), "MD5");
        assert_eq!(digest.to_string(), MD5_HEX);
    }

    #[test]
    fn digest_rejects_empty_inputs() {
        assert_eq!(Digest::new([], "MD5"), Err(DigestError::EmptyBytes));
        assert_eq!(
            Digest::new(MD5_BYTES, "   "),
            Err(DigestError::EmptyAlgorithm)
        );
    }

    #[test]
    fn digest_owns_bytes_and_compares_by_algorithm_then_bytes() {
        let mut bytes = MD5_BYTES.to_vec();
        let digest = Digest::new(&bytes, "MD5").unwrap();
        bytes[0] = 0xff;

        assert_eq!(bytes[0], 0xff);
        assert_eq!(digest.to_string(), MD5_HEX);
        assert_eq!(digest, Digest::new(MD5_BYTES, "MD5").unwrap());
        assert!(Digest::new(MD5_BYTES, "MD5").unwrap() < Digest::new(MD5_BYTES, "SHA-1").unwrap());
        assert!(
            Digest::new([0x10, 0x20, 0x30], "X").unwrap()
                > Digest::new([0x10, 0x10, 0x30], "X").unwrap()
        );
        assert!(
            Digest::new([0x10, 0x20], "X").unwrap() < Digest::new([0x10, 0x20, 0x30], "X").unwrap()
        );
    }

    #[test]
    fn expands_array_with_java_growth_factor() {
        let src = vec![1, 2, 3, 4];
        let expanded = expand_array(&src).unwrap();

        assert!(expanded.capacity() > src.len());
        assert_eq!(&expanded[..], &src[..]);
        assert_eq!(expand_array(&[0; 10]).unwrap().capacity(), 15);
        assert_eq!(expand_array::<u8>(&[]).unwrap().capacity(), 0);
    }

    #[test]
    fn null_message_digest_is_noop() {
        let mut digest = NullMessageDigest;

        assert_eq!(NullMessageDigest::ALGORITHM, "NO_DIGEST");
        assert_eq!(digest.algorithm(), "NO_DIGEST");
        digest.update([1, 2, 3]);
        digest.reset();
        assert_eq!(digest.digest(), [0]);
    }
}
