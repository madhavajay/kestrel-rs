use std::hash::{Hash, Hasher};

use thiserror::Error;

use super::Base;

const BASES_PER_WORD: usize = 16;
const BITS_PER_BASE: usize = 2;

/// Errors from k-mer encoding and decoding.
#[derive(Debug, Error, Eq, PartialEq)]
pub enum KmerError {
    /// K-mer size must be positive.
    #[error("k-mer size must be greater than 0: {0}")]
    InvalidKSize(usize),
    /// Input sequence or byte buffer length did not match the configured k-mer size.
    #[error("sequence length {actual} does not match k-mer size {expected}")]
    InvalidLength {
        /// Expected length.
        expected: usize,
        /// Actual length.
        actual: usize,
    },
    /// Input sequence contains a non-ACGT base.
    #[error("invalid base {base:?} at offset {offset}")]
    InvalidBase {
        /// Invalid byte.
        base: u8,
        /// Zero-based offset of the invalid byte.
        offset: usize,
    },
}

/// Encoded k-mer key in Java-compatible word layout.
#[derive(Clone, Debug, Eq, Ord, PartialEq, PartialOrd)]
pub struct KmerKey {
    words: Vec<u32>,
}

impl KmerKey {
    /// Creates a k-mer key from encoded words.
    #[must_use]
    pub fn from_words(words: Vec<u32>) -> Self {
        Self { words }
    }

    /// Returns the encoded words.
    #[must_use]
    pub fn words(&self) -> &[u32] {
        &self.words
    }
}

impl Hash for KmerKey {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.words.hash(state);
    }
}

/// Utility for encoding, decoding, and hashing fixed-size k-mers.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct KmerUtil {
    k_size: usize,
    word_size: usize,
    unused_bases: usize,
}

impl KmerUtil {
    /// Creates a k-mer utility for `k_size`.
    pub fn new(k_size: usize) -> Result<Self, KmerError> {
        if k_size == 0 {
            return Err(KmerError::InvalidKSize(k_size));
        }

        let word_size = ((k_size - 1) / BASES_PER_WORD) + 1;
        let unused_bases = word_size * BASES_PER_WORD - k_size;

        Ok(Self {
            k_size,
            word_size,
            unused_bases,
        })
    }

    /// Returns the configured k-mer size.
    #[must_use]
    pub fn k_size(&self) -> usize {
        self.k_size
    }

    /// Returns the number of encoded 32-bit words.
    #[must_use]
    pub fn word_size(&self) -> usize {
        self.word_size
    }

    /// Returns the number of bytes needed to encode one k-mer.
    #[must_use]
    pub fn word_size_bytes(&self) -> usize {
        (self.k_size * BITS_PER_BASE).div_ceil(8)
    }

    /// Returns the number of bytes used by the most-significant encoded word.
    #[must_use]
    pub fn msw_byte_size(&self) -> usize {
        self.word_size_bytes() - ((self.word_size - 1) * 4)
    }

    /// Encodes an ACGT sequence into a k-mer key.
    pub fn encode(&self, sequence: impl AsRef<[u8]>) -> Result<KmerKey, KmerError> {
        let sequence = sequence.as_ref();

        if sequence.len() != self.k_size {
            return Err(KmerError::InvalidLength {
                expected: self.k_size,
                actual: sequence.len(),
            });
        }

        let mut words = vec![0; self.word_size];

        for (offset, base) in sequence.iter().copied().enumerate() {
            let value = byte_to_base(base)
                .ok_or(KmerError::InvalidBase { base, offset })?
                .value() as u32;
            let (word_index, shift) = self.word_and_shift(offset);
            words[word_index] |= value << shift;
        }

        Ok(KmerKey { words })
    }

    /// Decodes a k-mer key to a string of bases.
    #[must_use]
    pub fn to_base_string(&self, kmer: &KmerKey) -> String {
        self.decode(kmer).into_iter().collect()
    }

    /// Decodes a k-mer key to base characters.
    #[must_use]
    pub fn decode(&self, kmer: &KmerKey) -> Vec<char> {
        (0..self.k_size)
            .map(|offset| {
                let (word_index, shift) = self.word_and_shift(offset);
                let value = ((kmer.words[word_index] >> shift) & 0x03) as u8;
                Base::from_value(value).expect("masked k-mer base value is always valid")
            })
            .map(Base::as_char)
            .collect()
    }

    /// Returns the reverse-complement k-mer key.
    pub fn reverse_complement(&self, kmer: &KmerKey) -> KmerKey {
        let mut words = vec![0; self.word_size];

        for offset in 0..self.k_size {
            let source_offset = self.k_size - 1 - offset;
            let (source_word, source_shift) = self.word_and_shift(source_offset);
            let value = ((kmer.words[source_word] >> source_shift) & 0x03) as u8;
            let complement = Base::from_value(value)
                .expect("masked k-mer base value is always valid")
                .complement()
                .value() as u32;
            let (word_index, shift) = self.word_and_shift(offset);
            words[word_index] |= complement << shift;
        }

        KmerKey { words }
    }

    /// Encodes a k-mer key as big-endian bytes.
    #[must_use]
    pub fn to_bytes(&self, kmer: &KmerKey) -> Vec<u8> {
        let mut bytes = Vec::with_capacity(self.word_size_bytes());
        let msw_size = self.msw_byte_size();
        let first = kmer.words[0].to_be_bytes();
        bytes.extend_from_slice(&first[4 - msw_size..]);
        for word in &kmer.words[1..] {
            bytes.extend_from_slice(&word.to_be_bytes());
        }
        bytes
    }

    /// Decodes a k-mer key from big-endian bytes.
    pub fn from_bytes(&self, bytes: &[u8]) -> Result<KmerKey, KmerError> {
        if bytes.len() != self.word_size_bytes() {
            return Err(KmerError::InvalidLength {
                expected: self.word_size_bytes(),
                actual: bytes.len(),
            });
        }

        let msw_size = self.msw_byte_size();
        let mut words = Vec::with_capacity(self.word_size);
        let mut first = 0_u32;
        for byte in &bytes[..msw_size] {
            first = (first << 8) | u32::from(*byte);
        }
        words.push(first);

        for chunk in bytes[msw_size..].chunks_exact(4) {
            words.push(u32::from_be_bytes([chunk[0], chunk[1], chunk[2], chunk[3]]));
        }

        Ok(KmerKey { words })
    }

    /// Computes the masked canonical minimizer over both strands.
    #[must_use]
    pub fn minimizer(&self, kmer: &KmerKey, min_size: usize, mask: u32) -> u32 {
        let sequence: Vec<u8> = self
            .decode(kmer)
            .into_iter()
            .map(|base| base as u8)
            .collect();
        let rev = self.reverse_complement(kmer);
        let rev_sequence: Vec<u8> = self
            .decode(&rev)
            .into_iter()
            .map(|base| base as u8)
            .collect();

        sequence
            .windows(min_size)
            .chain(rev_sequence.windows(min_size))
            .map(|window| minimizer_value(window) ^ mask)
            .min()
            .unwrap_or(0)
    }

    /// Compares two encoded k-mer keys.
    #[must_use]
    pub fn eq(&self, left: &KmerKey, right: &KmerKey) -> bool {
        left == right
    }

    /// Computes the Java-compatible k-mer hash code.
    #[must_use]
    pub fn hash_code(&self, kmer: &KmerKey) -> u64 {
        let mut hash = 1125899906842597_u64;

        for word in &kmer.words {
            hash = hash.wrapping_mul(31).wrapping_add(u64::from(*word));
        }

        hash
    }

    fn word_and_shift(&self, offset: usize) -> (usize, usize) {
        let slot = self.unused_bases + offset;
        let word_index = slot / BASES_PER_WORD;
        let base_index = slot % BASES_PER_WORD;
        let shift = (BASES_PER_WORD - 1 - base_index) * BITS_PER_BASE;

        (word_index, shift)
    }
}

fn minimizer_value(sequence: &[u8]) -> u32 {
    sequence.iter().fold(0_u32, |value, base| {
        let base = byte_to_base(*base)
            .expect("minimizer input is decoded from a valid k-mer")
            .value() as u32;
        (value << BITS_PER_BASE) | base
    })
}

fn byte_to_base(base: u8) -> Option<Base> {
    match base {
        b'A' | b'a' => Some(Base::A),
        b'C' | b'c' => Some(Base::C),
        b'G' | b'g' => Some(Base::G),
        b'T' | b't' | b'U' | b'u' => Some(Base::T),
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn encodes_to_java_word_layout() {
        let util = KmerUtil::new(17).unwrap();
        let kmer = util.encode("ACGTACGTACGTACGTA").unwrap();

        assert_eq!(kmer.words(), &[0x0000_0000, 0x6C6C_6C6C]);
        assert_eq!(util.to_base_string(&kmer), "ACGTACGTACGTACGTA");
    }

    #[test]
    fn accepts_lowercase_and_rna_u_as_t() {
        let util = KmerUtil::new(4).unwrap();
        let kmer = util.encode("acgu").unwrap();

        assert_eq!(util.to_base_string(&kmer), "ACGT");
    }

    #[test]
    fn rejects_invalid_input() {
        let util = KmerUtil::new(3).unwrap();

        assert_eq!(
            util.encode("AC"),
            Err(KmerError::InvalidLength {
                expected: 3,
                actual: 2
            })
        );
        assert_eq!(
            util.encode("ANC"),
            Err(KmerError::InvalidBase {
                base: b'N',
                offset: 1
            })
        );
    }

    #[test]
    fn reverse_complement_is_encoded() {
        let util = KmerUtil::new(5).unwrap();
        let kmer = util.encode("AACGT").unwrap();
        let rev = util.reverse_complement(&kmer);

        assert_eq!(util.to_base_string(&rev), "ACGTT");
    }

    #[test]
    fn hash_is_deterministic() {
        let util = KmerUtil::new(8).unwrap();
        let left = util.encode("ACGTACGT").unwrap();
        let right = util.encode("ACGTACGT").unwrap();
        let other = util.encode("ACGTACGA").unwrap();

        assert_eq!(util.hash_code(&left), util.hash_code(&right));
        assert_ne!(util.hash_code(&left), util.hash_code(&other));
        assert!(util.eq(&left, &right));
    }
}
