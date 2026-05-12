use std::collections::BTreeMap;
use std::fs;
use std::io::{self, Cursor, Read, Write};
use std::path::Path;

use thiserror::Error;

use crate::util::{KmerKey, KmerUtil};

const MAGIC: &[u8; 15] = b"Idx_Kmer_Count\0";
const VERSION: u8 = 1;
const HEADER_SIZE: usize = 80;
const INDEX_RECORD_SIZE: usize = 12;
const ID_SIZE: usize = 32;
const ID_STRING: &[u8] = b"KAnalyze Rust";

/// Header fields for an IKC version 1 file.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct IkcHeaderV1 {
    /// K-mer size stored in the file.
    pub k_size: usize,
    /// Minimizer size used to build the file index.
    pub k_min_size: usize,
    /// Minimizer mask used when indexing records.
    pub mask: u32,
    /// Byte offset where the index section begins.
    pub index_section_offset: u64,
    /// Byte offset where the metadata section begins.
    pub metadata_section_offset: u64,
    /// Writer identifier string from the file header.
    pub id_string: String,
}

/// Errors produced while reading or writing IKC data.
#[derive(Debug, Error)]
pub enum IkcError {
    /// Underlying I/O error.
    #[error(transparent)]
    Io(#[from] io::Error),
    /// File magic did not match the IKC signature.
    #[error("IKC magic mismatch at byte {index}: expected 0x{expected:02X}, got 0x{actual:02X}")]
    BadMagic {
        /// Header byte index.
        index: usize,
        /// Expected byte value.
        expected: u8,
        /// Actual byte value.
        actual: u8,
    },
    /// Unsupported IKC file version.
    #[error("unknown IKC file version: {0}")]
    UnknownVersion(u8),
    /// Reserved header byte was nonzero.
    #[error("reserved header byte {index} must be zero: 0x{value:02X}")]
    ReservedByte {
        /// Header byte index.
        index: usize,
        /// Unexpected byte value.
        value: u8,
    },
    /// Minimizer size is not supported.
    #[error("minimizer size is out of range [1, 15]: {0}")]
    InvalidMinimizerSize(usize),
    /// K-mer size is not supported.
    #[error("k-mer size is less than 1: {0}")]
    InvalidKmerSize(usize),
    /// Index section starts before the fixed header ends.
    #[error("index section offset is less than header size ({HEADER_SIZE}): {0}")]
    InvalidIndexOffset(u64),
    /// Metadata section starts before the index section.
    #[error("metadata section offset is less than index section offset ({index}): {metadata}")]
    InvalidMetadataOffset {
        /// Index section byte offset.
        index: u64,
        /// Metadata section byte offset.
        metadata: u64,
    },
    /// Index section length is malformed.
    #[error("index section length is not a multiple of {INDEX_RECORD_SIZE}: {0}")]
    InvalidIndexLength(u64),
    /// Data section length is not divisible by the encoded record size.
    #[error("data section length is not a multiple of the record size ({record_size}): {length}")]
    InvalidDataLength {
        /// Encoded record size in bytes.
        record_size: usize,
        /// Data section length in bytes.
        length: u64,
    },
    /// File k-mer size does not match the expected utility.
    #[error("expected k-mer size {expected}, but IKC file contains {actual}")]
    KmerSizeMismatch {
        /// Expected k-mer size.
        expected: usize,
        /// Actual file k-mer size.
        actual: usize,
    },
    /// Encoded count was negative.
    #[error("k-mer count is negative in IKC data record")]
    NegativeCount,
    /// Count is too large for the Java IKC signed integer field.
    #[error("k-mer count does not fit in a signed Java int: {0}")]
    CountTooLarge(u32),
    /// Minimizer size exceeds k-mer size.
    #[error("minimizer size {min_size} is greater than k-mer size {k_size}")]
    MinimizerLargerThanKmer {
        /// Requested minimizer size.
        min_size: usize,
        /// K-mer size.
        k_size: usize,
    },
}

/// Reader for indexed k-mer count data.
#[derive(Clone, Debug)]
pub struct IkcReader {
    header: IkcHeaderV1,
    kmer_util: KmerUtil,
    records: BTreeMap<KmerKey, u32>,
}

impl IkcReader {
    /// Opens an IKC file from disk.
    pub fn open(path: impl AsRef<Path>, expected: Option<&KmerUtil>) -> Result<Self, IkcError> {
        let bytes = fs::read(path)?;
        Self::from_bytes(&bytes, expected)
    }

    /// Reads IKC data from bytes.
    pub fn from_bytes(bytes: &[u8], expected: Option<&KmerUtil>) -> Result<Self, IkcError> {
        let header = read_header(bytes)?;
        if let Some(expected) = expected
            && expected.k_size() != header.k_size
        {
            return Err(IkcError::KmerSizeMismatch {
                expected: expected.k_size(),
                actual: header.k_size,
            });
        }

        let kmer_util =
            KmerUtil::new(header.k_size).map_err(|_| IkcError::InvalidKmerSize(header.k_size))?;
        let record_size = kmer_util.word_size_bytes() + 4;
        let data_len = header.index_section_offset - HEADER_SIZE as u64;
        if !data_len.is_multiple_of(record_size as u64) {
            return Err(IkcError::InvalidDataLength {
                record_size,
                length: data_len,
            });
        }

        let mut records = BTreeMap::new();
        let mut cursor = Cursor::new(&bytes[HEADER_SIZE..header.index_section_offset as usize]);
        while (cursor.position() as usize) < data_len as usize {
            let mut kmer_bytes = vec![0; kmer_util.word_size_bytes()];
            cursor.read_exact(&mut kmer_bytes)?;
            let kmer = kmer_util
                .from_bytes(&kmer_bytes)
                .map_err(|_| IkcError::InvalidKmerSize(header.k_size))?;
            let count = read_i32(&mut cursor)?;
            if count < 0 {
                return Err(IkcError::NegativeCount);
            }
            records.insert(kmer, count as u32);
        }

        Ok(Self {
            header,
            kmer_util,
            records,
        })
    }

    /// Returns the parsed file header.
    #[must_use]
    pub fn header(&self) -> &IkcHeaderV1 {
        &self.header
    }

    /// Returns the k-mer utility for the file.
    #[must_use]
    pub fn kmer_util(&self) -> &KmerUtil {
        &self.kmer_util
    }

    /// Returns the count for a k-mer, or zero if absent.
    #[must_use]
    pub fn get(&self, kmer: &KmerKey) -> u32 {
        self.records.get(kmer).copied().unwrap_or(0)
    }
}

/// Writer for indexed k-mer count data.
pub struct IkcWriter {
    kmer_util: KmerUtil,
    k_min_size: usize,
    mask: u32,
}

impl IkcWriter {
    /// Creates an IKC writer with k-mer and minimizer settings.
    pub fn new(kmer_util: KmerUtil, k_min_size: usize, mask: u32) -> Result<Self, IkcError> {
        if !(1..=15).contains(&k_min_size) {
            return Err(IkcError::InvalidMinimizerSize(k_min_size));
        }
        if k_min_size > kmer_util.k_size() {
            return Err(IkcError::MinimizerLargerThanKmer {
                min_size: k_min_size,
                k_size: kmer_util.k_size(),
            });
        }

        Ok(Self {
            kmer_util,
            k_min_size,
            mask,
        })
    }

    /// Writes IKC records to a file path.
    pub fn write_path(
        &self,
        path: impl AsRef<Path>,
        records: &[(KmerKey, u32)],
    ) -> Result<(), IkcError> {
        fs::write(path, self.to_bytes(records)?)?;
        Ok(())
    }

    /// Encodes IKC records into bytes.
    pub fn to_bytes(&self, records: &[(KmerKey, u32)]) -> Result<Vec<u8>, IkcError> {
        let mut sorted = records.to_vec();
        sorted.sort_by_key(|(kmer, _)| {
            (
                self.kmer_util.minimizer(kmer, self.k_min_size, self.mask),
                kmer.words().to_vec(),
            )
        });

        let record_size = self.kmer_util.word_size_bytes() + 4;
        let index_offset = HEADER_SIZE + (sorted.len() * record_size);
        let mut index_records = Vec::<(u32, u64)>::new();
        let mut current_minimizer = None;

        for (index, (kmer, _)) in sorted.iter().enumerate() {
            let minimizer = self.kmer_util.minimizer(kmer, self.k_min_size, self.mask);
            if current_minimizer != Some(minimizer) {
                current_minimizer = Some(minimizer);
                index_records.push((minimizer, (HEADER_SIZE + (index * record_size)) as u64));
            }
        }

        let metadata_offset = index_offset + (index_records.len() * INDEX_RECORD_SIZE);
        let mut out = Vec::with_capacity(metadata_offset);
        write_header(
            &mut out,
            self.kmer_util.k_size(),
            self.k_min_size,
            self.mask,
            index_offset as u64,
            metadata_offset as u64,
        )?;

        for (kmer, count) in &sorted {
            let count = i32::try_from(*count).map_err(|_| IkcError::CountTooLarge(*count))?;
            out.write_all(&self.kmer_util.to_bytes(kmer))?;
            out.write_all(&count.to_be_bytes())?;
        }

        for (minimizer, offset) in index_records {
            out.write_all(&minimizer.to_be_bytes())?;
            out.write_all(&offset.to_be_bytes())?;
        }

        Ok(out)
    }
}

fn read_header(bytes: &[u8]) -> Result<IkcHeaderV1, IkcError> {
    if bytes.len() < HEADER_SIZE {
        return Err(IkcError::Io(io::Error::new(
            io::ErrorKind::UnexpectedEof,
            "IKC file is shorter than the header",
        )));
    }
    for (index, (actual, expected)) in bytes[..MAGIC.len()].iter().zip(MAGIC.iter()).enumerate() {
        if actual != expected {
            return Err(IkcError::BadMagic {
                index,
                expected: *expected,
                actual: *actual,
            });
        }
    }

    let mut cursor = Cursor::new(bytes);
    cursor.set_position(MAGIC.len() as u64);
    let mut version = [0];
    cursor.read_exact(&mut version)?;
    if version[0] != VERSION {
        return Err(IkcError::UnknownVersion(version[0]));
    }

    let mut reserved = [0; 7];
    cursor.read_exact(&mut reserved)?;
    for (index, value) in reserved.into_iter().enumerate() {
        if value != 0 {
            return Err(IkcError::ReservedByte { index, value });
        }
    }

    let mut min_size = [0];
    cursor.read_exact(&mut min_size)?;
    let k_min_size = min_size[0] as usize;
    if !(1..=15).contains(&k_min_size) {
        return Err(IkcError::InvalidMinimizerSize(k_min_size));
    }

    let k_size = read_i32(&mut cursor)?;
    if k_size < 1 {
        return Err(IkcError::InvalidKmerSize(k_size as usize));
    }
    let mask = read_u32(&mut cursor)?;
    let index_section_offset = read_u64(&mut cursor)?;
    let metadata_section_offset = read_u64(&mut cursor)?;
    if index_section_offset < HEADER_SIZE as u64 {
        return Err(IkcError::InvalidIndexOffset(index_section_offset));
    }
    if metadata_section_offset < index_section_offset {
        return Err(IkcError::InvalidMetadataOffset {
            index: index_section_offset,
            metadata: metadata_section_offset,
        });
    }
    if (metadata_section_offset - index_section_offset) % INDEX_RECORD_SIZE as u64 != 0 {
        return Err(IkcError::InvalidIndexLength(
            metadata_section_offset - index_section_offset,
        ));
    }

    let mut id = [0; ID_SIZE];
    cursor.read_exact(&mut id)?;
    let id_end = id.iter().position(|byte| *byte == 0).unwrap_or(ID_SIZE);
    let id_string = String::from_utf8_lossy(&id[..id_end]).to_string();

    Ok(IkcHeaderV1 {
        k_size: k_size as usize,
        k_min_size,
        mask,
        index_section_offset,
        metadata_section_offset,
        id_string,
    })
}

fn write_header(
    mut out: impl Write,
    k_size: usize,
    k_min_size: usize,
    mask: u32,
    index_offset: u64,
    metadata_offset: u64,
) -> Result<(), IkcError> {
    out.write_all(MAGIC)?;
    out.write_all(&[VERSION])?;
    out.write_all(&[0; 7])?;
    out.write_all(&[k_min_size as u8])?;
    out.write_all(&(k_size as i32).to_be_bytes())?;
    out.write_all(&mask.to_be_bytes())?;
    out.write_all(&index_offset.to_be_bytes())?;
    out.write_all(&metadata_offset.to_be_bytes())?;
    let mut id = [0; ID_SIZE];
    let len = ID_STRING.len().min(ID_SIZE);
    id[..len].copy_from_slice(&ID_STRING[..len]);
    out.write_all(&id)?;
    Ok(())
}

fn read_i32(cursor: &mut Cursor<&[u8]>) -> Result<i32, IkcError> {
    let mut bytes = [0; 4];
    cursor.read_exact(&mut bytes)?;
    Ok(i32::from_be_bytes(bytes))
}

fn read_u32(cursor: &mut Cursor<&[u8]>) -> Result<u32, IkcError> {
    let mut bytes = [0; 4];
    cursor.read_exact(&mut bytes)?;
    Ok(u32::from_be_bytes(bytes))
}

fn read_u64(cursor: &mut Cursor<&[u8]>) -> Result<u64, IkcError> {
    let mut bytes = [0; 8];
    cursor.read_exact(&mut bytes)?;
    Ok(u64::from_be_bytes(bytes))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn round_trip_small_ikc_bytes() {
        let util = KmerUtil::new(5).unwrap();
        let acgta = util.encode("ACGTA").unwrap();
        let cgtac = util.encode("CGTAC").unwrap();
        let writer = IkcWriter::new(util.clone(), 3, 0).unwrap();

        let bytes = writer
            .to_bytes(&[(cgtac.clone(), 7), (acgta.clone(), 11)])
            .unwrap();
        let reader = IkcReader::from_bytes(&bytes, Some(&util)).unwrap();

        assert_eq!(reader.header().k_size, 5);
        assert_eq!(reader.header().k_min_size, 3);
        assert_eq!(reader.get(&acgta), 11);
        assert_eq!(reader.get(&cgtac), 7);
        assert_eq!(reader.get(&util.encode("TTTTT").unwrap()), 0);
    }

    #[test]
    fn rejects_bad_magic_and_expected_k_size_mismatch() {
        let util = KmerUtil::new(5).unwrap();
        let writer = IkcWriter::new(util.clone(), 2, 0).unwrap();
        let mut bytes = writer
            .to_bytes(&[(util.encode("ACGTA").unwrap(), 1)])
            .unwrap();
        bytes[0] = b'X';
        assert!(matches!(
            IkcReader::from_bytes(&bytes, Some(&util)),
            Err(IkcError::BadMagic { .. })
        ));

        let bytes = writer
            .to_bytes(&[(util.encode("ACGTA").unwrap(), 1)])
            .unwrap();
        let wrong_util = KmerUtil::new(4).unwrap();
        assert!(matches!(
            IkcReader::from_bytes(&bytes, Some(&wrong_util)),
            Err(IkcError::KmerSizeMismatch {
                expected: 4,
                actual: 5
            })
        ));
    }
}
