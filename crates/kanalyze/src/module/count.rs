use thiserror::Error;

use crate::comp::reader::{FileSequenceSource, ReaderError, SequenceReader};
use crate::util::{KmerCounter, KmerError, KmerUtil};

/// Errors from the k-mer count module.
#[derive(Debug, Error)]
pub enum CountError {
    /// Sequence reader error.
    #[error(transparent)]
    Reader(#[from] ReaderError),
    /// K-mer encoding error.
    #[error(transparent)]
    Kmer(#[from] KmerError),
}

/// Counts k-mers from one sequence source.
#[derive(Clone, Debug)]
pub struct CountModule {
    source: FileSequenceSource,
    kmer_util: KmerUtil,
}

impl CountModule {
    /// Creates a counter for a source and k-mer utility.
    #[must_use]
    pub fn new(source: FileSequenceSource, kmer_util: KmerUtil) -> Self {
        Self { source, kmer_util }
    }

    /// Counts all valid k-mers in the configured source.
    pub fn count(&self) -> Result<KmerCounter, CountError> {
        let records = SequenceReader::new(self.source.clone()).read_all()?;
        let mut counter = KmerCounter::new();
        let k = self.kmer_util.k_size();

        for record in records {
            if record.sequence.len() < k {
                continue;
            }

            for window in record.sequence.windows(k) {
                if let Ok(kmer) = self.kmer_util.encode(window) {
                    counter.increment(kmer);
                }
            }
        }

        Ok(counter)
    }
}

#[cfg(test)]
mod tests {
    use std::path::PathBuf;

    use super::*;
    use crate::comp::reader::FileSequenceSource;

    #[test]
    fn counts_kmers_from_fixture() {
        let source =
            FileSequenceSource::from_path(fixture_path("general.us-ascii.fasta"), 1).unwrap();
        let util = KmerUtil::new(3).unwrap();
        let counter = CountModule::new(source, util.clone()).count().unwrap();
        let acg = util.encode("ACG").unwrap();

        assert!(!counter.is_empty());
        assert!(counter.get(&acg) > 0);
    }

    fn fixture_path(name: &str) -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("../..")
            .join("fixtures/refreader")
            .join(name)
    }
}
