use std::path::PathBuf;
use std::sync::atomic::{AtomicBool, Ordering};

use kanalyze::io::ikc::{IkcError, IkcReader, IkcWriter};
use kanalyze::module::count::{CountError, CountModule};
use kanalyze::util::{KmerCounter, KmerKey, KmerUtil};
use thiserror::Error;

use crate::io::InputSample;

/// Common interface for k-mer count maps.
pub trait CountMap {
    /// Returns the count associated with a k-mer.
    fn get(&self, kmer: &KmerKey) -> u32;
    /// Counts all reads in a sample and stores the resulting map.
    fn set(&mut self, sample: InputSample) -> Result<(), CountMapError>;
    /// Requests that an in-progress count operation stop.
    fn abort(&self);
    /// Returns true if abort has been requested.
    fn is_aborted(&self) -> bool;
}

/// Errors returned by count-map implementations.
#[derive(Debug, Error)]
pub enum CountMapError {
    /// Error from the in-memory counter.
    #[error(transparent)]
    Count(#[from] CountError),
    /// Error from IKC reading or writing.
    #[error(transparent)]
    Ikc(#[from] IkcError),
    /// Counting was aborted.
    #[error("count map operation was aborted")]
    Aborted,
}

/// In-memory k-mer count map.
#[derive(Debug)]
pub struct MemoryCountMap {
    kmer_util: KmerUtil,
    counter: KmerCounter,
    sample: Option<InputSample>,
    aborted: AtomicBool,
    /// Drops k-mers with count strictly less than this threshold from the
    /// count map after counting completes. Mirrors Java
    /// `KestrelRunnerBase.getCountModule` which adds
    /// `kmercount:<min_kmer_count>` as a post-count filter whenever
    /// `min_kmer_count > 0`. Default `0` disables filtering.
    min_count: u32,
}

impl MemoryCountMap {
    /// Creates an empty in-memory count map.
    #[must_use]
    pub fn new(kmer_util: KmerUtil) -> Self {
        Self::with_min_count(kmer_util, 0)
    }

    /// Creates an empty in-memory count map that drops k-mers whose count
    /// is strictly less than `min_count` after counting (Java parity for
    /// `kmercount:<min_count>`).
    #[must_use]
    pub fn with_min_count(kmer_util: KmerUtil, min_count: u32) -> Self {
        Self {
            kmer_util,
            counter: KmerCounter::new(),
            sample: None,
            aborted: AtomicBool::new(false),
            min_count,
        }
    }

    /// Returns the k-mer utility used by this map.
    #[must_use]
    pub fn kmer_util(&self) -> &KmerUtil {
        &self.kmer_util
    }

    /// Returns the last sample loaded into this map.
    #[must_use]
    pub fn sample(&self) -> Option<&InputSample> {
        self.sample.as_ref()
    }

    /// Returns the number of counted k-mers.
    #[must_use]
    pub fn len(&self) -> usize {
        self.counter.len()
    }

    /// Returns true when no k-mers are stored.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.counter.is_empty()
    }

    fn pre_module_run(&mut self) {
        self.counter.clear();
    }
}

impl CountMap for MemoryCountMap {
    fn get(&self, kmer: &KmerKey) -> u32 {
        self.counter.get(kmer)
    }

    fn set(&mut self, sample: InputSample) -> Result<(), CountMapError> {
        self.aborted.store(false, Ordering::SeqCst);
        self.sample = Some(sample.clone());
        self.pre_module_run();

        for source in &sample.sources {
            if self.is_aborted() {
                self.sample = None;
                return Err(CountMapError::Aborted);
            }

            let source_counter =
                CountModule::new(source.clone(), self.kmer_util.clone()).count()?;
            for (kmer, count) in source_counter.iter() {
                self.counter.add(kmer.clone(), *count);
            }
        }

        if self.is_aborted() {
            self.sample = None;
            return Err(CountMapError::Aborted);
        }

        if self.min_count > 1 {
            let threshold = self.min_count;
            self.counter.retain(|_, count| *count >= threshold);
        }

        Ok(())
    }

    fn abort(&self) {
        self.aborted.store(true, Ordering::SeqCst);
    }

    fn is_aborted(&self) -> bool {
        self.aborted.load(Ordering::SeqCst)
    }
}

/// IKC-backed k-mer count map.
#[derive(Debug)]
pub struct IkcCountMap {
    kmer_util: KmerUtil,
    k_min_size: usize,
    mask: u32,
    reader: Option<IkcReader>,
    sample: Option<InputSample>,
    temp_file: Option<PathBuf>,
    aborted: AtomicBool,
    /// Drops k-mers with count strictly less than this threshold after
    /// counting (Java parity for `kmercount:<min_count>` post-count filter).
    /// Default `0` disables filtering.
    min_count: u32,
}

impl IkcCountMap {
    /// Creates an IKC-backed count map.
    pub fn new(kmer_util: KmerUtil, k_min_size: usize, mask: u32) -> Result<Self, CountMapError> {
        Self::with_min_count(kmer_util, k_min_size, mask, 0)
    }

    /// Creates an IKC-backed count map that drops k-mers below `min_count`
    /// after counting (Java parity).
    pub fn with_min_count(
        kmer_util: KmerUtil,
        k_min_size: usize,
        mask: u32,
        min_count: u32,
    ) -> Result<Self, CountMapError> {
        IkcWriter::new(kmer_util.clone(), k_min_size, mask)?;
        Ok(Self {
            kmer_util,
            k_min_size,
            mask,
            reader: None,
            sample: None,
            temp_file: None,
            aborted: AtomicBool::new(false),
            min_count,
        })
    }

    /// Returns the last sample loaded into this map.
    #[must_use]
    pub fn sample(&self) -> Option<&InputSample> {
        self.sample.as_ref()
    }

    /// Returns the temporary IKC file path, when counts have been written.
    #[must_use]
    pub fn temp_file(&self) -> Option<&PathBuf> {
        self.temp_file.as_ref()
    }

    fn temp_path(sample: &InputSample) -> PathBuf {
        let mut path = std::env::temp_dir();
        path.push(format!("IKC_{}_{}.ikc", sample.name, std::process::id()));
        path
    }
}

impl CountMap for IkcCountMap {
    fn get(&self, kmer: &KmerKey) -> u32 {
        self.reader.as_ref().map_or(0, |reader| reader.get(kmer))
    }

    fn set(&mut self, sample: InputSample) -> Result<(), CountMapError> {
        self.aborted.store(false, Ordering::SeqCst);
        self.reader = None;
        self.sample = Some(sample.clone());

        let mut counter = KmerCounter::new();
        for source in &sample.sources {
            if self.is_aborted() {
                self.sample = None;
                return Err(CountMapError::Aborted);
            }

            let source_counter =
                CountModule::new(source.clone(), self.kmer_util.clone()).count()?;
            for (kmer, count) in source_counter.iter() {
                counter.add(kmer.clone(), *count);
            }
        }

        if self.is_aborted() {
            self.sample = None;
            return Err(CountMapError::Aborted);
        }

        if self.min_count > 1 {
            let threshold = self.min_count;
            counter.retain(|_, count| *count >= threshold);
        }

        let records = counter
            .iter()
            .map(|(kmer, count)| (kmer.clone(), *count))
            .collect::<Vec<_>>();
        let path = Self::temp_path(&sample);
        IkcWriter::new(self.kmer_util.clone(), self.k_min_size, self.mask)?
            .write_path(&path, &records)?;
        self.reader = Some(IkcReader::open(&path, Some(&self.kmer_util))?);
        self.temp_file = Some(path);

        Ok(())
    }

    fn abort(&self) {
        self.aborted.store(true, Ordering::SeqCst);
    }

    fn is_aborted(&self) -> bool {
        self.aborted.load(Ordering::SeqCst)
    }
}

#[cfg(test)]
mod tests {
    use std::path::PathBuf;

    use kanalyze::comp::reader::FileSequenceSource;

    use super::*;
    use crate::io::InputSample;

    #[test]
    fn memory_count_map_counts_fixture_kmers() {
        let util = KmerUtil::new(3).unwrap();
        let source =
            FileSequenceSource::from_path(fixture_path("general.us-ascii.fasta"), 1).unwrap();
        let sample = InputSample::new(Some("sample"), vec![source]).unwrap();
        let mut map = MemoryCountMap::new(util.clone());

        map.set(sample).unwrap();

        let acg = util.encode("ACG").unwrap();
        assert!(map.get(&acg) > 0);
        assert!(!map.is_empty());
        assert_eq!(map.sample().unwrap().name, "sample");
    }

    #[test]
    fn memory_count_map_merges_multiple_sources_and_clears_between_samples() {
        let util = KmerUtil::new(3).unwrap();
        let source =
            FileSequenceSource::from_path(fixture_path("general.us-ascii.fasta"), 1).unwrap();
        let single = InputSample::new(Some("single"), vec![source.clone()]).unwrap();
        let double = InputSample::new(Some("double"), vec![source.clone(), source]).unwrap();
        let mut map = MemoryCountMap::new(util.clone());
        let acg = util.encode("ACG").unwrap();

        map.set(single).unwrap();
        let single_count = map.get(&acg);
        assert!(single_count > 0);

        map.set(double).unwrap();
        assert_eq!(map.get(&acg), single_count * 2);
    }

    #[test]
    fn memory_count_map_abort_flag_is_observable_and_reset_by_set() {
        let util = KmerUtil::new(3).unwrap();
        let source =
            FileSequenceSource::from_path(fixture_path("general.us-ascii.fasta"), 1).unwrap();
        let sample = InputSample::new(Some("sample"), vec![source]).unwrap();
        let mut map = MemoryCountMap::new(util);

        map.abort();
        assert!(map.is_aborted());

        map.set(sample).unwrap();
        assert!(!map.is_aborted());
    }

    #[test]
    fn ikc_count_map_writes_and_reads_fixture_counts() {
        let util = KmerUtil::new(3).unwrap();
        let source =
            FileSequenceSource::from_path(fixture_path("general.us-ascii.fasta"), 1).unwrap();
        let sample = InputSample::new(Some("ikc_sample"), vec![source]).unwrap();
        let mut map = IkcCountMap::new(util.clone(), 2, 0).unwrap();

        map.set(sample).unwrap();

        let acg = util.encode("ACG").unwrap();
        assert!(map.get(&acg) > 0);
        assert_eq!(map.sample().unwrap().name, "ikc_sample");
        assert!(map.temp_file().unwrap().is_file());
    }

    fn fixture_path(name: &str) -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("../..")
            .join("fixtures/refreader")
            .join(name)
    }
}
