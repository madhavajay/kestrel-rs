use std::collections::HashMap;

use super::KmerKey;

/// Counter map for encoded k-mers.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct KmerCounter {
    counts: HashMap<KmerKey, u32>,
}

impl KmerCounter {
    /// Creates an empty counter.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Increments a k-mer by one and returns the new count.
    pub fn increment(&mut self, kmer: KmerKey) -> u32 {
        let count = self.counts.entry(kmer).or_insert(0);
        *count = count.saturating_add(1);
        *count
    }

    /// Adds `amount` to a k-mer and returns the new count.
    pub fn add(&mut self, kmer: KmerKey, amount: u32) -> u32 {
        let count = self.counts.entry(kmer).or_insert(0);
        *count = count.saturating_add(amount);
        *count
    }

    /// Returns the count for a k-mer, or zero if absent.
    #[must_use]
    pub fn get(&self, kmer: &KmerKey) -> u32 {
        self.counts.get(kmer).copied().unwrap_or(0)
    }

    /// Removes all counts.
    pub fn clear(&mut self) {
        self.counts.clear();
    }

    /// Returns the number of distinct k-mers.
    #[must_use]
    pub fn len(&self) -> usize {
        self.counts.len()
    }

    /// Returns whether the counter contains no k-mers.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.counts.is_empty()
    }

    /// Iterates over k-mer/count pairs.
    pub fn iter(&self) -> impl Iterator<Item = (&KmerKey, &u32)> {
        self.counts.iter()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::util::KmerUtil;

    #[test]
    fn counts_kmers() {
        let util = KmerUtil::new(4).unwrap();
        let kmer = util.encode("ACGT").unwrap();
        let other = util.encode("TGCA").unwrap();
        let mut counter = KmerCounter::new();

        assert_eq!(counter.get(&kmer), 0);
        assert_eq!(counter.increment(kmer.clone()), 1);
        assert_eq!(counter.increment(kmer.clone()), 2);
        assert_eq!(counter.add(kmer.clone(), u32::MAX), u32::MAX);
        assert_eq!(counter.get(&kmer), u32::MAX);
        assert_eq!(counter.get(&other), 0);

        counter.clear();
        assert!(counter.is_empty());
    }
}
