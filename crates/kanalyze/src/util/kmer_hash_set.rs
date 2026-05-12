use std::collections::HashSet;

use super::KmerKey;

/// Hash set for encoded k-mer keys.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct KmerHashSet {
    kmers: HashSet<KmerKey>,
}

impl KmerHashSet {
    /// Creates an empty k-mer set.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Inserts a k-mer and returns whether it was newly inserted.
    pub fn insert(&mut self, kmer: KmerKey) -> bool {
        self.kmers.insert(kmer)
    }

    /// Returns whether the set contains the k-mer.
    #[must_use]
    pub fn contains(&self, kmer: &KmerKey) -> bool {
        self.kmers.contains(kmer)
    }

    /// Removes a k-mer and returns whether it was present.
    pub fn remove(&mut self, kmer: &KmerKey) -> bool {
        self.kmers.remove(kmer)
    }

    /// Removes all k-mers.
    pub fn clear(&mut self) {
        self.kmers.clear();
    }

    /// Returns the number of k-mers in the set.
    #[must_use]
    pub fn len(&self) -> usize {
        self.kmers.len()
    }

    /// Returns whether the set contains no k-mers.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.kmers.is_empty()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::util::KmerUtil;

    #[test]
    fn inserts_contains_removes_and_clones_independently() {
        let util = KmerUtil::new(4).unwrap();
        let kmer = util.encode("ACGT").unwrap();
        let other = util.encode("TGCA").unwrap();
        let mut set = KmerHashSet::new();

        assert!(set.insert(kmer.clone()));
        assert!(!set.insert(kmer.clone()));
        assert!(set.contains(&kmer));
        assert!(!set.contains(&other));

        let mut clone = set.clone();
        assert!(clone.remove(&kmer));
        assert!(!clone.contains(&kmer));
        assert!(set.contains(&kmer));

        assert!(set.remove(&kmer));
        assert!(set.is_empty());
    }

    #[test]
    fn grows_without_losing_entries() {
        let util = KmerUtil::new(6).unwrap();
        let mut set = KmerHashSet::new();
        let mut inserted = Vec::new();

        for value in 0..128_u32 {
            let sequence = format!(
                "{}{}{}{}{}{}",
                base_char((value >> 10) & 0x03),
                base_char((value >> 8) & 0x03),
                base_char((value >> 6) & 0x03),
                base_char((value >> 4) & 0x03),
                base_char((value >> 2) & 0x03),
                base_char(value & 0x03),
            );
            let kmer = util.encode(sequence).unwrap();
            assert!(set.insert(kmer.clone()));
            inserted.push(kmer);
        }

        assert_eq!(set.len(), inserted.len());
        for kmer in inserted {
            assert!(set.contains(&kmer));
        }
    }

    fn base_char(value: u32) -> char {
        match value {
            0 => 'A',
            1 => 'C',
            2 => 'G',
            3 => 'T',
            _ => unreachable!("value is masked to two bits"),
        }
    }
}
