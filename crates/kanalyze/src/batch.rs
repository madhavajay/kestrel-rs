use std::collections::VecDeque;

use crate::comp::reader::SequenceRecord;

/// Reusable collection of sequence records.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct SequenceBatch {
    records: Vec<SequenceRecord>,
}

impl SequenceBatch {
    /// Creates an empty sequence batch.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Creates an empty sequence batch with storage for `capacity` records.
    #[must_use]
    pub fn with_capacity(capacity: usize) -> Self {
        Self {
            records: Vec::with_capacity(capacity),
        }
    }

    /// Appends a sequence record.
    pub fn push(&mut self, record: SequenceRecord) {
        self.records.push(record);
    }

    /// Removes all records while retaining allocated storage.
    pub fn clear(&mut self) {
        self.records.clear();
    }

    /// Returns the number of records in the batch.
    #[must_use]
    pub fn len(&self) -> usize {
        self.records.len()
    }

    /// Returns whether the batch has no records.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.records.is_empty()
    }

    /// Returns the records currently held by the batch.
    #[must_use]
    pub fn records(&self) -> &[SequenceRecord] {
        &self.records
    }
}

/// Fixed-size cache of reusable sequence batches.
#[derive(Debug)]
pub struct BatchCache {
    batches: VecDeque<SequenceBatch>,
    capacity: usize,
    batch_capacity: usize,
}

impl BatchCache {
    /// Creates a cache that stores up to `capacity` batches.
    #[must_use]
    pub fn new(capacity: usize, batch_capacity: usize) -> Self {
        Self {
            batches: VecDeque::with_capacity(capacity),
            capacity,
            batch_capacity,
        }
    }

    /// Returns a cached batch or creates a new one.
    #[must_use]
    pub fn checkout(&mut self) -> SequenceBatch {
        self.batches
            .pop_front()
            .unwrap_or_else(|| SequenceBatch::with_capacity(self.batch_capacity))
    }

    /// Clears and stores a batch if the cache has room.
    pub fn recycle(&mut self, mut batch: SequenceBatch) {
        batch.clear();
        if self.batches.len() < self.capacity {
            self.batches.push_back(batch);
        }
    }

    /// Returns the number of batches currently cached.
    #[must_use]
    pub fn cached_len(&self) -> usize {
        self.batches.len()
    }
}

/// Compatibility alias for Java naming.
pub type SequenceBatchCache = BatchCache;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn batches_hold_records_and_are_recycled_empty() {
        let mut cache = BatchCache::new(1, 4);
        let mut batch = cache.checkout();
        batch.push(SequenceRecord {
            name: "seq".to_owned(),
            sequence: b"ACGT".to_vec(),
        });

        assert_eq!(batch.len(), 1);
        cache.recycle(batch);
        assert_eq!(cache.cached_len(), 1);

        let batch = cache.checkout();
        assert!(batch.is_empty());
    }

    #[test]
    fn cache_respects_capacity() {
        let mut cache = BatchCache::new(1, 1);
        cache.recycle(SequenceBatch::new());
        cache.recycle(SequenceBatch::new());

        assert_eq!(cache.cached_len(), 1);
    }
}
