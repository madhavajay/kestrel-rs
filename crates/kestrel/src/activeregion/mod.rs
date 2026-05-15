use std::fmt;

use kanalyze::util::{KmerKey, KmerUtil};
use thiserror::Error;

use crate::align::{AlignNode, AlignmentWeight, AlignmentWeightError, TraceMatrix};
use crate::counter::CountMap;
use crate::refreader::{ReferenceRegion, ReferenceSequenceError};

/// Errors from active-region summary statistics.
#[derive(Clone, Debug, Error, Eq, PartialEq)]
pub enum RegionStatsError {
    /// Start index was negative.
    #[error("Start index is negative: {0}")]
    NegativeStart(i32),
    /// End index did not define a nonempty range inside the count array.
    #[error(
        "End index is not less than the array length ({length}): the range from start to end is not at least 1 ({range}): {end}"
    )]
    InvalidEnd {
        /// Count array length.
        length: usize,
        /// Range length from start to end.
        range: i32,
        /// Requested end index.
        end: i32,
    },
}

/// Errors from active-region and haplotype construction.
#[derive(Clone, Debug, Error, PartialEq)]
pub enum ActiveRegionError {
    /// Start k-mer index was less than the left-open sentinel.
    #[error("start k-mer index is less than -1: {0}")]
    StartTooNegative(i32),
    /// End k-mer index was outside the count array.
    #[error(
        "end k-mer index is not between -1 and the last k-mer in the reference ({last}): {end}"
    )]
    EndOutOfRange {
        /// Requested end k-mer index.
        end: i32,
        /// Last valid k-mer index.
        last: i32,
    },
    /// End k-mer would extend past the reference sequence.
    #[error("end k-mer at index ({end}) extends past the right end of the reference sequence")]
    EndExtendsPastReference {
        /// Requested end k-mer index.
        end: i32,
    },
    /// Start k-mer index did not precede the end k-mer index.
    #[error("start k-mer index ({start}) must come before end k-mer index ({end})")]
    StartNotBeforeEnd {
        /// Requested start k-mer index.
        start: i32,
        /// Requested end k-mer index.
        end: i32,
    },
    /// Both ends were open, spanning the whole reference.
    #[error("active region may not span the entire reference sequence")]
    BothEndsOpen,
    /// Region endpoint k-mer contains ambiguous bases.
    #[error(
        "reference sequence contains an ambiguous base in the {side}-end k-mer at index {index}"
    )]
    AmbiguousEndKmer {
        /// Endpoint side.
        side: &'static str,
        /// K-mer index.
        index: i32,
    },
    /// Region statistics error.
    #[error(transparent)]
    Stats(#[from] RegionStatsError),
    /// Haplotype alignment list was empty.
    #[error("cannot create haplotype: alignment list is empty")]
    EmptyAlignmentList,
    /// Alternate alignment index was out of range.
    #[error(
        "altAlign is less than 0 or greater than the last alternate alignment index ({last}): {index}"
    )]
    AlignmentIndexOutOfRange {
        /// Requested alignment index.
        index: i32,
        /// Last valid alignment index.
        last: i32,
    },
    /// Region haplotype had no haplotypes.
    #[error("cannot create region haplotype: haplotype array is empty")]
    EmptyHaplotypeList,
}

/// Errors from active-region detection.
#[derive(Clone, Debug, Error, PartialEq)]
pub enum ActiveRegionDetectorError {
    /// Minimum difference must be positive.
    #[error("minimum difference must be at least 1: {0}")]
    InvalidMinimumDifference(i32),
    /// Difference quantile must be in `[0, 1)`.
    #[error("difference quantile must be in [0, 1): {0}")]
    InvalidDifferenceQuantile(f64),
    /// Decay minimum must be in `[0, 1]`.
    #[error("exponential decay minimum must be in [0, 1]: {0}")]
    InvalidDecayMinimum(f64),
    /// Decay alpha must be in `(0, 1)`.
    #[error("exponential decay alpha must be in (0, 1): {0}")]
    InvalidDecayAlpha(f64),
    /// Peak scan length cannot be negative.
    #[error("peak scan length may not be negative: {0}")]
    NegativePeakScanLength(i32),
    /// Scan limit factor cannot be negative.
    #[error("scan limit factor may not be negative: {0}")]
    NegativeScanLimitFactor(f64),
    /// Maximum repeat count cannot be negative.
    #[error("maximum repeat count may not be negative: {0}")]
    NegativeMaxRepeatCount(i32),
    /// Active-region construction error.
    #[error(transparent)]
    ActiveRegion(#[from] ActiveRegionError),
    /// Region statistics error.
    #[error(transparent)]
    RegionStats(#[from] RegionStatsError),
    /// Reference-region query error.
    #[error(transparent)]
    Reference(#[from] ReferenceSequenceError),
    /// Alignment-weight configuration error.
    #[error(transparent)]
    AlignmentWeight(#[from] AlignmentWeightError),
}

/// Error type produced by [`ActiveRegionDetector::detect_from_counts_with`].
/// Wraps either an internal detector error or an error returned from the
/// caller-supplied `accept` callback.
#[derive(Debug)]
pub enum AcceptError<E> {
    /// Error returned from the detector itself.
    Detector(ActiveRegionDetectorError),
    /// Error returned from the caller-supplied callback.
    Callback(E),
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum LeftScanResult {
    Candidate(Option<(i32, i32)>),
    SkipPeak { next_index: usize, count_l: i32 },
}

impl<E: fmt::Display> fmt::Display for AcceptError<E> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Detector(err) => fmt::Display::fmt(err, f),
            Self::Callback(err) => fmt::Display::fmt(err, f),
        }
    }
}

impl<E: std::error::Error + 'static> std::error::Error for AcceptError<E> {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Detector(err) => Some(err),
            Self::Callback(err) => Some(err),
        }
    }
}

/// Detects active regions from k-mer depth changes across a reference region.
#[derive(Clone, Debug, PartialEq)]
pub struct ActiveRegionDetector {
    kmer_util: KmerUtil,
    alignment_weight: AlignmentWeight,
    anchor_both_ends: bool,
    count_reverse_kmers: bool,
    minimum_difference: i32,
    difference_quantile: f64,
    peak_scan_length: i32,
    scan_limit_factor: f64,
    scan_limit: i32,
    decay_minimum: f64,
    decay_alpha: f64,
    decay_lambda: f64,
    max_repeat_count: i32,
    call_ambiguous_regions: bool,
    recover_right_anchor: bool,
    emit_wildtype_active_regions: bool,
    trace_haplotype_alignment: bool,
}

impl ActiveRegionDetector {
    /// Default minimum count difference required to trigger a region.
    pub const DEFAULT_MINIMUM_DIFFERENCE: i32 = 5;
    /// Default quantile used to derive the difference threshold.
    pub const DEFAULT_DIFFERENCE_QUANTILE: f64 = 0.90;
    /// Default setting for requiring both region anchors.
    pub const DEFAULT_ANCHOR_BOTH_ENDS: bool = true;
    /// Default setting for including reverse-complement k-mer counts.
    pub const DEFAULT_COUNT_REVERSE_KMERS: bool = true;
    /// Default setting for permitting ambiguous active regions.
    pub const DEFAULT_CALL_AMBIGUOUS_REGIONS: bool = true;
    /// Default peak scan length.
    pub const DEFAULT_PEAK_SCAN_LENGTH: i32 = 7;
    /// Default multiplier for the scan limit.
    pub const DEFAULT_SCAN_LIMIT_FACTOR: f64 = 7.0;
    /// Default setting for retaining trace matrices on haplotypes.
    pub const DEFAULT_TRACE_HAPLOTYPE_ALIGNMENT: bool = false;
    /// Default lower bound for exponential recovery.
    pub const DEFAULT_DECAY_MINIMUM: f64 = 0.55;
    /// Default exponential decay alpha.
    pub const DEFAULT_DECAY_ALPHA: f64 = 0.80;
    /// Default repeat limit.
    pub const DEFAULT_MAX_REPEAT_COUNT: i32 = 0;
    /// Default setting for recovering a right anchor.
    pub const DEFAULT_RECOVER_RIGHT_ANCHOR: bool = true;
    /// Default setting for emitting wildtype active regions.
    pub const DEFAULT_EMIT_WILDTYPE_ACTIVE_REGIONS: bool = false;

    /// Creates a detector for a k-mer size.
    pub fn new(kmer_util: KmerUtil) -> Result<Self, ActiveRegionDetectorError> {
        let alignment_weight = AlignmentWeight::defaults();
        let scan_limit = scan_limit(
            &alignment_weight,
            kmer_util.k_size() as i32,
            Self::DEFAULT_SCAN_LIMIT_FACTOR,
        )?;
        Ok(Self {
            decay_lambda: decay_lambda(kmer_util.k_size(), Self::DEFAULT_DECAY_ALPHA),
            kmer_util,
            alignment_weight,
            anchor_both_ends: Self::DEFAULT_ANCHOR_BOTH_ENDS,
            count_reverse_kmers: Self::DEFAULT_COUNT_REVERSE_KMERS,
            minimum_difference: Self::DEFAULT_MINIMUM_DIFFERENCE,
            difference_quantile: Self::DEFAULT_DIFFERENCE_QUANTILE,
            peak_scan_length: Self::DEFAULT_PEAK_SCAN_LENGTH,
            scan_limit_factor: Self::DEFAULT_SCAN_LIMIT_FACTOR,
            scan_limit,
            decay_minimum: Self::DEFAULT_DECAY_MINIMUM,
            decay_alpha: Self::DEFAULT_DECAY_ALPHA,
            max_repeat_count: Self::DEFAULT_MAX_REPEAT_COUNT,
            call_ambiguous_regions: Self::DEFAULT_CALL_AMBIGUOUS_REGIONS,
            recover_right_anchor: Self::DEFAULT_RECOVER_RIGHT_ANCHOR,
            emit_wildtype_active_regions: Self::DEFAULT_EMIT_WILDTYPE_ACTIVE_REGIONS,
            trace_haplotype_alignment: Self::DEFAULT_TRACE_HAPLOTYPE_ALIGNMENT,
        })
    }

    /// Returns the k-mer utility used by the detector.
    #[must_use]
    pub fn kmer_util(&self) -> &KmerUtil {
        &self.kmer_util
    }

    /// Builds an active-region container for a reference and count map.
    pub fn get_active_regions(
        &self,
        ref_region: ReferenceRegion,
        counter: &dyn CountMap,
    ) -> Result<ActiveRegionContainer, ActiveRegionDetectorError> {
        let counts = self.get_counts(&ref_region, counter);
        ActiveRegionContainer::new(ref_region, None, &counts).map_err(Into::into)
    }

    /// Detects active regions for a reference using counts from a count map.
    pub fn detect_active_regions(
        &self,
        ref_region: &ReferenceRegion,
        counter: &dyn CountMap,
    ) -> Result<Vec<ActiveRegion>, ActiveRegionDetectorError> {
        self.detect_from_counts(ref_region, &self.get_counts(ref_region, counter))
    }

    /// Returns per-reference-k-mer counts for a region.
    #[must_use]
    pub fn get_counts(&self, ref_region: &ReferenceRegion, counter: &dyn CountMap) -> Vec<i32> {
        let k_size = self.kmer_util.k_size();
        if ref_region.sequence.len() < k_size {
            return Vec::new();
        }

        ref_region
            .sequence
            .windows(k_size)
            .map(|window| {
                let Ok(kmer) = self.kmer_util.encode(window) else {
                    return 0;
                };
                let mut count = counter.get(&kmer) as i32;
                if self.count_reverse_kmers {
                    count += counter.get(&self.kmer_util.reverse_complement(&kmer)) as i32;
                }
                count
            })
            .collect()
    }

    /// Detects active regions from precomputed reference k-mer counts.
    ///
    /// All candidate regions are accepted unconditionally. Use
    /// [`Self::detect_from_counts_with`] to drive the scan with a callback
    /// that decides whether to accept each candidate based on downstream
    /// haplotype assembly (matching Java `KestrelRunner.exec`'s behaviour of
    /// retrying overlapping regions when haplotype building yields 0 or
    /// wildtype-only haplotypes).
    pub fn detect_from_counts(
        &self,
        ref_region: &ReferenceRegion,
        ref_count: &[i32],
    ) -> Result<Vec<ActiveRegion>, ActiveRegionDetectorError> {
        self.detect_from_counts_with(ref_region, ref_count, |_region| {
            Ok::<bool, ActiveRegionDetectorError>(true)
        })
        .map_err(|err| match err {
            AcceptError::Detector(err) => err,
            AcceptError::Callback(err) => err,
        })
    }

    /// Detects active regions from precomputed reference k-mer counts and
    /// calls `accept` for each candidate. When `accept` returns `Ok(true)`
    /// the candidate region is emitted and the scan advances past the
    /// region; when `accept` returns `Ok(false)` the region is discarded
    /// and the scan retries from `ref_count_index + 1`. This mirrors Java
    /// `KestrelRunner.exec`'s `REF_SEARCH` loop, which only advances past a
    /// region when haplotype assembly produced at least one non-wildtype
    /// haplotype.
    pub fn detect_from_counts_with<F, E>(
        &self,
        ref_region: &ReferenceRegion,
        ref_count: &[i32],
        mut accept: F,
    ) -> Result<Vec<ActiveRegion>, AcceptError<E>>
    where
        F: FnMut(&ActiveRegion) -> Result<bool, E>,
    {
        let ref_count_size = ref_count.len();
        if ref_count_size < 2 {
            return Ok(Vec::new());
        }

        let diff_threshold = self.difference_threshold(ref_count) - 1;
        let diff_threshold_l = -diff_threshold;
        let mut regions = Vec::new();
        let mut ref_count_index = 1_usize;
        let mut count_l = ref_count[0];
        let mut last_region_end = 0_i32;
        let k_size = self.kmer_util.k_size() as i32;

        while ref_count_index < ref_count_size {
            let count_r = ref_count[ref_count_index];
            let count_diff = count_l - count_r;

            if count_diff > diff_threshold {
                let scanned = self
                    .scan_right(
                        ref_region,
                        ref_count,
                        ref_count_index,
                        count_l,
                        count_r,
                        diff_threshold,
                    )
                    .map_err(AcceptError::Detector)?;
                if let Some((start, end, next_index, next_count)) = scanned {
                    let candidate = self
                        .make_region(ref_region, ref_count, start, end)
                        .map_err(AcceptError::Detector)?;
                    if let Some(region) = candidate {
                        let accepted = accept(&region).map_err(AcceptError::Callback)?;
                        if accepted {
                            last_region_end = region.end_kmer_index;
                            regions.push(region);
                            ref_count_index = next_index;
                            count_l = next_count;
                        } else {
                            count_l = count_r;
                            ref_count_index += 1;
                        }
                    } else {
                        count_l = count_r;
                        ref_count_index += 1;
                    }
                } else {
                    count_l = count_r;
                    ref_count_index += 1;
                }
            } else if count_diff < diff_threshold_l {
                let scanned = self
                    .scan_left(
                        ref_region,
                        ref_count,
                        ref_count_index,
                        diff_threshold,
                        last_region_end,
                    )
                    .map_err(AcceptError::Detector)?;
                let (start, end) = match scanned {
                    LeftScanResult::SkipPeak {
                        next_index,
                        count_l: next_count_l,
                    } => {
                        count_l = next_count_l;
                        ref_count_index = next_index;
                        continue;
                    }
                    LeftScanResult::Candidate(Some(region)) => region,
                    LeftScanResult::Candidate(None) => {
                        count_l = count_r;
                        ref_count_index += 1;
                        continue;
                    }
                };
                if start < last_region_end && last_region_end > 0 {
                    count_l = count_r;
                    ref_count_index += 1;
                    continue;
                }
                let candidate = self
                    .make_region(ref_region, ref_count, start, end)
                    .map_err(AcceptError::Detector)?;
                let mut accepted_left_region = false;
                if let Some(region) = candidate {
                    let accepted = accept(&region).map_err(AcceptError::Callback)?;
                    if accepted {
                        regions.push(region);
                        accepted_left_region = true;
                    }
                }
                count_l = count_r;
                ref_count_index += 1;
                if accepted_left_region {
                    last_region_end = ref_count_index as i32;
                }
            } else {
                count_l = count_r;
                ref_count_index += 1;
            }

            if regions.last().is_some_and(|region| {
                region.right_end || region.end_kmer_index + k_size >= ref_count_size as i32 + k_size
            }) {
                break;
            }
        }

        regions.sort();
        Ok(regions)
    }

    fn scan_right(
        &self,
        ref_region: &ReferenceRegion,
        ref_count: &[i32],
        ref_count_index: usize,
        count_l: i32,
        count_r: i32,
        diff_threshold: i32,
    ) -> Result<Option<(i32, i32, usize, i32)>, ActiveRegionDetectorError> {
        let ref_count_size = ref_count.len();
        let k_size = self.kmer_util.k_size() as i32;
        let last_scan_index = ref_count_index.saturating_add(self.scan_limit as usize);
        let mut scan_end_index = ref_count_index + 1;
        let mut n_peak = 0_usize;
        let mut peak_scan_index = 0_usize;
        let mut last_valley_index = 0_usize;

        'scan_loop: while scan_end_index <= last_scan_index {
            let recovery_value = if self.decay_minimum == 1.0 {
                let recovery_value = (count_l - diff_threshold).max(1);
                while scan_end_index < ref_count_size && ref_count[scan_end_index] < recovery_value
                {
                    scan_end_index += 1;
                }
                recovery_value
            } else {
                while scan_end_index < ref_count_size {
                    let distance = scan_end_index as i32 - ref_count_index as i32;
                    let recovery = self.recovery_value(count_l, distance, diff_threshold);
                    if ref_count[scan_end_index] >= recovery {
                        break;
                    }
                    scan_end_index += 1;
                }
                let distance = scan_end_index as i32 - ref_count_index as i32;
                self.recovery_value(count_l, distance, diff_threshold)
            };

            if self.peak_scan_length == 0 {
                break;
            }

            if peak_scan_index > 0
                && scan_end_index.saturating_sub(peak_scan_index) >= self.kmer_util.k_size()
            {
                last_valley_index = scan_end_index;
            } else if peak_scan_index == 0
                && scan_end_index.saturating_sub(ref_count_index) >= self.kmer_util.k_size()
            {
                last_valley_index = scan_end_index;
            }

            peak_scan_index = scan_end_index;
            let peak_scan_limit =
                (scan_end_index + self.peak_scan_length as usize).min(ref_count_size);

            while peak_scan_index < peak_scan_limit {
                if ref_count[peak_scan_index] < recovery_value {
                    n_peak += 1;
                    scan_end_index = peak_scan_index;

                    if n_peak > 3
                        && (scan_end_index - ref_count_index) / n_peak < self.kmer_util.k_size()
                    {
                        if last_valley_index > 0 {
                            scan_end_index = last_valley_index;
                            break 'scan_loop;
                        }
                        return Ok(None);
                    }

                    continue 'scan_loop;
                }

                peak_scan_index += 1;
            }

            if peak_scan_index == ref_count_size && last_valley_index > 0 {
                scan_end_index = last_valley_index;
                break;
            }

            break;
        }

        if scan_end_index > last_scan_index {
            return Ok(None);
        }

        if scan_end_index < ref_count_size {
            if scan_end_index as i32 - (ref_count_index as i32) < k_size - 1 {
                return Ok(None);
            }
            let end_base = scan_end_index as i32 + k_size - 1;
            if !self.call_ambiguous_regions
                && ref_region.contains_ambiguous_by_index(ref_count_index as i32, end_base)?
            {
                return Ok(None);
            }
            Ok(Some((
                ref_count_index as i32 - 1,
                scan_end_index as i32,
                scan_end_index + 1,
                ref_count[scan_end_index],
            )))
        } else {
            let mut end = -1;
            if self.recover_right_anchor {
                let mut recovery_index = ref_count_index + k_size as usize;
                while recovery_index < ref_count_size {
                    if ref_count[recovery_index] - ref_count[recovery_index - 1] > diff_threshold {
                        end = recovery_index as i32;
                        break;
                    }
                    recovery_index += 1;
                }
            }
            if end == -1 {
                if self.anchor_both_ends
                    || ref_count_size - ref_count_index > self.scan_limit as usize
                {
                    return Ok(None);
                }
                if !self.call_ambiguous_regions
                    && ref_region.contains_ambiguous_by_index(
                        ref_count_index as i32,
                        ref_count_size as i32 - 1,
                    )?
                {
                    return Ok(None);
                }
                Ok(Some((
                    ref_count_index as i32 - 1,
                    -1,
                    ref_count_size,
                    count_r,
                )))
            } else {
                Ok(Some((
                    ref_count_index as i32 - 1,
                    end,
                    end as usize + 1,
                    ref_count[end as usize],
                )))
            }
        }
    }

    fn scan_left(
        &self,
        ref_region: &ReferenceRegion,
        ref_count: &[i32],
        ref_count_index: usize,
        diff_threshold: i32,
        last_region_end: i32,
    ) -> Result<LeftScanResult, ActiveRegionDetectorError> {
        let count_l = ref_count[ref_count_index - 1];
        let count_r = ref_count[ref_count_index];

        if self.peak_scan_length > 0 {
            let recovery_value = count_l + diff_threshold;
            let last_scan_index =
                (ref_count_index + self.peak_scan_length as usize).min(ref_count.len());
            for scan_end_index in (ref_count_index + 1)..last_scan_index {
                if ref_count[scan_end_index] <= recovery_value
                    && ref_count[ref_count_index] - ref_count[scan_end_index] < diff_threshold
                {
                    return Ok(LeftScanResult::SkipPeak {
                        next_index: scan_end_index + 1,
                        count_l: ref_count[scan_end_index],
                    });
                }
            }
        }

        if ref_count_index > self.scan_limit as usize {
            return Ok(LeftScanResult::Candidate(None));
        }

        let mut scan_end_index = ref_count_index as i32 - 1;
        let last_scan_index = last_region_end.max(0);
        while scan_end_index >= last_scan_index {
            let distance = ref_count_index as i32 - scan_end_index;
            let recovery = self.left_recovery_value(count_r, distance, diff_threshold);
            if ref_count[scan_end_index as usize] >= recovery {
                break;
            }
            scan_end_index -= 1;
        }

        if scan_end_index > 0 {
            return Ok(LeftScanResult::Candidate(None));
        }

        let mut start = -1;
        if self.recover_right_anchor && ref_count_index as i32 - self.scan_limit < 0 {
            let mut recovery_index =
                (ref_count_index as i32 - self.kmer_util.k_size() as i32).max(0);
            while recovery_index > 0 {
                if ref_count[recovery_index as usize - 1] - ref_count[recovery_index as usize]
                    > diff_threshold
                {
                    start = recovery_index;
                    break;
                }
                recovery_index -= 1;
            }
        }

        if start < 0 && self.anchor_both_ends {
            return Ok(LeftScanResult::Candidate(None));
        }
        if !self.call_ambiguous_regions {
            let contains_ambiguous = if start < 0 {
                ref_region.contains_ambiguous_by_index(0, ref_count_index as i32)?
            } else {
                ref_region.contains_ambiguous_by_index(start, ref_count_index as i32)?
            };
            if contains_ambiguous {
                return Ok(LeftScanResult::Candidate(None));
            }
        }

        Ok(LeftScanResult::Candidate(Some((
            start,
            ref_count_index as i32,
        ))))
    }

    fn recovery_value(&self, anchor_count: i32, distance: i32, diff_threshold: i32) -> i32 {
        if self.decay_minimum == 1.0 {
            return (anchor_count - diff_threshold).max(1);
        }
        let decay_minimum = (anchor_count as f64 * self.decay_minimum).max(1.0);
        let decay_range = anchor_count as f64 - decay_minimum;
        (decay_range * (-(distance as f64) * self.decay_lambda).exp() + decay_minimum) as i32
    }

    fn left_recovery_value(&self, anchor_count: i32, distance: i32, diff_threshold: i32) -> i32 {
        if self.decay_minimum == 1.0 {
            return (anchor_count - diff_threshold).max(1);
        }
        let decay_minimum = (anchor_count as f64 * self.decay_minimum).max(1.0);
        let decay_range = anchor_count as f64 - decay_minimum;
        (decay_range * ((distance as f64) * self.decay_lambda).exp() + decay_minimum) as i32
    }

    fn make_region(
        &self,
        ref_region: &ReferenceRegion,
        ref_count: &[i32],
        start: i32,
        end: i32,
    ) -> Result<Option<ActiveRegion>, ActiveRegionDetectorError> {
        match ActiveRegion::new(ref_region.clone(), start, end, ref_count, &self.kmer_util) {
            Ok(region) => Ok(Some(region)),
            Err(ActiveRegionError::AmbiguousEndKmer { .. }) => Ok(None),
            Err(error) => Err(error.into()),
        }
    }

    /// Computes the count-difference threshold used to decide whether a depth change is active.
    #[must_use]
    pub fn difference_threshold(&self, count: &[i32]) -> i32 {
        let threshold = if self.difference_quantile > 0.0 && count.len() > 2 {
            active_region_count_diff_quantile(count, self.difference_quantile)
                .max(self.minimum_difference)
        } else {
            self.minimum_difference
        };
        threshold.max(1)
    }

    /// Sets whether detected active regions must have both left and right anchors.
    pub fn set_anchor_both_ends(&mut self, anchor_both_ends: bool) {
        self.anchor_both_ends = anchor_both_ends;
    }

    /// Returns whether detected active regions must have both left and right anchors.
    #[must_use]
    pub fn anchor_both_ends(&self) -> bool {
        self.anchor_both_ends
    }

    /// Sets whether reverse-complement k-mers contribute to reference k-mer depth.
    pub fn set_count_reverse_kmers(&mut self, count_reverse_kmers: bool) {
        self.count_reverse_kmers = count_reverse_kmers;
    }

    /// Returns whether reverse-complement k-mers contribute to reference k-mer depth.
    #[must_use]
    pub fn count_reverse_kmers(&self) -> bool {
        self.count_reverse_kmers
    }

    /// Sets the minimum depth difference required before a depth change can be active.
    pub fn set_minimum_difference(
        &mut self,
        minimum_difference: i32,
    ) -> Result<(), ActiveRegionDetectorError> {
        if minimum_difference < 1 {
            return Err(ActiveRegionDetectorError::InvalidMinimumDifference(
                minimum_difference,
            ));
        }
        self.minimum_difference = minimum_difference;
        Ok(())
    }

    /// Returns the minimum depth difference required before a depth change can be active.
    #[must_use]
    pub fn minimum_difference(&self) -> i32 {
        self.minimum_difference
    }

    /// Sets the count-difference quantile used to raise the active-region threshold.
    pub fn set_difference_quantile(
        &mut self,
        difference_quantile: f64,
    ) -> Result<(), ActiveRegionDetectorError> {
        if !(0.0..1.0).contains(&difference_quantile) {
            return Err(ActiveRegionDetectorError::InvalidDifferenceQuantile(
                difference_quantile,
            ));
        }
        self.difference_quantile = difference_quantile;
        Ok(())
    }

    /// Returns the count-difference quantile used to raise the active-region threshold.
    #[must_use]
    pub fn difference_quantile(&self) -> f64 {
        self.difference_quantile
    }

    /// Sets the minimum fractional recovery target for exponential depth recovery.
    pub fn set_decay_minimum(
        &mut self,
        decay_minimum: f64,
    ) -> Result<(), ActiveRegionDetectorError> {
        if !(0.0..=1.0).contains(&decay_minimum) {
            return Err(ActiveRegionDetectorError::InvalidDecayMinimum(
                decay_minimum,
            ));
        }
        self.decay_minimum = decay_minimum;
        Ok(())
    }

    /// Returns the minimum fractional recovery target for exponential depth recovery.
    #[must_use]
    pub fn decay_minimum(&self) -> f64 {
        self.decay_minimum
    }

    /// Sets the exponential decay alpha used while scanning for depth recovery.
    pub fn set_decay_alpha(&mut self, decay_alpha: f64) -> Result<(), ActiveRegionDetectorError> {
        if decay_alpha <= 0.0 || decay_alpha >= 1.0 {
            return Err(ActiveRegionDetectorError::InvalidDecayAlpha(decay_alpha));
        }
        self.decay_alpha = decay_alpha;
        self.decay_lambda = decay_lambda(self.kmer_util.k_size(), decay_alpha);
        Ok(())
    }

    /// Returns the exponential decay alpha used while scanning for depth recovery.
    #[must_use]
    pub fn decay_alpha(&self) -> f64 {
        self.decay_alpha
    }

    /// Returns the derived exponential decay lambda for the detector's k-mer size.
    #[must_use]
    pub fn decay_lambda(&self) -> f64 {
        self.decay_lambda
    }

    /// Sets how far the detector scans through local depth peaks after an apparent recovery.
    pub fn set_peak_scan_length(
        &mut self,
        peak_scan_length: i32,
    ) -> Result<(), ActiveRegionDetectorError> {
        if peak_scan_length < 0 {
            return Err(ActiveRegionDetectorError::NegativePeakScanLength(
                peak_scan_length,
            ));
        }
        self.peak_scan_length = peak_scan_length;
        Ok(())
    }

    /// Returns how far the detector scans through local depth peaks after an apparent recovery.
    #[must_use]
    pub fn peak_scan_length(&self) -> i32 {
        self.peak_scan_length
    }

    /// Sets the multiplier used to derive the maximum active-region scan length.
    pub fn set_scan_limit_factor(
        &mut self,
        scan_limit_factor: f64,
    ) -> Result<(), ActiveRegionDetectorError> {
        if scan_limit_factor < 0.0 {
            return Err(ActiveRegionDetectorError::NegativeScanLimitFactor(
                scan_limit_factor,
            ));
        }
        self.scan_limit_factor = scan_limit_factor;
        self.scan_limit = scan_limit(
            &self.alignment_weight,
            self.kmer_util.k_size() as i32,
            scan_limit_factor,
        )?;
        Ok(())
    }

    /// Returns the multiplier used to derive the maximum active-region scan length.
    #[must_use]
    pub fn scan_limit_factor(&self) -> f64 {
        self.scan_limit_factor
    }

    /// Returns the derived maximum active-region scan length in bases.
    #[must_use]
    pub fn scan_limit_length(&self) -> i32 {
        self.scan_limit
    }

    /// Sets alignment weights and refreshes the derived active-region scan limit.
    pub fn set_alignment_weight(
        &mut self,
        alignment_weight: AlignmentWeight,
    ) -> Result<(), ActiveRegionDetectorError> {
        self.alignment_weight = alignment_weight;
        self.scan_limit = scan_limit(
            &self.alignment_weight,
            self.kmer_util.k_size() as i32,
            self.scan_limit_factor,
        )?;
        Ok(())
    }

    /// Returns the alignment weights used to derive scan limits.
    #[must_use]
    pub fn alignment_weight(&self) -> &AlignmentWeight {
        &self.alignment_weight
    }

    /// Sets the maximum repeat count allowed while assembling region haplotypes.
    pub fn set_max_repeat_count(
        &mut self,
        max_repeat_count: i32,
    ) -> Result<(), ActiveRegionDetectorError> {
        if max_repeat_count < 0 {
            return Err(ActiveRegionDetectorError::NegativeMaxRepeatCount(
                max_repeat_count,
            ));
        }
        self.max_repeat_count = max_repeat_count;
        Ok(())
    }

    /// Returns the maximum repeat count allowed while assembling region haplotypes.
    #[must_use]
    pub fn max_repeat_count(&self) -> i32 {
        self.max_repeat_count
    }

    /// Sets whether active regions containing ambiguous bases may be emitted.
    pub fn set_call_ambiguous_regions(&mut self, call_ambiguous_regions: bool) {
        self.call_ambiguous_regions = call_ambiguous_regions;
    }

    /// Returns whether active regions containing ambiguous bases may be emitted.
    #[must_use]
    pub fn call_ambiguous_regions(&self) -> bool {
        self.call_ambiguous_regions
    }

    /// Sets whether the detector attempts to recover a missing right anchor.
    pub fn set_recover_right_anchor(&mut self, recover_right_anchor: bool) {
        self.recover_right_anchor = recover_right_anchor;
    }

    /// Returns whether the detector attempts to recover a missing right anchor.
    #[must_use]
    pub fn recover_right_anchor(&self) -> bool {
        self.recover_right_anchor
    }

    /// Sets whether wildtype active regions are emitted.
    pub fn set_emit_wildtype_active_regions(&mut self, emit_wildtype_active_regions: bool) {
        self.emit_wildtype_active_regions = emit_wildtype_active_regions;
    }

    /// Returns whether wildtype active regions are emitted.
    #[must_use]
    pub fn emit_wildtype_active_regions(&self) -> bool {
        self.emit_wildtype_active_regions
    }

    /// Sets whether haplotype trace matrices are retained.
    pub fn set_trace_haplotype_alignment(&mut self, trace_haplotype_alignment: bool) {
        self.trace_haplotype_alignment = trace_haplotype_alignment;
    }

    /// Returns whether haplotype trace matrices are retained.
    #[must_use]
    pub fn trace_haplotype_alignment(&self) -> bool {
        self.trace_haplotype_alignment
    }
}

fn active_region_count_diff_quantile(count: &[i32], quantile: f64) -> i32 {
    match count.len() {
        0 | 1 => 0,
        2 => (count[1] - count[0]).abs(),
        len => {
            let mut diffs = Vec::with_capacity(len - 1);
            let mut last_count = count[0];
            for this_count in count.iter().take(len - 1) {
                diffs.push((last_count - *this_count).abs());
                last_count = *this_count;
            }
            diffs.sort_unstable();

            let n_count = len - 2;
            let raw = n_count as f64 * quantile;
            let loc = raw as usize;
            let offset = raw - loc as f64;
            let value = diffs[loc] as f64 * (1.0 - offset) + diffs[loc + 1] as f64 * offset;

            value as i32
        }
    }
}

fn scan_limit(
    alignment_weight: &AlignmentWeight,
    k_size: i32,
    scan_limit_factor: f64,
) -> Result<i32, AlignmentWeightError> {
    let max_gap_size = alignment_weight.max_exclusive_gap_size(k_size)?;
    let scan_limit = max_gap_size as f64 + scan_limit_factor * k_size as f64;
    let scan_limit = if scan_limit > i32::MAX as f64 {
        i32::MAX
    } else {
        scan_limit as i32
    };
    Ok(scan_limit.max(k_size))
}

fn decay_lambda(k_size: usize, decay_alpha: f64) -> f64 {
    -decay_alpha.ln() / k_size as f64
}

/// A reference interval where k-mer depth suggests a possible variant haplotype.
#[derive(Clone, Debug, PartialEq)]
pub struct ActiveRegion {
    /// Reference region containing this active interval.
    pub ref_region: ReferenceRegion,
    /// Inclusive start base index within `ref_region`.
    pub start_index: i32,
    /// Inclusive end base index within `ref_region`.
    pub end_index: i32,
    /// Left anchor k-mer index, or `0` when the left end is open.
    pub start_kmer_index: i32,
    /// Right anchor k-mer index, or the final k-mer index when the right end is open.
    pub end_kmer_index: i32,
    /// True when this region is open at the left end of the reference.
    pub left_end: bool,
    /// True when this region is open at the right end of the reference.
    pub right_end: bool,
    /// Summary depth statistics across the region.
    pub stats: RegionStats,
    left_end_kmer: Option<KmerKey>,
    right_end_kmer: Option<KmerKey>,
}

impl ActiveRegion {
    /// Creates an active region from left and right anchor k-mer indexes.
    pub fn new(
        ref_region: ReferenceRegion,
        start_kmer_index: i32,
        end_kmer_index: i32,
        count: &[i32],
        kmer_util: &KmerUtil,
    ) -> Result<Self, ActiveRegionError> {
        if start_kmer_index < -1 {
            return Err(ActiveRegionError::StartTooNegative(start_kmer_index));
        }
        if end_kmer_index < -1 || (end_kmer_index >= 0 && end_kmer_index as usize >= count.len()) {
            return Err(ActiveRegionError::EndOutOfRange {
                end: end_kmer_index,
                last: count.len() as i32 - 1,
            });
        }
        if end_kmer_index > ref_region.size - kmer_util.k_size() as i32 {
            return Err(ActiveRegionError::EndExtendsPastReference {
                end: end_kmer_index,
            });
        }
        if start_kmer_index >= 0 && end_kmer_index >= 0 {
            if start_kmer_index >= end_kmer_index {
                return Err(ActiveRegionError::StartNotBeforeEnd {
                    start: start_kmer_index,
                    end: end_kmer_index,
                });
            }
        } else if start_kmer_index == -1 && end_kmer_index == -1 {
            return Err(ActiveRegionError::BothEndsOpen);
        }

        let (left_end, start_index, start_kmer_index, left_end_kmer) = if start_kmer_index == -1 {
            (true, 0, 0, None)
        } else {
            let kmer = end_kmer(&ref_region, start_kmer_index, kmer_util).ok_or(
                ActiveRegionError::AmbiguousEndKmer {
                    side: "left",
                    index: start_kmer_index,
                },
            )?;
            (false, start_kmer_index, start_kmer_index, Some(kmer))
        };

        let (right_end, end_index, end_kmer_index, right_end_kmer) = if end_kmer_index == -1 {
            (true, ref_region.size - 1, count.len() as i32 - 1, None)
        } else {
            let kmer = end_kmer(&ref_region, end_kmer_index, kmer_util).ok_or(
                ActiveRegionError::AmbiguousEndKmer {
                    side: "right",
                    index: end_kmer_index,
                },
            )?;
            (
                false,
                end_kmer_index + kmer_util.k_size() as i32 - 1,
                end_kmer_index,
                Some(kmer),
            )
        };

        let stats = RegionStats::from_counts(count, start_kmer_index, end_kmer_index)?;

        Ok(Self {
            ref_region,
            start_index,
            end_index,
            start_kmer_index,
            end_kmer_index,
            left_end,
            right_end,
            stats,
            left_end_kmer,
            right_end_kmer,
        })
    }

    /// Returns the encoded left-end anchor k-mer, when the left end is anchored.
    #[must_use]
    pub fn left_end_kmer(&self) -> Option<KmerKey> {
        self.left_end_kmer.clone()
    }

    /// Returns the encoded right-end anchor k-mer, when the right end is anchored.
    #[must_use]
    pub fn right_end_kmer(&self) -> Option<KmerKey> {
        self.right_end_kmer.clone()
    }

    /// Returns true if `kmer` matches this region's left-end anchor.
    #[must_use]
    pub fn match_left_end(&self, kmer: Option<&KmerKey>) -> bool {
        matches!((kmer, &self.left_end_kmer), (Some(a), Some(b)) if a == b)
    }

    /// Returns true if `kmer` matches this region's right-end anchor.
    #[must_use]
    pub fn match_right_end(&self, kmer: Option<&KmerKey>) -> bool {
        matches!((kmer, &self.right_end_kmer), (Some(a), Some(b)) if a == b)
    }
}

impl Eq for ActiveRegion {}

impl Ord for ActiveRegion {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.ref_region
            .cmp(&other.ref_region)
            .then_with(|| self.start_index.cmp(&other.start_index))
            .then_with(|| self.end_index.cmp(&other.end_index))
    }
}

impl PartialOrd for ActiveRegion {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl fmt::Display for ActiveRegion {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "ActiveRegion[name={}, start={}, end={}]",
            self.ref_region.name, self.start_index, self.end_index
        )
    }
}

fn end_kmer(ref_region: &ReferenceRegion, index: i32, kmer_util: &KmerUtil) -> Option<KmerKey> {
    let index = usize::try_from(index).ok()?;
    let end = index.checked_add(kmer_util.k_size())?;
    let sequence = ref_region.sequence.get(index..end)?;
    kmer_util.encode(sequence).ok()
}

/// A candidate haplotype sequence aligned to an active region.
#[derive(Clone, Debug, PartialEq)]
pub struct Haplotype {
    /// Active region this haplotype covers.
    pub active_region: ActiveRegion,
    /// Haplotype consensus sequence.
    pub sequence: Vec<u8>,
    /// Haplotype sequence length in bases.
    pub length: usize,
    /// Depth statistics for this haplotype.
    pub stats: RegionStats,
    /// Best alignment for this haplotype.
    pub alignment: AlignNode,
    alignment_list: Vec<AlignNode>,
    /// Number of alternate alignments retained.
    pub n_align: usize,
    /// Alignment score assigned to this haplotype.
    pub alignment_score: f32,
    /// Optional dynamic-programming trace matrix for diagnostics.
    pub trace_matrix: Option<TraceMatrix>,
}

impl Haplotype {
    /// Creates a haplotype and selects the best alignment from `alignment_list`.
    pub fn new(
        sequence: impl Into<Vec<u8>>,
        active_region: ActiveRegion,
        alignment_list: Vec<AlignNode>,
        alignment_score: f32,
        trace_matrix: Option<TraceMatrix>,
        stats: RegionStats,
    ) -> Result<Self, ActiveRegionError> {
        if alignment_list.is_empty() {
            return Err(ActiveRegionError::EmptyAlignmentList);
        }

        let sequence = sequence.into();
        let length = sequence.len();
        let mut alignment_list = alignment_list;
        let alignment = alignment_list[0].clone();
        alignment_list.sort_by(|left, right| left.compare_to(right).cmp(&0));
        let n_align = alignment_list.len();

        Ok(Self {
            active_region,
            sequence,
            length,
            stats,
            alignment,
            alignment_list,
            n_align,
            alignment_score,
            trace_matrix,
        })
    }

    /// Renders the requested alignment as reference, match-bar, and haplotype strings.
    pub fn alignment_string(&self, alignment_index: i32) -> Result<[String; 3], ActiveRegionError> {
        if alignment_index < 0 || alignment_index as usize >= self.n_align {
            return Err(ActiveRegionError::AlignmentIndexOutOfRange {
                index: alignment_index,
                last: self.n_align as i32 - 1,
            });
        }

        let mut ref_builder = String::new();
        let mut bar_builder = String::new();
        let mut con_builder = String::new();
        let mut ref_index = self.active_region.start_index as usize;
        let mut con_index = 0_usize;
        let mut node = Some(&self.alignment_list[alignment_index as usize]);

        while let Some(current) = node {
            for _ in 0..current.n {
                match current.node_type {
                    AlignNode::MATCH => {
                        ref_builder.push(char::from(
                            self.active_region.ref_region.sequence[ref_index],
                        ));
                        con_builder.push(char::from(self.sequence[con_index]));
                        bar_builder.push('|');
                        ref_index += 1;
                        con_index += 1;
                    }
                    AlignNode::MISMATCH => {
                        ref_builder.push(char::from(
                            self.active_region.ref_region.sequence[ref_index],
                        ));
                        con_builder.push(char::from(self.sequence[con_index]));
                        bar_builder.push(' ');
                        ref_index += 1;
                        con_index += 1;
                    }
                    AlignNode::INS => {
                        ref_builder.push('-');
                        con_builder.push(char::from(self.sequence[con_index]));
                        bar_builder.push(' ');
                        con_index += 1;
                    }
                    AlignNode::DEL => {
                        ref_builder.push(char::from(
                            self.active_region.ref_region.sequence[ref_index],
                        ));
                        con_builder.push('-');
                        bar_builder.push(' ');
                        ref_index += 1;
                    }
                    _ => unreachable!("AlignNode constructor validates node type"),
                }
            }
            node = current.next.as_deref();
        }

        Ok([ref_builder, bar_builder, con_builder])
    }

    /// Returns true when the best alignment is a full-length reference match.
    #[must_use]
    pub fn is_wildtype(&self) -> bool {
        self.alignment.node_type == AlignNode::MATCH && self.alignment.next.is_none()
    }

    /// Returns all alignments retained for this haplotype, sorted best first.
    #[must_use]
    pub fn alignment_list(&self) -> &[AlignNode] {
        &self.alignment_list
    }
}

impl fmt::Display for Haplotype {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "Haplotype[length={}, min_depth={}, region={}]",
            self.length, self.stats.min, self.active_region
        )
    }
}

/// Haplotype calls associated with one active region.
#[derive(Clone, Debug, PartialEq)]
pub struct RegionHaplotype {
    /// Active region covered by the haplotypes.
    pub active_region: ActiveRegion,
    /// Candidate haplotypes for the region.
    pub haplotypes: Vec<Haplotype>,
    /// Combined minimum depth for the region and its haplotypes.
    pub min_depth: i32,
    /// True when all haplotypes are wildtype.
    pub is_wildtype: bool,
}

impl RegionHaplotype {
    /// Creates a region haplotype bundle.
    pub fn new(
        active_region: ActiveRegion,
        haplotypes: Vec<Haplotype>,
    ) -> Result<Self, ActiveRegionError> {
        if haplotypes.is_empty() {
            return Err(ActiveRegionError::EmptyHaplotypeList);
        }

        let mut min_depth = active_region.stats.min;
        let mut is_wildtype = true;

        for haplotype in &haplotypes {
            min_depth += haplotype.stats.min;
            if !haplotype.is_wildtype() {
                is_wildtype = false;
            }
        }

        Ok(Self {
            active_region,
            haplotypes,
            min_depth,
            is_wildtype,
        })
    }
}

impl Eq for RegionHaplotype {}

impl Ord for RegionHaplotype {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.active_region.cmp(&other.active_region)
    }
}

impl PartialOrd for RegionHaplotype {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl fmt::Display for RegionHaplotype {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "RegionHaplotype[region={}, haplotypes={}, wildtype={}, mindepth={}]",
            self.active_region,
            self.haplotypes.len(),
            self.is_wildtype,
            self.min_depth
        )
    }
}

/// Active-region output for one reference region.
#[derive(Clone, Debug, PartialEq)]
pub struct ActiveRegionContainer {
    /// Reference region that was scanned.
    pub ref_region: ReferenceRegion,
    /// Haplotype calls discovered in the reference region.
    pub haplotypes: Vec<RegionHaplotype>,
    /// Depth statistics across the scanned reference k-mers.
    pub stats: RegionStats,
}

impl ActiveRegionContainer {
    /// Creates a container from an optional list of region haplotypes and reference counts.
    pub fn new(
        ref_region: ReferenceRegion,
        haplotypes: Option<&[RegionHaplotype]>,
        count: &[i32],
    ) -> Result<Self, RegionStatsError> {
        let stats = RegionStats::from_counts(count, 0, count.len() as i32)?;

        Ok(Self {
            ref_region,
            haplotypes: haplotypes.map_or_else(Vec::new, <[RegionHaplotype]>::to_vec),
            stats,
        })
    }
}

/// Summary statistics for k-mer depths in a region.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RegionStats {
    /// Minimum observed depth.
    pub min: i32,
    /// 25th percentile depth.
    pub pct25: f32,
    /// Median depth.
    pub pct50: f32,
    /// 75th percentile depth.
    pub pct75: f32,
    /// Maximum observed depth.
    pub max: i32,
    /// Number of depth observations.
    pub n: usize,
}

impl RegionStats {
    /// Computes summary statistics over `count[start..end]`.
    pub fn from_counts(count: &[i32], start: i32, end: i32) -> Result<Self, RegionStatsError> {
        if start < 0 {
            return Err(RegionStatsError::NegativeStart(start));
        }

        let count_slice_size = end - start;
        if end as usize > count.len() || count_slice_size < 1 {
            return Err(RegionStatsError::InvalidEnd {
                length: count.len(),
                range: count_slice_size,
                end,
            });
        }

        let start = start as usize;
        let end = end as usize;
        if count_slice_size == 1 {
            let value = count[start];
            return Ok(Self {
                min: value,
                pct25: value as f32,
                pct50: value as f32,
                pct75: value as f32,
                max: value,
                n: 1,
            });
        }

        let mut count_slice = count[start..end].to_vec();
        count_slice.sort_unstable();
        let n_less_one = count_slice.len() - 1;

        Ok(Self {
            min: count_slice[0],
            pct25: percentile(&count_slice, n_less_one, 0.25),
            pct50: percentile(&count_slice, n_less_one, 0.50),
            pct75: percentile(&count_slice, n_less_one, 0.75),
            max: count_slice[count_slice.len() - 1],
            n: count_slice.len(),
        })
    }
}

impl fmt::Display for RegionStats {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "Stats[min={}, 25th={:.6}, 50th={:.6}, 75th={:.6}, max={}, n={}]",
            self.min, self.pct25, self.pct50, self.pct75, self.max, self.n
        )
    }
}

fn percentile(count_slice: &[i32], n_less_one: usize, q: f32) -> f32 {
    let scaled = n_less_one as f32 * q;
    let loc = scaled as usize;
    let offset = scaled - loc as f32;

    count_slice[loc] as f32 * (1.0 - offset) + count_slice[loc + 1] as f32 * offset
}

#[cfg(test)]
mod tests {
    use std::collections::HashMap;

    use crate::counter::{CountMap, CountMapError};
    use crate::io::InputSample;
    use crate::refreader::{ReferenceRegion, ReferenceSequence};
    use crate::util::digest::Digest;

    use super::*;

    const EPS: f32 = 0.01;

    #[test]
    fn single_element_has_all_stats_equal_to_value() {
        let stats = RegionStats::from_counts(&[42], 0, 1).unwrap();
        assert_eq!(stats.min, 42);
        assert_eq!(stats.max, 42);
        assert_eq!(stats.pct25, 42.0);
        assert_eq!(stats.pct50, 42.0);
        assert_eq!(stats.pct75, 42.0);
        assert_eq!(stats.n, 1);
    }

    #[test]
    fn two_elements_min_max() {
        let stats = RegionStats::from_counts(&[1, 10], 0, 2).unwrap();
        assert_eq!(stats.min, 1);
        assert_eq!(stats.max, 10);
        assert_eq!(stats.n, 2);
    }

    #[test]
    fn quartile_interpolation() {
        let stats = RegionStats::from_counts(&[1, 2, 3, 4, 5], 0, 5).unwrap();
        assert_eq!(stats.min, 1);
        assert_eq!(stats.max, 5);
        assert_eq!(stats.n, 5);
        assert!((stats.pct25 - 2.0).abs() < EPS);
        assert!((stats.pct50 - 3.0).abs() < EPS);
        assert!((stats.pct75 - 4.0).abs() < EPS);
    }

    #[test]
    fn sorts_before_stats_and_respects_slice() {
        let stats = RegionStats::from_counts(&[5, 1, 4, 2, 3], 0, 5).unwrap();
        assert_eq!(stats.min, 1);
        assert_eq!(stats.max, 5);
        assert!((stats.pct50 - 3.0).abs() < EPS);

        let sliced = RegionStats::from_counts(&[99, 99, 1, 2, 3, 99], 2, 5).unwrap();
        assert_eq!(sliced.min, 1);
        assert_eq!(sliced.max, 3);
        assert_eq!(sliced.n, 3);
    }

    #[test]
    fn invalid_ranges_error_like_java() {
        assert_eq!(
            RegionStats::from_counts(&[1, 2], -1, 1),
            Err(RegionStatsError::NegativeStart(-1))
        );
        assert_eq!(
            RegionStats::from_counts(&[1, 2], 0, 10),
            Err(RegionStatsError::InvalidEnd {
                length: 2,
                range: 10,
                end: 10
            })
        );
        assert_eq!(
            RegionStats::from_counts(&[1, 2, 3], 1, 1),
            Err(RegionStatsError::InvalidEnd {
                length: 3,
                range: 0,
                end: 1
            })
        );
    }

    #[test]
    fn display_includes_all_fields() {
        let display = RegionStats::from_counts(&[1, 2, 3, 4, 5], 0, 5)
            .unwrap()
            .to_string();
        assert!(display.contains("min="));
        assert!(display.contains("max="));
        assert!(display.contains("n="));
    }

    #[test]
    fn active_region_constructs_interior_and_end_regions() {
        let count = vec![10; 13];
        let util = KmerUtil::new(4).unwrap();
        let active = ActiveRegion::new(ref16(), 2, 10, &count, &util).unwrap();
        assert!(!active.left_end);
        assert!(!active.right_end);
        assert_eq!(active.start_kmer_index, 2);
        assert_eq!(active.end_kmer_index, 10);
        assert_eq!(active.start_index, 2);
        assert_eq!(active.end_index, 13);
        assert_eq!(active.stats.n, 8);

        let left = ActiveRegion::new(ref16(), -1, 8, &count, &util).unwrap();
        assert!(left.left_end);
        assert!(!left.right_end);
        assert_eq!(left.start_index, 0);
        assert_eq!(left.start_kmer_index, 0);
        assert_eq!(left.left_end_kmer(), None);

        let right = ActiveRegion::new(ref16(), 2, -1, &count, &util).unwrap();
        assert!(!right.left_end);
        assert!(right.right_end);
        assert_eq!(right.end_kmer_index, 12);
        assert_eq!(right.end_index, 15);
        assert_eq!(right.right_end_kmer(), None);
    }

    #[test]
    fn active_region_kmers_are_cloned_and_match() {
        let count = vec![1; 13];
        let util = KmerUtil::new(4).unwrap();
        let active = ActiveRegion::new(ref16(), 2, 10, &count, &util).unwrap();

        let left_a = active.left_end_kmer().unwrap();
        let left_b = active.left_end_kmer().unwrap();
        assert_eq!(left_a, left_b);
        assert!(active.match_left_end(Some(&left_a)));
        assert!(active.match_right_end(active.right_end_kmer().as_ref()));
        assert!(!active.match_left_end(None));
    }

    #[test]
    fn active_region_validates_indices_and_ambiguous_end_kmers() {
        let count = vec![1; 13];
        let util = KmerUtil::new(4).unwrap();
        assert!(matches!(
            ActiveRegion::new(ref16(), -2, 5, &count, &util),
            Err(ActiveRegionError::StartTooNegative(-2))
        ));
        assert!(matches!(
            ActiveRegion::new(ref16(), 0, 100, &count, &util),
            Err(ActiveRegionError::EndOutOfRange { .. })
        ));
        assert!(matches!(
            ActiveRegion::new(ref16(), 5, 5, &count, &util),
            Err(ActiveRegionError::StartNotBeforeEnd { .. })
        ));
        assert!(matches!(
            ActiveRegion::new(ref16(), 8, 5, &count, &util),
            Err(ActiveRegionError::StartNotBeforeEnd { .. })
        ));
        assert_eq!(
            ActiveRegion::new(ref16(), -1, -1, &count, &util),
            Err(ActiveRegionError::BothEndsOpen)
        );

        let ambiguous = ref_region(b"ANAACCCCGGGGTTTT");
        assert!(matches!(
            ActiveRegion::new(ambiguous, 0, 8, &count, &util),
            Err(ActiveRegionError::AmbiguousEndKmer { side: "left", .. })
        ));
    }

    #[test]
    fn active_region_orders_and_displays() {
        let count = vec![1; 13];
        let util = KmerUtil::new(4).unwrap();
        let early = ActiveRegion::new(ref16(), 1, 5, &count, &util).unwrap();
        let late = ActiveRegion::new(ref16(), 6, 10, &count, &util).unwrap();
        assert!(early < late);
        assert!(late > early);
        assert!(!early.to_string().is_empty());
    }

    #[test]
    fn haplotype_stores_fields_and_detects_wildtype() {
        let active_region = ar_default();
        let alignment = AlignNode::new(AlignNode::MATCH, 16, None).unwrap();
        let haplotype = Haplotype::new(
            b"AAAACCCCGGGGTTTT".to_vec(),
            active_region.clone(),
            vec![alignment.clone()],
            75.0,
            None,
            stats(),
        )
        .unwrap();

        assert_eq!(haplotype.active_region, active_region);
        assert_eq!(haplotype.sequence, b"AAAACCCCGGGGTTTT");
        assert_eq!(haplotype.length, 16);
        assert_eq!(haplotype.alignment_score, 75.0);
        assert_eq!(haplotype.trace_matrix, None);
        assert_eq!(haplotype.alignment, alignment);
        assert_eq!(haplotype.n_align, 1);
        assert!(haplotype.is_wildtype());
    }

    #[test]
    fn haplotype_wildtype_requires_single_match_node() {
        let tail = AlignNode::new(AlignNode::MATCH, 8, None).unwrap();
        let head = AlignNode::new(AlignNode::MATCH, 8, Some(Box::new(tail))).unwrap();
        let haplotype = Haplotype::new(
            b"AAAACCCCGGGGTTTT".to_vec(),
            ar_default(),
            vec![head],
            75.0,
            None,
            stats(),
        )
        .unwrap();

        assert!(!haplotype.is_wildtype());
    }

    #[test]
    fn haplotype_rejects_empty_alignment_list() {
        assert_eq!(
            Haplotype::new(
                b"ACGT".to_vec(),
                ar_default(),
                Vec::new(),
                0.0,
                None,
                stats()
            ),
            Err(ActiveRegionError::EmptyAlignmentList)
        );
    }

    #[test]
    fn haplotype_primary_alignment_preserves_java_unsorted_constructor_quirk() {
        let earlier = AlignNode::new(
            AlignNode::MATCH,
            1,
            Some(Box::new(AlignNode::new(AlignNode::INS, 1, None).unwrap())),
        )
        .unwrap();
        let later = AlignNode::new(
            AlignNode::MATCH,
            2,
            Some(Box::new(AlignNode::new(AlignNode::INS, 1, None).unwrap())),
        )
        .unwrap();

        assert!(earlier.compare_to(&later) < 0);

        let haplotype = Haplotype::new(
            b"AAAACCCCGGGGTTTT".to_vec(),
            ar_default(),
            vec![later.clone(), earlier.clone()],
            75.0,
            None,
            stats(),
        )
        .unwrap();

        assert_eq!(haplotype.alignment, later);
        assert_eq!(haplotype.alignment_list()[0], earlier);
    }

    #[test]
    fn haplotype_alignment_string_and_bounds_match_java_tests() {
        let active_region = ar_default();
        let n = active_region.end_index - active_region.start_index + 1;
        let alignment = AlignNode::new(AlignNode::MATCH, n, None).unwrap();
        let haplotype = Haplotype::new(
            b"AAAACCCCGGGGTTTT".to_vec(),
            active_region,
            vec![alignment],
            0.0,
            None,
            stats(),
        )
        .unwrap();

        let rows = haplotype.alignment_string(0).unwrap();
        assert_eq!(rows.len(), 3);
        assert!(rows[1].chars().all(|base| base == '|'));
        assert_eq!(
            haplotype.alignment_string(-1),
            Err(ActiveRegionError::AlignmentIndexOutOfRange { index: -1, last: 0 })
        );
        assert_eq!(
            haplotype.alignment_string(99),
            Err(ActiveRegionError::AlignmentIndexOutOfRange { index: 99, last: 0 })
        );
    }

    #[test]
    fn haplotype_alignment_string_renders_mismatch_insertion_and_deletion() {
        let del = AlignNode::new(AlignNode::DEL, 1, None).unwrap();
        let ins = AlignNode::new(AlignNode::INS, 1, Some(Box::new(del))).unwrap();
        let mismatch = AlignNode::new(AlignNode::MISMATCH, 1, Some(Box::new(ins))).unwrap();
        let match_node = AlignNode::new(AlignNode::MATCH, 1, Some(Box::new(mismatch))).unwrap();
        let haplotype = Haplotype::new(
            b"ATG".to_vec(),
            ar_default(),
            vec![match_node],
            0.0,
            None,
            stats(),
        )
        .unwrap();

        let rows = haplotype.alignment_string(0).unwrap();
        assert_eq!(rows[0], "AA-C");
        assert_eq!(rows[1], "|   ");
        assert_eq!(rows[2], "ATG-");
        assert!(!haplotype.to_string().is_empty());
    }

    #[test]
    fn region_haplotype_stores_fields_and_sums_min_depth() {
        let active_region = ar_default();
        let haplotype = make_haplotype(
            &active_region,
            AlignNode::new(AlignNode::MATCH, 16, None).unwrap(),
        );
        let region_haplotype =
            RegionHaplotype::new(active_region.clone(), vec![haplotype.clone()]).unwrap();

        assert_eq!(region_haplotype.active_region, active_region);
        assert_eq!(region_haplotype.haplotypes, vec![haplotype.clone()]);
        assert_eq!(
            region_haplotype.min_depth,
            region_haplotype.active_region.stats.min + haplotype.stats.min
        );
        assert!(region_haplotype.is_wildtype);
    }

    #[test]
    fn region_haplotype_wildtype_requires_all_haplotypes_wildtype() {
        let active_region = ar_default();
        let mismatch = AlignNode::new(AlignNode::MISMATCH, 1, None).unwrap();
        let match_node = AlignNode::new(AlignNode::MATCH, 15, Some(Box::new(mismatch))).unwrap();
        let haplotype = make_haplotype(&active_region, match_node);
        let region_haplotype = RegionHaplotype::new(active_region, vec![haplotype]).unwrap();

        assert!(!region_haplotype.is_wildtype);
    }

    #[test]
    fn region_haplotype_rejects_empty_haplotype_list() {
        assert_eq!(
            RegionHaplotype::new(ar_default(), Vec::new()),
            Err(ActiveRegionError::EmptyHaplotypeList)
        );
    }

    #[test]
    fn region_haplotype_orders_by_active_region_and_displays() {
        let count = vec![5; 13];
        let util = KmerUtil::new(4).unwrap();
        let early = ActiveRegion::new(ref16(), 0, 4, &count, &util).unwrap();
        let late = ActiveRegion::new(ref16(), 5, 10, &count, &util).unwrap();
        let early_region = RegionHaplotype::new(
            early.clone(),
            vec![make_haplotype(
                &early,
                AlignNode::new(AlignNode::MATCH, 16, None).unwrap(),
            )],
        )
        .unwrap();
        let late_region = RegionHaplotype::new(
            late.clone(),
            vec![make_haplotype(
                &late,
                AlignNode::new(AlignNode::MATCH, 16, None).unwrap(),
            )],
        )
        .unwrap();

        assert!(early_region < late_region);
        assert!(late_region > early_region);
        let display = early_region.to_string();
        assert!(display.contains("wildtype"));
        assert!(display.contains("mindepth"));
    }

    #[test]
    fn active_region_container_allows_empty_and_none_haplotypes() {
        let empty = ActiveRegionContainer::new(ref16(), Some(&[]), &[1, 2, 3]).unwrap();
        assert_eq!(empty.haplotypes.len(), 0);
        assert_eq!(empty.stats.n, 3);

        let none = ActiveRegionContainer::new(ref16(), None, &[1, 2, 3]).unwrap();
        assert_eq!(none.haplotypes.len(), 0);
    }

    #[test]
    fn active_region_container_stores_reference_and_clones_haplotypes() {
        let ref_region = ref16();
        let region_haplotype = make_region_haplotype();
        let input = vec![region_haplotype.clone()];
        let container =
            ActiveRegionContainer::new(ref_region.clone(), Some(input.as_slice()), &[1, 2, 3])
                .unwrap();

        assert_eq!(container.ref_region, ref_region);
        assert_eq!(container.haplotypes, input);
    }

    #[test]
    fn active_region_container_stats_cover_full_count_array() {
        let container =
            ActiveRegionContainer::new(ref16(), Some(&[]), &[10, 20, 30, 40, 50]).unwrap();

        assert_eq!(container.stats.min, 10);
        assert_eq!(container.stats.max, 50);
        assert_eq!(container.stats.n, 5);
    }

    #[test]
    fn active_region_detector_defaults_match_java() {
        let detector = ActiveRegionDetector::new(KmerUtil::new(4).unwrap()).unwrap();
        assert_eq!(
            detector.minimum_difference(),
            ActiveRegionDetector::DEFAULT_MINIMUM_DIFFERENCE
        );
        assert_eq!(
            detector.difference_quantile(),
            ActiveRegionDetector::DEFAULT_DIFFERENCE_QUANTILE
        );
        assert!(detector.anchor_both_ends());
        assert!(detector.count_reverse_kmers());
        assert!(detector.call_ambiguous_regions());
        assert_eq!(detector.peak_scan_length(), 7);
        assert_eq!(detector.scan_limit_factor(), 7.0);
        assert_eq!(detector.decay_minimum(), 0.55);
        assert_eq!(detector.decay_alpha(), 0.80);
        assert_eq!(detector.max_repeat_count(), 0);
        assert!(detector.recover_right_anchor());
        assert!(!detector.emit_wildtype_active_regions());
        assert!(!detector.trace_haplotype_alignment());
        assert_eq!(detector.scan_limit_length(), 28);
    }

    #[test]
    fn active_region_detector_difference_threshold_matches_java_detector_quantile_quirk() {
        let mut detector = ActiveRegionDetector::new(KmerUtil::new(20).unwrap()).unwrap();
        detector.set_minimum_difference(5).unwrap();
        detector.set_difference_quantile(0.90).unwrap();
        let counts = [
            20320, 21214, 23717, 24751, 24555, 21382, 21499, 26513, 25154, 26661, 26536, 26633,
            21662, 20471, 20483, 21048, 21226, 21403, 21503, 21805, 21694, 21648, 21646, 21419,
            21460, 23762, 24142, 23891, 22801, 22787, 22938, 23009, 22823, 23764, 29079, 28929,
            28820, 29199, 29139, 29036, 28896, 28766, 26133, 6331, 6325, 6347, 6331, 6249, 5879,
            5860, 5912, 5848, 5866, 5869, 5887, 5712, 5802, 5708, 5691, 5689, 5669, 5755, 5715,
            5741, 5684, 4056, 4060, 26513, 25154, 26661, 26536, 26633, 5849, 5737, 5732, 5871,
            5911, 5944, 5985, 6002, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            28896,
        ];

        assert_eq!(detector.difference_threshold(&counts), 2322);
    }

    #[test]
    fn active_region_detector_left_exponential_recovery_matches_java_direction() {
        let mut detector = ActiveRegionDetector::new(KmerUtil::new(20).unwrap()).unwrap();
        detector.set_minimum_difference(5).unwrap();
        detector.set_difference_quantile(0.90).unwrap();
        detector.set_anchor_both_ends(true);
        detector.set_decay_minimum(0.55).unwrap();
        detector.set_decay_alpha(0.80).unwrap();
        detector.set_peak_scan_length(7).unwrap();
        detector.set_scan_limit_factor(7.0).unwrap();
        let sequence = vec![b'A'; 120];
        let ref_region = ref_region(&sequence);
        let counts = [
            0, 0, 0, 0, 0, 7, 21499, 26513, 25154, 26661, 26536, 26633, 21662, 20471, 20483, 21048,
            21226, 21403, 21503, 21805, 21694, 21648, 21646, 21, 22, 26, 24, 25, 24, 35, 34, 34,
            1576, 23, 29, 29, 28, 29, 29, 29, 29, 28, 24, 6331, 6325, 6347, 6331, 6249, 5879, 5860,
            5912, 5848, 5866, 5869, 5887, 5712, 5802, 5708, 5691, 5689, 5669, 5755, 5715, 5741,
            5684, 4056, 4060, 26513, 25154, 26661, 26536, 26633, 5849, 5737, 5732, 5871, 5911,
            5944, 5985, 6002, 6004, 5993, 5969, 5938, 5944, 5853, 5936, 5906, 5714, 5704, 5741,
            5787, 5604, 5698, 29079, 28929, 28820, 29199, 29139, 29036, 28896,
        ];
        let mut candidates = Vec::new();

        let regions = detector
            .detect_from_counts_with(&ref_region, &counts, |region| {
                candidates.push((region.start_index, region.end_index));
                Ok::<bool, std::convert::Infallible>(matches!(
                    (region.start_index, region.end_index),
                    (33, 86) | (71, 113)
                ))
            })
            .unwrap();

        assert!(candidates.contains(&(33, 86)), "{candidates:?}");
        assert!(
            regions
                .iter()
                .any(|region| { (region.start_index, region.end_index) == (33, 86) })
        );
    }

    #[test]
    fn active_region_detector_validates_settings() {
        let mut detector = ActiveRegionDetector::new(KmerUtil::new(4).unwrap()).unwrap();
        assert!(matches!(
            detector.set_minimum_difference(0),
            Err(ActiveRegionDetectorError::InvalidMinimumDifference(0))
        ));
        assert!(matches!(
            detector.set_difference_quantile(1.0),
            Err(ActiveRegionDetectorError::InvalidDifferenceQuantile(1.0))
        ));
        assert!(matches!(
            detector.set_decay_minimum(-0.1),
            Err(ActiveRegionDetectorError::InvalidDecayMinimum(_))
        ));
        assert!(matches!(
            detector.set_decay_alpha(0.0),
            Err(ActiveRegionDetectorError::InvalidDecayAlpha(0.0))
        ));
        assert!(matches!(
            detector.set_peak_scan_length(-1),
            Err(ActiveRegionDetectorError::NegativePeakScanLength(-1))
        ));
        assert!(matches!(
            detector.set_scan_limit_factor(-1.0),
            Err(ActiveRegionDetectorError::NegativeScanLimitFactor(_))
        ));
        assert!(matches!(
            detector.set_max_repeat_count(-1),
            Err(ActiveRegionDetectorError::NegativeMaxRepeatCount(-1))
        ));
    }

    #[test]
    fn active_region_detector_counts_reference_kmers_and_reverse_complements() {
        let util = KmerUtil::new(4).unwrap();
        let mut map = StaticCountMap::default();
        map.insert(&util, "AAAA", 7);
        map.insert(&util, "TTTT", 11);
        map.insert(&util, "AAAC", 3);

        let detector = ActiveRegionDetector::new(util).unwrap();
        let counts = detector.get_counts(&ref16(), &map);
        assert_eq!(counts[0], 18);
        assert_eq!(counts[1], 3);

        let ambiguous = ref_region(b"AAAANCCCGGGGTTTT");
        let ambiguous_counts = detector.get_counts(&ambiguous, &map);
        assert_eq!(ambiguous_counts[1], 0);
        assert_eq!(ambiguous_counts[2], 0);
        assert_eq!(ambiguous_counts[3], 0);
        assert_eq!(ambiguous_counts[4], 0);
    }

    #[test]
    fn active_region_detector_finds_right_anchored_drop_and_recovery() {
        let mut detector = ActiveRegionDetector::new(KmerUtil::new(4).unwrap()).unwrap();
        detector.set_difference_quantile(0.0).unwrap();
        detector.set_minimum_difference(5).unwrap();
        detector.set_decay_minimum(1.0).unwrap();
        detector.set_peak_scan_length(0).unwrap();

        let regions = detector
            .detect_from_counts(&ref16(), &[30, 30, 5, 5, 5, 5, 30, 30, 30, 30, 30, 30, 30])
            .unwrap();

        assert_eq!(regions.len(), 1);
        assert_eq!(regions[0].start_kmer_index, 1);
        assert_eq!(regions[0].end_kmer_index, 6);
        assert!(!regions[0].left_end);
        assert!(!regions[0].right_end);
    }

    #[test]
    fn active_region_detector_splits_repetitive_peaks_at_last_stable_valley() {
        let mut detector = ActiveRegionDetector::new(KmerUtil::new(20).unwrap()).unwrap();
        detector.set_minimum_difference(5).unwrap();
        detector.set_difference_quantile(0.90).unwrap();
        detector.set_decay_minimum(0.55).unwrap();
        detector.set_decay_alpha(0.80).unwrap();
        detector.set_peak_scan_length(7).unwrap();
        detector.set_scan_limit_factor(7.0).unwrap();
        detector.set_anchor_both_ends(true);

        let counts = [
            114, 149, 153, 59202, 59100, 52224, 53222, 66378, 64023, 66797, 66598, 67415, 16599,
            16308, 16240, 16992, 17180, 17448, 17327, 17719, 16701, 16696, 16617, 9, 5, 5, 5, 5, 5,
            6, 6, 6, 9, 6, 925, 964, 962, 976, 944, 926, 916, 903, 903, 50498, 50138, 48915, 47040,
            44388, 44729, 44481, 44683, 44318, 44677, 44692, 45030, 45854, 47400, 47460, 47328,
            47362, 47909, 18, 18, 19, 22, 15, 16, 20, 18, 22, 21, 23, 7, 7, 6, 7, 7, 4, 4, 6, 7,
            16696, 16617, 16494, 16438, 16333, 16546, 16184, 16001, 15990, 15958, 16076, 15435,
            15844, 77353, 77541, 76438, 77276, 75778, 74633, 75227,
        ];

        let regions = detector
            .detect_from_counts(&vntyper_ns_region(), &counts)
            .unwrap();

        let region_bounds = regions
            .iter()
            .map(|region| (region.start_kmer_index, region.end_kmer_index))
            .collect::<Vec<_>>();
        assert_eq!(region_bounds, [(4, 43), (60, 94)]);
    }

    #[test]
    fn active_region_detector_respects_anchor_both_ends_for_right_end_regions() {
        let mut detector = ActiveRegionDetector::new(KmerUtil::new(4).unwrap()).unwrap();
        detector.set_difference_quantile(0.0).unwrap();
        detector.set_decay_minimum(1.0).unwrap();
        detector.set_peak_scan_length(0).unwrap();
        detector.set_recover_right_anchor(false);

        let counts = [30, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5];
        assert_eq!(detector.detect_from_counts(&ref16(), &counts).unwrap(), []);

        detector.set_anchor_both_ends(false);
        let regions = detector.detect_from_counts(&ref16(), &counts).unwrap();
        assert_eq!(regions.len(), 1);
        assert_eq!(regions[0].start_kmer_index, 0);
        assert!(regions[0].right_end);
    }

    #[test]
    fn active_region_detector_finds_left_end_region_when_allowed() {
        let mut detector = ActiveRegionDetector::new(KmerUtil::new(4).unwrap()).unwrap();
        detector.set_difference_quantile(0.0).unwrap();
        detector.set_decay_minimum(1.0).unwrap();
        detector.set_peak_scan_length(0).unwrap();
        detector.set_anchor_both_ends(false);

        let regions = detector
            .detect_from_counts(&ref16(), &[5, 5, 5, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30])
            .unwrap();

        assert_eq!(regions.len(), 1);
        assert!(regions[0].left_end);
        assert_eq!(regions[0].end_kmer_index, 3);
    }

    #[test]
    fn active_region_detector_filters_ambiguous_regions_when_disabled() {
        let mut detector = ActiveRegionDetector::new(KmerUtil::new(4).unwrap()).unwrap();
        detector.set_difference_quantile(0.0).unwrap();
        detector.set_decay_minimum(1.0).unwrap();
        detector.set_peak_scan_length(0).unwrap();
        detector.set_call_ambiguous_regions(false);

        let ambiguous = ref_region(b"AAANCCCGGGGTTTT");
        let counts = [30, 30, 5, 5, 5, 5, 30, 30, 30, 30, 30, 30];
        let regions = detector.detect_from_counts(&ambiguous, &counts).unwrap();
        assert!(regions.is_empty());
    }

    fn ref16() -> ReferenceRegion {
        ref_region(b"AAAACCCCGGGGTTTT")
    }

    fn vntyper_ns_region() -> ReferenceRegion {
        ref_region(
            b"TGCGGGGGCGGTGGAGCCCGGGGCCGGCCTGCTCTCCGGGGCTGAGGTGACACCGTGGGCTGGGGGGGCGGTGGAGCCCGTGGCCGGCCTGCTCTCCGGGGCCGAGGTGACACCGTGGGC",
        )
    }

    fn ar_default() -> ActiveRegion {
        let count = vec![5; 13];
        let util = KmerUtil::new(4).unwrap();
        ActiveRegion::new(ref16(), 2, 10, &count, &util).unwrap()
    }

    fn stats() -> RegionStats {
        RegionStats::from_counts(&[1, 2, 3, 4, 5], 0, 5).unwrap()
    }

    fn make_haplotype(active_region: &ActiveRegion, alignment: AlignNode) -> Haplotype {
        let stats = RegionStats::from_counts(
            &[5; 13],
            active_region.start_kmer_index,
            active_region.end_kmer_index,
        )
        .unwrap();
        Haplotype::new(
            b"AAAACCCCGGGGTTTT".to_vec(),
            active_region.clone(),
            vec![alignment],
            100.0,
            None,
            stats,
        )
        .unwrap()
    }

    fn make_region_haplotype() -> RegionHaplotype {
        let active_region = ar_default();
        RegionHaplotype::new(
            active_region.clone(),
            vec![make_haplotype(
                &active_region,
                AlignNode::new(AlignNode::MATCH, 16, None).unwrap(),
            )],
        )
        .unwrap()
    }

    fn ref_region(sequence: &[u8]) -> ReferenceRegion {
        let reference =
            ReferenceSequence::new("chr1", sequence.len() as i32, Some(digest()), Some("test"))
                .unwrap();
        ReferenceRegion::whole(reference, sequence, 0).unwrap()
    }

    fn digest() -> Digest {
        Digest::new((0..16).collect::<Vec<_>>(), "MD5").unwrap()
    }

    #[derive(Default)]
    struct StaticCountMap {
        counts: HashMap<KmerKey, u32>,
    }

    impl StaticCountMap {
        fn insert(&mut self, util: &KmerUtil, sequence: &str, count: u32) {
            self.counts.insert(util.encode(sequence).unwrap(), count);
        }
    }

    impl CountMap for StaticCountMap {
        fn get(&self, kmer: &KmerKey) -> u32 {
            self.counts.get(kmer).copied().unwrap_or(0)
        }

        fn set(&mut self, _sample: InputSample) -> Result<(), CountMapError> {
            Ok(())
        }

        fn abort(&self) {}

        fn is_aborted(&self) -> bool {
            false
        }
    }
}
