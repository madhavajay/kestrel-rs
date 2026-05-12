use thiserror::Error;

use crate::constants::{ARRAY_EXPAND_FACTOR, MAX_ARRAY_SIZE};
use crate::variant::{VariantCall, VariantType};

/// Errors returned while parsing or applying variant filters.
#[derive(Clone, Debug, Error, Eq, PartialEq)]
pub enum VariantFilterError {
    /// Filter specification was empty.
    #[error("cannot get variant filter with an empty specification")]
    EmptySpec,
    /// Filter name was not recognized.
    #[error("unknown variant filter: {0}")]
    UnknownFilter(String),
    /// Variant type name was not recognized.
    #[error("unrecognized type: {0}")]
    UnknownType(String),
    /// Filter arguments were empty.
    #[error("empty arguments for {0} variant filter")]
    EmptyArgs(&'static str),
    /// Filter argument was invalid.
    #[error("invalid argument for {filter} variant filter: {argument}")]
    InvalidArgument {
        /// Filter name.
        filter: &'static str,
        /// Invalid argument text.
        argument: String,
    },
    /// Filter runner capacity was negative.
    #[error("cannot set a negative filter capacity: {0}")]
    NegativeCapacity(i32),
    /// Filter runner reached maximum capacity.
    #[error("maximum number of filters reached: {0}")]
    MaximumCapacity(usize),
}

/// A configured variant filter.
#[derive(Clone, Debug, PartialEq)]
pub enum VariantFilterKind {
    /// Type-based filter.
    Type(TypeVariantFilter),
    /// Location-based filter.
    Location(LocationVariantFilter),
    /// Coverage-based filter.
    Coverage(CoverageVariantFilter),
}

impl VariantFilterKind {
    /// Parses a filter specification into a variant filter.
    pub fn get_filter(spec: &str) -> Result<Self, VariantFilterError> {
        let spec = spec.trim();
        if spec.is_empty() {
            return Err(VariantFilterError::EmptySpec);
        }

        let mut parts = spec.splitn(2, ':');
        let name = parts.next().unwrap_or_default().trim();
        let args = parts.next().unwrap_or_default().trim();

        match name.to_ascii_lowercase().as_str() {
            "type" => Ok(Self::Type(TypeVariantFilter::new(args)?)),
            "location" => Ok(Self::Location(LocationVariantFilter::new(args)?)),
            "coverage" => Ok(Self::Coverage(CoverageVariantFilter::new(args)?)),
            _ => Err(VariantFilterError::UnknownFilter(name.to_owned())),
        }
    }

    /// Returns the description for a named filter.
    #[must_use]
    pub fn description(filter_name: &str) -> Option<&'static str> {
        match filter_name.to_ascii_lowercase().as_str() {
            "type" => Some(TypeVariantFilter::DESCRIPTION),
            "location" => Some(LocationVariantFilter::DESCRIPTION),
            "coverage" => Some(CoverageVariantFilter::DESCRIPTION),
            _ => None,
        }
    }

    /// Applies this filter to a variant.
    #[must_use]
    pub fn filter<'a, V: VariantCall>(&self, variant: Option<&'a V>) -> Option<&'a V> {
        match self {
            Self::Type(filter) => filter.filter(variant),
            Self::Location(filter) => filter.filter(variant),
            Self::Coverage(filter) => filter.filter(variant),
        }
    }
}

/// Variant filter that keeps selected variant types.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct TypeVariantFilter {
    type_snp: bool,
    type_ins: bool,
    type_del: bool,
}

impl Default for TypeVariantFilter {
    fn default() -> Self {
        Self {
            type_snp: true,
            type_ins: true,
            type_del: true,
        }
    }
}

impl TypeVariantFilter {
    /// Filter description.
    pub const DESCRIPTION: &'static str = "Filter variants by type.";

    /// Creates a type filter from comma-separated type names.
    pub fn new(args: &str) -> Result<Self, VariantFilterError> {
        let args = args.trim();
        if args.is_empty() {
            return Err(VariantFilterError::EmptyArgs("type"));
        }

        let mut filter = Self {
            type_snp: false,
            type_ins: false,
            type_del: false,
        };

        for token in args.split(',').map(str::trim) {
            match token {
                "snp" | "snv" => filter.type_snp = true,
                "ins" | "insertion" => filter.type_ins = true,
                "del" | "deletion" => filter.type_del = true,
                "indel" | "insdel" | "insertiondeletion" => {
                    filter.type_ins = true;
                    filter.type_del = true;
                }
                _ => return Err(VariantFilterError::UnknownType(token.to_owned())),
            }
        }

        Ok(filter)
    }

    /// Applies this type filter to a variant.
    #[must_use]
    pub fn filter<'a, V: VariantCall>(&self, variant: Option<&'a V>) -> Option<&'a V> {
        let variant = variant?;
        let keep = match variant.data().variant_type {
            VariantType::Snp => self.type_snp,
            VariantType::Insertion => self.type_ins,
            VariantType::Deletion => self.type_del,
        };
        keep.then_some(variant)
    }
}

/// Variant filter that removes calls near reference-region ends.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct LocationVariantFilter {
    left_end_length: i32,
    right_end_length: i32,
}

impl LocationVariantFilter {
    /// Filter description.
    pub const DESCRIPTION: &'static str = "Filter variants by distance from reference ends.";

    /// Creates a location filter from length arguments.
    pub fn new(args: &str) -> Result<Self, VariantFilterError> {
        let args = args.trim();
        if args.is_empty() {
            return Err(VariantFilterError::EmptyArgs("location"));
        }

        let mut filter = Self::default();
        for token in args.split(',').map(str::trim) {
            let (attribute, value) =
                parse_assignment(token).ok_or_else(|| VariantFilterError::InvalidArgument {
                    filter: "location",
                    argument: token.to_owned(),
                })?;
            let value = parse_nonnegative_i32("location", token, value)?;

            match attribute {
                "length" | "len" | "l" => {
                    filter.left_end_length = value;
                    filter.right_end_length = value;
                }
                "leftlength" | "leftlen" | "leftl" | "llength" | "llen" | "ll" => {
                    filter.left_end_length = value;
                }
                "rightlength" | "rightlen" | "rightl" | "rlength" | "rlen" | "rl" => {
                    filter.right_end_length = value;
                }
                _ => {
                    return Err(VariantFilterError::InvalidArgument {
                        filter: "location",
                        argument: token.to_owned(),
                    });
                }
            }
        }

        if filter.left_end_length == 0 && filter.right_end_length == 0 {
            return Err(VariantFilterError::InvalidArgument {
                filter: "location",
                argument: args.to_owned(),
            });
        }

        Ok(filter)
    }

    /// Applies this location filter to a variant.
    #[must_use]
    pub fn filter<'a, V: VariantCall>(&self, variant: Option<&'a V>) -> Option<&'a V> {
        let variant = variant?;
        if variant.data().start <= self.left_end_length {
            return None;
        }
        if variant.reference_end()
            >= variant.data().active_region.ref_region.size - self.right_end_length
        {
            return None;
        }
        Some(variant)
    }
}

/// Variant filter that removes low coverage or low depth calls.
#[derive(Clone, Debug, PartialEq)]
pub struct CoverageVariantFilter {
    min_coverage: f64,
    min_depth: i32,
}

impl Default for CoverageVariantFilter {
    fn default() -> Self {
        Self {
            min_coverage: 0.0,
            min_depth: i32::MAX,
        }
    }
}

impl CoverageVariantFilter {
    /// Filter description.
    pub const DESCRIPTION: &'static str = "Filter variants by coverage and depth.";

    /// Creates a coverage filter from coverage and depth arguments.
    pub fn new(args: &str) -> Result<Self, VariantFilterError> {
        let args = args.trim();
        if args.is_empty() {
            return Err(VariantFilterError::EmptyArgs("coverage"));
        }

        let mut filter = Self {
            min_coverage: 0.0,
            min_depth: 0,
        };
        let mut arg0_is_bare = false;

        for (index, token) in args.split(',').map(str::trim).enumerate() {
            if let Some((attribute, _value)) = parse_assignment(token) {
                match attribute {
                    "coverage" | "cov" | "c" => {
                        // Preserve Java bug: parse the attribute name, not the value.
                        filter.min_coverage = parse_coverage("coverage", token, attribute)?;
                    }
                    "depth" | "dep" | "d" => {
                        // Preserve Java bug: parse the attribute name, not the value.
                        filter.min_depth = parse_nonnegative_i32("coverage", token, attribute)?;
                    }
                    _ => {
                        return Err(VariantFilterError::InvalidArgument {
                            filter: "coverage",
                            argument: token.to_owned(),
                        });
                    }
                }
            } else if index == 0 {
                filter.min_coverage = parse_coverage("coverage", token, token)?;
                arg0_is_bare = true;
            } else if index == 1 && arg0_is_bare {
                filter.min_depth = parse_nonnegative_i32("coverage", token, token)?;
            } else {
                return Err(VariantFilterError::InvalidArgument {
                    filter: "coverage",
                    argument: token.to_owned(),
                });
            }
        }

        Ok(filter)
    }

    /// Applies this coverage filter to a variant.
    #[must_use]
    pub fn filter<'a, V: VariantCall>(&self, variant: Option<&'a V>) -> Option<&'a V> {
        let variant = variant?;
        if variant.data().variant_depth < self.min_depth {
            return None;
        }
        if f64::from(variant.data().variant_depth) / f64::from(variant.data().locus_depth)
            < self.min_coverage
        {
            return None;
        }
        Some(variant)
    }
}

/// Ordered collection of variant filters.
#[derive(Clone, Debug)]
pub struct VariantFilterRunner {
    filters: Vec<VariantFilterKind>,
    filter_capacity: usize,
}

impl Default for VariantFilterRunner {
    fn default() -> Self {
        Self::new(Self::DEFAULT_FILTER_CAPACITY).expect("default filter capacity is valid")
    }
}

impl VariantFilterRunner {
    /// Default initial filter capacity.
    pub const DEFAULT_FILTER_CAPACITY: i32 = 5;

    /// Creates an empty filter runner.
    pub fn new(filter_capacity: i32) -> Result<Self, VariantFilterError> {
        if filter_capacity < 0 {
            return Err(VariantFilterError::NegativeCapacity(filter_capacity));
        }
        Ok(Self {
            filters: Vec::with_capacity(filter_capacity as usize),
            filter_capacity: filter_capacity as usize,
        })
    }

    /// Adds one filter, expanding capacity as needed.
    pub fn add_filter(
        &mut self,
        variant_filter: Option<VariantFilterKind>,
    ) -> Result<(), VariantFilterError> {
        let Some(variant_filter) = variant_filter else {
            return Ok(());
        };

        if self.filters.len() == self.filter_capacity {
            let new_capacity = ((self.filter_capacity as f32) * ARRAY_EXPAND_FACTOR) as usize;
            if new_capacity == 0 {
                panic!("index out of bounds for zero-capacity VariantFilterRunner");
            }
            if new_capacity > MAX_ARRAY_SIZE {
                if self.filter_capacity == MAX_ARRAY_SIZE {
                    return Err(VariantFilterError::MaximumCapacity(MAX_ARRAY_SIZE));
                }
                self.filter_capacity = MAX_ARRAY_SIZE;
            } else {
                self.filter_capacity = new_capacity;
            }
            self.filters
                .reserve(self.filter_capacity.saturating_sub(self.filters.capacity()));
        }

        self.filters.push(variant_filter);
        Ok(())
    }

    /// Adds several filters.
    pub fn add_filters(
        &mut self,
        variant_filters: impl IntoIterator<Item = VariantFilterKind>,
    ) -> Result<(), VariantFilterError> {
        for variant_filter in variant_filters {
            self.add_filter(Some(variant_filter))?;
        }
        Ok(())
    }

    /// Applies all filters in order to a variant.
    #[must_use]
    pub fn filter<'a, V: VariantCall>(&self, variant: Option<&'a V>) -> Option<&'a V> {
        let mut variant = variant;
        for filter in &self.filters {
            variant = filter.filter(variant);
            if variant.is_none() {
                break;
            }
        }
        variant
    }
}

fn parse_assignment(token: &str) -> Option<(&str, &str)> {
    token
        .split_once('=')
        .map(|(attribute, value)| (attribute.trim(), value.trim()))
}

fn parse_nonnegative_i32(
    filter: &'static str,
    token: &str,
    value: &str,
) -> Result<i32, VariantFilterError> {
    let parsed = value
        .parse::<i32>()
        .map_err(|_| VariantFilterError::InvalidArgument {
            filter,
            argument: token.to_owned(),
        })?;
    if parsed < 0 {
        return Err(VariantFilterError::InvalidArgument {
            filter,
            argument: token.to_owned(),
        });
    }
    Ok(parsed)
}

fn parse_coverage(
    filter: &'static str,
    token: &str,
    value: &str,
) -> Result<f64, VariantFilterError> {
    let parsed = value
        .parse::<f64>()
        .map_err(|_| VariantFilterError::InvalidArgument {
            filter,
            argument: token.to_owned(),
        })?;
    if !(0.0..=1.0).contains(&parsed) {
        return Err(VariantFilterError::InvalidArgument {
            filter,
            argument: token.to_owned(),
        });
    }
    Ok(parsed)
}

#[cfg(test)]
mod tests {
    use kanalyze::util::KmerUtil;

    use crate::activeregion::{ActiveRegion, Haplotype, RegionStats};
    use crate::align::AlignNode;
    use crate::refreader::{ReferenceRegion, ReferenceSequence};
    use crate::util::digest::Digest;
    use crate::variant::{VariantInsertion, VariantSnp};

    use super::*;

    #[test]
    fn factory_returns_known_filters_and_descriptions() {
        assert!(matches!(
            VariantFilterKind::get_filter("TYPE:snp").unwrap(),
            VariantFilterKind::Type(_)
        ));
        assert!(matches!(
            VariantFilterKind::get_filter("LOCATION:length=3").unwrap(),
            VariantFilterKind::Location(_)
        ));
        assert!(matches!(
            VariantFilterKind::get_filter("COVERAGE:0.5").unwrap(),
            VariantFilterKind::Coverage(_)
        ));
        assert!(
            !VariantFilterKind::description("coverage")
                .unwrap()
                .is_empty()
        );
        assert_eq!(VariantFilterKind::description("nope"), None);
        assert_eq!(
            VariantFilterKind::get_filter("NOPE:x").unwrap_err(),
            VariantFilterError::UnknownFilter("NOPE".to_owned())
        );
    }

    #[test]
    fn type_filter_keeps_and_drops_by_type() {
        let snp_filter = TypeVariantFilter::new("snp").unwrap();
        let ins_filter = TypeVariantFilter::new("ins").unwrap();
        let indel_filter = TypeVariantFilter::new("indel").unwrap();
        let snp = snp_at(5, 5, 10);
        let insertion = insertion();

        assert!(snp_filter.filter(Some(&snp)).is_some());
        assert!(snp_filter.filter(Some(&insertion)).is_none());
        assert!(ins_filter.filter(Some(&insertion)).is_some());
        assert!(ins_filter.filter(Some(&snp)).is_none());
        assert!(indel_filter.filter(Some(&insertion)).is_some());
        assert!(indel_filter.filter(Some(&snp)).is_none());
        assert_eq!(
            TypeVariantFilter::new("not_a_type").unwrap_err(),
            VariantFilterError::UnknownType("not_a_type".to_owned())
        );
    }

    #[test]
    fn type_filter_default_allows_all_types() {
        let filter = TypeVariantFilter::default();
        assert!(filter.filter(Some(&snp_at(5, 5, 10))).is_some());
        assert!(filter.filter(Some(&insertion())).is_some());
        assert!(filter.filter::<VariantSnp>(None).is_none());
    }

    #[test]
    fn location_filter_matches_java_edges() {
        let left = LocationVariantFilter::new("ll=5").unwrap();
        assert!(left.filter(Some(&snp_at(1, 5, 10))).is_none());
        assert!(left.filter(Some(&snp_at(5, 5, 10))).is_none());
        assert!(left.filter(Some(&snp_at(6, 5, 10))).is_some());

        let right = LocationVariantFilter::new("rl=5").unwrap();
        assert!(right.filter(Some(&snp_at(10, 5, 10))).is_some());
        assert!(right.filter(Some(&snp_at(12, 5, 10))).is_none());

        let both = LocationVariantFilter::new("length=3").unwrap();
        assert!(both.filter(Some(&snp_at(2, 5, 10))).is_none());
        assert!(both.filter(Some(&snp_at(8, 5, 10))).is_some());
        assert!(both.filter(Some(&snp_at(14, 5, 10))).is_none());
    }

    #[test]
    fn location_filter_rejects_bad_args() {
        assert!(LocationVariantFilter::new("").is_err());
        assert!(LocationVariantFilter::new("ll=0,rl=0").is_err());
        assert!(LocationVariantFilter::new("nonsense=5").is_err());
        assert!(LocationVariantFilter::new("ll=-3").is_err());
        assert!(LocationVariantFilter::new("ll=abc").is_err());
    }

    #[test]
    fn coverage_filter_matches_bare_arg_behavior() {
        let coverage = CoverageVariantFilter::new("0.5").unwrap();
        assert!(coverage.filter(Some(&snp_at(5, 60, 100))).is_some());
        assert!(coverage.filter(Some(&snp_at(5, 40, 100))).is_none());

        let coverage_and_depth = CoverageVariantFilter::new("0.5,5").unwrap();
        assert!(
            coverage_and_depth
                .filter(Some(&snp_at(5, 60, 100)))
                .is_some()
        );
        assert!(
            coverage_and_depth
                .filter(Some(&snp_at(5, 60, 100_000)))
                .is_none()
        );
        assert!(coverage_and_depth.filter(Some(&snp_at(5, 3, 4))).is_none());
    }

    #[test]
    fn coverage_filter_preserves_attribute_value_parse_bug() {
        assert!(CoverageVariantFilter::new("cov=0.5").is_err());
        assert!(CoverageVariantFilter::new("depth=10").is_err());
        assert!(CoverageVariantFilter::new("0.5,depth=10").is_err());
    }

    #[test]
    fn coverage_filter_default_drops_real_variants() {
        assert!(
            CoverageVariantFilter::default()
                .filter(Some(&snp_at(5, 100, 100)))
                .is_none()
        );
        assert!(
            CoverageVariantFilter::default()
                .filter::<VariantSnp>(None)
                .is_none()
        );
    }

    #[test]
    fn runner_chains_filters_and_ignores_none() {
        let snp = snp_at(5, 5, 10);
        let mut runner = VariantFilterRunner::default();
        assert!(runner.filter(Some(&snp)).is_some());
        assert!(runner.filter::<VariantSnp>(None).is_none());

        runner.add_filter(None).unwrap();
        runner
            .add_filter(Some(VariantFilterKind::get_filter("TYPE:snp").unwrap()))
            .unwrap();
        runner
            .add_filter(Some(
                VariantFilterKind::get_filter("LOCATION:length=2").unwrap(),
            ))
            .unwrap();
        assert!(runner.filter(Some(&snp)).is_some());

        runner
            .add_filter(Some(VariantFilterKind::get_filter("TYPE:ins").unwrap()))
            .unwrap();
        assert!(runner.filter(Some(&snp)).is_none());
    }

    #[test]
    fn runner_validates_capacity_and_preserves_zero_capacity_bug() {
        assert_eq!(
            VariantFilterRunner::new(-1).unwrap_err(),
            VariantFilterError::NegativeCapacity(-1)
        );

        let mut runner = VariantFilterRunner::new(0).unwrap();
        let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            runner
                .add_filter(Some(VariantFilterKind::get_filter("TYPE:snp").unwrap()))
                .unwrap();
        }));
        assert!(result.is_err());
    }

    fn snp_at(start: i32, variant_depth: i32, locus_depth: i32) -> VariantSnp {
        VariantSnp::new(
            start,
            variant_depth,
            locus_depth,
            "C",
            "A",
            &[haplotype()],
            1,
            true,
        )
        .unwrap()
    }

    fn insertion() -> VariantInsertion {
        VariantInsertion::new(5, 5, 10, "ACG", &[haplotype()], 1, true).unwrap()
    }

    fn haplotype() -> Haplotype {
        let reference = ReferenceSequence::new("chr1", 16, Some(digest()), Some("test")).unwrap();
        let ref_region = ReferenceRegion::whole(reference, b"AAAACCCCGGGGTTTT", 0).unwrap();
        let count = vec![10; 13];
        let kmer_util = KmerUtil::new(4).unwrap();
        let active_region = ActiveRegion::new(ref_region, 0, 12, &count, &kmer_util).unwrap();
        let stats = RegionStats::from_counts(
            &count,
            active_region.start_kmer_index,
            active_region.end_kmer_index,
        )
        .unwrap();
        let alignment = AlignNode::new(AlignNode::MATCH, 16, None).unwrap();
        Haplotype::new(
            b"AAAACCCCGGGGTTTT".to_vec(),
            active_region,
            vec![alignment],
            100.0,
            None,
            stats,
        )
        .unwrap()
    }

    fn digest() -> Digest {
        Digest::new(vec![0; 16], "MD5").unwrap()
    }
}
