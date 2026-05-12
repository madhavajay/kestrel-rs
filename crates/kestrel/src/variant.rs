use std::cmp::Ordering;
use std::fmt;

use thiserror::Error;

use crate::activeregion::{ActiveRegion, Haplotype};
use crate::align::AlignNode;
use crate::refreader::ReferenceSequenceError;

/// Type of small variant call.
#[derive(Clone, Copy, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub enum VariantType {
    /// Single-nucleotide polymorphism.
    Snp,
    /// Insertion relative to the reference.
    Insertion,
    /// Deletion relative to the reference.
    Deletion,
}

impl VariantType {
    /// All variant types in Java Kestrel order.
    pub const ALL: [Self; 3] = [Self::Snp, Self::Insertion, Self::Deletion];

    /// Returns the Java-compatible hash code for this variant type.
    #[must_use]
    pub const fn java_hash_code(self) -> i32 {
        match self {
            Self::Snp => 0o11111111111,
            Self::Insertion => 0o22222222222_u32 as i32,
            Self::Deletion => 0o33333333333_u32 as i32,
        }
    }
}

/// Errors returned while building or calling variants.
#[derive(Clone, Debug, Error, Eq, PartialEq)]
pub enum VariantError {
    /// Haplotype list was empty.
    #[error("haplotype array is empty")]
    EmptyHaplotypeList,
    /// Requested haplotype prefix size was invalid.
    #[error(
        "haplotype size (haplo_size) is less than 1 or greater than the haplotype length: {haplo_size} (haplotype length = {len})"
    )]
    InvalidHaplotypeSize {
        /// Requested haplotype count.
        haplo_size: usize,
        /// Available haplotype count.
        len: usize,
    },
    /// Variant start coordinate was invalid.
    #[error("start position is less than 1: {0}")]
    InvalidStart(i32),
    /// Variant depth was negative.
    #[error("variant depth is negative: {0}")]
    NegativeVariantDepth(i32),
    /// Locus depth was smaller than variant depth.
    #[error("locus depth ({locus_depth}) is less than variant depth ({variant_depth})")]
    LocusDepthLessThanVariantDepth {
        /// Locus depth.
        locus_depth: i32,
        /// Variant-supporting depth.
        variant_depth: i32,
    },
    /// Haplotype belonged to a different active region.
    #[error("haplotype array contains haplotypes from different active regions at index {0}")]
    DifferentActiveRegion(usize),
    /// SNP reference allele was not one base.
    #[error("reference allele must contain exactly one base for a SNP variant: {0}")]
    InvalidSnpRef(String),
    /// SNP alternate allele was not one base.
    #[error("alternate allele must contain exactly one base for a SNP variant: {0}")]
    InvalidSnpAlt(String),
    /// Insertion alternate allele was empty.
    #[error("alternate allele (alt) is empty on an insertion variant")]
    EmptyInsertionAlt,
    /// Deletion reference allele was empty.
    #[error("reference allele (ref) is empty on a deletion variant")]
    EmptyDeletionRef,
    /// Variant caller was used before initialization.
    #[error("variant caller has not been initialized")]
    CallerNotInitialized,
    /// Haplotype active region did not match the caller's active region.
    #[error("haplotype active region does not match the caller active region")]
    CallerRegionMismatch,
    /// Alignment node type was not recognized.
    #[error("unrecognized alignment node type: {0}")]
    UnknownAlignmentNode(i8),
}

/// Common interface implemented by variant call records.
pub trait VariantCall {
    /// Returns shared call data.
    fn data(&self) -> &VariantCallData;
    /// Returns the one-based inclusive reference end coordinate.
    fn reference_end(&self) -> i32;
    /// Returns the VCF position.
    fn vcf_pos(&self) -> i32;
    /// Returns the VCF reference allele.
    fn vcf_ref(&self) -> Result<String, ReferenceSequenceError>;
    /// Returns the VCF alternate allele.
    fn vcf_alt(&self) -> Result<String, ReferenceSequenceError>;
    /// Returns an HGVS-like variant string.
    fn hgvs(&self) -> String;

    /// Returns a fully qualified display string including the reference name.
    fn to_string_full(&self) -> String {
        format!(
            "{}:{}",
            self.data().active_region.ref_region.reference_sequence.name,
            self.hgvs()
        )
    }
}

/// Shared data stored by all variant call types.
#[derive(Clone, Debug, PartialEq)]
pub struct VariantCallData {
    /// Active region that produced this variant.
    pub active_region: ActiveRegion,
    haplotypes: Vec<Haplotype>,
    /// Variant type.
    pub variant_type: VariantType,
    /// One-based variant start coordinate.
    pub start: i32,
    /// Depth supporting this variant.
    pub variant_depth: i32,
    /// Total locus depth.
    pub locus_depth: i32,
    /// Reference allele.
    pub ref_allele: String,
    /// Alternate allele.
    pub alt_allele: String,
    /// True when the reference allele contains ambiguous bases.
    pub is_ambiguous: bool,
    /// True when the call was aligned against the reference.
    pub is_reference_aligned: bool,
}

impl VariantCallData {
    #[allow(clippy::too_many_arguments)]
    fn new(
        variant_type: VariantType,
        start: i32,
        variant_depth: i32,
        locus_depth: i32,
        ref_allele: impl Into<String>,
        alt_allele: impl Into<String>,
        haplotypes: &[Haplotype],
        haplo_size: usize,
        is_reference_aligned: bool,
    ) -> Result<Self, VariantError> {
        if haplotypes.is_empty() {
            return Err(VariantError::EmptyHaplotypeList);
        }
        if haplo_size < 1 || haplo_size > haplotypes.len() {
            return Err(VariantError::InvalidHaplotypeSize {
                haplo_size,
                len: haplotypes.len(),
            });
        }
        if start < 1 {
            return Err(VariantError::InvalidStart(start));
        }
        if variant_depth < 0 {
            return Err(VariantError::NegativeVariantDepth(variant_depth));
        }
        if locus_depth < variant_depth {
            return Err(VariantError::LocusDepthLessThanVariantDepth {
                locus_depth,
                variant_depth,
            });
        }

        let haplotypes = haplotypes[..haplo_size].to_vec();
        let active_region = haplotypes[0].active_region.clone();
        for (index, haplotype) in haplotypes.iter().enumerate() {
            if haplotype.active_region != active_region {
                return Err(VariantError::DifferentActiveRegion(index));
            }
        }

        let ref_allele = ref_allele.into().trim().to_owned();
        let alt_allele = alt_allele.into().trim().to_owned();
        let is_ambiguous = ref_allele.bytes().any(|base| {
            !matches!(
                base,
                b'A' | b'C' | b'G' | b'T' | b'U' | b'a' | b'c' | b'g' | b't' | b'u'
            )
        });

        Ok(Self {
            active_region,
            haplotypes,
            variant_type,
            start,
            variant_depth,
            locus_depth,
            ref_allele,
            alt_allele,
            is_ambiguous,
            is_reference_aligned,
        })
    }

    /// Returns cloned haplotypes supporting this variant.
    #[must_use]
    pub fn haplotypes(&self) -> Vec<Haplotype> {
        self.haplotypes.clone()
    }
}

impl Eq for VariantCallData {}

impl Ord for VariantCallData {
    fn cmp(&self, other: &Self) -> Ordering {
        self.active_region
            .ref_region
            .reference_sequence
            .name
            .cmp(&other.active_region.ref_region.reference_sequence.name)
            .then_with(|| self.start.cmp(&other.start))
            .then_with(|| self.variant_type.cmp(&other.variant_type))
            .then_with(|| self.ref_allele.cmp(&other.ref_allele))
            .then_with(|| self.alt_allele.cmp(&other.alt_allele))
    }
}

impl PartialOrd for VariantCallData {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

/// Single-nucleotide variant call.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct VariantSnp {
    /// Shared variant call data.
    pub data: VariantCallData,
}

impl VariantSnp {
    /// Creates a SNP variant.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        start: i32,
        variant_depth: i32,
        locus_depth: i32,
        ref_allele: impl Into<String>,
        alt_allele: impl Into<String>,
        haplotypes: &[Haplotype],
        haplo_size: usize,
        is_reference_aligned: bool,
    ) -> Result<Self, VariantError> {
        let data = VariantCallData::new(
            VariantType::Snp,
            start,
            variant_depth,
            locus_depth,
            ref_allele,
            alt_allele,
            haplotypes,
            haplo_size,
            is_reference_aligned,
        )?;
        if data.ref_allele.len() != 1 {
            return Err(VariantError::InvalidSnpRef(data.ref_allele));
        }
        if data.alt_allele.len() != 1 {
            return Err(VariantError::InvalidSnpAlt(data.alt_allele));
        }
        Ok(Self { data })
    }
}

impl VariantCall for VariantSnp {
    fn data(&self) -> &VariantCallData {
        &self.data
    }

    fn reference_end(&self) -> i32 {
        self.data.start
    }

    fn vcf_pos(&self) -> i32 {
        self.data.start
    }

    fn vcf_ref(&self) -> Result<String, ReferenceSequenceError> {
        Ok(self.data.ref_allele.clone())
    }

    fn vcf_alt(&self) -> Result<String, ReferenceSequenceError> {
        Ok(self.data.alt_allele.clone())
    }

    fn hgvs(&self) -> String {
        format!(
            "{}{}>{}",
            self.data.start, self.data.ref_allele, self.data.alt_allele
        )
    }
}

impl fmt::Display for VariantSnp {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.hgvs())
    }
}

impl Ord for VariantSnp {
    fn cmp(&self, other: &Self) -> Ordering {
        self.data.cmp(&other.data)
    }
}

impl PartialOrd for VariantSnp {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

/// Insertion variant call.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct VariantInsertion {
    /// Shared variant call data.
    pub data: VariantCallData,
    /// Inserted sequence length.
    pub length: usize,
}

impl VariantInsertion {
    /// Creates an insertion variant.
    pub fn new(
        start: i32,
        variant_depth: i32,
        locus_depth: i32,
        alt_allele: impl Into<String>,
        haplotypes: &[Haplotype],
        haplo_size: usize,
        is_reference_aligned: bool,
    ) -> Result<Self, VariantError> {
        let data = VariantCallData::new(
            VariantType::Insertion,
            start,
            variant_depth,
            locus_depth,
            "",
            alt_allele,
            haplotypes,
            haplo_size,
            is_reference_aligned,
        )?;
        let length = data.alt_allele.len();
        if length == 0 {
            return Err(VariantError::EmptyInsertionAlt);
        }
        Ok(Self { data, length })
    }
}

impl VariantCall for VariantInsertion {
    fn data(&self) -> &VariantCallData {
        &self.data
    }

    fn reference_end(&self) -> i32 {
        self.data.start
    }

    fn vcf_pos(&self) -> i32 {
        if self.data.start == 1 {
            self.data.start
        } else {
            self.data.start - 1
        }
    }

    fn vcf_ref(&self) -> Result<String, ReferenceSequenceError> {
        Ok(char::from(insertion_anchor(&self.data)?).to_string())
    }

    fn vcf_alt(&self) -> Result<String, ReferenceSequenceError> {
        if self.data.start == 1 {
            return Ok(format!(
                "{}{}",
                self.data.alt_allele,
                char::from(
                    self.data
                        .active_region
                        .ref_region
                        .get_base(self.data.start)?
                )
            ));
        }

        Ok(format!(
            "{}{}",
            char::from(insertion_anchor(&self.data)?),
            self.data.alt_allele
        ))
    }

    fn hgvs(&self) -> String {
        format!(
            "{}_{}ins{}",
            self.data.start,
            self.data.start + 1,
            self.data.alt_allele
        )
    }
}

impl fmt::Display for VariantInsertion {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.hgvs())
    }
}

/// Deletion variant call.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct VariantDeletion {
    /// Shared variant call data.
    pub data: VariantCallData,
    /// Deleted sequence length.
    pub length: usize,
}

impl VariantDeletion {
    /// Creates a deletion variant.
    pub fn new(
        start: i32,
        variant_depth: i32,
        locus_depth: i32,
        ref_allele: impl Into<String>,
        haplotypes: &[Haplotype],
        haplo_size: usize,
        is_reference_aligned: bool,
    ) -> Result<Self, VariantError> {
        let data = VariantCallData::new(
            VariantType::Insertion,
            start,
            variant_depth,
            locus_depth,
            ref_allele,
            "",
            haplotypes,
            haplo_size,
            is_reference_aligned,
        )?;
        let length = data.ref_allele.len();
        if length == 0 {
            return Err(VariantError::EmptyDeletionRef);
        }
        Ok(Self { data, length })
    }
}

impl VariantCall for VariantDeletion {
    fn data(&self) -> &VariantCallData {
        &self.data
    }

    fn reference_end(&self) -> i32 {
        self.data.start + self.data.ref_allele.len() as i32 - 1
    }

    fn vcf_pos(&self) -> i32 {
        if self.data.start > 1 {
            self.data.start - 1
        } else {
            self.data.start
        }
    }

    fn vcf_ref(&self) -> Result<String, ReferenceSequenceError> {
        if self.data.start == 1 {
            return Ok(format!(
                "{}{}",
                self.data.ref_allele,
                char::from(
                    self.data
                        .active_region
                        .ref_region
                        .get_base(self.reference_end() + 1)?
                )
            ));
        }

        Ok(format!(
            "{}{}",
            char::from(deletion_anchor(&self.data, self.reference_end())?),
            self.data.ref_allele
        ))
    }

    fn vcf_alt(&self) -> Result<String, ReferenceSequenceError> {
        if self.data.start == 1 {
            return Ok(char::from(
                self.data
                    .active_region
                    .ref_region
                    .get_base(self.reference_end() + 1)?,
            )
            .to_string());
        }

        Ok(char::from(deletion_anchor(&self.data, self.reference_end())?).to_string())
    }

    fn hgvs(&self) -> String {
        format!(
            "{}_{}del{}",
            self.data.start,
            self.data.start + self.length as i32 - 1,
            self.data.alt_allele
        )
    }
}

impl fmt::Display for VariantDeletion {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.hgvs())
    }
}

/// Any supported variant call.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum Variant {
    /// SNP call.
    Snp(VariantSnp),
    /// Insertion call.
    Insertion(VariantInsertion),
    /// Deletion call.
    Deletion(VariantDeletion),
}

impl VariantCall for Variant {
    fn data(&self) -> &VariantCallData {
        match self {
            Self::Snp(variant) => variant.data(),
            Self::Insertion(variant) => variant.data(),
            Self::Deletion(variant) => variant.data(),
        }
    }

    fn reference_end(&self) -> i32 {
        match self {
            Self::Snp(variant) => variant.reference_end(),
            Self::Insertion(variant) => variant.reference_end(),
            Self::Deletion(variant) => variant.reference_end(),
        }
    }

    fn vcf_pos(&self) -> i32 {
        match self {
            Self::Snp(variant) => variant.vcf_pos(),
            Self::Insertion(variant) => variant.vcf_pos(),
            Self::Deletion(variant) => variant.vcf_pos(),
        }
    }

    fn vcf_ref(&self) -> Result<String, ReferenceSequenceError> {
        match self {
            Self::Snp(variant) => variant.vcf_ref(),
            Self::Insertion(variant) => variant.vcf_ref(),
            Self::Deletion(variant) => variant.vcf_ref(),
        }
    }

    fn vcf_alt(&self) -> Result<String, ReferenceSequenceError> {
        match self {
            Self::Snp(variant) => variant.vcf_alt(),
            Self::Insertion(variant) => variant.vcf_alt(),
            Self::Deletion(variant) => variant.vcf_alt(),
        }
    }

    fn hgvs(&self) -> String {
        match self {
            Self::Snp(variant) => variant.hgvs(),
            Self::Insertion(variant) => variant.hgvs(),
            Self::Deletion(variant) => variant.hgvs(),
        }
    }
}

/// Calls variants from aligned haplotypes in one active region.
#[derive(Clone, Debug, Default)]
pub struct VariantCaller {
    active_region: Option<ActiveRegion>,
    variants: std::collections::BTreeMap<CallerVarKey, CallerVarNode>,
    total_depth: i32,
    variant_call_by_reference: bool,
    call_ambiguous_variant: bool,
}

impl VariantCaller {
    /// Default setting for emitting ambiguous variants.
    pub const DEFAULT_CALL_AMBIGUOUS_VARIANT: bool = true;

    /// Creates an uninitialized variant caller.
    #[must_use]
    pub fn new() -> Self {
        Self {
            active_region: None,
            variants: std::collections::BTreeMap::new(),
            total_depth: 0,
            variant_call_by_reference: true,
            call_ambiguous_variant: Self::DEFAULT_CALL_AMBIGUOUS_VARIANT,
        }
    }

    /// Initializes the caller for an active region and clears previous calls.
    pub fn init(&mut self, active_region: ActiveRegion) {
        self.active_region = Some(active_region);
        self.variants.clear();
        self.total_depth = 0;
    }

    /// Groups equivalent variants by reference coordinates.
    pub fn set_variant_call_by_reference(&mut self) {
        self.variant_call_by_reference = true;
    }

    /// Groups equivalent variants by active-region coordinates.
    pub fn set_variant_call_by_region(&mut self) {
        self.variant_call_by_reference = false;
    }

    /// Returns true when variants are grouped by reference coordinates.
    #[must_use]
    pub fn is_variant_call_by_reference(&self) -> bool {
        self.variant_call_by_reference
    }

    /// Sets whether ambiguous variants are emitted.
    pub fn set_call_ambiguous_variant(&mut self, call_ambiguous_variant: bool) {
        self.call_ambiguous_variant = call_ambiguous_variant;
    }

    /// Returns whether ambiguous variants are emitted.
    #[must_use]
    pub fn call_ambiguous_variant(&self) -> bool {
        self.call_ambiguous_variant
    }

    /// Adds one aligned haplotype and records its variant events.
    pub fn add(&mut self, haplotype: Haplotype) -> Result<(), VariantError> {
        let active_region = self
            .active_region
            .as_ref()
            .ok_or(VariantError::CallerNotInitialized)?
            .clone();
        if haplotype.active_region != active_region {
            return Err(VariantError::CallerRegionMismatch);
        }

        self.total_depth += haplotype.stats.min;
        let mut ref_position = active_region.start_index as usize;
        let mut alt_position = 0_usize;
        let mut align_node = Some(&haplotype.alignment);

        while let Some(node) = align_node {
            match node.node_type {
                AlignNode::MATCH => {
                    ref_position += node.n as usize;
                    alt_position += node.n as usize;
                }
                AlignNode::MISMATCH => {
                    for _ in 0..node.n {
                        let ref_base = char::from(active_region.ref_region.sequence[ref_position]);
                        let alt_base = char::from(haplotype.sequence[alt_position]);
                        self.add_node(CallerVarNode::new(
                            VariantType::Snp,
                            ref_position as i32 + 1,
                            ref_base.to_string(),
                            alt_base.to_string(),
                            haplotype.clone(),
                        ));
                        ref_position += 1;
                        alt_position += 1;
                    }
                }
                AlignNode::INS => {
                    let start = alt_position;
                    alt_position += node.n as usize;
                    let alt = String::from_utf8_lossy(&haplotype.sequence[start..alt_position])
                        .to_string();
                    self.add_node(CallerVarNode::new(
                        VariantType::Insertion,
                        ref_position as i32,
                        String::new(),
                        alt,
                        haplotype.clone(),
                    ));
                }
                AlignNode::DEL => {
                    let start_position = ref_position;
                    ref_position += node.n as usize;
                    let ref_allele = String::from_utf8_lossy(
                        &active_region.ref_region.sequence[start_position..ref_position],
                    )
                    .to_string();
                    self.add_node(CallerVarNode::new(
                        VariantType::Deletion,
                        start_position as i32 + 1,
                        ref_allele,
                        String::new(),
                        haplotype.clone(),
                    ));
                }
                other => return Err(VariantError::UnknownAlignmentNode(other)),
            }

            align_node = node.next.as_deref();
        }

        Ok(())
    }

    /// Returns called variants sorted by genomic order.
    #[must_use]
    pub fn variants(&self) -> Vec<Variant> {
        self.variants
            .values()
            .map(|node| node.to_variant(self.total_depth, self.variant_call_by_reference))
            .collect::<Result<Vec<_>, _>>()
            .expect("caller nodes are built from validated variant inputs")
    }

    fn add_node(&mut self, node: CallerVarNode) {
        if node.is_flank {
            return;
        }
        if node.is_ambiguous && !self.call_ambiguous_variant {
            return;
        }

        let key = node.key.clone();
        self.variants
            .entry(key)
            .and_modify(|existing| existing.add_variant(&node))
            .or_insert(node);
    }
}

#[derive(Clone, Debug, Eq, Ord, PartialEq, PartialOrd)]
struct CallerVarKey {
    start: i32,
    variant_type: VariantType,
    ref_allele: String,
    alt_allele: String,
}

#[derive(Clone, Debug)]
struct CallerVarNode {
    key: CallerVarKey,
    is_ambiguous: bool,
    is_flank: bool,
    var_depth: i32,
    haplotypes: Vec<Haplotype>,
}

impl CallerVarNode {
    fn new(
        variant_type: VariantType,
        start: i32,
        ref_allele: String,
        alt_allele: String,
        haplotype: Haplotype,
    ) -> Self {
        let is_ambiguous = ref_allele.bytes().any(|base| {
            !matches!(
                base,
                b'A' | b'C' | b'G' | b'T' | b'U' | b'a' | b'c' | b'g' | b't' | b'u'
            )
        });
        let end = if alt_allele.is_empty() {
            start
        } else {
            start + alt_allele.len() as i32 - 1
        };
        let is_flank = haplotype
            .active_region
            .ref_region
            .is_flank_by_coordinate(start, end)
            .unwrap_or(true);
        let var_depth = haplotype.stats.min;

        Self {
            key: CallerVarKey {
                start,
                variant_type,
                ref_allele,
                alt_allele,
            },
            is_ambiguous,
            is_flank,
            var_depth,
            haplotypes: vec![haplotype],
        }
    }

    fn add_variant(&mut self, other: &Self) {
        self.var_depth += other.var_depth;
        self.haplotypes.extend(other.haplotypes.iter().cloned());
    }

    fn to_variant(
        &self,
        total_depth: i32,
        variant_call_by_reference: bool,
    ) -> Result<Variant, VariantError> {
        let ref_region = &self.haplotypes[0].active_region.ref_region;
        let mut start = self.key.start - ref_region.left_flank_length;
        if variant_call_by_reference {
            start += ref_region.interval.start - 1;
        }

        match self.key.variant_type {
            VariantType::Snp => Ok(Variant::Snp(VariantSnp::new(
                start,
                self.var_depth,
                total_depth,
                self.key.ref_allele.clone(),
                self.key.alt_allele.clone(),
                &self.haplotypes,
                self.haplotypes.len(),
                variant_call_by_reference,
            )?)),
            VariantType::Insertion => Ok(Variant::Insertion(VariantInsertion::new(
                start,
                self.var_depth,
                total_depth,
                self.key.alt_allele.clone(),
                &self.haplotypes,
                self.haplotypes.len(),
                variant_call_by_reference,
            )?)),
            VariantType::Deletion => Ok(Variant::Deletion(VariantDeletion::new(
                start,
                self.var_depth,
                total_depth,
                self.key.ref_allele.clone(),
                &self.haplotypes,
                self.haplotypes.len(),
                variant_call_by_reference,
            )?)),
        }
    }
}

fn insertion_anchor(data: &VariantCallData) -> Result<u8, ReferenceSequenceError> {
    data.active_region
        .ref_region
        .get_base(data.start - 1)
        .or_else(|_| data.active_region.ref_region.get_base(data.start))
}

fn deletion_anchor(
    data: &VariantCallData,
    reference_end: i32,
) -> Result<u8, ReferenceSequenceError> {
    data.active_region
        .ref_region
        .get_base(data.start - 1)
        .or_else(|_| data.active_region.ref_region.get_base(reference_end + 1))
}

#[cfg(test)]
mod tests {
    use kanalyze::util::KmerUtil;

    use crate::activeregion::{ActiveRegion, RegionStats};
    use crate::align::AlignNode;
    use crate::refreader::{ReferenceRegion, ReferenceSequence};
    use crate::util::digest::Digest;

    use super::*;

    #[test]
    fn enum_has_three_values_in_java_order() {
        assert_eq!(
            VariantType::ALL,
            [
                VariantType::Snp,
                VariantType::Insertion,
                VariantType::Deletion
            ]
        );
        assert!(VariantType::Snp < VariantType::Insertion);
        assert!(VariantType::Insertion < VariantType::Deletion);
    }

    #[test]
    fn hash_codes_match_java_octal_literals() {
        assert_eq!(VariantType::Snp.java_hash_code(), 0o11111111111);
        assert_eq!(
            VariantType::Insertion.java_hash_code(),
            0o22222222222_u32 as i32
        );
        assert_eq!(
            VariantType::Deletion.java_hash_code(),
            0o33333333333_u32 as i32
        );
        assert_ne!(
            VariantType::Snp.java_hash_code(),
            VariantType::Insertion.java_hash_code()
        );
        assert_ne!(
            VariantType::Insertion.java_hash_code(),
            VariantType::Deletion.java_hash_code()
        );
    }

    #[test]
    fn variant_base_validation_matches_java_boundaries() {
        let haplotype = haplotype("AAAACCCCGGGGTTTT", 0, 12, b"AAAACCCCGGGGTTTT");
        assert_eq!(
            VariantSnp::new(5, 1, 5, "C", "A", &[], 0, true).unwrap_err(),
            VariantError::EmptyHaplotypeList
        );
        assert_eq!(
            VariantSnp::new(5, 1, 5, "C", "A", std::slice::from_ref(&haplotype), 0, true)
                .unwrap_err(),
            VariantError::InvalidHaplotypeSize {
                haplo_size: 0,
                len: 1
            }
        );
        assert_eq!(
            VariantSnp::new(0, 1, 5, "C", "A", std::slice::from_ref(&haplotype), 1, true)
                .unwrap_err(),
            VariantError::InvalidStart(0)
        );
        assert_eq!(
            VariantSnp::new(
                5,
                -1,
                5,
                "C",
                "A",
                std::slice::from_ref(&haplotype),
                1,
                true
            )
            .unwrap_err(),
            VariantError::NegativeVariantDepth(-1)
        );
        assert_eq!(
            VariantSnp::new(5, 10, 5, "C", "A", &[haplotype], 1, true).unwrap_err(),
            VariantError::LocusDepthLessThanVariantDepth {
                locus_depth: 5,
                variant_depth: 10
            }
        );
    }

    #[test]
    fn snp_fields_vcf_hgvs_and_ambiguity_match_java_tests() {
        let haplotype = haplotype("AAAACCCCGGGGTTTT", 0, 12, b"AAAACCCCGGGGTTTT");
        let snp = VariantSnp::new(5, 5, 10, "N", "A", &[haplotype], 1, true).unwrap();

        assert_eq!(snp.data.variant_type, VariantType::Snp);
        assert_eq!(snp.vcf_pos(), 5);
        assert_eq!(snp.vcf_ref().unwrap(), "N");
        assert_eq!(snp.vcf_alt().unwrap(), "A");
        assert_eq!(snp.reference_end(), 5);
        assert_eq!(snp.to_string(), "5N>A");
        assert!(snp.to_string_full().contains("chr1"));
        assert!(snp.data.is_ambiguous);
        assert!(snp.data.is_reference_aligned);
        assert_eq!(snp.data.haplotypes().len(), 1);
    }

    #[test]
    fn snp_rejects_non_single_base_alleles() {
        let haplotype = haplotype("AAAACCCCGGGGTTTT", 0, 12, b"AAAACCCCGGGGTTTT");
        assert!(matches!(
            VariantSnp::new(
                5,
                5,
                10,
                "CC",
                "A",
                std::slice::from_ref(&haplotype),
                1,
                true
            ),
            Err(VariantError::InvalidSnpRef(_))
        ));
        assert!(matches!(
            VariantSnp::new(5, 5, 10, "C", "AA", &[haplotype], 1, true),
            Err(VariantError::InvalidSnpAlt(_))
        ));
    }

    #[test]
    fn variants_compare_by_reference_start_type_ref_alt() {
        let haplotype = haplotype("AAAACCCCGGGGTTTT", 0, 12, b"AAAACCCCGGGGTTTT");
        let early =
            VariantSnp::new(3, 1, 5, "G", "A", std::slice::from_ref(&haplotype), 1, true).unwrap();
        let late = VariantSnp::new(
            10,
            1,
            5,
            "G",
            "A",
            std::slice::from_ref(&haplotype),
            1,
            true,
        )
        .unwrap();
        let alt_a =
            VariantSnp::new(5, 1, 5, "A", "G", std::slice::from_ref(&haplotype), 1, true).unwrap();
        let alt_t = VariantSnp::new(5, 1, 5, "A", "T", &[haplotype], 1, true).unwrap();

        assert!(early < late);
        assert!(late > early);
        assert!(alt_a < alt_t);
    }

    #[test]
    fn insertion_vcf_and_hgvs_match_java_tests() {
        let haplotype = haplotype("AAAACCCCGGGGTTTT", 0, 12, b"AAAACCCCGGGGTTTT");
        let insertion =
            VariantInsertion::new(5, 5, 10, "GGG", std::slice::from_ref(&haplotype), 1, true)
                .unwrap();
        let at_start = VariantInsertion::new(1, 5, 10, "ACG", &[haplotype], 1, true).unwrap();

        assert_eq!(insertion.data.variant_type, VariantType::Insertion);
        assert_eq!(insertion.length, 3);
        assert_eq!(insertion.reference_end(), 5);
        assert_eq!(insertion.vcf_pos(), 4);
        assert_eq!(insertion.vcf_ref().unwrap(), "A");
        assert_eq!(insertion.vcf_alt().unwrap(), "AGGG");
        assert_eq!(insertion.to_string(), "5_6insGGG");
        assert_eq!(at_start.vcf_pos(), 1);
        assert_eq!(at_start.vcf_alt().unwrap(), "ACGA");
    }

    #[test]
    fn insertion_rejects_empty_alt() {
        let haplotype = haplotype("AAAACCCCGGGGTTTT", 0, 12, b"AAAACCCCGGGGTTTT");
        assert_eq!(
            VariantInsertion::new(5, 5, 10, "", &[haplotype], 1, true).unwrap_err(),
            VariantError::EmptyInsertionAlt
        );
    }

    #[test]
    fn deletion_preserves_java_type_bug_and_vcf_behavior() {
        let haplotype = haplotype("AAAACCCCGGGGTTTT", 0, 12, b"AAAACCCCGGGGTTTT");
        let deletion =
            VariantDeletion::new(5, 3, 10, "CC", std::slice::from_ref(&haplotype), 1, true)
                .unwrap();
        let at_start = VariantDeletion::new(1, 3, 10, "AA", &[haplotype], 1, true).unwrap();

        assert_eq!(deletion.data.variant_type, VariantType::Insertion);
        assert_ne!(deletion.data.variant_type, VariantType::Deletion);
        assert_eq!(deletion.length, 2);
        assert_eq!(deletion.reference_end(), 6);
        assert_eq!(deletion.vcf_pos(), 4);
        assert_eq!(deletion.vcf_ref().unwrap(), "ACC");
        assert_eq!(deletion.vcf_alt().unwrap(), "A");
        assert_eq!(deletion.to_string(), "5_6del");
        assert_eq!(at_start.vcf_pos(), 1);
        assert_eq!(at_start.vcf_ref().unwrap(), "AAA");
        assert_eq!(at_start.vcf_alt().unwrap(), "A");
    }

    #[test]
    fn deletion_rejects_empty_ref() {
        let haplotype = haplotype("AAAACCCCGGGGTTTT", 0, 12, b"AAAACCCCGGGGTTTT");
        assert_eq!(
            VariantDeletion::new(5, 3, 10, "", &[haplotype], 1, true).unwrap_err(),
            VariantError::EmptyDeletionRef
        );
    }

    #[test]
    fn variant_caller_calls_snp_from_mismatch_node() {
        let mismatch = AlignNode::new(AlignNode::MISMATCH, 1, None).unwrap();
        let head = AlignNode::new(AlignNode::MATCH, 4, Some(Box::new(mismatch))).unwrap();
        let haplotype =
            haplotype_with_alignment("AAAACCCCGGGGTTTT", 0, 12, b"AAAATCCCGGGGTTTT", head);
        let mut caller = VariantCaller::new();
        caller.init(haplotype.active_region.clone());
        caller.add(haplotype).unwrap();

        let variants = caller.variants();
        assert_eq!(variants.len(), 1);
        match &variants[0] {
            Variant::Snp(snp) => {
                assert_eq!(snp.data.start, 5);
                assert_eq!(snp.data.ref_allele, "C");
                assert_eq!(snp.data.alt_allele, "T");
                assert_eq!(snp.data.variant_depth, 10);
                assert_eq!(snp.data.locus_depth, 10);
            }
            other => panic!("expected SNP, got {other:?}"),
        }
    }

    #[test]
    fn variant_caller_calls_insertion_and_deletion() {
        let ins_tail = AlignNode::new(AlignNode::MATCH, 11, None).unwrap();
        let ins = AlignNode::new(AlignNode::INS, 2, Some(Box::new(ins_tail))).unwrap();
        let ins_head = AlignNode::new(AlignNode::MATCH, 5, Some(Box::new(ins))).unwrap();
        let insertion_haplotype =
            haplotype_with_alignment("AAAACCCCGGGGTTTT", 0, 12, b"AAAACGGCCCGGGGTTTT", ins_head);
        let mut caller = VariantCaller::new();
        caller.init(insertion_haplotype.active_region.clone());
        caller.add(insertion_haplotype).unwrap();
        let variants = caller.variants();
        assert_eq!(variants.len(), 1);
        assert!(matches!(&variants[0], Variant::Insertion(_)));
        assert_eq!(variants[0].data().start, 5);
        assert_eq!(variants[0].data().alt_allele, "GG");

        let del_tail = AlignNode::new(AlignNode::MATCH, 10, None).unwrap();
        let del = AlignNode::new(AlignNode::DEL, 2, Some(Box::new(del_tail))).unwrap();
        let del_head = AlignNode::new(AlignNode::MATCH, 4, Some(Box::new(del))).unwrap();
        let deletion_haplotype =
            haplotype_with_alignment("AAAACCCCGGGGTTTT", 0, 12, b"AAAACCGGGGTTTT", del_head);
        let mut caller = VariantCaller::new();
        caller.init(deletion_haplotype.active_region.clone());
        caller.add(deletion_haplotype).unwrap();
        let variants = caller.variants();
        assert_eq!(variants.len(), 1);
        assert!(matches!(&variants[0], Variant::Deletion(_)));
        assert_eq!(variants[0].data().start, 5);
        assert_eq!(variants[0].data().ref_allele, "CC");
    }

    #[test]
    fn variant_caller_merges_equivalent_variants_and_tracks_depth() {
        let mismatch = AlignNode::new(AlignNode::MISMATCH, 1, None).unwrap();
        let head = AlignNode::new(AlignNode::MATCH, 4, Some(Box::new(mismatch))).unwrap();
        let haplotype_a =
            haplotype_with_alignment("AAAACCCCGGGGTTTT", 0, 12, b"AAAATCCCGGGGTTTT", head.clone());
        let haplotype_b =
            haplotype_with_alignment("AAAACCCCGGGGTTTT", 0, 12, b"AAAATCCCGGGGTTTT", head);
        let mut caller = VariantCaller::new();
        caller.init(haplotype_a.active_region.clone());
        caller.add(haplotype_a).unwrap();
        caller.add(haplotype_b).unwrap();

        let variants = caller.variants();
        assert_eq!(variants.len(), 1);
        assert_eq!(variants[0].data().variant_depth, 20);
        assert_eq!(variants[0].data().locus_depth, 20);
        assert_eq!(variants[0].data().haplotypes().len(), 2);
    }

    #[test]
    fn variant_caller_requires_init_and_matching_region() {
        let matching_haplotype = haplotype("AAAACCCCGGGGTTTT", 0, 12, b"AAAACCCCGGGGTTTT");
        let mut caller = VariantCaller::new();
        assert_eq!(
            caller.add(matching_haplotype.clone()),
            Err(VariantError::CallerNotInitialized)
        );

        caller.init(haplotype("AAAACCCCGGGGTTTT", 1, 12, b"AAAACCCCGGGGTTTT").active_region);
        assert_eq!(
            caller.add(matching_haplotype),
            Err(VariantError::CallerRegionMismatch)
        );
    }

    fn haplotype(ref_seq: &str, start_kmer: i32, end_kmer: i32, consensus: &[u8]) -> Haplotype {
        let alignment = AlignNode::new(AlignNode::MATCH, consensus.len() as i32, None).unwrap();
        haplotype_with_alignment(ref_seq, start_kmer, end_kmer, consensus, alignment)
    }

    fn haplotype_with_alignment(
        ref_seq: &str,
        start_kmer: i32,
        end_kmer: i32,
        consensus: &[u8],
        alignment: AlignNode,
    ) -> Haplotype {
        let ref_region = reference_region(ref_seq);
        let count = vec![10; ref_seq.len() - 4 + 1];
        let kmer_util = KmerUtil::new(4).unwrap();
        let active_region =
            ActiveRegion::new(ref_region, start_kmer, end_kmer, &count, &kmer_util).unwrap();
        let stats = RegionStats::from_counts(
            &count,
            active_region.start_kmer_index,
            active_region.end_kmer_index,
        )
        .unwrap();

        Haplotype::new(
            consensus.to_vec(),
            active_region,
            vec![alignment],
            100.0,
            None,
            stats,
        )
        .unwrap()
    }

    fn reference_region(sequence: &str) -> ReferenceRegion {
        let digest = Digest::new(
            (0..16)
                .map(|index| sequence.len() as u8 + index)
                .collect::<Vec<_>>(),
            "MD5",
        )
        .unwrap();
        let reference_sequence =
            ReferenceSequence::new("chr1", sequence.len() as i32, Some(digest), Some("test"))
                .unwrap();
        ReferenceRegion::whole(reference_sequence, sequence.as_bytes(), 0).unwrap()
    }
}
