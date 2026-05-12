use std::fmt;

use kanalyze::Base;
use kanalyze::util::{KmerHashSet, KmerKey, KmerUtil};
use thiserror::Error;

use crate::activeregion::{ActiveRegion, ActiveRegionError, Haplotype, RegionStats};
use crate::constants::{ARRAY_EXPAND_FACTOR, MAX_ARRAY_SIZE, MIN_KMER_SIZE};
use crate::counter::CountMap;
use crate::util::number::is_zero;

/// Errors returned while constructing trace nodes.
#[derive(Clone, Debug, Error, PartialEq)]
pub enum TraceNodeError {
    /// Alignment score was negative.
    #[error("score is negative")]
    NegativeScore,
    /// Trace node type was outside the valid range.
    #[error("type is out of range [{min}, {max}]: {value}")]
    TypeOutOfRange {
        /// Provided node type.
        value: i8,
        /// Minimum valid node type.
        min: i8,
        /// Maximum valid node type.
        max: i8,
    },
}

/// A traceback node in the k-mer alignment dynamic-programming graph.
#[derive(Clone, Debug, PartialEq)]
pub struct TraceNode {
    /// Alignment score at this node.
    pub score: f32,
    /// Trace node type.
    pub node_type: i8,
    /// Previous node in the selected trace path.
    pub next_node: Option<Box<TraceNode>>,
    /// Alternate predecessor with the same score.
    pub branch_node: Option<Box<TraceNode>>,
}

impl TraceNode {
    /// Sentinel node type.
    pub const TYPE_NONE: i8 = 0;
    /// Match node type.
    pub const TYPE_MATCH: i8 = 1;
    /// Mismatch node type.
    pub const TYPE_MISMATCH: i8 = 2;
    /// Gap in the reference sequence.
    pub const TYPE_GAP_REF: i8 = 3;
    /// Gap in the consensus sequence.
    pub const TYPE_GAP_CON: i8 = 4;
    /// CIGAR characters indexed by trace node type.
    pub const CIGAR_CHARS: [char; 5] = ['*', '=', 'X', 'I', 'D'];
    /// Zero-valued sentinel trace node.
    pub const ZERO_NODE: Self = Self {
        score: 0.0,
        node_type: Self::TYPE_NONE,
        next_node: None,
        branch_node: None,
    };

    /// Creates a validated trace node.
    pub fn new(
        score: f32,
        node_type: i8,
        next_node: Option<Box<Self>>,
        branch_node: Option<Box<Self>>,
    ) -> Result<Self, TraceNodeError> {
        if score < 0.0 {
            return Err(TraceNodeError::NegativeScore);
        }
        if !(Self::TYPE_NONE..=Self::TYPE_GAP_CON).contains(&node_type) {
            return Err(TraceNodeError::TypeOutOfRange {
                value: node_type,
                min: Self::TYPE_NONE,
                max: Self::TYPE_GAP_CON,
            });
        }

        Ok(Self {
            score,
            node_type,
            next_node,
            branch_node,
        })
    }

    /// Creates a trace node without a branch predecessor.
    #[must_use]
    pub fn with_next(score: f32, node_type: i8, next_node: Option<Box<Self>>) -> Self {
        Self {
            score,
            node_type,
            next_node,
            branch_node: None,
        }
    }

    /// Returns the CIGAR character table.
    #[must_use]
    pub fn cigar_array() -> Vec<char> {
        Self::CIGAR_CHARS.to_vec()
    }

    /// Returns a display name for this node's trace type.
    #[must_use]
    pub fn type_string(&self) -> &'static str {
        match self.node_type {
            Self::TYPE_NONE => "NONE",
            Self::TYPE_MATCH => "ALIGN_MATCH",
            Self::TYPE_MISMATCH => "ALIGN_MISMATCH",
            Self::TYPE_GAP_REF => "GAP_REFERENCE",
            Self::TYPE_GAP_CON => "GAP_CONSENSUS",
            _ => "UNKNOWN",
        }
    }
}

impl fmt::Display for TraceNode {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "TraceNode[score={:.6}, type={}]",
            self.score,
            self.type_string()
        )
    }
}

/// Errors returned while constructing alignment nodes.
#[derive(Clone, Debug, Error, Eq, PartialEq)]
pub enum AlignNodeError {
    /// Alignment node type was outside the valid range.
    #[error("type is out of range [1, 4]: {0}")]
    TypeOutOfRange(i8),
}

/// A run-length encoded alignment operation.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct AlignNode {
    /// Alignment node type.
    pub node_type: i8,
    /// CIGAR operation character.
    pub cigar_char: char,
    /// Run length for this operation.
    pub n: i32,
    /// Next alignment operation.
    pub next: Option<Box<AlignNode>>,
}

impl AlignNode {
    /// Match operation.
    pub const MATCH: i8 = TraceNode::TYPE_MATCH;
    /// Mismatch operation.
    pub const MISMATCH: i8 = TraceNode::TYPE_MISMATCH;
    /// Insertion relative to the reference.
    pub const INS: i8 = TraceNode::TYPE_GAP_REF;
    /// Deletion relative to the reference.
    pub const DEL: i8 = TraceNode::TYPE_GAP_CON;

    /// Creates a run-length encoded alignment operation.
    pub fn new(node_type: i8, n: i32, next: Option<Box<Self>>) -> Result<Self, AlignNodeError> {
        if !(Self::MATCH..=Self::DEL).contains(&node_type) {
            return Err(AlignNodeError::TypeOutOfRange(node_type));
        }

        Ok(Self {
            node_type,
            cigar_char: TraceNode::CIGAR_CHARS[node_type as usize],
            n,
            next,
        })
    }

    /// Renders this alignment chain as a CIGAR string.
    #[must_use]
    pub fn cigar_string(&self) -> String {
        let mut cigar = String::new();
        let mut node = Some(self);
        while let Some(current) = node {
            cigar.push_str(&current.n.to_string());
            cigar.push(current.cigar_char);
            node = current.next.as_deref();
        }
        cigar
    }

    /// Compares two alignment chains using Kestrel's tie-breaking order.
    #[must_use]
    pub fn compare_to(&self, other: &Self) -> i32 {
        let mut t_node = Some(self);
        let mut o_node = Some(other);

        while let (Some(t_current), Some(o_current)) = (t_node, o_node) {
            let mut t_type = t_current.node_type;
            let mut o_type = o_current.node_type;
            let n_diff = t_current.n - o_current.n;

            if t_type != o_type || n_diff != 0 {
                if t_type == o_type {
                    if n_diff > 0 {
                        o_type = o_current
                            .next
                            .as_deref()
                            .map_or(TraceNode::TYPE_NONE, |node| node.node_type);
                    } else {
                        t_type = t_current
                            .next
                            .as_deref()
                            .map_or(TraceNode::TYPE_NONE, |node| node.node_type);
                    }
                }

                return compare_rank(t_type) - compare_rank(o_type);
            }

            t_node = t_current.next.as_deref();
            o_node = o_current.next.as_deref();
        }

        if t_node.is_none() {
            return -1;
        }
        if o_node.is_none() {
            return 1;
        }
        0
    }
}

fn compare_rank(node_type: i8) -> i32 {
    if node_type == AlignNode::MATCH {
        5
    } else {
        i32::from(node_type)
    }
}

/// Errors returned while tracking maximum-scoring alignment nodes.
#[derive(Clone, Debug, Error, Eq, PartialEq)]
pub enum MaxAlignmentScoreNodeError {
    /// No trace node was provided.
    #[error("traceNode is null")]
    NullTraceNode,
    /// Consensus base count must be positive.
    #[error("nConsensusBases is less than 1: {0}")]
    InvalidConsensusBases(i32),
}

/// Errors returned by [`HaplotypeContainer`].
#[derive(Clone, Debug, Error, Eq, PartialEq)]
pub enum HaplotypeContainerError {
    /// Haplotype limit must be positive.
    #[error("Haplotype container limit is less than 1: {0}")]
    InvalidLimit(i32),
}

/// Bounded collection that keeps the deepest haplotypes.
#[derive(Clone, Debug, PartialEq)]
pub struct HaplotypeContainer {
    haplotypes: Vec<Haplotype>,
    /// Maximum number of haplotypes retained.
    pub limit: usize,
}

impl HaplotypeContainer {
    /// Creates an empty haplotype container with a positive limit.
    pub fn new(limit: i32) -> Result<Self, HaplotypeContainerError> {
        if limit < 1 {
            return Err(HaplotypeContainerError::InvalidLimit(limit));
        }

        Ok(Self {
            haplotypes: Vec::new(),
            limit: limit as usize,
        })
    }

    /// Adds a haplotype, replacing the current shallowest haplotype when full.
    pub fn add(&mut self, haplotype: Haplotype) {
        if self.haplotypes.len() >= self.limit && !self.remove_min_haplotype(haplotype.stats.min) {
            return;
        }

        self.haplotypes.insert(0, haplotype);
    }

    /// Returns the number of retained haplotypes.
    #[must_use]
    pub fn size(&self) -> usize {
        self.haplotypes.len()
    }

    /// Returns retained haplotypes as a cloned vector.
    #[must_use]
    pub fn to_array(&self) -> Vec<Haplotype> {
        self.haplotypes.clone()
    }

    fn remove_min_haplotype(&mut self, limit: i32) -> bool {
        let Some((min_index, _)) = self
            .haplotypes
            .iter()
            .enumerate()
            .filter(|(_, haplotype)| haplotype.stats.min < limit)
            .min_by_key(|(_, haplotype)| haplotype.stats.min)
        else {
            return false;
        };

        self.haplotypes.remove(min_index);
        true
    }
}

/// A maximum-scoring trace node and the consensus length that produced it.
#[derive(Clone, Debug, PartialEq)]
pub struct MaxAlignmentScoreNode {
    /// Trace node at the maximum score.
    pub trace_node: Box<TraceNode>,
    /// Number of consensus bases consumed at this score.
    pub n_consensus_bases: i32,
    /// Next maximum-score node with the same score.
    pub next: Option<Box<MaxAlignmentScoreNode>>,
    /// True once the node has been converted into a haplotype.
    pub haplotype_built: bool,
}

impl MaxAlignmentScoreNode {
    /// Creates a maximum-score node.
    pub fn new(
        trace_node: Option<TraceNode>,
        n_consensus_bases: i32,
        next: Option<Box<Self>>,
    ) -> Result<Self, MaxAlignmentScoreNodeError> {
        let trace_node = trace_node.ok_or(MaxAlignmentScoreNodeError::NullTraceNode)?;
        if n_consensus_bases < 1 {
            return Err(MaxAlignmentScoreNodeError::InvalidConsensusBases(
                n_consensus_bases,
            ));
        }

        Ok(Self {
            trace_node: Box::new(trace_node),
            n_consensus_bases,
            next,
            haplotype_built: false,
        })
    }
}

impl fmt::Display for MaxAlignmentScoreNode {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "MaxAlignment[len={}, last={}, trace={}]",
            self.n_consensus_bases, self.trace_node.node_type, self.trace_node
        )
    }
}

/// Errors returned by [`TraceMatrix`].
#[derive(Clone, Debug, Error, Eq, PartialEq)]
pub enum TraceMatrixError {
    /// Trace matrix row count must be positive.
    #[error("Cannot create trace matrix with less than 1 row: {0}")]
    InvalidRowCount(i32),
    /// Requested row was out of range.
    #[error("Row is out of range [0, {max}]: {row}")]
    RowOutOfRange {
        /// Requested row.
        row: i32,
        /// Maximum valid row.
        max: i32,
    },
    /// No column has been allocated yet.
    #[error("No columns in matrix: must call nextCol() before set()")]
    NoColumns,
    /// Matrix cannot expand beyond the maximum capacity.
    #[error("Cannot expand matrix: Column capacity is at its maximum size: {0}")]
    MaximumCapacity(usize),
}

/// Compact diagnostic matrix encoding traceback transitions.
#[derive(Clone, Debug, PartialEq)]
pub struct TraceMatrix {
    /// Matrix storage, indexed by row then column.
    pub matrix: Vec<Vec<i16>>,
    /// Number of rows in the matrix.
    pub n_row: usize,
    column: i32,
    col_capacity: usize,
}

impl TraceMatrix {
    /// Initial column capacity.
    pub const DEFAULT_COLUMN_CAPACITY: usize = 250;

    /// Creates a trace matrix with `n_row` rows.
    pub fn new(n_row: i32) -> Result<Self, TraceMatrixError> {
        if n_row < 1 {
            return Err(TraceMatrixError::InvalidRowCount(n_row));
        }

        let n_row = n_row as usize;
        Ok(Self {
            matrix: vec![vec![0; Self::DEFAULT_COLUMN_CAPACITY]; n_row],
            n_row,
            column: -1,
            col_capacity: Self::DEFAULT_COLUMN_CAPACITY,
        })
    }

    /// Sets a transition bit in the current column for a row.
    pub fn set(
        &mut self,
        row: i32,
        mut type_from: i32,
        mut type_to: i32,
    ) -> Result<(), TraceMatrixError> {
        if row < 0 || row as usize >= self.n_row {
            return Err(TraceMatrixError::RowOutOfRange {
                row,
                max: self.n_row as i32 - 1,
            });
        }
        if self.column < 0 {
            return Err(TraceMatrixError::NoColumns);
        }

        if type_from > 1 {
            type_from -= 1;
        }
        if type_to > 1 {
            type_to -= 1;
        }
        let bit = 0x1_i16 << ((3 - type_from) + 3 * (3 - type_to));
        self.matrix[row as usize][self.column as usize] |= bit;
        Ok(())
    }

    /// Advances to the next column, expanding storage if needed.
    pub fn next_col(&mut self) -> Result<i32, TraceMatrixError> {
        self.column += 1;
        if self.column as usize >= self.col_capacity
            && let Err(err) = self.expand_matrix()
        {
            self.column -= 1;
            return Err(err);
        }
        Ok(self.column)
    }

    /// Returns the current column index.
    #[must_use]
    pub fn column(&self) -> i32 {
        self.column
    }

    fn expand_matrix(&mut self) -> Result<(), TraceMatrixError> {
        if self.col_capacity == MAX_ARRAY_SIZE {
            return Err(TraceMatrixError::MaximumCapacity(MAX_ARRAY_SIZE));
        }

        let mut new_capacity = ((self.col_capacity as f32) * ARRAY_EXPAND_FACTOR) as usize;
        if new_capacity <= self.col_capacity {
            new_capacity = self.col_capacity + 1;
        }
        new_capacity = new_capacity.min(MAX_ARRAY_SIZE);

        for row in &mut self.matrix {
            row.resize(new_capacity, 0);
        }
        self.col_capacity = new_capacity;
        Ok(())
    }

    /// Renders the trace matrix, optionally annotating rows and columns with sequences.
    #[must_use]
    pub fn matrix_string(
        &self,
        ref_sequence: Option<&[u8]>,
        con_sequence: Option<&[u8]>,
        ref_start: i32,
    ) -> String {
        let ref_start = ref_start.max(0) as usize;
        let write_ref_seq =
            ref_sequence.is_some_and(|sequence| sequence.len() >= self.n_row + ref_start);
        let write_con_seq = con_sequence
            .is_some_and(|sequence| self.column >= 0 && sequence.len() > self.column as usize);
        let row_space = digit_width(self.n_row.saturating_sub(1));
        let mut builder = String::new();

        if write_con_seq {
            if write_ref_seq {
                builder.push_str("  ");
            }
            builder.push_str(&" ".repeat(row_space));
            if let Some(sequence) = con_sequence {
                for base in sequence.iter().take(self.column as usize + 1) {
                    builder.push_str("      ");
                    builder.push(char::from(*base));
                    builder.push_str("     ");
                }
            }
            builder.push('\n');
        }

        if write_ref_seq {
            builder.push_str("  ");
        }
        builder.push_str(&" ".repeat(row_space));
        for col in 0..=self.column {
            let n_width = digit_width(col as usize);
            let left = (11 - n_width) / 2;
            let right = 11 - left - n_width;
            builder.push(' ');
            builder.push_str(&" ".repeat(left));
            builder.push_str(&col.to_string());
            builder.push_str(&" ".repeat(right));
        }
        builder.push_str("\n\n");

        for row in 0..self.n_row {
            if let Some(sequence) = ref_sequence.filter(|_| write_ref_seq) {
                builder.push(char::from(sequence[ref_start + row]));
                builder.push(' ');
            }
            builder.push_str(&format!("{row:>row_space$}"));

            for col in 0..=self.column {
                builder.push(' ');
                push_trace_bits(&mut builder, self.matrix[row][col as usize]);
            }
            builder.push('\n');
        }

        builder
    }
}

impl fmt::Display for TraceMatrix {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.matrix_string(None, None, 0))
    }
}

fn digit_width(value: usize) -> usize {
    value.to_string().len()
}

fn push_trace_bits(builder: &mut String, value: i16) {
    builder.push(if ((value >> 8) & 0x1) != 0 { 'a' } else { '-' });
    builder.push(if ((value >> 7) & 0x1) != 0 { 'r' } else { '-' });
    builder.push(if ((value >> 6) & 0x1) != 0 { 'c' } else { '-' });
    builder.push(',');
    builder.push(if ((value >> 5) & 0x1) != 0 { 'a' } else { '-' });
    builder.push(if ((value >> 4) & 0x1) != 0 { 'r' } else { '-' });
    builder.push(if ((value >> 3) & 0x1) != 0 { 'c' } else { '-' });
    builder.push(',');
    builder.push(if ((value >> 2) & 0x1) != 0 { 'a' } else { '-' });
    builder.push(if ((value >> 1) & 0x1) != 0 { 'r' } else { '-' });
    builder.push(if (value & 0x1) != 0 { 'c' } else { '-' });
}

/// Errors returned by [`TypeList`].
#[derive(Clone, Debug, Error, Eq, PartialEq)]
pub enum TypeListError {
    /// Alignment type was outside the valid range.
    #[error("type is out of range [1, 4]: {0}")]
    TypeOutOfRange(i8),
    /// Type list cannot expand beyond the maximum capacity.
    #[error("Cannot expand type list: Reached maximum capacity: {0}")]
    MaximumCapacity(usize),
}

/// Run-length encoded list of alignment operation types.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct TypeList {
    types: Vec<i8>,
    counts: Vec<i32>,
    capacity: usize,
}

impl TypeList {
    /// Initial operation capacity.
    pub const INITIAL_CAPACITY: usize = 100;

    /// Creates an empty type list.
    #[must_use]
    pub fn new() -> Self {
        let mut types = Vec::with_capacity(Self::INITIAL_CAPACITY);
        let mut counts = Vec::with_capacity(Self::INITIAL_CAPACITY);
        types.push(0);
        counts.push(0);

        Self {
            types,
            counts,
            capacity: Self::INITIAL_CAPACITY,
        }
    }

    /// Clones an existing type list, or creates an empty list when absent.
    #[must_use]
    pub fn new_from(other: Option<&Self>) -> Self {
        other.cloned().unwrap_or_default()
    }

    /// Adds one operation type, coalescing adjacent identical operations.
    pub fn add(&mut self, type_byte: i8) -> Result<(), TypeListError> {
        if !(AlignNode::MATCH..=AlignNode::DEL).contains(&type_byte) {
            return Err(TypeListError::TypeOutOfRange(type_byte));
        }

        let last_index = self.types.len() - 1;
        if self.types[last_index] == type_byte {
            self.counts[last_index] += 1;
            return Ok(());
        }

        if self.types.len() == self.capacity {
            self.expand()?;
        }

        self.types.push(type_byte);
        self.counts.push(1);
        Ok(())
    }

    /// Converts this run-length list into an alignment chain.
    #[must_use]
    pub fn to_alignment(&self, reverse: bool) -> Option<Box<AlignNode>> {
        let mut align_head = None;
        if reverse {
            for index in 1..self.types.len() {
                align_head = Some(Box::new(
                    AlignNode::new(self.types[index], self.counts[index], align_head)
                        .expect("TypeList only stores valid alignment node types"),
                ));
            }
        } else {
            for index in (1..self.types.len()).rev() {
                align_head = Some(Box::new(
                    AlignNode::new(self.types[index], self.counts[index], align_head)
                        .expect("TypeList only stores valid alignment node types"),
                ));
            }
        }
        align_head
    }

    /// Clears all stored operation runs.
    pub fn clear(&mut self) {
        self.types.truncate(1);
        self.counts.truncate(1);
    }

    /// Replaces this list with another list, or clears it when absent.
    pub fn clone_from_option(&mut self, other: Option<&Self>) {
        if let Some(other) = other {
            self.types = other.types.clone();
            self.counts = other.counts.clone();
            self.capacity = other.capacity;
        } else {
            self.clear();
        }
    }

    fn expand(&mut self) -> Result<(), TypeListError> {
        if self.capacity == MAX_ARRAY_SIZE {
            return Err(TypeListError::MaximumCapacity(MAX_ARRAY_SIZE));
        }

        let mut new_capacity = ((self.capacity as f32) * ARRAY_EXPAND_FACTOR) as usize;
        if new_capacity <= self.capacity {
            new_capacity = self.capacity + 1;
        }
        new_capacity = new_capacity.min(MAX_ARRAY_SIZE);
        self.types.reserve_exact(new_capacity - self.capacity);
        self.counts.reserve_exact(new_capacity - self.capacity);
        self.capacity = new_capacity;
        Ok(())
    }
}

impl Default for TypeList {
    fn default() -> Self {
        Self::new()
    }
}

/// Errors returned by [`KmerAlignmentBuilder`].
#[derive(Clone, Debug, Error, Eq, PartialEq)]
pub enum KmerAlignmentBuilderError {
    /// Alignment type-list error.
    #[error(transparent)]
    TypeList(#[from] TypeListError),
}

/// Errors returned by [`KmerAligner`].
#[derive(Clone, Debug, Error, PartialEq)]
pub enum KmerAlignerError {
    /// K-mer size was smaller than supported.
    #[error("k-mer size is less than {min}: {actual}")]
    KmerSizeTooSmall {
        /// Minimum supported k-mer size.
        min: usize,
        /// Actual k-mer size.
        actual: usize,
    },
    /// No active region was provided.
    #[error("Cannot initialize alignment with active region: null")]
    NullActiveRegion,
    /// Aligner has not been initialized.
    #[error("Aligner was not initialized (must call init())")]
    NotInitialized,
    /// Base argument was absent.
    #[error("Cannot add base to alignment: null")]
    NullBase,
    /// Consensus sequence cannot grow further.
    #[error("Consensus sequence has reached maximum length: {0}")]
    ConsensusAtMaximum(usize),
    /// Maximum saved-state count must be positive.
    #[error("Maximum number of states must not be less than 1: {0}")]
    InvalidMaxState(i32),
    /// K-mer argument was absent.
    #[error("Cannot save state for k-mer: null")]
    NullKmer,
    /// Next base argument was absent.
    #[error("Cannot save state with next base: null")]
    NullNextBase,
    /// Minimum depth must be positive.
    #[error("Minimum depth may not be less than 1: {0}")]
    InvalidMinDepth(i32),
    /// K-mer hash argument was absent.
    #[error("Cannot save state with k-mer hash: null")]
    NullKmerHash,
    /// Repeat count cannot be negative.
    #[error("Repeat count must be non-negative: {0}")]
    NegativeRepeatCount(i32),
    /// Alignment-weight error.
    #[error(transparent)]
    AlignmentWeight(#[from] AlignmentWeightError),
    /// Trace-matrix error.
    #[error(transparent)]
    TraceMatrix(#[from] TraceMatrixError),
    /// Alignment-builder error.
    #[error(transparent)]
    Builder(#[from] KmerAlignmentBuilderError),
    /// Active-region error.
    #[error(transparent)]
    ActiveRegion(#[from] ActiveRegionError),
}

/// Builds alignment chains from traceback branches.
#[derive(Clone, Debug, Default)]
pub struct KmerAlignmentBuilder {
    type_list_cache: Vec<TypeList>,
}

impl KmerAlignmentBuilder {
    /// Default maximum number of haplotypes retained per region.
    pub const DEFAULT_MAX_HAPLOTYPES: i32 = 15;

    /// Creates an empty alignment builder.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Converts a trace graph into candidate alignment chains.
    pub fn alignments_from_trace(
        &mut self,
        trace_node: Option<&TraceNode>,
        aligner_reverse: bool,
    ) -> Result<Vec<AlignNode>, KmerAlignmentBuilderError> {
        let Some(trace_node) = trace_node else {
            return Ok(Vec::new());
        };

        let mut starts = vec![AlignStart::new(trace_node, self.get_type_list(None))];
        let mut start_index = 0;

        while start_index < starts.len() {
            let mut current_trace = Some(starts[start_index].trace_node);
            let mut type_list = starts[start_index].type_list.clone();

            while let Some(trace) = current_trace {
                if let Some(branch) = trace.branch_node.as_deref() {
                    starts.push(AlignStart::new(
                        branch,
                        self.get_type_list(Some(&type_list)),
                    ));
                }

                type_list.add(trace.node_type)?;
                current_trace = trace.next_node.as_deref();
            }

            self.return_type_list(std::mem::replace(
                &mut starts[start_index].type_list,
                type_list,
            ));
            start_index += 1;
        }

        let reverse_alignment = !aligner_reverse;
        let mut alignments = Vec::with_capacity(starts.len());
        for start in starts {
            if let Some(alignment) = start.type_list.to_alignment(reverse_alignment) {
                alignments.push(*alignment);
            }
        }

        Ok(alignments)
    }

    fn get_type_list(&mut self, source: Option<&TypeList>) -> TypeList {
        if let Some(mut type_list) = self.type_list_cache.pop() {
            type_list.clone_from_option(source);
            type_list
        } else {
            TypeList::new_from(source)
        }
    }

    fn return_type_list(&mut self, mut type_list: TypeList) {
        type_list.clear();
        self.type_list_cache.push(type_list);
    }

    /// Returns the number of reusable type lists currently cached.
    #[must_use]
    pub fn cached_type_lists(&self) -> usize {
        self.type_list_cache.len()
    }
}

#[derive(Clone, Debug)]
struct AlignStart<'a> {
    trace_node: &'a TraceNode,
    type_list: TypeList,
}

impl<'a> AlignStart<'a> {
    fn new(trace_node: &'a TraceNode, type_list: TypeList) -> Self {
        Self {
            trace_node,
            type_list,
        }
    }
}

/// Incremental k-mer aligner for assembling haplotypes through an active region.
#[derive(Clone, Debug)]
pub struct KmerAligner {
    /// K-mer utility used for encoding consensus and reference k-mers.
    pub kmer_util: KmerUtil,
    /// Alignment scoring weights.
    pub alignment_weight: AlignmentWeight,
    /// Whether trace matrices are retained.
    pub trace: bool,
    active_region: Option<ActiveRegion>,
    ref_base_start: usize,
    ref_base_end: usize,
    ref_length: usize,
    reverse: bool,
    allow_end_deletion: bool,
    matrix_col_align: Vec<Option<TraceNode>>,
    matrix_col_gap_ref: Vec<Option<TraceNode>>,
    matrix_col_gap_con: Vec<Option<TraceNode>>,
    matrix_col_align_next: Vec<Option<TraceNode>>,
    matrix_col_gap_ref_next: Vec<Option<TraceNode>>,
    matrix_col_gap_con_next: Vec<Option<TraceNode>>,
    consensus: Vec<u8>,
    consensus_capacity: usize,
    max_alignment_score: f32,
    max_alignment_score_node: Option<Box<MaxAlignmentScoreNode>>,
    trace_matrix: Option<TraceMatrix>,
    saved_states: Vec<SavedAlignmentState>,
    max_state: i32,
    alignment_builder: KmerAlignmentBuilder,
    initialized: bool,
}

impl KmerAligner {
    /// Default maximum number of saved traversal states.
    pub const DEFAULT_MAX_STATE: i32 = 10;
    const CONS_SIZE_MULTIPLIER: f32 = 1.2;

    /// Creates an uninitialized aligner.
    pub fn new(
        kmer_util: KmerUtil,
        alignment_weight: AlignmentWeight,
        trace: bool,
    ) -> Result<Self, KmerAlignerError> {
        if kmer_util.k_size() < MIN_KMER_SIZE {
            return Err(KmerAlignerError::KmerSizeTooSmall {
                min: MIN_KMER_SIZE,
                actual: kmer_util.k_size(),
            });
        }

        Ok(Self {
            kmer_util,
            alignment_weight,
            trace,
            active_region: None,
            ref_base_start: 0,
            ref_base_end: 0,
            ref_length: 0,
            reverse: false,
            allow_end_deletion: false,
            matrix_col_align: Vec::new(),
            matrix_col_gap_ref: Vec::new(),
            matrix_col_gap_con: Vec::new(),
            matrix_col_align_next: Vec::new(),
            matrix_col_gap_ref_next: Vec::new(),
            matrix_col_gap_con_next: Vec::new(),
            consensus: Vec::new(),
            consensus_capacity: 0,
            max_alignment_score: 0.0,
            max_alignment_score_node: None,
            trace_matrix: None,
            saved_states: Vec::new(),
            max_state: Self::DEFAULT_MAX_STATE,
            alignment_builder: KmerAlignmentBuilder::new(),
            initialized: false,
        })
    }

    /// Initializes the aligner for an active region.
    pub fn init(&mut self, active_region: ActiveRegion) -> Result<(), KmerAlignerError> {
        self.reverse = active_region.left_end;
        self.allow_end_deletion = active_region.left_end || active_region.right_end;
        self.ref_base_start = active_region.start_index as usize;
        self.ref_base_end = active_region.end_index as usize;
        self.ref_length = self.ref_base_end - self.ref_base_start + 1;

        self.matrix_col_align = vec![None; self.ref_length];
        self.matrix_col_gap_ref = vec![None; self.ref_length];
        self.matrix_col_gap_con = vec![None; self.ref_length];
        self.matrix_col_align_next = vec![None; self.ref_length];
        self.matrix_col_gap_ref_next = vec![None; self.ref_length];
        self.matrix_col_gap_con_next = vec![None; self.ref_length];

        let new_capacity = ((self.ref_length as f32) * Self::CONS_SIZE_MULTIPLIER) as usize;
        self.consensus_capacity = new_capacity
            .max(self.kmer_util.k_size())
            .min(MAX_ARRAY_SIZE);
        self.consensus.clear();
        self.consensus.reserve(self.consensus_capacity);

        self.max_alignment_score = 0.0;
        self.max_alignment_score_node = None;
        self.saved_states.clear();
        self.active_region = Some(active_region);
        self.init_alignment()?;
        self.initialized = true;
        Ok(())
    }

    fn init_alignment(&mut self) -> Result<(), KmerAlignerError> {
        let active_region = self.active_region.as_ref().unwrap();
        let init_score = self
            .alignment_weight
            .initial_score(self.kmer_util.k_size() as i32)?;
        let mut last_node = None;

        self.trace_matrix = if self.trace {
            Some(TraceMatrix::new(self.ref_length as i32)?)
        } else {
            None
        };

        for count in 0..self.kmer_util.k_size() {
            last_node = Some(TraceNode::with_next(
                init_score,
                TraceNode::TYPE_MATCH,
                last_node.map(Box::new),
            ));
            if let Some(trace_matrix) = &mut self.trace_matrix {
                trace_matrix.next_col()?;
                trace_matrix.set(
                    count as i32,
                    TraceNode::TYPE_MATCH.into(),
                    TraceNode::TYPE_MATCH.into(),
                )?;
            }
        }

        self.matrix_col_align[self.kmer_util.k_size() - 1] = last_node.clone();

        if self.kmer_util.k_size() < self.ref_length
            && init_score + self.alignment_weight.new_gap > 0.0
        {
            let index = self.kmer_util.k_size();
            self.matrix_col_gap_con[index] = Some(TraceNode::with_next(
                init_score + self.alignment_weight.new_gap,
                TraceNode::TYPE_GAP_CON,
                last_node.map(Box::new),
            ));
            let mut previous = self.matrix_col_gap_con[index].clone();
            for index in self.kmer_util.k_size() + 1..self.ref_length {
                let Some(previous_node) = previous else {
                    break;
                };
                let score = previous_node.score + self.alignment_weight.gap_extend;
                if score <= 0.0 {
                    break;
                }
                self.matrix_col_gap_con[index] = Some(TraceNode::with_next(
                    score,
                    TraceNode::TYPE_GAP_CON,
                    Some(Box::new(previous_node.clone())),
                ));
                previous = self.matrix_col_gap_con[index].clone();
            }
        }

        let sequence = &active_region.ref_region.sequence;
        if self.reverse {
            for index in 0..self.kmer_util.k_size() {
                self.consensus.push(sequence[self.ref_base_end - index]);
            }
        } else {
            for index in 0..self.kmer_util.k_size() {
                self.consensus.push(sequence[self.ref_base_start + index]);
            }
        }
        Ok(())
    }

    /// Adds one consensus base and advances the alignment matrix.
    pub fn add_base(&mut self, base: Base) -> Result<bool, KmerAlignerError> {
        if !self.initialized {
            return Err(KmerAlignerError::NotInitialized);
        }
        if self.consensus.len() == MAX_ARRAY_SIZE {
            return Err(KmerAlignerError::ConsensusAtMaximum(MAX_ARRAY_SIZE));
        }
        if self.consensus.len() == self.consensus_capacity {
            self.expand_consensus()?;
        }

        self.consensus.push(base.as_byte());
        self.matrix_col_align_next.fill(None);
        self.matrix_col_gap_ref_next.fill(None);
        self.matrix_col_gap_con_next.fill(None);

        let active_region = self.active_region.as_ref().unwrap();
        let mut max_potential_score: f32 = 0.0;

        for index in 1..self.ref_length {
            let ref_index = if self.reverse {
                self.ref_base_end - index
            } else {
                self.ref_base_start + index
            };
            let (add_align_score, align_type) =
                if active_region.ref_region.sequence[ref_index] == base.as_byte() {
                    (self.alignment_weight.match_score, TraceNode::TYPE_MATCH)
                } else {
                    (self.alignment_weight.mismatch, TraceNode::TYPE_MISMATCH)
                };

            let candidates = [
                (
                    score_from(&self.matrix_col_align[index - 1], add_align_score),
                    &self.matrix_col_align[index - 1],
                ),
                (
                    score_from(&self.matrix_col_gap_ref[index - 1], add_align_score),
                    &self.matrix_col_gap_ref[index - 1],
                ),
                (
                    score_from(&self.matrix_col_gap_con[index - 1], add_align_score),
                    &self.matrix_col_gap_con[index - 1],
                ),
            ];
            let max_score = candidates
                .iter()
                .map(|(score, _)| *score)
                .fold(0.0, f32::max);

            if max_score > 0.0 {
                self.matrix_col_align_next[index] =
                    trace_branch(max_score, align_type, &candidates);
                max_potential_score = max_potential_score.max(
                    max_score
                        + (self.ref_length - index - 1) as f32 * self.alignment_weight.match_score,
                );
            }
        }

        if let Some(node) = &self.matrix_col_align_next[self.ref_length - 1] {
            self.record_max_node(node.clone());
        }

        for index in 0..self.ref_length {
            let candidates = [
                (
                    score_from(&self.matrix_col_align[index], self.alignment_weight.new_gap),
                    &self.matrix_col_align[index],
                ),
                (
                    score_from(
                        &self.matrix_col_gap_ref[index],
                        self.alignment_weight.gap_extend,
                    ),
                    &self.matrix_col_gap_ref[index],
                ),
                (
                    score_from(
                        &self.matrix_col_gap_con[index],
                        self.alignment_weight.new_gap,
                    ),
                    &self.matrix_col_gap_con[index],
                ),
            ];
            let max_score = candidates
                .iter()
                .map(|(score, _)| *score)
                .fold(0.0, f32::max);
            if max_score > 0.0 {
                self.matrix_col_gap_ref_next[index] =
                    trace_branch(max_score, TraceNode::TYPE_GAP_REF, &candidates);
            }
        }

        for index in 1..self.ref_length {
            let candidates = [
                (
                    score_from(
                        &self.matrix_col_align_next[index - 1],
                        self.alignment_weight.new_gap,
                    ),
                    &self.matrix_col_align_next[index - 1],
                ),
                (
                    score_from(
                        &self.matrix_col_gap_ref_next[index - 1],
                        self.alignment_weight.new_gap,
                    ),
                    &self.matrix_col_gap_ref_next[index - 1],
                ),
                (
                    score_from(
                        &self.matrix_col_gap_con_next[index - 1],
                        self.alignment_weight.gap_extend,
                    ),
                    &self.matrix_col_gap_con_next[index - 1],
                ),
            ];
            let max_score = candidates
                .iter()
                .map(|(score, _)| *score)
                .fold(0.0, f32::max);
            if max_score > 0.0 {
                self.matrix_col_gap_con_next[index] =
                    trace_branch(max_score, TraceNode::TYPE_GAP_CON, &candidates);
                if self.allow_end_deletion {
                    max_potential_score = max_potential_score.max(
                        max_score
                            + (self.ref_length - index - 1) as f32
                                * self.alignment_weight.match_score,
                    );
                }
            }
        }

        if self.allow_end_deletion
            && let Some(node) = &self.matrix_col_gap_con_next[self.ref_length - 1]
        {
            self.record_max_node(node.clone());
        }

        if let Some(trace_matrix) = &mut self.trace_matrix {
            trace_matrix.next_col()?;
            set_trace_matrix(trace_matrix, &self.matrix_col_align_next)?;
            set_trace_matrix(trace_matrix, &self.matrix_col_gap_ref_next)?;
            set_trace_matrix(trace_matrix, &self.matrix_col_gap_con_next)?;
        }

        std::mem::swap(&mut self.matrix_col_align, &mut self.matrix_col_align_next);
        std::mem::swap(
            &mut self.matrix_col_gap_ref,
            &mut self.matrix_col_gap_ref_next,
        );
        std::mem::swap(
            &mut self.matrix_col_gap_con,
            &mut self.matrix_col_gap_con_next,
        );

        Ok(max_potential_score >= self.max_alignment_score && max_potential_score > 0.0)
    }

    fn record_max_node(&mut self, node: TraceNode) {
        let max_score = node.score;
        if max_score >= self.max_alignment_score && max_score > 0.0 {
            let next = if max_score > self.max_alignment_score {
                None
            } else {
                self.max_alignment_score_node.take()
            };
            self.max_alignment_score_node = Some(Box::new(MaxAlignmentScoreNode {
                trace_node: Box::new(node),
                n_consensus_bases: self.consensus.len() as i32,
                next,
                haplotype_built: false,
            }));
            self.max_alignment_score = max_score;
        }
    }

    fn expand_consensus(&mut self) -> Result<(), KmerAlignerError> {
        if self.consensus_capacity == MAX_ARRAY_SIZE {
            return Err(KmerAlignerError::ConsensusAtMaximum(MAX_ARRAY_SIZE));
        }
        let new_capacity = ((self.consensus_capacity as f32) * ARRAY_EXPAND_FACTOR) as usize + 1;
        self.consensus_capacity = new_capacity.min(MAX_ARRAY_SIZE);
        self.consensus.reserve(
            self.consensus_capacity
                .saturating_sub(self.consensus.capacity()),
        );
        Ok(())
    }

    /// Saves the current aligner state for later backtracking.
    pub fn save_state(
        &mut self,
        kmer: Option<KmerKey>,
        next_base: Option<Base>,
        min_depth: i32,
        kmer_hash: Option<KmerHashSet>,
        repeat_count: i32,
    ) -> Result<(), KmerAlignerError> {
        let kmer = kmer.ok_or(KmerAlignerError::NullKmer)?;
        let next_base = next_base.ok_or(KmerAlignerError::NullNextBase)?;
        if min_depth < 1 {
            return Err(KmerAlignerError::InvalidMinDepth(min_depth));
        }
        let kmer_hash = kmer_hash.ok_or(KmerAlignerError::NullKmerHash)?;
        if repeat_count < 0 {
            return Err(KmerAlignerError::NegativeRepeatCount(repeat_count));
        }

        if self.saved_states.len() == self.max_state as usize && !self.remove_min_state(min_depth) {
            return Ok(());
        }

        self.saved_states.push(SavedAlignmentState {
            kmer,
            next_base,
            consensus_size: self.consensus.len(),
            matrix_col_align: self.matrix_col_align.clone(),
            matrix_col_gap_ref: self.matrix_col_gap_ref.clone(),
            matrix_col_gap_con: self.matrix_col_gap_con.clone(),
            max_alignment_score: self.max_alignment_score,
            max_alignment_score_node: self.max_alignment_score_node.clone(),
            min_depth,
            kmer_hash,
            repeat_count,
        });
        Ok(())
    }

    /// Restores the most recently saved state, if one exists.
    pub fn restore_state(&mut self) -> Result<Option<state::RestoredState>, KmerAlignerError> {
        let Some(saved) = self.saved_states.pop() else {
            return Ok(None);
        };

        self.consensus.truncate(saved.consensus_size);
        self.matrix_col_align = saved.matrix_col_align.clone();
        self.matrix_col_gap_ref = saved.matrix_col_gap_ref.clone();
        self.matrix_col_gap_con = saved.matrix_col_gap_con.clone();
        self.max_alignment_score = saved.max_alignment_score;
        self.max_alignment_score_node = saved.max_alignment_score_node.clone();
        self.add_base(saved.next_base)?;

        Ok(Some(state::RestoredState {
            kmer: saved.kmer.words().to_vec(),
            consensus_size: saved.consensus_size as i32,
            min_depth: saved.min_depth,
            kmer_hash: saved.kmer_hash,
            repeat_count: saved.repeat_count,
        }))
    }

    fn remove_min_state(&mut self, min_depth_limit: i32) -> bool {
        let Some((index, _)) = self
            .saved_states
            .iter()
            .enumerate()
            .filter(|(_, state)| state.min_depth < min_depth_limit)
            .min_by_key(|(_, state)| state.min_depth)
        else {
            return false;
        };
        self.saved_states.remove(index);
        true
    }

    /// Returns true when backtracking states are cached.
    #[must_use]
    pub fn has_cached_states(&self) -> bool {
        !self.saved_states.is_empty()
    }

    /// Builds haplotypes from maximum-scoring alignment traces.
    pub fn get_haplotypes(
        &mut self,
        counter: &dyn CountMap,
        count_reverse_kmers: bool,
    ) -> Result<Vec<Haplotype>, KmerAlignerError> {
        self.trim_haplotypes();
        let active_region = self
            .active_region
            .clone()
            .ok_or(KmerAlignerError::NotInitialized)?;
        let mut haplotypes = Vec::new();
        let mut current = self.max_alignment_score_node.as_mut();

        while let Some(score_node) = current {
            if !score_node.haplotype_built {
                let consensus_size = score_node.n_consensus_bases as usize;
                let mut sequence = self.consensus[..consensus_size].to_vec();
                if self.reverse {
                    sequence.reverse();
                }
                let stats =
                    haplotype_stats(&self.kmer_util, &sequence, counter, count_reverse_kmers)?;
                let alignments = self
                    .alignment_builder
                    .alignments_from_trace(Some(&score_node.trace_node), self.reverse)?;
                haplotypes.push(Haplotype::new(
                    sequence,
                    active_region.clone(),
                    alignments,
                    score_node.trace_node.score,
                    self.trace_matrix.clone(),
                    stats,
                )?);
                score_node.haplotype_built = true;
            }
            current = score_node.next.as_mut();
        }

        Ok(haplotypes)
    }

    fn trim_haplotypes(&mut self) {
        let Some(active_region) = &self.active_region else {
            return;
        };
        if active_region.left_end || active_region.right_end {
            return;
        }

        let k_size = self.kmer_util.k_size();
        let ref_seq = &active_region.ref_region.sequence;
        let mut kept = Vec::new();
        let mut current = self.max_alignment_score_node.take();
        while let Some(mut node) = current {
            current = node.next.take();
            let con_index = node.n_consensus_bases as usize - k_size;
            let remove = if self.reverse {
                (0..k_size).any(|offset| {
                    ref_seq[active_region.start_index as usize + k_size - 1 - offset]
                        != self.consensus[con_index + offset]
                })
            } else {
                (0..k_size).any(|offset| {
                    ref_seq[active_region.end_index as usize - k_size + 1 + offset]
                        != self.consensus[con_index + offset]
                })
            };
            if !remove {
                kept.push(node);
            }
        }

        self.max_alignment_score_node = kept.into_iter().rev().fold(None, |next, mut node| {
            node.next = next;
            Some(node)
        });
    }

    /// Returns true when this aligner traverses the active region in reverse.
    #[must_use]
    pub fn is_reverse(&self) -> bool {
        self.reverse
    }

    /// Returns true when terminal deletions are allowed for open-ended regions.
    #[must_use]
    pub fn is_allow_end_deletion(&self) -> bool {
        self.allow_end_deletion
    }

    /// Sets the maximum number of saved backtracking states.
    pub fn set_max_state(&mut self, max_state: i32) -> Result<(), KmerAlignerError> {
        if max_state < 1 {
            return Err(KmerAlignerError::InvalidMaxState(max_state));
        }
        self.max_state = max_state;
        while self.saved_states.len() > max_state as usize {
            self.saved_states.remove(0);
        }
        Ok(())
    }

    /// Returns the maximum number of saved backtracking states.
    #[must_use]
    pub fn max_state(&self) -> i32 {
        self.max_state
    }

    /// Returns the current consensus sequence.
    #[must_use]
    pub fn consensus(&self) -> &[u8] {
        &self.consensus
    }

    /// Returns the best alignment score seen so far.
    #[must_use]
    pub fn max_alignment_score(&self) -> f32 {
        self.max_alignment_score
    }
}

#[derive(Clone, Debug)]
struct SavedAlignmentState {
    kmer: KmerKey,
    next_base: Base,
    consensus_size: usize,
    matrix_col_align: Vec<Option<TraceNode>>,
    matrix_col_gap_ref: Vec<Option<TraceNode>>,
    matrix_col_gap_con: Vec<Option<TraceNode>>,
    max_alignment_score: f32,
    max_alignment_score_node: Option<Box<MaxAlignmentScoreNode>>,
    min_depth: i32,
    kmer_hash: KmerHashSet,
    repeat_count: i32,
}

fn score_from(node: &Option<TraceNode>, add: f32) -> f32 {
    node.as_ref().map_or(0.0, |node| node.score + add)
}

fn trace_branch(
    max_score: f32,
    node_type: i8,
    candidates: &[(f32, &Option<TraceNode>); 3],
) -> Option<TraceNode> {
    let mut last = None;
    for (score, previous) in candidates {
        if *score == max_score
            && let Some(previous) = previous
        {
            last = Some(
                TraceNode::new(
                    max_score,
                    node_type,
                    Some(Box::new(previous.clone())),
                    last.map(Box::new),
                )
                .expect("KmerAligner only creates valid trace nodes"),
            );
        }
    }
    last
}

fn set_trace_matrix(
    trace_matrix: &mut TraceMatrix,
    column: &[Option<TraceNode>],
) -> Result<(), TraceMatrixError> {
    for (index, node) in column.iter().enumerate() {
        let mut node = node.as_ref();
        while let Some(current) = node {
            if let Some(next) = current.next_node.as_deref() {
                trace_matrix.set(
                    index as i32,
                    next.node_type.into(),
                    current.node_type.into(),
                )?;
            }
            node = current.branch_node.as_deref();
        }
    }
    Ok(())
}

fn haplotype_stats(
    kmer_util: &KmerUtil,
    sequence: &[u8],
    counter: &dyn CountMap,
    count_reverse_kmers: bool,
) -> Result<RegionStats, ActiveRegionError> {
    if sequence.len() < kmer_util.k_size() {
        return RegionStats::from_counts(&[0], 0, 1).map_err(ActiveRegionError::Stats);
    }

    let mut counts = Vec::with_capacity(sequence.len() - kmer_util.k_size() + 1);
    for window in sequence.windows(kmer_util.k_size()) {
        let kmer = kmer_util
            .encode(window)
            .expect("haplotype sequence only contains encoded DNA bases");
        let mut count = counter.get(&kmer) as i32;
        if count_reverse_kmers {
            count += counter.get(&kmer_util.reverse_complement(&kmer)) as i32;
        }
        counts.push(count);
    }
    RegionStats::from_counts(&counts, 0, counts.len() as i32).map_err(ActiveRegionError::Stats)
}

/// Saved aligner state types used during graph backtracking.
pub mod state {
    use kanalyze::Base;
    use kanalyze::util::KmerHashSet;
    use thiserror::Error;

    use super::{MaxAlignmentScoreNode, TraceNode};

    /// Errors returned by [`TraceNodeContainer`].
    #[derive(Clone, Debug, Error, Eq, PartialEq)]
    pub enum TraceNodeContainerError {
        /// Container index cannot be negative.
        #[error("TraceNodeContainer(): Negative index: {0}")]
        NegativeIndex(i32),
    }

    /// Linked container for trace nodes at one matrix column.
    #[derive(Clone, Debug, PartialEq)]
    pub struct TraceNodeContainer {
        /// Matrix row or column index represented by this container.
        pub index: i32,
        /// Trace node stored at this index.
        pub node: Box<TraceNode>,
        /// Next container in the linked list.
        pub next: Option<Box<TraceNodeContainer>>,
    }

    impl TraceNodeContainer {
        /// Creates a trace-node container.
        pub fn new(
            index: i32,
            node: TraceNode,
            next: Option<Box<Self>>,
        ) -> Result<Self, TraceNodeContainerError> {
            if index < 0 {
                return Err(TraceNodeContainerError::NegativeIndex(index));
            }

            Ok(Self {
                index,
                node: Box::new(node),
                next,
            })
        }
    }

    /// Public representation of a restored aligner state.
    #[derive(Clone, Debug, Eq, PartialEq)]
    pub struct RestoredState {
        /// Encoded k-mer words restored with the state.
        pub kmer: Vec<u32>,
        /// Consensus size before replaying the next base.
        pub consensus_size: i32,
        /// Minimum depth associated with the restored path.
        pub min_depth: i32,
        /// K-mer hash set associated with the restored path.
        pub kmer_hash: KmerHashSet,
        /// Repeat count associated with the restored path.
        pub repeat_count: i32,
    }

    /// Full saved aligner state node used by compatibility APIs.
    #[derive(Clone, Debug, PartialEq)]
    pub struct StateStackNode {
        /// Encoded k-mer words.
        pub kmer: Vec<u32>,
        /// Next base to replay after restoration.
        pub next_base: Base,
        /// Consensus size at the saved point.
        pub consensus_size: i32,
        /// K-mer hash set at the saved point.
        pub kmer_hash: KmerHashSet,
        /// Repeat count at the saved point.
        pub repeat_count: i32,
        /// Saved alignment-column trace nodes.
        pub align_container: Option<Box<TraceNodeContainer>>,
        /// Saved reference-gap-column trace nodes.
        pub gap_ref_container: Option<Box<TraceNodeContainer>>,
        /// Saved consensus-gap-column trace nodes.
        pub gap_con_container: Option<Box<TraceNodeContainer>>,
        /// Maximum alignment score at the saved point.
        pub max_alignment_score: f32,
        /// Maximum-score trace nodes at the saved point.
        pub max_alignment_score_node: Option<Box<MaxAlignmentScoreNode>>,
        /// Minimum depth associated with the saved path.
        pub min_depth: i32,
        /// Next node deeper in the saved-state stack.
        pub next_node_down: Option<Box<StateStackNode>>,
        /// Next node toward the top of the saved-state stack.
        pub next_node_up: Option<Box<StateStackNode>>,
    }

    /// Constructor fields for [`StateStackNode`].
    #[derive(Clone, Debug, PartialEq)]
    pub struct StateStackNodeFields {
        /// Encoded k-mer words.
        pub kmer: Vec<u32>,
        /// Next base to replay after restoration.
        pub next_base: Base,
        /// Consensus size at the saved point.
        pub consensus_size: i32,
        /// Saved alignment-column trace nodes.
        pub align_container: Option<Box<TraceNodeContainer>>,
        /// Saved reference-gap-column trace nodes.
        pub gap_ref_container: Option<Box<TraceNodeContainer>>,
        /// Saved consensus-gap-column trace nodes.
        pub gap_con_container: Option<Box<TraceNodeContainer>>,
        /// Maximum alignment score at the saved point.
        pub max_alignment_score: f32,
        /// Maximum-score trace nodes at the saved point.
        pub max_alignment_score_node: Option<Box<MaxAlignmentScoreNode>>,
        /// Minimum depth associated with the saved path.
        pub min_depth: i32,
        /// Next node deeper in the saved-state stack.
        pub next_node_down: Option<Box<StateStackNode>>,
        /// K-mer hash set at the saved point.
        pub kmer_hash: KmerHashSet,
        /// Repeat count at the saved point.
        pub repeat_count: i32,
    }

    impl StateStackNode {
        /// Creates a saved-state stack node from grouped fields.
        #[must_use]
        pub fn new(fields: StateStackNodeFields) -> Self {
            Self {
                kmer: fields.kmer,
                next_base: fields.next_base,
                consensus_size: fields.consensus_size,
                kmer_hash: fields.kmer_hash,
                repeat_count: fields.repeat_count,
                align_container: fields.align_container,
                gap_ref_container: fields.gap_ref_container,
                gap_con_container: fields.gap_con_container,
                max_alignment_score: fields.max_alignment_score,
                max_alignment_score_node: fields.max_alignment_score_node,
                min_depth: fields.min_depth,
                next_node_down: fields.next_node_down,
                next_node_up: None,
            }
        }

        /// Returns the subset of this node exposed after state restoration.
        #[must_use]
        pub fn restored_state(&self) -> RestoredState {
            RestoredState {
                kmer: self.kmer.clone(),
                consensus_size: self.consensus_size,
                min_depth: self.min_depth,
                kmer_hash: self.kmer_hash.clone(),
                repeat_count: self.repeat_count,
            }
        }
    }
}

/// Alignment scoring weights for k-mer dynamic programming.
#[derive(Clone, Copy, Debug)]
pub struct AlignmentWeight {
    /// Score added for a matching base.
    pub match_score: f32,
    /// Score added for a mismatching base.
    pub mismatch: f32,
    /// Score added when opening a gap.
    pub gap_open: f32,
    /// Score added when extending a gap.
    pub gap_extend: f32,
    /// Initial score, or zero to derive it from the k-mer size.
    pub init_score: f32,
    /// Precomputed gap-open plus gap-extend score.
    pub new_gap: f32,
}

/// Errors returned while parsing or using alignment weights.
#[derive(Clone, Debug, Error, Eq, PartialEq)]
pub enum AlignmentWeightError {
    /// Required weight was zero or too close to zero.
    #[error("weight is zero or too close to zero")]
    ZeroWeight,
    /// Weight string delimiters were mismatched.
    #[error("mismatched delimiters in alignment weights: {0}")]
    MismatchedDelimiters(String),
    /// Weight string contained a closing delimiter with no opener.
    #[error("closing delimiter with no opening delimiter in alignment weights: {0}")]
    UnopenedDelimiter(String),
    /// Weight string contained too many comma-separated values.
    #[error("weight vector has more than 5 comma-separated values: {0}")]
    TooManyValues(usize),
    /// Weight token could not be parsed as a number.
    #[error("cannot convert weight to a numeric value: {0}")]
    InvalidNumber(String),
    /// Formatting precision was outside the supported range.
    #[error("precision out of range [0, 10]: {0}")]
    InvalidPrecision(i32),
    /// K-mer size must be positive.
    #[error("k-mer size is less than 1: {0}")]
    InvalidKSize(i32),
}

impl AlignmentWeight {
    /// Default match score.
    pub const DEFAULT_MATCH: f32 = 10.0;
    /// Default mismatch score.
    pub const DEFAULT_MISMATCH: f32 = -10.0;
    /// Default gap-open score.
    pub const DEFAULT_GAP_OPEN: f32 = -40.0;
    /// Default gap-extend score.
    pub const DEFAULT_GAP_EXTEND: f32 = -4.0;
    /// Default initial score.
    pub const DEFAULT_INIT_SCORE: f32 = 0.0;

    /// Returns the default Kestrel alignment weights.
    #[must_use]
    pub fn defaults() -> Self {
        Self::raw(
            Self::DEFAULT_MATCH,
            Self::DEFAULT_MISMATCH,
            Self::DEFAULT_GAP_OPEN,
            Self::DEFAULT_GAP_EXTEND,
            Self::DEFAULT_INIT_SCORE,
        )
    }

    /// Creates normalized alignment weights.
    pub fn new(
        match_score: f32,
        mismatch: f32,
        gap_open: f32,
        gap_extend: f32,
        init_score: f32,
    ) -> Result<Self, AlignmentWeightError> {
        Ok(Self::raw(
            normalize_match(match_score)?,
            normalize_mismatch(mismatch)?,
            normalize_gap_open(gap_open),
            normalize_gap_extend(gap_extend)?,
            normalize_init_score(init_score),
        ))
    }

    /// Parses optional comma-separated alignment weights.
    pub fn parse(value: Option<&str>) -> Result<Self, AlignmentWeightError> {
        let Some(value) = value.map(str::trim).filter(|value| !value.is_empty()) else {
            return Ok(Self::defaults());
        };

        let value = strip_delimiters(value)?;
        let tokens = value.split(',').map(str::trim).collect::<Vec<_>>();
        if tokens.len() > 5 {
            return Err(AlignmentWeightError::TooManyValues(tokens.len()));
        }

        let mut weights = Self::defaults();
        if let Some(token) = tokens.first().filter(|token| !token.is_empty()) {
            weights.match_score = normalize_match(parse_weight(token)?)?;
        }
        if let Some(token) = tokens.get(1).filter(|token| !token.is_empty()) {
            weights.mismatch = normalize_mismatch(parse_weight(token)?)?;
        }
        if let Some(token) = tokens.get(2).filter(|token| !token.is_empty()) {
            weights.gap_open = normalize_gap_open(parse_weight(token)?);
        }
        if let Some(token) = tokens.get(3).filter(|token| !token.is_empty()) {
            weights.gap_extend = normalize_gap_extend(parse_weight(token)?)?;
        }
        if let Some(token) = tokens.get(4).filter(|token| !token.is_empty()) {
            weights.init_score = normalize_init_score(parse_weight(token)?);
        }
        weights.new_gap = weights.gap_open + weights.gap_extend;

        Ok(weights)
    }

    #[must_use]
    fn raw(
        match_score: f32,
        mismatch: f32,
        gap_open: f32,
        gap_extend: f32,
        init_score: f32,
    ) -> Self {
        Self {
            match_score,
            mismatch,
            gap_open,
            gap_extend,
            init_score,
            new_gap: gap_open + gap_extend,
        }
    }

    /// Returns a copy with a different match score.
    pub fn with_match(self, match_score: f32) -> Result<Self, AlignmentWeightError> {
        Self::new(
            match_score,
            self.mismatch,
            self.gap_open,
            self.gap_extend,
            self.init_score,
        )
    }

    /// Returns a copy with a different mismatch score.
    pub fn with_mismatch(self, mismatch: f32) -> Result<Self, AlignmentWeightError> {
        Self::new(
            self.match_score,
            mismatch,
            self.gap_open,
            self.gap_extend,
            self.init_score,
        )
    }

    /// Returns a copy with a different gap-open score.
    pub fn with_gap_open(self, gap_open: f32) -> Result<Self, AlignmentWeightError> {
        Self::new(
            self.match_score,
            self.mismatch,
            gap_open,
            self.gap_extend,
            self.init_score,
        )
    }

    /// Returns a copy with a different gap-extend score.
    pub fn with_gap_extend(self, gap_extend: f32) -> Result<Self, AlignmentWeightError> {
        Self::new(
            self.match_score,
            self.mismatch,
            self.gap_open,
            gap_extend,
            self.init_score,
        )
    }

    /// Returns a copy with a different initial score.
    pub fn with_initial_score(self, init_score: f32) -> Result<Self, AlignmentWeightError> {
        Self::new(
            self.match_score,
            self.mismatch,
            self.gap_open,
            self.gap_extend,
            init_score,
        )
    }

    /// Returns the alignment initial score for a k-mer size.
    pub fn initial_score(&self, k_size: i32) -> Result<f32, AlignmentWeightError> {
        if k_size < 1 {
            return Err(AlignmentWeightError::InvalidKSize(k_size));
        }

        if is_zero(self.init_score) {
            Ok(self.match_score * k_size as f32)
        } else {
            Ok(self.init_score)
        }
    }

    /// Returns the largest gap size that can still beat zero score after initialization.
    pub fn max_exclusive_gap_size(&self, k_size: i32) -> Result<i32, AlignmentWeightError> {
        let init_score = self.initial_score(k_size)? as i32 as f32;
        if init_score > self.gap_open {
            Ok(((init_score + self.gap_open) / -self.gap_extend) as i32)
        } else {
            Ok(0)
        }
    }

    /// Returns true when these weights match the defaults.
    #[must_use]
    pub fn is_default(&self) -> bool {
        *self == Self::defaults()
    }

    /// Formats this weight vector with a fixed decimal precision.
    pub fn format_with_precision(&self, precision: i32) -> Result<String, AlignmentWeightError> {
        if !(0..=10).contains(&precision) {
            return Err(AlignmentWeightError::InvalidPrecision(precision));
        }

        let precision = precision as usize;
        Ok(format!(
            "({:.precision$}, {:.precision$}, {:.precision$}, {:.precision$}, {:.precision$})",
            self.match_score, self.mismatch, self.gap_open, self.gap_extend, self.init_score
        ))
    }
}

impl Default for AlignmentWeight {
    fn default() -> Self {
        Self::defaults()
    }
}

impl PartialEq for AlignmentWeight {
    fn eq(&self, other: &Self) -> bool {
        is_zero(self.match_score - other.match_score)
            && is_zero(self.mismatch - other.mismatch)
            && is_zero(self.gap_open - other.gap_open)
            && is_zero(self.gap_extend - other.gap_extend)
            && is_zero(self.init_score - other.init_score)
    }
}

impl Eq for AlignmentWeight {}

impl std::hash::Hash for AlignmentWeight {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        state.write_u32(self.match_score.to_bits());
        state.write_u32(self.mismatch.to_bits());
        state.write_u32(self.gap_open.to_bits());
        state.write_u32(self.gap_extend.to_bits());
        state.write_u32(self.init_score.to_bits());
    }
}

impl fmt::Display for AlignmentWeight {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(
            &self
                .format_with_precision(4)
                .expect("default display precision is valid"),
        )
    }
}

fn normalize_match(value: f32) -> Result<f32, AlignmentWeightError> {
    if is_zero(value) {
        return Err(AlignmentWeightError::ZeroWeight);
    }
    Ok(value.abs())
}

fn normalize_mismatch(value: f32) -> Result<f32, AlignmentWeightError> {
    if is_zero(value) {
        return Err(AlignmentWeightError::ZeroWeight);
    }
    Ok(-value.abs())
}

fn normalize_gap_open(value: f32) -> f32 {
    -value.abs()
}

fn normalize_gap_extend(value: f32) -> Result<f32, AlignmentWeightError> {
    if is_zero(value) {
        return Err(AlignmentWeightError::ZeroWeight);
    }
    Ok(-value.abs())
}

fn normalize_init_score(value: f32) -> f32 {
    if is_zero(value) { 0.0 } else { value.abs() }
}

fn strip_delimiters(value: &str) -> Result<&str, AlignmentWeightError> {
    let Some(first) = value.chars().next() else {
        return Ok(value);
    };
    let Some(last) = value.chars().last() else {
        return Ok(value);
    };

    let expected = match first {
        '(' => Some(')'),
        '<' => Some('>'),
        '[' => Some(']'),
        '{' => Some('}'),
        _ => None,
    };

    if let Some(expected) = expected {
        if last == expected {
            return Ok(&value[1..value.len() - 1]);
        }
        return Err(AlignmentWeightError::MismatchedDelimiters(value.to_owned()));
    }

    if matches!(last, ')' | '>' | ']' | '}') {
        return Err(AlignmentWeightError::UnopenedDelimiter(value.to_owned()));
    }

    Ok(value)
}

fn parse_weight(value: &str) -> Result<f32, AlignmentWeightError> {
    let value = value.trim();
    if let Ok(parsed) = value.parse::<f32>() {
        return Ok(parsed);
    }
    parse_java_integer(value)
        .map(|value| value as f32)
        .ok_or_else(|| AlignmentWeightError::InvalidNumber(value.to_owned()))
}

fn parse_java_integer(value: &str) -> Option<i32> {
    let (sign, rest) = if let Some(rest) = value.strip_prefix('-') {
        (-1_i64, rest)
    } else if let Some(rest) = value.strip_prefix('+') {
        (1_i64, rest)
    } else {
        (1_i64, value)
    };

    let (radix, digits) =
        if let Some(rest) = rest.strip_prefix("0x").or_else(|| rest.strip_prefix("0X")) {
            (16, rest)
        } else if let Some(rest) = rest.strip_prefix('#') {
            (16, rest)
        } else if rest.len() > 1 && rest.starts_with('0') {
            (8, rest)
        } else {
            (10, rest)
        };

    i64::from_str_radix(digits, radix)
        .ok()
        .and_then(|value| i32::try_from(sign * value).ok())
}

#[cfg(test)]
mod tests {
    use std::collections::hash_map::DefaultHasher;
    use std::hash::{Hash, Hasher};

    use kanalyze::Base;
    use kanalyze::util::{KmerHashSet, KmerUtil};

    use crate::activeregion::ActiveRegion;
    use crate::counter::{CountMap, CountMapError};
    use crate::io::InputSample;
    use crate::refreader::{ReferenceRegion, ReferenceSequence};
    use crate::util::digest::Digest;

    use super::*;

    const EPS: f32 = 0.0001;

    #[test]
    fn trace_node_constants_and_zero_node_match_java() {
        assert_eq!(TraceNode::TYPE_NONE, 0);
        assert_eq!(TraceNode::TYPE_MATCH, 1);
        assert_eq!(TraceNode::TYPE_MISMATCH, 2);
        assert_eq!(TraceNode::TYPE_GAP_REF, 3);
        assert_eq!(TraceNode::TYPE_GAP_CON, 4);
        assert_eq!(TraceNode::ZERO_NODE.node_type, TraceNode::TYPE_NONE);
        assert_eq!(TraceNode::ZERO_NODE.score, 0.0);
    }

    #[test]
    fn trace_node_constructors_validate_like_java() {
        let prev = TraceNode::with_next(5.0, TraceNode::TYPE_MATCH, None);
        let branch = TraceNode::with_next(3.0, TraceNode::TYPE_MISMATCH, None);
        let node = TraceNode::new(
            10.0,
            TraceNode::TYPE_MATCH,
            Some(Box::new(prev.clone())),
            Some(Box::new(branch.clone())),
        )
        .unwrap();

        assert_eq!(node.score, 10.0);
        assert_eq!(node.node_type, TraceNode::TYPE_MATCH);
        assert_eq!(node.next_node.as_deref(), Some(&prev));
        assert_eq!(node.branch_node.as_deref(), Some(&branch));

        let no_branch = TraceNode::with_next(5.0, TraceNode::TYPE_MISMATCH, None);
        assert_eq!(no_branch.branch_node, None);

        assert_eq!(
            TraceNode::new(-1.0, TraceNode::TYPE_MATCH, None, None),
            Err(TraceNodeError::NegativeScore)
        );
        assert!(matches!(
            TraceNode::new(1.0, 5, None, None),
            Err(TraceNodeError::TypeOutOfRange { .. })
        ));
        assert!(matches!(
            TraceNode::new(1.0, -1, None, None),
            Err(TraceNodeError::TypeOutOfRange { .. })
        ));
    }

    #[test]
    fn trace_node_cigar_array_and_display_match_java() {
        let mut cigar = TraceNode::cigar_array();
        assert_eq!(cigar[TraceNode::TYPE_NONE as usize], '*');
        assert_eq!(cigar[TraceNode::TYPE_MATCH as usize], '=');
        assert_eq!(cigar[TraceNode::TYPE_MISMATCH as usize], 'X');
        assert_eq!(cigar[TraceNode::TYPE_GAP_REF as usize], 'I');
        assert_eq!(cigar[TraceNode::TYPE_GAP_CON as usize], 'D');

        cigar[TraceNode::TYPE_MATCH as usize] = '?';
        assert_eq!(
            TraceNode::cigar_array()[TraceNode::TYPE_MATCH as usize],
            '='
        );

        let display = TraceNode::with_next(5.5, TraceNode::TYPE_GAP_REF, None).to_string();
        assert!(display.contains("5.5"));
        assert!(display.contains("GAP_REFERENCE"));
    }

    #[test]
    fn align_node_constants_cigar_and_chain_match_java() {
        assert_eq!(AlignNode::MATCH, 1);
        assert_eq!(AlignNode::MISMATCH, 2);
        assert_eq!(AlignNode::INS, 3);
        assert_eq!(AlignNode::DEL, 4);

        assert_eq!(
            AlignNode::new(AlignNode::MATCH, 10, None)
                .unwrap()
                .cigar_char,
            '='
        );
        assert_eq!(
            AlignNode::new(AlignNode::MISMATCH, 1, None)
                .unwrap()
                .cigar_char,
            'X'
        );
        assert_eq!(
            AlignNode::new(AlignNode::INS, 2, None).unwrap().cigar_char,
            'I'
        );
        assert_eq!(
            AlignNode::new(AlignNode::DEL, 3, None).unwrap().cigar_char,
            'D'
        );

        let last = AlignNode::new(AlignNode::MATCH, 5, None).unwrap();
        let mid = AlignNode::new(AlignNode::MISMATCH, 1, Some(Box::new(last))).unwrap();
        let first = AlignNode::new(AlignNode::MATCH, 10, Some(Box::new(mid))).unwrap();
        assert_eq!(first.cigar_string(), "10=1X5=");

        assert_eq!(
            AlignNode::new(AlignNode::DEL, 7, None)
                .unwrap()
                .cigar_string(),
            "7D"
        );
        assert_eq!(
            AlignNode::new(5, 1, None),
            Err(AlignNodeError::TypeOutOfRange(5))
        );
        assert_eq!(
            AlignNode::new(0, 1, None),
            Err(AlignNodeError::TypeOutOfRange(0))
        );
    }

    #[test]
    fn align_node_compare_to_preserves_java_ordering() {
        let a = AlignNode::new(AlignNode::MATCH, 10, None).unwrap();
        let b = AlignNode::new(AlignNode::MATCH, 10, None).unwrap();
        assert_eq!(a.compare_to(&b), -1);
        assert_eq!(b.compare_to(&a), -1);

        let shorter = AlignNode::new(AlignNode::MATCH, 10, None).unwrap();
        let longer_tail = AlignNode::new(AlignNode::MATCH, 5, None).unwrap();
        let longer = AlignNode::new(AlignNode::MATCH, 10, Some(Box::new(longer_tail))).unwrap();
        assert!(shorter.compare_to(&longer) < 0);
        assert!(longer.compare_to(&shorter) > 0);

        let mismatch = AlignNode::new(AlignNode::MISMATCH, 1, None).unwrap();
        let matched = AlignNode::new(AlignNode::MATCH, 1, None).unwrap();
        assert!(mismatch.compare_to(&matched) < 0);
    }

    #[test]
    fn max_alignment_score_node_stores_fields() {
        let trace = TraceNode::with_next(5.0, TraceNode::TYPE_MATCH, None);
        let node = MaxAlignmentScoreNode::new(Some(trace.clone()), 10, None).unwrap();
        assert_eq!(node.trace_node.as_ref(), &trace);
        assert_eq!(node.n_consensus_bases, 10);
        assert_eq!(node.next, None);
        assert!(!node.haplotype_built);

        let tail = MaxAlignmentScoreNode::new(
            Some(TraceNode::with_next(1.0, TraceNode::TYPE_MATCH, None)),
            1,
            None,
        )
        .unwrap();
        let head = MaxAlignmentScoreNode::new(
            Some(TraceNode::with_next(1.0, TraceNode::TYPE_MATCH, None)),
            2,
            Some(Box::new(tail.clone())),
        )
        .unwrap();
        assert_eq!(head.next.as_deref(), Some(&tail));
    }

    #[test]
    fn max_alignment_score_node_validates_and_displays() {
        let mut node = MaxAlignmentScoreNode::new(
            Some(TraceNode::with_next(1.0, TraceNode::TYPE_MATCH, None)),
            1,
            None,
        )
        .unwrap();
        node.haplotype_built = true;
        assert!(node.haplotype_built);

        assert_eq!(
            MaxAlignmentScoreNode::new(None, 1, None),
            Err(MaxAlignmentScoreNodeError::NullTraceNode)
        );
        assert_eq!(
            MaxAlignmentScoreNode::new(
                Some(TraceNode::with_next(1.0, TraceNode::TYPE_MATCH, None)),
                0,
                None
            ),
            Err(MaxAlignmentScoreNodeError::InvalidConsensusBases(0))
        );
        assert_eq!(
            MaxAlignmentScoreNode::new(
                Some(TraceNode::with_next(1.0, TraceNode::TYPE_MATCH, None)),
                -1,
                None
            ),
            Err(MaxAlignmentScoreNodeError::InvalidConsensusBases(-1))
        );
        let display = MaxAlignmentScoreNode::new(
            Some(TraceNode::with_next(1.0, TraceNode::TYPE_MATCH, None)),
            42,
            None,
        )
        .unwrap()
        .to_string();
        assert!(display.contains("len=42"));
    }

    #[test]
    fn haplotype_container_stores_limit_and_starts_empty() {
        let container = HaplotypeContainer::new(5).unwrap();
        assert_eq!(container.size(), 0);
        assert_eq!(container.limit, 5);
        assert!(container.to_array().is_empty());
    }

    #[test]
    fn haplotype_container_rejects_invalid_limits() {
        assert_eq!(
            HaplotypeContainer::new(0),
            Err(HaplotypeContainerError::InvalidLimit(0))
        );
        assert_eq!(
            HaplotypeContainer::new(-1),
            Err(HaplotypeContainerError::InvalidLimit(-1))
        );
    }

    #[test]
    fn haplotype_container_adds_to_head_like_java_list() {
        let mut container = HaplotypeContainer::new(5).unwrap();
        let low = make_haplotype(5);
        let high = make_haplotype(7);

        container.add(low.clone());
        container.add(high.clone());

        assert_eq!(container.size(), 2);
        assert_eq!(container.to_array(), vec![high, low]);
    }

    #[test]
    fn haplotype_container_rejects_lower_depth_when_at_capacity() {
        let mut container = HaplotypeContainer::new(2).unwrap();
        let first = make_haplotype(10);
        let second = make_haplotype(10);

        container.add(first.clone());
        container.add(second.clone());
        container.add(make_haplotype(1));

        assert_eq!(container.size(), 2);
        assert_eq!(container.to_array(), vec![second, first]);
    }

    #[test]
    fn haplotype_container_evicts_lower_depth_for_incoming_haplotype() {
        let mut container = HaplotypeContainer::new(2).unwrap();
        let low_a = make_haplotype(1);
        let low_b = make_haplotype(1);
        let high = make_haplotype(10);

        container.add(low_a);
        container.add(low_b.clone());
        container.add(high.clone());

        assert_eq!(container.size(), 2);
        assert_eq!(container.to_array(), vec![high, low_b]);
    }

    #[test]
    fn trace_matrix_dimensions_and_validation_match_java() {
        let matrix = TraceMatrix::new(5).unwrap();
        assert_eq!(matrix.n_row, 5);
        assert_eq!(matrix.matrix.len(), 5);
        assert_eq!(matrix.matrix[0].len(), TraceMatrix::DEFAULT_COLUMN_CAPACITY);
        assert_eq!(
            TraceMatrix::new(0),
            Err(TraceMatrixError::InvalidRowCount(0))
        );
        assert_eq!(
            TraceMatrix::new(-3),
            Err(TraceMatrixError::InvalidRowCount(-3))
        );
    }

    #[test]
    fn trace_matrix_next_col_and_set_match_java() {
        let mut matrix = TraceMatrix::new(5).unwrap();
        assert_eq!(matrix.next_col().unwrap(), 0);
        assert_eq!(matrix.next_col().unwrap(), 1);
        assert_eq!(matrix.next_col().unwrap(), 2);

        let mut empty = TraceMatrix::new(5).unwrap();
        assert_eq!(empty.set(0, 1, 1), Err(TraceMatrixError::NoColumns));

        let mut bounded = TraceMatrix::new(3).unwrap();
        bounded.next_col().unwrap();
        assert_eq!(
            bounded.set(5, 1, 1),
            Err(TraceMatrixError::RowOutOfRange { row: 5, max: 2 })
        );
        assert_eq!(
            bounded.set(-1, 1, 1),
            Err(TraceMatrixError::RowOutOfRange { row: -1, max: 2 })
        );

        bounded.set(0, 1, 1).unwrap();
        assert_ne!(bounded.matrix[0][0], 0);
    }

    #[test]
    fn trace_matrix_expands_and_prints() {
        let mut matrix = TraceMatrix::new(3).unwrap();
        for _ in 0..TraceMatrix::DEFAULT_COLUMN_CAPACITY + 10 {
            matrix.next_col().unwrap();
        }
        assert!(matrix.matrix[0].len() >= TraceMatrix::DEFAULT_COLUMN_CAPACITY + 10);
        assert!(!matrix.to_string().is_empty());

        let mut labelled = TraceMatrix::new(2).unwrap();
        labelled.next_col().unwrap();
        labelled
            .set(
                0,
                TraceNode::TYPE_MATCH.into(),
                TraceNode::TYPE_MATCH.into(),
            )
            .unwrap();
        let display = labelled.matrix_string(Some(b"AC"), Some(b"G"), 0);
        assert!(display.contains('A'));
        assert!(display.contains('G'));
        assert!(display.contains("a--"));
    }

    #[test]
    fn type_list_builds_and_merges_alignments() {
        assert_eq!(TypeList::new().to_alignment(false), None);

        let mut list = TypeList::new();
        list.add(AlignNode::MATCH).unwrap();
        let head = list.to_alignment(false).unwrap();
        assert_eq!(head.node_type, AlignNode::MATCH);
        assert_eq!(head.n, 1);
        assert_eq!(head.next, None);

        let mut merged = TypeList::new();
        for _ in 0..5 {
            merged.add(AlignNode::MATCH).unwrap();
        }
        let head = merged.to_alignment(false).unwrap();
        assert_eq!(head.n, 5);
        assert_eq!(head.next, None);

        let mut chain = TypeList::new();
        for _ in 0..3 {
            chain.add(AlignNode::MATCH).unwrap();
        }
        chain.add(AlignNode::MISMATCH).unwrap();
        for _ in 0..2 {
            chain.add(AlignNode::MATCH).unwrap();
        }
        let head = chain.to_alignment(false).unwrap();
        assert_eq!(head.node_type, AlignNode::MATCH);
        assert_eq!(head.n, 3);
        let middle = head.next.as_deref().unwrap();
        assert_eq!(middle.node_type, AlignNode::MISMATCH);
        assert_eq!(middle.n, 1);
        let tail = middle.next.as_deref().unwrap();
        assert_eq!(tail.node_type, AlignNode::MATCH);
        assert_eq!(tail.n, 2);
        assert_eq!(tail.next, None);
    }

    #[test]
    fn type_list_reverse_clear_and_expand_match_java() {
        let mut list = TypeList::new();
        list.add(AlignNode::MATCH).unwrap();
        list.add(AlignNode::MISMATCH).unwrap();
        list.add(AlignNode::INS).unwrap();
        let head = list.to_alignment(true).unwrap();
        assert_eq!(head.node_type, AlignNode::INS);
        assert_eq!(head.next.as_deref().unwrap().node_type, AlignNode::MISMATCH);
        assert_eq!(
            head.next
                .as_deref()
                .unwrap()
                .next
                .as_deref()
                .unwrap()
                .node_type,
            AlignNode::MATCH
        );

        list.clear();
        assert_eq!(list.to_alignment(false), None);

        let mut expanded = TypeList::new();
        for index in 0..200 {
            expanded
                .add(if index % 2 == 0 {
                    AlignNode::MATCH
                } else {
                    AlignNode::MISMATCH
                })
                .unwrap();
        }
        let mut count = 0;
        let mut node = expanded.to_alignment(false);
        while let Some(current) = node {
            count += 1;
            node = current.next;
        }
        assert_eq!(count, 200);

        assert_eq!(
            TypeList::new().add(99),
            Err(TypeListError::TypeOutOfRange(99))
        );
    }

    #[test]
    fn type_list_clone_variants_are_independent() {
        let mut original = TypeList::new();
        original.add(AlignNode::MATCH).unwrap();
        original.add(AlignNode::MATCH).unwrap();

        let mut cloned = original.clone();
        cloned.add(AlignNode::MISMATCH).unwrap();
        let original_head = original.to_alignment(false).unwrap();
        assert_eq!(original_head.n, 2);
        assert_eq!(original_head.next, None);

        let mut source = TypeList::new();
        source.add(AlignNode::DEL).unwrap();
        source.add(AlignNode::INS).unwrap();
        let copied = TypeList::new_from(Some(&source));
        let head = copied.to_alignment(false).unwrap();
        assert_eq!(head.node_type, AlignNode::DEL);
        assert_eq!(head.next.as_deref().unwrap().node_type, AlignNode::INS);
        assert_eq!(TypeList::new_from(None).to_alignment(false), None);

        original.clone_from_option(Some(&source));
        let head = original.to_alignment(false).unwrap();
        assert_eq!(head.node_type, AlignNode::DEL);
        assert_eq!(head.next.as_deref().unwrap().node_type, AlignNode::INS);
        original.clone_from_option(None);
        assert_eq!(original.to_alignment(false), None);
    }

    #[test]
    fn kmer_alignment_builder_expands_linear_trace_like_java_get_alignment() {
        let tail = TraceNode::with_next(3.0, TraceNode::TYPE_MATCH, None);
        let middle = TraceNode::with_next(4.0, TraceNode::TYPE_MISMATCH, Some(Box::new(tail)));
        let root = TraceNode::with_next(5.0, TraceNode::TYPE_MATCH, Some(Box::new(middle)));
        let mut builder = KmerAlignmentBuilder::new();

        let fwd = builder.alignments_from_trace(Some(&root), false).unwrap();
        assert_eq!(fwd.len(), 1);
        assert_eq!(fwd[0].cigar_string(), "1=1X1=");

        let rev = builder.alignments_from_trace(Some(&root), true).unwrap();
        assert_eq!(rev.len(), 1);
        assert_eq!(rev[0].cigar_string(), "1=1X1=");
        assert!(builder.cached_type_lists() > 0);
    }

    #[test]
    fn kmer_alignment_builder_preserves_trace_order_difference_for_asymmetric_paths() {
        let tail = TraceNode::with_next(3.0, TraceNode::TYPE_GAP_REF, None);
        let root = TraceNode::with_next(5.0, TraceNode::TYPE_MATCH, Some(Box::new(tail)));
        let mut builder = KmerAlignmentBuilder::new();

        let fwd = builder.alignments_from_trace(Some(&root), false).unwrap();
        assert_eq!(fwd[0].cigar_string(), "1I1=");

        let rev = builder.alignments_from_trace(Some(&root), true).unwrap();
        assert_eq!(rev[0].cigar_string(), "1=1I");
    }

    #[test]
    fn kmer_alignment_builder_expands_branches_into_distinct_alignments() {
        let main_tail = TraceNode::with_next(3.0, TraceNode::TYPE_MATCH, None);
        let branch_tail = TraceNode::with_next(2.0, TraceNode::TYPE_GAP_CON, None);
        let branch =
            TraceNode::with_next(2.0, TraceNode::TYPE_GAP_REF, Some(Box::new(branch_tail)));
        let root = TraceNode::new(
            5.0,
            TraceNode::TYPE_MATCH,
            Some(Box::new(main_tail)),
            Some(Box::new(branch)),
        )
        .unwrap();
        let mut builder = KmerAlignmentBuilder::new();
        let mut cigars = builder
            .alignments_from_trace(Some(&root), true)
            .unwrap()
            .into_iter()
            .map(|alignment| alignment.cigar_string())
            .collect::<Vec<_>>();
        cigars.sort();

        assert_eq!(cigars, vec!["1I1D", "2="]);
    }

    #[test]
    fn kmer_alignment_builder_handles_empty_and_invalid_trace() {
        let mut builder = KmerAlignmentBuilder::new();
        assert!(
            builder
                .alignments_from_trace(None, false)
                .unwrap()
                .is_empty()
        );

        let invalid = TraceNode::with_next(0.0, TraceNode::TYPE_NONE, None);
        assert_eq!(
            builder.alignments_from_trace(Some(&invalid), false),
            Err(KmerAlignmentBuilderError::TypeList(
                TypeListError::TypeOutOfRange(TraceNode::TYPE_NONE)
            ))
        );
    }

    #[test]
    fn kmer_aligner_initializes_from_active_region_and_extends_matching_reference() {
        let mut aligner =
            KmerAligner::new(KmerUtil::new(4).unwrap(), AlignmentWeight::defaults(), true).unwrap();
        let active_region = make_active_region(10, 0, 12);
        aligner.init(active_region).unwrap();

        assert!(!aligner.is_reverse());
        assert!(!aligner.is_allow_end_deletion());
        assert_eq!(aligner.consensus(), b"AAAA");
        for base in b"CCCCGGGGTTTT" {
            let keep_going = aligner
                .add_base(Base::from_char(char::from(*base)).unwrap())
                .unwrap();
            if *base != b'T' {
                assert!(keep_going);
            }
        }

        assert_eq!(aligner.consensus(), b"AAAACCCCGGGGTTTT");
        assert!(aligner.max_alignment_score() > 0.0);

        let mut haplotypes = aligner.get_haplotypes(&StaticCountMap, false).unwrap();
        assert_eq!(haplotypes.len(), 1);
        let haplotype = haplotypes.pop().unwrap();
        assert_eq!(haplotype.sequence, b"AAAACCCCGGGGTTTT");
        assert_eq!(haplotype.stats.min, 10);
        assert_eq!(haplotype.alignment.cigar_string(), "16=");
        assert!(haplotype.trace_matrix.is_some());
    }

    #[test]
    fn kmer_aligner_builds_mismatch_haplotype_from_crafted_consensus() {
        let mut aligner = KmerAligner::new(
            KmerUtil::new(4).unwrap(),
            AlignmentWeight::defaults(),
            false,
        )
        .unwrap();
        aligner.init(make_active_region(10, 0, 12)).unwrap();
        for base in b"TCCCGGGGTTTT" {
            aligner
                .add_base(Base::from_char(char::from(*base)).unwrap())
                .unwrap();
        }

        let haplotypes = aligner.get_haplotypes(&StaticCountMap, false).unwrap();
        assert_eq!(haplotypes.len(), 1);
        assert_eq!(haplotypes[0].sequence, b"AAAATCCCGGGGTTTT");
        assert_eq!(haplotypes[0].alignment.cigar_string(), "4=1X11=");
    }

    #[test]
    fn kmer_aligner_chooses_reverse_for_left_end_regions() {
        let mut aligner = KmerAligner::new(
            KmerUtil::new(4).unwrap(),
            AlignmentWeight::defaults(),
            false,
        )
        .unwrap();
        let active_region = make_active_region(10, -1, 12);
        aligner.init(active_region).unwrap();

        assert!(aligner.is_reverse());
        assert!(aligner.is_allow_end_deletion());
        assert_eq!(aligner.consensus(), b"TTTT");
    }

    #[test]
    fn kmer_aligner_state_save_restore_and_capacity_match_java_shape() {
        let kmer_util = KmerUtil::new(4).unwrap();
        let mut aligner =
            KmerAligner::new(kmer_util.clone(), AlignmentWeight::defaults(), false).unwrap();
        aligner.init(make_active_region(10, 0, 12)).unwrap();
        aligner.set_max_state(1).unwrap();

        let first_kmer = kmer_util.encode(b"AAAC").unwrap();
        let second_kmer = kmer_util.encode(b"AAAG").unwrap();
        aligner
            .save_state(
                Some(first_kmer),
                Some(Base::C),
                1,
                Some(KmerHashSet::new()),
                0,
            )
            .unwrap();
        aligner
            .save_state(
                Some(second_kmer.clone()),
                Some(Base::G),
                5,
                Some(KmerHashSet::new()),
                0,
            )
            .unwrap();

        assert!(aligner.has_cached_states());
        let restored = aligner.restore_state().unwrap().unwrap();
        assert_eq!(restored.kmer, second_kmer.words());
        assert_eq!(restored.min_depth, 5);
        assert_eq!(aligner.consensus(), b"AAAAG");
        assert!(!aligner.has_cached_states());
    }

    #[test]
    fn kmer_aligner_validates_constructor_and_mutators() {
        assert!(matches!(
            KmerAligner::new(
                KmerUtil::new(3).unwrap(),
                AlignmentWeight::defaults(),
                false
            ),
            Err(KmerAlignerError::KmerSizeTooSmall { .. })
        ));
        let mut aligner = KmerAligner::new(
            KmerUtil::new(4).unwrap(),
            AlignmentWeight::defaults(),
            false,
        )
        .unwrap();
        assert_eq!(aligner.max_state(), KmerAligner::DEFAULT_MAX_STATE);
        assert_eq!(
            aligner.set_max_state(0),
            Err(KmerAlignerError::InvalidMaxState(0))
        );
        assert_eq!(
            aligner.add_base(Base::A),
            Err(KmerAlignerError::NotInitialized)
        );
    }

    #[test]
    fn trace_node_container_stores_and_chains() {
        let trace = TraceNode::with_next(5.0, TraceNode::TYPE_MATCH, None);
        let container = state::TraceNodeContainer::new(3, trace.clone(), None).unwrap();
        assert_eq!(container.index, 3);
        assert_eq!(container.node.as_ref(), &trace);
        assert_eq!(container.next, None);

        let tail = state::TraceNodeContainer::new(
            0,
            TraceNode::with_next(1.0, TraceNode::TYPE_MATCH, None),
            None,
        )
        .unwrap();
        let head = state::TraceNodeContainer::new(
            1,
            TraceNode::with_next(1.0, TraceNode::TYPE_MATCH, None),
            Some(Box::new(tail.clone())),
        )
        .unwrap();
        assert_eq!(head.next.as_deref(), Some(&tail));
        assert_eq!(tail.next, None);
        assert_eq!(
            state::TraceNodeContainer::new(-1, trace, None),
            Err(state::TraceNodeContainerError::NegativeIndex(-1))
        );
    }

    #[test]
    fn state_stack_node_stores_fields_and_restores_subset() {
        let max = MaxAlignmentScoreNode::new(
            Some(TraceNode::with_next(1.0, TraceNode::TYPE_MATCH, None)),
            10,
            None,
        )
        .unwrap();
        let node = make_state_node(TestStateNodeArgs {
            kmer: vec![0x123, 0x456],
            next_base: Base::A,
            consensus_size: 5,
            max_alignment_score: 50.0,
            max_alignment_score_node: max,
            min_depth: 3,
            next_node_down: None,
            repeat_count: 1,
        });

        assert_eq!(node.consensus_size, 5);
        assert_eq!(node.next_base, Base::A);
        assert_eq!(node.max_alignment_score, 50.0);
        assert_eq!(node.min_depth, 3);
        assert_eq!(node.repeat_count, 1);
        assert_eq!(node.next_node_down, None);
        assert_eq!(node.next_node_up, None);
        assert!(node.max_alignment_score_node.is_some());

        let restored = node.restored_state();
        assert_eq!(restored.kmer, node.kmer);
        assert_eq!(restored.kmer_hash, node.kmer_hash);
        assert_eq!(restored.consensus_size, node.consensus_size);
        assert_eq!(restored.min_depth, node.min_depth);
        assert_eq!(restored.repeat_count, node.repeat_count);
    }

    #[test]
    fn state_stack_nodes_chain_downward() {
        let bottom = make_state_node(TestStateNodeArgs {
            kmer: vec![0x123, 0x456],
            next_base: Base::A,
            consensus_size: 5,
            max_alignment_score: 50.0,
            max_alignment_score_node: MaxAlignmentScoreNode::new(
                Some(TraceNode::with_next(1.0, TraceNode::TYPE_MATCH, None)),
                10,
                None,
            )
            .unwrap(),
            min_depth: 3,
            next_node_down: None,
            repeat_count: 1,
        });
        let top = make_state_node(TestStateNodeArgs {
            kmer: vec![1, 2],
            next_base: Base::C,
            consensus_size: 6,
            max_alignment_score: 60.0,
            max_alignment_score_node: MaxAlignmentScoreNode::new(
                Some(TraceNode::with_next(2.0, TraceNode::TYPE_MATCH, None)),
                1,
                None,
            )
            .unwrap(),
            min_depth: 4,
            next_node_down: Some(Box::new(bottom.clone())),
            repeat_count: 0,
        });
        assert_eq!(top.next_node_down.as_deref(), Some(&bottom));
    }

    struct TestStateNodeArgs {
        kmer: Vec<u32>,
        next_base: Base,
        consensus_size: i32,
        max_alignment_score: f32,
        max_alignment_score_node: MaxAlignmentScoreNode,
        min_depth: i32,
        next_node_down: Option<Box<state::StateStackNode>>,
        repeat_count: i32,
    }

    fn make_state_node(args: TestStateNodeArgs) -> state::StateStackNode {
        state::StateStackNode::new(state::StateStackNodeFields {
            kmer: args.kmer,
            next_base: args.next_base,
            consensus_size: args.consensus_size,
            align_container: None,
            gap_ref_container: None,
            gap_con_container: None,
            max_alignment_score: args.max_alignment_score,
            max_alignment_score_node: Some(Box::new(args.max_alignment_score_node)),
            min_depth: args.min_depth,
            next_node_down: args.next_node_down,
            kmer_hash: KmerHashSet::new(),
            repeat_count: args.repeat_count,
        })
    }

    #[test]
    fn defaults_match_published_constants() {
        let weight = AlignmentWeight::defaults();
        assert_eq!(AlignmentWeight::DEFAULT_MATCH, 10.0);
        assert_eq!(AlignmentWeight::DEFAULT_MISMATCH, -10.0);
        assert_eq!(AlignmentWeight::DEFAULT_GAP_OPEN, -40.0);
        assert_eq!(AlignmentWeight::DEFAULT_GAP_EXTEND, -4.0);
        assert_eq!(AlignmentWeight::DEFAULT_INIT_SCORE, 0.0);
        assert!((weight.new_gap - (weight.gap_open + weight.gap_extend)).abs() < EPS);
        assert!(weight.is_default());
    }

    #[test]
    fn factory_normalizes_and_validates_signs() {
        let weight = AlignmentWeight::new(-5.0, 5.0, 30.0, 3.0, -2.0).unwrap();
        assert!(weight.match_score > 0.0);
        assert!(weight.mismatch < 0.0);
        assert!(weight.gap_open < 0.0);
        assert!(weight.gap_extend < 0.0);
        assert!(weight.init_score >= 0.0);

        assert_eq!(
            AlignmentWeight::new(0.0, -10.0, -40.0, -4.0, 0.0),
            Err(AlignmentWeightError::ZeroWeight)
        );
        assert_eq!(
            AlignmentWeight::new(10.0, 0.0, -40.0, -4.0, 0.0),
            Err(AlignmentWeightError::ZeroWeight)
        );
        assert_eq!(
            AlignmentWeight::new(10.0, -10.0, -40.0, 0.0, 0.0),
            Err(AlignmentWeightError::ZeroWeight)
        );
        assert_eq!(
            AlignmentWeight::new(10.0, -10.0, 0.0, -4.0, 0.0)
                .unwrap()
                .gap_open,
            0.0
        );
    }

    #[test]
    fn parses_strings_with_defaults_and_delimiters() {
        assert!(AlignmentWeight::parse(None).unwrap().is_default());
        assert!(AlignmentWeight::parse(Some("   ")).unwrap().is_default());
        assert_eq!(
            AlignmentWeight::parse(Some("5,5,40,4"))
                .unwrap()
                .match_score,
            5.0
        );
        assert_eq!(
            AlignmentWeight::parse(Some("5,5,40,4,100"))
                .unwrap()
                .init_score,
            100.0
        );

        for value in ["(5,5,40,4)", "<5,5,40,4>", "[5,5,40,4]", "{5,5,40,4}"] {
            assert_eq!(
                AlignmentWeight::parse(Some(value)).unwrap().match_score,
                5.0
            );
        }

        assert!(AlignmentWeight::parse(Some(",,,")).unwrap().is_default());
        let partial = AlignmentWeight::parse(Some("20")).unwrap();
        assert_eq!(partial.match_score, 20.0);
        assert_eq!(partial.mismatch, -10.0);
    }

    #[test]
    fn parses_exponential_hex_and_whitespace() {
        let exp = AlignmentWeight::parse(Some("1.0e1,1.0e1,4.0e1,4.0e0")).unwrap();
        assert_eq!(exp.match_score, 10.0);
        assert_eq!(exp.mismatch, -10.0);
        assert_eq!(exp.gap_open, -40.0);
        assert_eq!(exp.gap_extend, -4.0);

        let hex = AlignmentWeight::parse(Some("0xa,0xa,0x28,0x4")).unwrap();
        assert_eq!(hex.match_score, 10.0);

        let spaced = AlignmentWeight::parse(Some("5 ,  5 , 40 , 4")).unwrap();
        assert_eq!(spaced.match_score, 5.0);
        assert_eq!(spaced.mismatch, -5.0);
    }

    #[test]
    fn parse_errors_match_java_cases() {
        assert!(matches!(
            AlignmentWeight::parse(Some("(5,5,40,4")),
            Err(AlignmentWeightError::MismatchedDelimiters(_))
        ));
        assert!(matches!(
            AlignmentWeight::parse(Some("[5,5,40,4")),
            Err(AlignmentWeightError::MismatchedDelimiters(_))
        ));
        assert!(matches!(
            AlignmentWeight::parse(Some("5,5,40,4)")),
            Err(AlignmentWeightError::UnopenedDelimiter(_))
        ));
        assert_eq!(
            AlignmentWeight::parse(Some("1,2,3,4,5,6")),
            Err(AlignmentWeightError::TooManyValues(6))
        );
        assert_eq!(
            AlignmentWeight::parse(Some("5,5,xxx,4")),
            Err(AlignmentWeightError::InvalidNumber("xxx".to_owned()))
        );
    }

    #[test]
    fn mutators_preserve_other_fields() {
        let base = AlignmentWeight::defaults();
        let changed = base.with_match(20.0).unwrap();
        assert_eq!(changed.match_score, 20.0);
        assert_eq!(changed.mismatch, base.mismatch);
        assert_eq!(base.with_mismatch(20.0).unwrap().mismatch, -20.0);
        assert_eq!(base.with_gap_open(50.0).unwrap().gap_open, -50.0);
        assert_eq!(base.with_gap_extend(5.0).unwrap().gap_extend, -5.0);
        assert_eq!(base.with_initial_score(50.0).unwrap().init_score, 50.0);
    }

    #[test]
    fn score_helpers_match_java() {
        let base = AlignmentWeight::defaults();
        assert_eq!(base.initial_score(25).unwrap(), 250.0);
        assert_eq!(base.initial_score(1).unwrap(), 10.0);
        assert_eq!(
            base.with_initial_score(75.0)
                .unwrap()
                .initial_score(25)
                .unwrap(),
            75.0
        );
        assert_eq!(
            base.initial_score(0),
            Err(AlignmentWeightError::InvalidKSize(0))
        );
        assert_eq!(base.max_exclusive_gap_size(5).unwrap(), 2);
        assert!(base.max_exclusive_gap_size(1).unwrap() < i32::MAX);
    }

    #[test]
    fn equality_hash_and_formatting() {
        let a = AlignmentWeight::defaults();
        let b = AlignmentWeight::new(10.00001, -10.0, -40.0, -4.0, 0.0).unwrap();
        assert_eq!(a, b);
        assert_ne!(a, a.with_match(20.0).unwrap());
        assert_eq!(hash(&a), hash(&AlignmentWeight::defaults()));

        assert_eq!(a.to_string().matches(',').count(), 4);
        assert!(a.format_with_precision(2).unwrap().contains("10.00"));
        assert_eq!(
            a.format_with_precision(-1),
            Err(AlignmentWeightError::InvalidPrecision(-1))
        );
        assert_eq!(
            a.format_with_precision(11),
            Err(AlignmentWeightError::InvalidPrecision(11))
        );
        assert!(!a.with_match(20.0).unwrap().is_default());
    }

    fn hash(value: &AlignmentWeight) -> u64 {
        let mut hasher = DefaultHasher::new();
        value.hash(&mut hasher);
        hasher.finish()
    }

    #[derive(Debug)]
    struct StaticCountMap;

    impl CountMap for StaticCountMap {
        fn get(&self, _kmer: &KmerKey) -> u32 {
            10
        }

        fn set(&mut self, _sample: InputSample) -> Result<(), CountMapError> {
            Ok(())
        }

        fn abort(&self) {}

        fn is_aborted(&self) -> bool {
            false
        }
    }

    fn make_haplotype(stat_value: i32) -> Haplotype {
        let active_region = make_active_region(stat_value, 2, 10);
        let stats = active_region.stats;
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

    fn make_active_region(stat_value: i32, start_kmer: i32, end_kmer: i32) -> ActiveRegion {
        let reference = ReferenceSequence::new("chr1", 16, Some(digest()), Some("test")).unwrap();
        let ref_region = ReferenceRegion::whole(reference, b"AAAACCCCGGGGTTTT", 0).unwrap();
        let count = vec![stat_value; 13];
        let kmer_util = KmerUtil::new(4).unwrap();
        ActiveRegion::new(ref_region, start_kmer, end_kmer, &count, &kmer_util).unwrap()
    }

    fn digest() -> Digest {
        Digest::new(vec![0; 16], "MD5").unwrap()
    }
}
