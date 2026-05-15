//! Diagnostic test: reproduce the J-R:4-119 active region traversal using the
//! real post-`kmercount:5` k-mer count map from kanalyze and the actual J-R
//! reference sequence. Run with `KESTREL_RUN_JR_DIAGNOSTIC=1` (off by default
//! so it stays out of the normal test rotation).

use std::collections::HashMap;
use std::fs;
use std::path::PathBuf;

use kanalyze::util::{KmerKey, KmerUtil};
use kestrel::activeregion::ActiveRegion;
use kestrel::counter::{CountMap, CountMapError};
use kestrel::interval::RegionInterval;
use kestrel::io::InputSample;
use kestrel::refreader::{ReferenceRegion, ReferenceSequence};
use kestrel::runner::{RunConfig, graph_haplotypes_for_test};
use kestrel::variant::{VariantCall, VariantCaller};

const J_R_REFERENCE: &str = "TGGGGGGGCGGTGGAGCCCGGGGCCGGCCTGGTGTCCGTGCCCGAGGTGACACCGTGGGCTGCGGGCGCGGTGGAGCCCGGGGCCGGCCTGCTCTCCGGTGCCGAGGTGACACCGTGGGC";
const FIVE_C_N_REFERENCE: &str = "TGGGGGGGCGGTGGAGCCCGTGGCCGGCCTGCTCTCCGGGGCCGAGGTGACACCGTGGGCTTGGGGGGCGGTGGAGCCCGGGGCCGGCCTGGTGTCCGGGGCTGAGGTGACATCGTGGGC";

#[test]
fn jr_traversal_chain_growth() {
    if std::env::var_os("KESTREL_RUN_JR_DIAGNOSTIC").is_none() {
        eprintln!("set KESTREL_RUN_JR_DIAGNOSTIC=1 to run J-R traversal diagnostic");
        return;
    }

    let kmer_util = KmerUtil::new(20).unwrap();
    let counts = load_counts(&kmer_util);
    let count_map = TsvCountMap { counts };

    let ref_seq =
        ReferenceSequence::new("J-R", J_R_REFERENCE.len() as i32, None, Some("test")).unwrap();
    let ref_region = ReferenceRegion::whole(ref_seq, J_R_REFERENCE.as_bytes(), 0).unwrap();

    let counts_array = collect_counts(&kmer_util, &ref_region, &count_map);
    let region = ActiveRegion::new(ref_region.clone(), 4, 100, &counts_array, &kmer_util).unwrap();

    eprintln!(
        "[JR] region start_index={} end_index={} end_kmer_index={} ref_len={}",
        region.start_index,
        region.end_index,
        region.end_kmer_index,
        region.end_index - region.start_index + 1
    );

    let mut config = RunConfig::default();
    config.k_size = 20;
    config.min_kmer_count = 5;
    config.minimum_difference = 5;
    config.count_reverse_kmers = true;
    config.max_haplotypes = 15;
    config.max_aligner_state = 10;
    config.max_repeat_count = 0;
    config.peak_scan_length = 7;

    let haplotypes = graph_haplotypes_for_test(&config, &kmer_util, &count_map, &region).unwrap();

    eprintln!("[JR] haplotypes produced: {}", haplotypes.len());
    for (idx, hap) in haplotypes.iter().enumerate() {
        eprintln!(
            "[JR] hap {}: len={} min={} cigar={}",
            idx,
            hap.sequence.len(),
            hap.stats.min,
            hap.alignment.cigar_string()
        );
    }

    // Java reports 0 haplotypes for J-R:4-119. We expect Rust to also produce
    // 0 after the algorithmic divergence is fixed.
    assert!(
        haplotypes.is_empty() || (haplotypes.len() == 1 && haplotypes[0].is_wildtype()),
        "Expected 0 non-wildtype haplotypes for J-R:4-119 to match Java, got {}",
        haplotypes.len()
    );
}

#[test]
fn five_c_n_insertion_region_diagnostic() {
    if std::env::var_os("KESTREL_RUN_MUC1_INS_DIAGNOSTIC").is_none() {
        eprintln!("set KESTREL_RUN_MUC1_INS_DIAGNOSTIC=1 to run MUC1 INS diagnostic");
        return;
    }

    let kmer_util = KmerUtil::new(20).unwrap();
    let counts = load_counts(&kmer_util);
    let count_map = TsvCountMap { counts };

    let ref_region = reference_region("5C-N", FIVE_C_N_REFERENCE);
    let counts_array = collect_counts(&kmer_util, &ref_region, &count_map);
    let region = ActiveRegion::new(ref_region.clone(), 0, 34, &counts_array, &kmer_util).unwrap();

    let mut config = RunConfig::default();
    config.k_size = 20;
    config.min_kmer_count = 5;
    config.minimum_difference = 5;
    config.count_reverse_kmers = true;
    config.max_haplotypes = 15;
    config.max_aligner_state = 10;
    config.max_repeat_count = 0;
    config.peak_scan_length = 7;

    let haplotypes = graph_haplotypes_for_test(&config, &kmer_util, &count_map, &region).unwrap();

    eprintln!("[5C-N] haplotypes produced: {}", haplotypes.len());
    for (idx, hap) in haplotypes.iter().enumerate() {
        eprintln!(
            "[5C-N] hap {}: len={} min={} cigar={} seq={}",
            idx,
            hap.sequence.len(),
            hap.stats.min,
            hap.alignment.cigar_string(),
            String::from_utf8_lossy(&hap.sequence)
        );
        for row in hap.alignment_string(0).unwrap() {
            eprintln!("[5C-N]   {row}");
        }
    }

    let mut caller = VariantCaller::new();
    caller.init(region);
    caller.set_variant_call_by_reference();
    for haplotype in haplotypes {
        caller.add(haplotype).unwrap();
    }
    let variants = caller.variants();
    for variant in &variants {
        eprintln!(
            "[5C-N] variant pos={} ref={} alt={} gdp={} dp={}",
            variant.vcf_pos(),
            variant.vcf_ref().unwrap(),
            variant.vcf_alt().unwrap(),
            variant.data().variant_depth,
            variant.data().locus_depth
        );
    }

    assert!(
        variants.iter().any(|variant| {
            variant.vcf_pos() == 26
                && variant.vcf_ref().unwrap() == "G"
                && variant.vcf_alt().unwrap() == "GGGTGGAGCCCGGGGCCGG"
        }),
        "expected Java's 5C-N:26 18-base insertion in diagnostic output"
    );
}

fn load_counts(util: &KmerUtil) -> HashMap<KmerKey, u32> {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/jr_counts.tsv");
    let contents = fs::read_to_string(&path)
        .unwrap_or_else(|err| panic!("failed to read {}: {err}", path.display()));
    let mut counts = HashMap::new();
    for line in contents.lines() {
        let mut parts = line.split('\t');
        let Some(seq) = parts.next() else { continue };
        let Some(count_str) = parts.next() else {
            continue;
        };
        let count: u32 = count_str.parse().expect("invalid count");
        if seq.len() != util.k_size() {
            continue;
        }
        let kmer = util.encode(seq).expect("invalid kmer");
        counts.insert(kmer, count);
    }
    counts
}

fn reference_region(name: &str, sequence: &str) -> ReferenceRegion {
    let ref_seq = ReferenceSequence::new(name, sequence.len() as i32, None, Some("test")).unwrap();
    ReferenceRegion::whole(ref_seq, sequence.as_bytes(), 0).unwrap()
}

fn collect_counts(
    util: &KmerUtil,
    ref_region: &ReferenceRegion,
    counter: &TsvCountMap,
) -> Vec<i32> {
    let k_size = util.k_size();
    ref_region
        .sequence
        .windows(k_size)
        .map(|window| {
            let Ok(kmer) = util.encode(window) else {
                return 0;
            };
            let mut count = counter.get(&kmer) as i32;
            count += counter.get(&util.reverse_complement(&kmer)) as i32;
            count
        })
        .collect()
}

struct TsvCountMap {
    counts: HashMap<KmerKey, u32>,
}

impl CountMap for TsvCountMap {
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

#[allow(dead_code)]
fn _unused_region_interval(_interval: RegionInterval) {}
