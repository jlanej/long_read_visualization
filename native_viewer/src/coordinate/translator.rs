//! Coordinate translation between reference and assembly haplotype space.
//!
//! Ported from the Python `coordinate_mapper.py` module. Uses CIGAR-aware
//! projection for base-level precision when translating genomic coordinates
//! between the reference genome and phased haplotype assemblies.

use std::collections::HashMap;
use std::fs::File;
use std::io::Read;
use std::path::Path;

use anyhow::{Context, Result};
use serde::{Deserialize, Serialize};

use crate::genome::region::GenomicRegion;

/// A single alignment block from the mapping index.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AlignmentBlock {
    /// Reference start (0-based).
    #[serde(rename = "rs")]
    pub ref_start: u64,
    /// Reference end (0-based, exclusive).
    #[serde(rename = "re")]
    pub ref_end: u64,
    /// Assembly contig name.
    #[serde(rename = "ac")]
    pub asm_chrom: String,
    /// Assembly start (0-based).
    #[serde(rename = "as")]
    pub asm_start: u64,
    /// Assembly end (0-based, exclusive).
    #[serde(rename = "ae")]
    pub asm_end: u64,
    /// Strand: "+" or "-".
    #[serde(rename = "st")]
    pub strand: String,
    /// Mapping quality.
    #[serde(rename = "mq", default)]
    pub mapq: u8,
    /// CIGAR string (optional, for fine-grained projection).
    #[serde(rename = "cg", default)]
    pub cigar: Option<String>,
    /// Alignment type (P=primary, S=supplementary).
    #[serde(rename = "tp", default)]
    pub aln_type: Option<String>,
}

/// Loaded mapping index for one haplotype, keyed by reference chromosome.
#[derive(Debug, Clone)]
pub struct MappingIndex {
    /// Alignment blocks per reference chromosome, sorted by ref_start.
    pub blocks: HashMap<String, Vec<AlignmentBlock>>,
}

/// Result of a coordinate translation query for one haplotype.
#[derive(Debug, Clone)]
pub struct TranslationHit {
    /// Assembly contig name.
    pub asm_chrom: String,
    /// Assembly start (0-based).
    pub asm_start: u64,
    /// Assembly end (0-based, exclusive).
    pub asm_end: u64,
    /// Strand.
    pub strand: String,
    /// Event type: "alignment" or SV gap type.
    pub event_type: String,
}

/// Bidirectional coordinate translator for reference ↔ haplotype mapping.
pub struct CoordinateTranslator {
    /// Mapping indices: (sample_id, "hap1"|"hap2") → MappingIndex.
    indices: HashMap<(String, String), MappingIndex>,
}

impl CoordinateTranslator {
    pub fn new() -> Self {
        Self {
            indices: HashMap::new(),
        }
    }

    /// Load a mapping index from a compressed JSON file.
    pub fn load_index(&mut self, sample_id: &str, haplotype: &str, path: &Path) -> Result<()> {
        let index = load_mapping_index(path)?;
        self.indices
            .insert((sample_id.to_string(), haplotype.to_string()), index);
        Ok(())
    }

    /// Translate a reference region to assembly coordinates for both haplotypes.
    pub fn translate(
        &self,
        sample_id: &str,
        region: &GenomicRegion,
        min_mapq: u8,
    ) -> TranslationResult {
        let mut result = TranslationResult {
            hap1: Vec::new(),
            hap2: Vec::new(),
        };

        for (hap, hits) in [
            ("hap1", &mut result.hap1),
            ("hap2", &mut result.hap2),
        ] {
            let key = (sample_id.to_string(), hap.to_string());
            if let Some(index) = self.indices.get(&key) {
                *hits = query_index(index, region, min_mapq);
            }
        }

        result
    }

    /// Check if index data is loaded for a given sample/haplotype.
    pub fn has_index(&self, sample_id: &str, haplotype: &str) -> bool {
        self.indices
            .contains_key(&(sample_id.to_string(), haplotype.to_string()))
    }
}

impl Default for CoordinateTranslator {
    fn default() -> Self {
        Self::new()
    }
}

/// Result of translating a reference region.
#[derive(Debug, Clone)]
pub struct TranslationResult {
    pub hap1: Vec<TranslationHit>,
    pub hap2: Vec<TranslationHit>,
}

/// Load a mapping index from a .mapping.json.gz file.
pub fn load_mapping_index(path: &Path) -> Result<MappingIndex> {
    let file = File::open(path)
        .with_context(|| format!("Failed to open mapping index: {}", path.display()))?;

    let path_str = path.to_string_lossy();
    let data: HashMap<String, Vec<AlignmentBlock>> = if path_str.ends_with(".gz") {
        let decoder = flate2::read::GzDecoder::new(file);
        let mut reader = std::io::BufReader::new(decoder);
        let mut content = String::new();
        reader.read_to_string(&mut content)?;
        serde_json::from_str(&content)?
    } else {
        let reader = std::io::BufReader::new(file);
        serde_json::from_reader(reader)?
    };

    // Ensure blocks are sorted by ref_start within each chromosome
    let mut blocks = data;
    for block_list in blocks.values_mut() {
        block_list.sort_by_key(|b| b.ref_start);
    }

    Ok(MappingIndex { blocks })
}

/// Query a mapping index for alignment blocks overlapping a reference region.
pub fn query_index(
    index: &MappingIndex,
    region: &GenomicRegion,
    min_mapq: u8,
) -> Vec<TranslationHit> {
    let Some(blocks) = index.blocks.get(&region.chrom) else {
        return Vec::new();
    };

    let mut hits = Vec::new();

    for block in blocks {
        // Filter by mapping quality
        if block.mapq < min_mapq {
            continue;
        }

        // Check overlap with query region
        if block.ref_start >= region.end || block.ref_end <= region.start {
            continue;
        }

        // Project coordinates
        let (asm_start, asm_end) = if let Some(ref cigar) = block.cigar {
            // CIGAR-based precise projection
            project_cigar(block, region, cigar)
        } else {
            // Linear interpolation fallback
            project_linear(block, region)
        };

        hits.push(TranslationHit {
            asm_chrom: block.asm_chrom.clone(),
            asm_start,
            asm_end,
            strand: block.strand.clone(),
            event_type: "alignment".to_string(),
        });
    }

    // Detect gaps between consecutive alignment blocks (SV evidence)
    detect_sv_gaps(&mut hits, blocks, region);

    hits
}

/// CIGAR-based coordinate projection for precise reference→assembly mapping.
///
/// Walks the CIGAR string to translate reference positions to assembly
/// positions with base-level accuracy.
fn project_cigar(block: &AlignmentBlock, region: &GenomicRegion, cigar: &str) -> (u64, u64) {
    let ops = parse_cigar(cigar);
    if ops.is_empty() {
        return project_linear(block, region);
    }

    let is_reverse = block.strand == "-";

    // Clamp query range to the block's reference extent
    let qstart = region.start.max(block.ref_start);
    let qend = region.end.min(block.ref_end);

    let mut ref_pos = block.ref_start;
    let mut asm_pos = if is_reverse {
        block.asm_end
    } else {
        block.asm_start
    };

    let mut asm_start: Option<u64> = None;
    let mut asm_end: Option<u64> = None;

    for (len, op) in &ops {
        let len = *len as u64;
        match op {
            CigarOp::Match | CigarOp::SeqMatch | CigarOp::SeqMismatch => {
                // Consumes both ref and asm
                let ref_end_op = ref_pos + len;
                let overlap_start = qstart.max(ref_pos);
                let overlap_end = qend.min(ref_end_op);

                if overlap_start < overlap_end {
                    let offset_start = overlap_start - ref_pos;
                    let offset_end = overlap_end - ref_pos;

                    let (a_start, a_end) = if is_reverse {
                        (asm_pos - offset_end, asm_pos - offset_start)
                    } else {
                        (asm_pos + offset_start, asm_pos + offset_end)
                    };

                    asm_start = Some(asm_start.map_or(a_start, |s: u64| s.min(a_start)));
                    asm_end = Some(asm_end.map_or(a_end, |e: u64| e.max(a_end)));
                }

                ref_pos += len;
                if is_reverse {
                    asm_pos -= len;
                } else {
                    asm_pos += len;
                }
            }
            CigarOp::Insertion => {
                // Consumes only asm
                if is_reverse {
                    asm_pos -= len;
                } else {
                    asm_pos += len;
                }
            }
            CigarOp::Deletion | CigarOp::Skip => {
                // Consumes only ref
                ref_pos += len;
            }
            CigarOp::SoftClip | CigarOp::HardClip | CigarOp::Pad => {
                // Doesn't consume ref
                if matches!(op, CigarOp::SoftClip) {
                    if is_reverse {
                        asm_pos -= len;
                    } else {
                        asm_pos += len;
                    }
                }
            }
        }
    }

    match (asm_start, asm_end) {
        (Some(s), Some(e)) => (s, e),
        _ => project_linear(block, region),
    }
}

/// Linear interpolation fallback when no CIGAR is available.
fn project_linear(block: &AlignmentBlock, region: &GenomicRegion) -> (u64, u64) {
    let ref_span = (block.ref_end - block.ref_start).max(1) as f64;
    let asm_span = (block.asm_end - block.asm_start).max(1) as f64;
    let ratio = asm_span / ref_span;

    let qstart = region.start.max(block.ref_start);
    let qend = region.end.min(block.ref_end);

    let offset_start = (qstart - block.ref_start) as f64;
    let offset_end = (qend - block.ref_start) as f64;

    if block.strand == "-" {
        let asm_end = block.asm_end - (offset_start * ratio) as u64;
        let asm_start = block.asm_end - (offset_end * ratio) as u64;
        (asm_start, asm_end)
    } else {
        let asm_start = block.asm_start + (offset_start * ratio) as u64;
        let asm_end = block.asm_start + (offset_end * ratio) as u64;
        (asm_start, asm_end)
    }
}

/// Detect structural variant gaps between consecutive alignment blocks.
fn detect_sv_gaps(
    hits: &mut Vec<TranslationHit>,
    blocks: &[AlignmentBlock],
    region: &GenomicRegion,
) {
    // Find blocks overlapping the region (alignment-type hits only)
    let overlapping: Vec<&AlignmentBlock> = blocks
        .iter()
        .filter(|b| b.ref_start < region.end && b.ref_end > region.start)
        .collect();

    // Look at gaps between consecutive overlapping blocks
    for pair in overlapping.windows(2) {
        let b1 = pair[0];
        let b2 = pair[1];

        if b2.ref_start <= b1.ref_end {
            continue; // No gap
        }

        let gap_type = classify_gap(b1, b2);
        if gap_type != "none" {
            hits.push(TranslationHit {
                asm_chrom: b1.asm_chrom.clone(),
                asm_start: b1.asm_end,
                asm_end: b2.asm_start,
                strand: b1.strand.clone(),
                event_type: gap_type,
            });
        }
    }
}

/// Classify the type of structural variant between two alignment blocks.
fn classify_gap(b1: &AlignmentBlock, b2: &AlignmentBlock) -> String {
    if b1.asm_chrom != b2.asm_chrom {
        return "translocation".to_string();
    }

    if b1.strand != b2.strand {
        return "inversion".to_string();
    }

    let ref_gap = b2.ref_start.saturating_sub(b1.ref_end);
    let asm_gap = if b1.strand == "+" {
        b2.asm_start.saturating_sub(b1.asm_end)
    } else {
        b1.asm_start.saturating_sub(b2.asm_end)
    };

    if ref_gap > asm_gap + 50 {
        "deletion".to_string()
    } else if asm_gap > ref_gap + 50 {
        "insertion".to_string()
    } else {
        "none".to_string()
    }
}

/// CIGAR operation types.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum CigarOp {
    Match,        // M
    Insertion,    // I
    Deletion,     // D
    Skip,         // N
    SoftClip,     // S
    HardClip,     // H
    Pad,          // P
    SeqMatch,     // =
    SeqMismatch,  // X
}

/// Parse a CIGAR string into (length, operation) pairs.
fn parse_cigar(cigar: &str) -> Vec<(u32, CigarOp)> {
    let mut ops = Vec::new();
    let mut num = 0u32;

    for ch in cigar.chars() {
        if ch.is_ascii_digit() {
            num = num * 10 + ch.to_digit(10).unwrap();
        } else {
            let op = match ch {
                'M' => CigarOp::Match,
                'I' => CigarOp::Insertion,
                'D' => CigarOp::Deletion,
                'N' => CigarOp::Skip,
                'S' => CigarOp::SoftClip,
                'H' => CigarOp::HardClip,
                'P' => CigarOp::Pad,
                '=' => CigarOp::SeqMatch,
                'X' => CigarOp::SeqMismatch,
                _ => continue,
            };
            if num > 0 {
                ops.push((num, op));
            }
            num = 0;
        }
    }

    ops
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_cigar() {
        let ops = parse_cigar("100M5I50M3D20M");
        assert_eq!(ops.len(), 5);
        assert_eq!(ops[0], (100, CigarOp::Match));
        assert_eq!(ops[1], (5, CigarOp::Insertion));
        assert_eq!(ops[2], (50, CigarOp::Match));
        assert_eq!(ops[3], (3, CigarOp::Deletion));
        assert_eq!(ops[4], (20, CigarOp::Match));
    }

    #[test]
    fn test_parse_cigar_extended() {
        let ops = parse_cigar("50=3X20=5I10D");
        assert_eq!(ops.len(), 5);
        assert_eq!(ops[0], (50, CigarOp::SeqMatch));
        assert_eq!(ops[1], (3, CigarOp::SeqMismatch));
        assert_eq!(ops[2], (20, CigarOp::SeqMatch));
    }

    #[test]
    fn test_parse_cigar_empty() {
        let ops = parse_cigar("");
        assert!(ops.is_empty());
    }

    #[test]
    fn test_project_linear_forward() {
        let block = AlignmentBlock {
            ref_start: 1000,
            ref_end: 2000,
            asm_chrom: "ctg1".to_string(),
            asm_start: 5000,
            asm_end: 6000,
            strand: "+".to_string(),
            mapq: 60,
            cigar: None,
            aln_type: None,
        };

        let region = GenomicRegion::new("chr1", 1200, 1800);
        let (start, end) = project_linear(&block, &region);
        assert_eq!(start, 5200);
        assert_eq!(end, 5800);
    }

    #[test]
    fn test_project_linear_reverse() {
        let block = AlignmentBlock {
            ref_start: 1000,
            ref_end: 2000,
            asm_chrom: "ctg1".to_string(),
            asm_start: 5000,
            asm_end: 6000,
            strand: "-".to_string(),
            mapq: 60,
            cigar: None,
            aln_type: None,
        };

        let region = GenomicRegion::new("chr1", 1200, 1800);
        let (start, end) = project_linear(&block, &region);
        assert_eq!(start, 5200);
        assert_eq!(end, 5800);
    }

    #[test]
    fn test_project_cigar_simple_match() {
        let block = AlignmentBlock {
            ref_start: 1000,
            ref_end: 1100,
            asm_chrom: "ctg1".to_string(),
            asm_start: 5000,
            asm_end: 5100,
            strand: "+".to_string(),
            mapq: 60,
            cigar: Some("100M".to_string()),
            aln_type: None,
        };

        let region = GenomicRegion::new("chr1", 1020, 1080);
        let (start, end) = project_cigar(&block, &region, "100M");
        assert_eq!(start, 5020);
        assert_eq!(end, 5080);
    }

    #[test]
    fn test_project_cigar_with_insertion() {
        // 50M10I50M → ref consumes 100bp, asm consumes 110bp
        let block = AlignmentBlock {
            ref_start: 1000,
            ref_end: 1100,
            asm_chrom: "ctg1".to_string(),
            asm_start: 5000,
            asm_end: 5110,
            strand: "+".to_string(),
            mapq: 60,
            cigar: Some("50M10I50M".to_string()),
            aln_type: None,
        };

        // Query region that spans the insertion
        let region = GenomicRegion::new("chr1", 1020, 1080);
        let (start, end) = project_cigar(&block, &region, "50M10I50M");
        // First 50 M: ref 1000-1050 → asm 5000-5050
        // 10 I: asm 5050-5060
        // Second 50 M: ref 1050-1100 → asm 5060-5110
        // Query 1020: within first match → asm 5020
        // Query 1080: within second match (ref_pos 1080 is at offset 30 in second match)
        //   → asm 5060 + 30 = 5090
        assert_eq!(start, 5020);
        assert_eq!(end, 5090);
    }

    #[test]
    fn test_classify_gap_deletion() {
        let b1 = AlignmentBlock {
            ref_start: 1000,
            ref_end: 2000,
            asm_chrom: "ctg1".to_string(),
            asm_start: 5000,
            asm_end: 6000,
            strand: "+".to_string(),
            mapq: 60,
            cigar: None,
            aln_type: None,
        };
        let b2 = AlignmentBlock {
            ref_start: 3000,
            ref_end: 4000,
            asm_chrom: "ctg1".to_string(),
            asm_start: 6100,
            asm_end: 7100,
            strand: "+".to_string(),
            mapq: 60,
            cigar: None,
            aln_type: None,
        };

        // ref gap = 1000, asm gap = 100 → deletion
        assert_eq!(classify_gap(&b1, &b2), "deletion");
    }

    #[test]
    fn test_classify_gap_insertion() {
        let b1 = AlignmentBlock {
            ref_start: 1000,
            ref_end: 2000,
            asm_chrom: "ctg1".to_string(),
            asm_start: 5000,
            asm_end: 6000,
            strand: "+".to_string(),
            mapq: 60,
            cigar: None,
            aln_type: None,
        };
        let b2 = AlignmentBlock {
            ref_start: 2100,
            ref_end: 3100,
            asm_chrom: "ctg1".to_string(),
            asm_start: 7000,
            asm_end: 8000,
            strand: "+".to_string(),
            mapq: 60,
            cigar: None,
            aln_type: None,
        };

        // ref gap = 100, asm gap = 1000 → insertion
        assert_eq!(classify_gap(&b1, &b2), "insertion");
    }

    #[test]
    fn test_classify_gap_inversion() {
        let b1 = AlignmentBlock {
            ref_start: 1000,
            ref_end: 2000,
            asm_chrom: "ctg1".to_string(),
            asm_start: 5000,
            asm_end: 6000,
            strand: "+".to_string(),
            mapq: 60,
            cigar: None,
            aln_type: None,
        };
        let b2 = AlignmentBlock {
            ref_start: 3000,
            ref_end: 4000,
            asm_chrom: "ctg1".to_string(),
            asm_start: 7000,
            asm_end: 8000,
            strand: "-".to_string(), // Different strand
            mapq: 60,
            cigar: None,
            aln_type: None,
        };

        assert_eq!(classify_gap(&b1, &b2), "inversion");
    }

    #[test]
    fn test_classify_gap_translocation() {
        let b1 = AlignmentBlock {
            ref_start: 1000,
            ref_end: 2000,
            asm_chrom: "ctg1".to_string(),
            asm_start: 5000,
            asm_end: 6000,
            strand: "+".to_string(),
            mapq: 60,
            cigar: None,
            aln_type: None,
        };
        let b2 = AlignmentBlock {
            ref_start: 3000,
            ref_end: 4000,
            asm_chrom: "ctg2".to_string(), // Different contig
            asm_start: 1000,
            asm_end: 2000,
            strand: "+".to_string(),
            mapq: 60,
            cigar: None,
            aln_type: None,
        };

        assert_eq!(classify_gap(&b1, &b2), "translocation");
    }

    #[test]
    fn test_query_index_no_overlap() {
        let mut blocks = HashMap::new();
        blocks.insert(
            "chr1".to_string(),
            vec![AlignmentBlock {
                ref_start: 1000,
                ref_end: 2000,
                asm_chrom: "ctg1".to_string(),
                asm_start: 5000,
                asm_end: 6000,
                strand: "+".to_string(),
                mapq: 60,
                cigar: None,
                aln_type: None,
            }],
        );
        let index = MappingIndex { blocks };

        // Query region that doesn't overlap
        let region = GenomicRegion::new("chr1", 3000, 4000);
        let hits = query_index(&index, &region, 0);
        assert!(hits.is_empty());
    }

    #[test]
    fn test_query_index_overlap() {
        let mut blocks = HashMap::new();
        blocks.insert(
            "chr1".to_string(),
            vec![AlignmentBlock {
                ref_start: 1000,
                ref_end: 2000,
                asm_chrom: "ctg1".to_string(),
                asm_start: 5000,
                asm_end: 6000,
                strand: "+".to_string(),
                mapq: 60,
                cigar: None,
                aln_type: None,
            }],
        );
        let index = MappingIndex { blocks };

        let region = GenomicRegion::new("chr1", 1500, 1800);
        let hits = query_index(&index, &region, 0);
        assert_eq!(hits.len(), 1);
        assert_eq!(hits[0].event_type, "alignment");
    }

    #[test]
    fn test_query_index_mapq_filter() {
        let mut blocks = HashMap::new();
        blocks.insert(
            "chr1".to_string(),
            vec![AlignmentBlock {
                ref_start: 1000,
                ref_end: 2000,
                asm_chrom: "ctg1".to_string(),
                asm_start: 5000,
                asm_end: 6000,
                strand: "+".to_string(),
                mapq: 10,
                cigar: None,
                aln_type: None,
            }],
        );
        let index = MappingIndex { blocks };

        let region = GenomicRegion::new("chr1", 1500, 1800);
        // Should be filtered out by mapq threshold
        let hits = query_index(&index, &region, 20);
        assert!(hits.is_empty());
    }

    #[test]
    fn test_coordinator_translator_basic() {
        let translator = CoordinateTranslator::new();
        assert!(!translator.has_index("sample1", "hap1"));

        let result = translator.translate(
            "sample1",
            &GenomicRegion::new("chr1", 1000, 2000),
            0,
        );
        assert!(result.hap1.is_empty());
        assert!(result.hap2.is_empty());
    }
}
