//! CIGAR-aware coordinate translation between reference genome and haplotype
//! assembly coordinates.
//!
//! This module loads a gzipped JSON mapping index (produced by the Python
//! `coordinate_mapper.py build` command) and provides efficient interval-based
//! lookup that translates reference coordinates to assembly coordinates.
//! Structural-variation aware: when a CIGAR string is present, coordinate
//! projection walks the CIGAR base-by-base rather than using linear
//! interpolation.  Gaps between alignment blocks are detected and classified
//! as deletions, insertions, inversions, or translocations.

use flate2::read::GzDecoder;
use serde::Deserialize;
use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufReader, Read};

// ---------------------------------------------------------------------------
// Types
// ---------------------------------------------------------------------------

/// A single CIGAR operation (length + operation character).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CigarOp {
    pub len: u64,
    pub op: char,
}

/// Pre-computed cumulative reference and assembly offsets for binary search
/// into a CIGAR string.
#[derive(Debug, Clone)]
pub struct CigarIndex {
    pub ref_cumul: Vec<u64>,
    pub asm_cumul: Vec<u64>,
}

/// An alignment block as stored in the JSON index.
#[derive(Debug, Clone, Deserialize)]
pub struct AlignmentBlock {
    /// Reference start (0-based, inclusive).
    pub rs: u64,
    /// Reference end (0-based, exclusive).
    pub re: u64,
    /// Assembly contig name.
    pub ac: String,
    /// Assembly start (0-based, inclusive).  JSON key is `"as"`.
    #[serde(rename = "as")]
    pub as_: u64,
    /// Assembly end (0-based, exclusive).
    pub ae: u64,
    /// Strand (`+` or `-`).
    pub st: String,
    /// Mapping quality.
    pub mq: u64,
    /// Optional CIGAR string.
    #[serde(default)]
    pub cg: Option<String>,
}

/// Per-chromosome data: sorted blocks with auxiliary arrays for binary search.
#[derive(Debug, Clone)]
pub struct ChromData {
    pub blocks: Vec<AlignmentBlock>,
    pub starts: Vec<u64>,
    pub ends: Vec<u64>,
    pub max_block_len: u64,
    /// Pre-parsed CIGAR ops per block (empty vec if no CIGAR).
    cigar_ops: Vec<Vec<CigarOp>>,
    /// Pre-built CIGAR index per block (None for short CIGARs).
    cigar_indices: Vec<Option<CigarIndex>>,
}

/// The full mapping index keyed by reference chromosome.
pub type MappingIndex = HashMap<String, ChromData>;

/// Classification of a structural-variation event.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum EventType {
    Alignment,
    Deletion,
    Insertion,
    Inversion,
    Translocation,
    Complex,
}

impl std::fmt::Display for EventType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            EventType::Alignment => write!(f, "alignment"),
            EventType::Deletion => write!(f, "deletion"),
            EventType::Insertion => write!(f, "insertion"),
            EventType::Inversion => write!(f, "inversion"),
            EventType::Translocation => write!(f, "translocation"),
            EventType::Complex => write!(f, "complex"),
        }
    }
}

/// Result of a coordinate query.
#[derive(Debug, Clone)]
pub struct QueryResult {
    pub ref_chrom: String,
    pub ref_start: u64,
    pub ref_end: u64,
    pub asm_chrom: String,
    pub asm_start: u64,
    pub asm_end: u64,
    pub strand: String,
    pub mapq: u64,
    pub event_type: EventType,
    /// Present only for gap events.
    pub ref_gap_size: Option<u64>,
    /// Present only for gap events. `None` for translocations and inversions.
    pub asm_gap_size: Option<i64>,
}

// ---------------------------------------------------------------------------
// Constants
// ---------------------------------------------------------------------------

/// Minimum number of CIGAR operations before building a binary-search index.
const CIGAR_INDEX_THRESHOLD: usize = 64;

// ---------------------------------------------------------------------------
// CIGAR helpers
// ---------------------------------------------------------------------------

#[inline]
fn is_ref_consuming(op: char) -> bool {
    matches!(op, 'M' | 'D' | 'N' | 'X' | '=')
}

#[inline]
fn is_asm_consuming(op: char) -> bool {
    matches!(op, 'M' | 'I' | 'S' | 'X' | '=')
}

/// Parse a CIGAR string into a list of [`CigarOp`]s.
pub fn parse_cigar(cigar_str: &str) -> Vec<CigarOp> {
    if cigar_str.is_empty() {
        return Vec::new();
    }
    let mut ops = Vec::new();
    let mut num_start: Option<usize> = None;
    for (i, c) in cigar_str.char_indices() {
        if c.is_ascii_digit() {
            if num_start.is_none() {
                num_start = Some(i);
            }
        } else if matches!(c, 'M' | 'I' | 'D' | 'N' | 'S' | 'H' | 'P' | 'X' | '=') {
            if let Some(start) = num_start {
                if let Ok(len) = cigar_str[start..i].parse::<u64>() {
                    ops.push(CigarOp { len, op: c });
                }
            }
            num_start = None;
        }
    }
    ops
}

/// Build a cumulative-offset index for binary search into a CIGAR.
pub fn build_cigar_index(ops: &[CigarOp]) -> CigarIndex {
    let mut ref_cumul = Vec::with_capacity(ops.len() + 1);
    let mut asm_cumul = Vec::with_capacity(ops.len() + 1);
    let mut ref_cur: u64 = 0;
    let mut asm_cur: u64 = 0;
    ref_cumul.push(0);
    asm_cumul.push(0);
    for op in ops {
        if is_ref_consuming(op.op) {
            ref_cur += op.len;
        }
        if is_asm_consuming(op.op) {
            asm_cur += op.len;
        }
        ref_cumul.push(ref_cur);
        asm_cumul.push(asm_cur);
    }
    CigarIndex {
        ref_cumul,
        asm_cumul,
    }
}

/// Walk a CIGAR to project a reference sub-interval to assembly coordinates.
///
/// All coordinates are 0-based, half-open `[start, end)`.
pub fn project_cigar(
    ops: &[CigarOp],
    ref_offset_start: u64,
    ref_offset_end: u64,
    strand: &str,
    asm_block_start: u64,
    asm_block_end: u64,
    cigar_index: Option<&CigarIndex>,
) -> (u64, u64) {
    let mut start_op: usize = 0;
    let mut ref_cur: u64 = 0;
    let mut asm_cur: u64 = 0;

    // Fast path: binary-search into pre-computed cumulative offsets.
    if let Some(idx) = cigar_index {
        if ref_offset_start > 0 {
            let pos = idx
                .ref_cumul
                .partition_point(|&v| v <= ref_offset_start);
            let pos = if pos > 0 { pos - 1 } else { 0 };
            let pos = pos.min(ops.len() - 1);
            start_op = pos;
            ref_cur = idx.ref_cumul[pos];
            asm_cur = idx.asm_cumul[pos];
        }
    }

    let mut asm_q_start: Option<u64> = None;
    let mut asm_q_end: Option<u64> = None;

    for i in start_op..ops.len() {
        let length = ops[i].len;
        let op = ops[i].op;
        let consumes_ref = is_ref_consuming(op);
        let consumes_asm = is_asm_consuming(op);

        if consumes_ref {
            // Does ref_offset_start fall inside this operation?
            if asm_q_start.is_none() && ref_cur + length > ref_offset_start {
                let inner = ref_offset_start - ref_cur;
                asm_q_start = Some(asm_cur + if consumes_asm { inner } else { 0 });
            }

            // Does ref_offset_end fall inside (or at the end of) this op?
            if asm_q_start.is_some() && asm_q_end.is_none() && ref_cur + length >= ref_offset_end {
                let inner = ref_offset_end - ref_cur;
                asm_q_end = Some(asm_cur + if consumes_asm { inner } else { 0 });
                break;
            }
        }

        if consumes_ref {
            ref_cur += length;
        }
        if consumes_asm {
            asm_cur += length;
        }
    }

    // Fallback: clamp to last position reached.
    let asm_q_start = asm_q_start.unwrap_or(asm_cur);
    let asm_q_end = asm_q_end.unwrap_or(asm_cur);

    if strand == "+" {
        (
            asm_block_start + asm_q_start,
            asm_block_start + asm_q_end,
        )
    } else {
        (
            asm_block_end - asm_q_end,
            asm_block_end - asm_q_start,
        )
    }
}

// ---------------------------------------------------------------------------
// SV gap classification
// ---------------------------------------------------------------------------

/// Classify the structural-variation event represented by a gap between blocks.
pub fn classify_sv_gap(
    prev_block: &AlignmentBlock,
    next_block: &AlignmentBlock,
    gap_ref_start: u64,
    gap_ref_end: u64,
    chrom: &str,
) -> QueryResult {
    let ref_gap_size = gap_ref_end - gap_ref_start;
    let same_contig = prev_block.ac == next_block.ac;
    let same_strand = prev_block.st == next_block.st;

    if !same_contig {
        let asm_pos = if prev_block.st == "+" {
            prev_block.ae
        } else {
            prev_block.as_
        };
        return QueryResult {
            ref_chrom: chrom.to_string(),
            ref_start: gap_ref_start,
            ref_end: gap_ref_end,
            asm_chrom: prev_block.ac.clone(),
            asm_start: asm_pos,
            asm_end: asm_pos,
            strand: prev_block.st.clone(),
            mapq: prev_block.mq.min(next_block.mq),
            event_type: EventType::Translocation,
            ref_gap_size: Some(ref_gap_size),
            asm_gap_size: None,
        };
    }

    if !same_strand {
        let asm_pos = if prev_block.st == "+" {
            prev_block.ae
        } else {
            prev_block.as_
        };
        return QueryResult {
            ref_chrom: chrom.to_string(),
            ref_start: gap_ref_start,
            ref_end: gap_ref_end,
            asm_chrom: prev_block.ac.clone(),
            asm_start: asm_pos,
            asm_end: asm_pos,
            strand: prev_block.st.clone(),
            mapq: prev_block.mq.min(next_block.mq),
            event_type: EventType::Inversion,
            ref_gap_size: Some(ref_gap_size),
            asm_gap_size: None,
        };
    }

    // Same contig, same strand.
    let (asm_gap_size, asm_start, asm_end) = if prev_block.st == "+" {
        let gap = next_block.as_ as i64 - prev_block.ae as i64;
        (gap, prev_block.ae, next_block.as_)
    } else {
        let gap = prev_block.as_ as i64 - next_block.ae as i64;
        (gap, next_block.ae, prev_block.as_)
    };

    let event_type = if asm_gap_size < 0 {
        EventType::Complex
    } else if (asm_gap_size as u64) < ref_gap_size {
        EventType::Deletion
    } else {
        EventType::Insertion
    };

    QueryResult {
        ref_chrom: chrom.to_string(),
        ref_start: gap_ref_start,
        ref_end: gap_ref_end,
        asm_chrom: prev_block.ac.clone(),
        asm_start,
        asm_end,
        strand: prev_block.st.clone(),
        mapq: prev_block.mq.min(next_block.mq),
        event_type,
        ref_gap_size: Some(ref_gap_size),
        asm_gap_size: Some(asm_gap_size),
    }
}

// ---------------------------------------------------------------------------
// Index I/O
// ---------------------------------------------------------------------------

/// Load a gzipped JSON mapping index from disk.
///
/// Pre-parses CIGARs and builds cumulative-offset indices for large CIGARs.
pub fn load_index(path: &str) -> Result<MappingIndex, io::Error> {
    let file = File::open(path)?;
    let reader = BufReader::new(GzDecoder::new(file));
    let mut buf = String::new();
    { let mut r = reader; r.read_to_string(&mut buf)?; }
    let data: HashMap<String, Vec<AlignmentBlock>> =
        serde_json::from_str(&buf).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

    let mut index = MappingIndex::new();
    for (chrom, blocks) in data {
        let starts: Vec<u64> = blocks.iter().map(|b| b.rs).collect();
        let ends: Vec<u64> = blocks.iter().map(|b| b.re).collect();
        let max_block_len = starts
            .iter()
            .zip(ends.iter())
            .map(|(s, e)| e - s)
            .max()
            .unwrap_or(0);

        let mut cigar_ops_vec = Vec::with_capacity(blocks.len());
        let mut cigar_idx_vec = Vec::with_capacity(blocks.len());

        for block in &blocks {
            if let Some(cg) = &block.cg {
                let ops = parse_cigar(cg);
                let idx = if ops.len() >= CIGAR_INDEX_THRESHOLD {
                    Some(build_cigar_index(&ops))
                } else {
                    None
                };
                cigar_ops_vec.push(ops);
                cigar_idx_vec.push(idx);
            } else {
                cigar_ops_vec.push(Vec::new());
                cigar_idx_vec.push(None);
            }
        }

        index.insert(
            chrom,
            ChromData {
                blocks,
                starts,
                ends,
                max_block_len,
                cigar_ops: cigar_ops_vec,
                cigar_indices: cigar_idx_vec,
            },
        );
    }
    Ok(index)
}

// ---------------------------------------------------------------------------
// Query
// ---------------------------------------------------------------------------

/// Return `true` if there is a structural event at the junction between
/// ref-adjacent blocks (different contigs, strands, or non-zero asm gap).
fn has_sv_at_junction(prev_block: &AlignmentBlock, block: &AlignmentBlock) -> bool {
    if prev_block.ac != block.ac {
        return true;
    }
    if prev_block.st != block.st {
        return true;
    }
    let asm_gap = if prev_block.st == "+" {
        block.as_ as i64 - prev_block.ae as i64
    } else {
        prev_block.as_ as i64 - block.ae as i64
    };
    asm_gap != 0
}

/// Find assembly regions and SV events overlapping `[start, end)` on `chrom`.
pub fn query(
    index: &MappingIndex,
    chrom: &str,
    start: u64,
    end: u64,
    min_mapq: u64,
) -> Vec<QueryResult> {
    let chrom_data = match index.get(chrom) {
        Some(cd) => cd,
        None => return Vec::new(),
    };

    let starts = &chrom_data.starts;
    let ends = &chrom_data.ends;
    let blocks = &chrom_data.blocks;

    // Binary search: block_start < end
    let right_idx = starts.partition_point(|&s| s < end);

    // Use max_block_len to skip blocks that end too far left.
    let max_len = chrom_data.max_block_len;
    let left_idx = if max_len > 0 && start >= max_len {
        starts.partition_point(|&s| s < start - max_len)
    } else {
        0
    };

    // Collect overlapping blocks (indices).
    let mut overlapping: Vec<usize> = Vec::new();
    for i in left_idx..right_idx {
        if ends[i] > start && blocks[i].mq >= min_mapq {
            overlapping.push(i);
        }
    }

    if overlapping.is_empty() {
        return Vec::new();
    }

    let mut results = Vec::new();
    let mut prev_block_idx: Option<usize> = None;

    for &block_i in &overlapping {
        let block = &blocks[block_i];

        // Detect and classify gap between previous block and this one.
        if let Some(prev_i) = prev_block_idx {
            let prev = &blocks[prev_i];
            let junction = prev.re;
            let next_start = block.rs;

            let gap_start = start.max(junction);
            let gap_end = end.min(next_start);

            let has_ref_gap = gap_start < gap_end;
            let has_junction_event = start <= junction
                && junction <= end
                && junction == next_start
                && has_sv_at_junction(prev, block);

            if has_ref_gap || has_junction_event {
                results.push(classify_sv_gap(prev, block, gap_start, gap_end, chrom));
            }
        }

        // Intersect query with this block.
        let overlap_start = start.max(block.rs);
        let overlap_end = end.min(block.re);

        if overlap_start >= overlap_end {
            // No actual overlap; still track high-water mark.
            if prev_block_idx.is_none()
                || block.re > blocks[prev_block_idx.unwrap()].re
            {
                prev_block_idx = Some(block_i);
            }
            continue;
        }

        let ref_offset_start = overlap_start - block.rs;
        let ref_offset_end = overlap_end - block.rs;

        // Prefer CIGAR projection; fall back to linear interpolation.
        let ops = &chrom_data.cigar_ops[block_i];
        let (asm_start, asm_end) = if !ops.is_empty() {
            let cigar_idx = chrom_data.cigar_indices[block_i].as_ref();
            project_cigar(
                ops,
                ref_offset_start,
                ref_offset_end,
                &block.st,
                block.as_,
                block.ae,
                cigar_idx,
            )
        } else {
            if block.st == "+" {
                (block.as_ + ref_offset_start, block.as_ + ref_offset_end)
            } else {
                (block.ae - ref_offset_end, block.ae - ref_offset_start)
            }
        };

        results.push(QueryResult {
            ref_chrom: chrom.to_string(),
            ref_start: overlap_start,
            ref_end: overlap_end,
            asm_chrom: block.ac.clone(),
            asm_start,
            asm_end,
            strand: block.st.clone(),
            mapq: block.mq,
            event_type: EventType::Alignment,
            ref_gap_size: None,
            asm_gap_size: None,
        });

        // High-water mark.
        if prev_block_idx.is_none()
            || block.re > blocks[prev_block_idx.unwrap()].re
        {
            prev_block_idx = Some(block_i);
        }
    }

    results
}

// ===========================================================================
// Tests
// ===========================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use flate2::write::GzEncoder;
    use flate2::Compression;
    use std::io::Write;

    // Helper: write a gzipped JSON index to a temp file, load it.
    fn write_gz_json(dir: &tempfile::TempDir, data: &serde_json::Value) -> String {
        let path = dir.path().join("test.mapping.json.gz");
        let file = File::create(&path).unwrap();
        let mut gz = GzEncoder::new(file, Compression::default());
        gz.write_all(serde_json::to_string(data).unwrap().as_bytes())
            .unwrap();
        gz.finish().unwrap();
        path.to_str().unwrap().to_string()
    }

    // Helper: build a JSON value from block specs.
    fn make_index_json(
        blocks: &[(&str, u64, u64, &str, u64, u64, &str, u64, Option<&str>)],
    ) -> serde_json::Value {
        let mut map: HashMap<String, Vec<serde_json::Value>> = HashMap::new();
        for &(chrom, rs, re, ac, as_, ae, st, mq, cg) in blocks {
            let mut obj = serde_json::json!({
                "rs": rs, "re": re, "ac": ac, "as": as_, "ae": ae,
                "st": st, "mq": mq,
            });
            if let Some(c) = cg {
                obj["cg"] = serde_json::json!(c);
            }
            map.entry(chrom.to_string()).or_default().push(obj);
        }
        serde_json::json!(map)
    }

    // -----------------------------------------------------------------------
    // CIGAR parsing
    // -----------------------------------------------------------------------

    #[test]
    fn test_parse_cigar_empty() {
        assert!(parse_cigar("").is_empty());
    }

    #[test]
    fn test_parse_cigar_simple_match() {
        let ops = parse_cigar("100M");
        assert_eq!(ops, vec![CigarOp { len: 100, op: 'M' }]);
    }

    #[test]
    fn test_parse_cigar_complex() {
        let ops = parse_cigar("50M5D30M10I20M");
        assert_eq!(
            ops,
            vec![
                CigarOp { len: 50, op: 'M' },
                CigarOp { len: 5, op: 'D' },
                CigarOp { len: 30, op: 'M' },
                CigarOp { len: 10, op: 'I' },
                CigarOp { len: 20, op: 'M' },
            ]
        );
    }

    #[test]
    fn test_parse_cigar_extended() {
        let ops = parse_cigar("10=5X3=");
        assert_eq!(
            ops,
            vec![
                CigarOp { len: 10, op: '=' },
                CigarOp { len: 5, op: 'X' },
                CigarOp { len: 3, op: '=' },
            ]
        );
    }

    #[test]
    fn test_parse_cigar_soft_hard_clip() {
        let ops = parse_cigar("5H10S90M5S3H");
        assert_eq!(ops.len(), 5);
        assert_eq!(ops[0], CigarOp { len: 5, op: 'H' });
        assert_eq!(ops[1], CigarOp { len: 10, op: 'S' });
        assert_eq!(ops[4], CigarOp { len: 3, op: 'H' });
    }

    // -----------------------------------------------------------------------
    // build_cigar_index
    // -----------------------------------------------------------------------

    #[test]
    fn test_build_cigar_index_match_only() {
        let ops = parse_cigar("100M");
        let idx = build_cigar_index(&ops);
        assert_eq!(idx.ref_cumul, vec![0, 100]);
        assert_eq!(idx.asm_cumul, vec![0, 100]);
    }

    #[test]
    fn test_build_cigar_index_with_deletion() {
        let ops = parse_cigar("50M10D40M");
        let idx = build_cigar_index(&ops);
        assert_eq!(idx.ref_cumul, vec![0, 50, 60, 100]);
        assert_eq!(idx.asm_cumul, vec![0, 50, 50, 90]);
    }

    #[test]
    fn test_build_cigar_index_with_insertion() {
        let ops = parse_cigar("50M5I40M");
        let idx = build_cigar_index(&ops);
        assert_eq!(idx.ref_cumul, vec![0, 50, 50, 90]);
        assert_eq!(idx.asm_cumul, vec![0, 50, 55, 95]);
    }

    #[test]
    fn test_build_cigar_index_complex() {
        let ops = parse_cigar("100M5D200M3I50M");
        let idx = build_cigar_index(&ops);
        assert_eq!(idx.ref_cumul, vec![0, 100, 105, 305, 305, 355]);
        assert_eq!(idx.asm_cumul, vec![0, 100, 100, 300, 303, 353]);
    }

    // -----------------------------------------------------------------------
    // project_cigar: pure match + strand
    // -----------------------------------------------------------------------

    #[test]
    fn test_project_pure_match_plus_full() {
        let ops = parse_cigar("100M");
        let (s, e) = project_cigar(&ops, 0, 100, "+", 1000, 1100, None);
        assert_eq!((s, e), (1000, 1100));
    }

    #[test]
    fn test_project_pure_match_plus_partial() {
        let ops = parse_cigar("100M");
        let (s, e) = project_cigar(&ops, 10, 40, "+", 1000, 1100, None);
        assert_eq!((s, e), (1010, 1040));
    }

    // -----------------------------------------------------------------------
    // project_cigar: deletion on + strand
    // -----------------------------------------------------------------------

    #[test]
    fn test_deletion_query_before() {
        let ops = parse_cigar("100M5D95M");
        let (s, e) = project_cigar(&ops, 0, 50, "+", 1000, 1195, None);
        assert_eq!((s, e), (1000, 1050));
    }

    #[test]
    fn test_deletion_query_spanning() {
        let ops = parse_cigar("100M5D95M");
        let (s, e) = project_cigar(&ops, 90, 110, "+", 1000, 1195, None);
        assert_eq!((s, e), (1090, 1105));
    }

    #[test]
    fn test_deletion_query_after() {
        let ops = parse_cigar("100M5D95M");
        let (s, e) = project_cigar(&ops, 110, 150, "+", 1000, 1195, None);
        assert_eq!((s, e), (1105, 1145));
    }

    // -----------------------------------------------------------------------
    // project_cigar: insertion on + strand
    // -----------------------------------------------------------------------

    #[test]
    fn test_insertion_query_before() {
        let ops = parse_cigar("100M10I100M");
        let (s, e) = project_cigar(&ops, 0, 50, "+", 2000, 2210, None);
        assert_eq!((s, e), (2000, 2050));
    }

    #[test]
    fn test_insertion_query_straddling() {
        let ops = parse_cigar("100M10I100M");
        let (s, e) = project_cigar(&ops, 90, 110, "+", 2000, 2210, None);
        assert_eq!((s, e), (2090, 2120));
    }

    #[test]
    fn test_insertion_query_after() {
        let ops = parse_cigar("100M10I100M");
        let (s, e) = project_cigar(&ops, 110, 150, "+", 2000, 2210, None);
        assert_eq!((s, e), (2120, 2160));
    }

    // -----------------------------------------------------------------------
    // project_cigar: minus strand
    // -----------------------------------------------------------------------

    #[test]
    fn test_minus_strand_full() {
        let ops = parse_cigar("100M5D95M");
        let (s, e) = project_cigar(&ops, 0, 200, "-", 3000, 3195, None);
        assert_eq!((s, e), (3000, 3195));
    }

    #[test]
    fn test_minus_strand_before_deletion() {
        let ops = parse_cigar("100M5D95M");
        let (s, e) = project_cigar(&ops, 0, 50, "-", 3000, 3195, None);
        assert_eq!((s, e), (3145, 3195));
    }

    #[test]
    fn test_minus_strand_spanning_deletion() {
        let ops = parse_cigar("100M5D95M");
        let (s, e) = project_cigar(&ops, 90, 110, "-", 3000, 3195, None);
        assert_eq!((s, e), (3090, 3105));
    }

    // -----------------------------------------------------------------------
    // SV gap classification (unit)
    // -----------------------------------------------------------------------

    fn make_block(
        ac: &str,
        as_: u64,
        ae: u64,
        st: &str,
        mq: u64,
        rs: u64,
        re: u64,
    ) -> AlignmentBlock {
        AlignmentBlock {
            rs, re, ac: ac.to_string(), as_, ae,
            st: st.to_string(), mq, cg: None,
        }
    }

    #[test]
    fn test_classify_deletion_plus() {
        let prev = make_block("ctg1", 0, 1000, "+", 60, 0, 1000);
        let nxt = make_block("ctg1", 1000, 2000, "+", 60, 0, 1000);
        let r = classify_sv_gap(&prev, &nxt, 1000, 6000, "chr1");
        assert_eq!(r.event_type, EventType::Deletion);
        assert_eq!(r.ref_gap_size, Some(5000));
        assert_eq!(r.asm_gap_size, Some(0));
    }

    #[test]
    fn test_classify_insertion_plus() {
        let prev = make_block("ctg1", 0, 1000, "+", 60, 0, 1000);
        let nxt = make_block("ctg1", 6000, 7000, "+", 60, 0, 1000);
        let r = classify_sv_gap(&prev, &nxt, 1000, 1000, "chr1");
        assert_eq!(r.event_type, EventType::Insertion);
        assert_eq!(r.ref_gap_size, Some(0));
        assert_eq!(r.asm_gap_size, Some(5000));
    }

    #[test]
    fn test_classify_complex_plus() {
        let prev = make_block("ctg1", 0, 5000, "+", 60, 0, 1000);
        let nxt = make_block("ctg1", 3000, 6000, "+", 60, 0, 1000);
        let r = classify_sv_gap(&prev, &nxt, 1000, 2000, "chr1");
        assert_eq!(r.event_type, EventType::Complex);
        assert_eq!(r.asm_gap_size, Some(-2000));
    }

    #[test]
    fn test_classify_inversion() {
        let prev = make_block("ctg1", 0, 1000, "+", 60, 0, 1000);
        let nxt = make_block("ctg1", 0, 1000, "-", 60, 0, 1000);
        let r = classify_sv_gap(&prev, &nxt, 1000, 2000, "chr1");
        assert_eq!(r.event_type, EventType::Inversion);
        assert!(r.asm_gap_size.is_none());
    }

    #[test]
    fn test_classify_translocation() {
        let prev = make_block("ctgA", 0, 1000, "+", 60, 0, 1000);
        let nxt = make_block("ctgB", 0, 1000, "+", 60, 0, 1000);
        let r = classify_sv_gap(&prev, &nxt, 1000, 2000, "chr1");
        assert_eq!(r.event_type, EventType::Translocation);
        assert!(r.asm_gap_size.is_none());
    }

    #[test]
    fn test_classify_deletion_minus() {
        let prev = make_block("ctg1", 5000, 8000, "-", 60, 0, 8000);
        let nxt = make_block("ctg1", 3000, 4000, "-", 60, 0, 8000);
        let r = classify_sv_gap(&prev, &nxt, 8000, 13000, "chr1");
        assert_eq!(r.event_type, EventType::Deletion);
        assert_eq!(r.ref_gap_size, Some(5000));
        assert_eq!(r.asm_gap_size, Some(1000));
    }

    #[test]
    fn test_classify_insertion_minus() {
        let prev = make_block("ctg1", 10000, 13000, "-", 60, 0, 13000);
        let nxt = make_block("ctg1", 3000, 6000, "-", 60, 0, 13000);
        let r = classify_sv_gap(&prev, &nxt, 13000, 13000, "chr1");
        assert_eq!(r.event_type, EventType::Insertion);
        assert_eq!(r.asm_gap_size, Some(4000));
    }

    #[test]
    fn test_classify_mapq_is_minimum() {
        let prev = make_block("ctg1", 0, 1000, "+", 60, 0, 1000);
        let nxt = make_block("ctg1", 1000, 2000, "+", 20, 0, 1000);
        let r = classify_sv_gap(&prev, &nxt, 1000, 2000, "chr1");
        assert_eq!(r.mapq, 20);
    }

    // -----------------------------------------------------------------------
    // Full query round-trip with CIGAR index loaded from JSON
    // -----------------------------------------------------------------------

    #[test]
    fn test_round_trip_plus_strand_deletion_cigar() {
        // Block: ref chr1 [10000, 10200), asm_chr1 [1000, 1195), +
        // CIGAR: 100M5D95M
        let dir = tempfile::tempdir().unwrap();
        let data = make_index_json(&[(
            "chr1", 10000, 10200, "asm_chr1", 1000, 1195, "+", 60,
            Some("100M5D95M"),
        )]);
        let path = write_gz_json(&dir, &data);
        let index = load_index(&path).unwrap();

        // Query spanning deletion: ref [10090, 10110)
        let results = query(&index, "chr1", 10090, 10110, 0);
        assert_eq!(results.len(), 1);
        assert_eq!(results[0].event_type, EventType::Alignment);
        assert_eq!(results[0].asm_start, 1090);
        assert_eq!(results[0].asm_end, 1105);
    }

    #[test]
    fn test_round_trip_plus_strand_insertion_cigar() {
        // Block: ref chr1 [20000, 20200), asm_chr1 [2000, 2210), +
        // CIGAR: 100M10I100M
        let dir = tempfile::tempdir().unwrap();
        let data = make_index_json(&[(
            "chr1", 20000, 20200, "asm_chr1", 2000, 2210, "+", 60,
            Some("100M10I100M"),
        )]);
        let path = write_gz_json(&dir, &data);
        let index = load_index(&path).unwrap();

        let results = query(&index, "chr1", 20090, 20110, 0);
        assert_eq!(results.len(), 1);
        assert_eq!(results[0].asm_start, 2090);
        assert_eq!(results[0].asm_end, 2120);
    }

    #[test]
    fn test_round_trip_minus_strand_deletion_cigar() {
        // Block: ref chr2 [50000, 50200), asm_chr2 [3000, 3195), -
        // CIGAR: 100M5D95M
        let dir = tempfile::tempdir().unwrap();
        let data = make_index_json(&[(
            "chr2", 50000, 50200, "asm_chr2", 3000, 3195, "-", 50,
            Some("100M5D95M"),
        )]);
        let path = write_gz_json(&dir, &data);
        let index = load_index(&path).unwrap();

        let results = query(&index, "chr2", 50090, 50110, 0);
        assert_eq!(results.len(), 1);
        assert_eq!(results[0].asm_start, 3090);
        assert_eq!(results[0].asm_end, 3105);
    }

    // -----------------------------------------------------------------------
    // Biologically realistic SV scenarios
    // -----------------------------------------------------------------------

    fn bio_index() -> (tempfile::TempDir, MappingIndex) {
        let dir = tempfile::tempdir().unwrap();
        let data = make_index_json(&[
            // Scenario 1: het deletion (5kb)
            ("chr1", 1000000, 1005000, "ctg1_hap1", 500000, 505000, "+", 60, Some("5000=")),
            ("chr1", 1010000, 1015000, "ctg1_hap1", 505000, 510000, "+", 60, Some("5000=")),
            // Scenario 2: het insertion (3kb)
            ("chr1", 2000000, 2005000, "ctg1_hap1", 600000, 605000, "+", 60, Some("5000=")),
            ("chr1", 2005000, 2010000, "ctg1_hap1", 608000, 613000, "+", 60, Some("5000=")),
            // Scenario 3: inversion
            ("chr2", 3000000, 3005000, "ctg2_hap1", 700000, 705000, "+", 60, Some("5000=")),
            ("chr2", 3005000, 3015000, "ctg2_hap1", 705000, 715000, "-", 60, Some("10000=")),
            ("chr2", 3015000, 3020000, "ctg2_hap1", 715000, 720000, "+", 60, Some("5000=")),
            // Scenario 4: tandem dup
            ("chr3", 4000000, 4002000, "ctg3_hap1", 800000, 802000, "+", 60, Some("2000=")),
            ("chr3", 4000000, 4002000, "ctg3_hap1", 802000, 804000, "+", 60, Some("2000=")),
            // Scenario 5: minus-strand deletion
            ("chr5", 5000000, 5003000, "ctg5_hap2", 900000, 903000, "-", 60, Some("3000=")),
            ("chr5", 5008000, 5011000, "ctg5_hap2", 895000, 898000, "-", 60, Some("3000=")),
            // Scenario 6: minus-strand insertion
            ("chr5", 6000000, 6003000, "ctg5_hap2", 1000000, 1003000, "-", 60, Some("3000=")),
            ("chr5", 6003000, 6006000, "ctg5_hap2", 993000, 996000, "-", 60, Some("3000=")),
            // Scenario 7: translocation
            ("chr6", 7000000, 7005000, "ctg6a", 100000, 105000, "+", 60, Some("5000=")),
            ("chr6", 7006000, 7011000, "ctg6b", 200000, 205000, "+", 60, Some("5000=")),
        ]);
        let path = write_gz_json(&dir, &data);
        let index = load_index(&path).unwrap();
        (dir, index)
    }

    #[test]
    fn test_het_deletion_event_detected() {
        let (_dir, index) = bio_index();
        let results = query(&index, "chr1", 999000, 1016000, 0);
        let types: Vec<_> = results.iter().map(|r| &r.event_type).collect();
        assert_eq!(
            types,
            vec![&EventType::Alignment, &EventType::Deletion, &EventType::Alignment]
        );
    }

    #[test]
    fn test_het_deletion_gap_sizes() {
        let (_dir, index) = bio_index();
        let results = query(&index, "chr1", 999000, 1016000, 0);
        let gap = results.iter().find(|r| r.event_type == EventType::Deletion).unwrap();
        assert_eq!(gap.ref_gap_size, Some(5000));
        assert_eq!(gap.asm_gap_size, Some(0));
    }

    #[test]
    fn test_het_deletion_asm_coords() {
        let (_dir, index) = bio_index();
        let results = query(&index, "chr1", 999000, 1016000, 0);
        let gap = results.iter().find(|r| r.event_type == EventType::Deletion).unwrap();
        assert_eq!(gap.asm_start, 505000);
        assert_eq!(gap.asm_end, 505000);
    }

    #[test]
    fn test_het_deletion_flanking_alignments() {
        let (_dir, index) = bio_index();
        let results = query(&index, "chr1", 999000, 1016000, 0);
        let alns: Vec<_> = results.iter().filter(|r| r.event_type == EventType::Alignment).collect();
        assert_eq!(alns.len(), 2);
        assert_eq!(alns[0].ref_start, 1000000);
        assert_eq!(alns[0].ref_end, 1005000);
        assert_eq!(alns[0].asm_start, 500000);
        assert_eq!(alns[0].asm_end, 505000);
        assert_eq!(alns[1].ref_start, 1010000);
        assert_eq!(alns[1].ref_end, 1015000);
        assert_eq!(alns[1].asm_start, 505000);
        assert_eq!(alns[1].asm_end, 510000);
    }

    #[test]
    fn test_het_insertion_event_detected() {
        let (_dir, index) = bio_index();
        let results = query(&index, "chr1", 1999000, 2011000, 0);
        let types: Vec<_> = results.iter().map(|r| &r.event_type).collect();
        assert_eq!(
            types,
            vec![&EventType::Alignment, &EventType::Insertion, &EventType::Alignment]
        );
    }

    #[test]
    fn test_het_insertion_gap_sizes() {
        let (_dir, index) = bio_index();
        let results = query(&index, "chr1", 1999000, 2011000, 0);
        let gap = results.iter().find(|r| r.event_type == EventType::Insertion).unwrap();
        assert_eq!(gap.ref_gap_size, Some(0));
        assert_eq!(gap.asm_gap_size, Some(3000));
    }

    #[test]
    fn test_het_insertion_asm_coords() {
        let (_dir, index) = bio_index();
        let results = query(&index, "chr1", 1999000, 2011000, 0);
        let gap = results.iter().find(|r| r.event_type == EventType::Insertion).unwrap();
        assert_eq!(gap.asm_start, 605000);
        assert_eq!(gap.asm_end, 608000);
    }

    #[test]
    fn test_inversion_three_blocks() {
        let (_dir, index) = bio_index();
        let results = query(&index, "chr2", 2999000, 3021000, 0);
        let types: Vec<_> = results.iter().map(|r| &r.event_type).collect();
        assert_eq!(
            types,
            vec![
                &EventType::Alignment,
                &EventType::Inversion,
                &EventType::Alignment,
                &EventType::Inversion,
                &EventType::Alignment,
            ]
        );
    }

    #[test]
    fn test_inversion_flanking_strands() {
        let (_dir, index) = bio_index();
        let results = query(&index, "chr2", 2999000, 3021000, 0);
        let alns: Vec<_> = results.iter().filter(|r| r.event_type == EventType::Alignment).collect();
        assert_eq!(alns.len(), 3);
        assert_eq!(alns[0].strand, "+");
        assert_eq!(alns[1].strand, "-");
        assert_eq!(alns[2].strand, "+");
    }

    #[test]
    fn test_inversion_inverted_block_asm_coords() {
        let (_dir, index) = bio_index();
        let results = query(&index, "chr2", 3005000, 3015000, 0);
        assert_eq!(results.len(), 1);
        assert_eq!(results[0].strand, "-");
        assert_eq!(results[0].asm_start, 705000);
        assert_eq!(results[0].asm_end, 715000);
    }

    #[test]
    fn test_inversion_partial_query_of_inverted_block() {
        let (_dir, index) = bio_index();
        let results = query(&index, "chr2", 3007000, 3012000, 0);
        assert_eq!(results.len(), 1);
        assert_eq!(results[0].asm_start, 708000);
        assert_eq!(results[0].asm_end, 713000);
    }

    #[test]
    fn test_tandem_dup_two_hits() {
        let (_dir, index) = bio_index();
        let results = query(&index, "chr3", 4000000, 4002000, 0);
        assert_eq!(results.len(), 2);
        let mut asm_starts: Vec<u64> = results.iter().map(|r| r.asm_start).collect();
        asm_starts.sort();
        assert_eq!(asm_starts, vec![800000, 802000]);
    }

    #[test]
    fn test_tandem_dup_partial_query() {
        let (_dir, index) = bio_index();
        let results = query(&index, "chr3", 4000500, 4001500, 0);
        assert_eq!(results.len(), 2);
    }

    #[test]
    fn test_minus_strand_deletion_event() {
        let (_dir, index) = bio_index();
        let results = query(&index, "chr5", 4999000, 5012000, 0);
        let types: Vec<_> = results.iter().map(|r| &r.event_type).collect();
        assert!(types.contains(&&EventType::Deletion));
    }

    #[test]
    fn test_minus_strand_deletion_gap_sizes() {
        let (_dir, index) = bio_index();
        let results = query(&index, "chr5", 4999000, 5012000, 0);
        let gap = results.iter().find(|r| r.event_type == EventType::Deletion).unwrap();
        assert_eq!(gap.ref_gap_size, Some(5000));
        assert_eq!(gap.asm_gap_size, Some(2000));
    }

    #[test]
    fn test_minus_strand_deletion_asm_coords() {
        let (_dir, index) = bio_index();
        let results = query(&index, "chr5", 4999000, 5012000, 0);
        let gap = results.iter().find(|r| r.event_type == EventType::Deletion).unwrap();
        assert_eq!(gap.asm_start, 898000);
        assert_eq!(gap.asm_end, 900000);
    }

    #[test]
    fn test_minus_strand_insertion_event() {
        let (_dir, index) = bio_index();
        let results = query(&index, "chr5", 5999000, 6007000, 0);
        let types: Vec<_> = results.iter().map(|r| &r.event_type).collect();
        assert!(types.contains(&&EventType::Insertion));
    }

    #[test]
    fn test_minus_strand_insertion_gap_sizes() {
        let (_dir, index) = bio_index();
        let results = query(&index, "chr5", 5999000, 6007000, 0);
        let gap = results.iter().find(|r| r.event_type == EventType::Insertion).unwrap();
        assert_eq!(gap.ref_gap_size, Some(0));
        assert_eq!(gap.asm_gap_size, Some(4000));
    }

    #[test]
    fn test_minus_strand_insertion_asm_coords() {
        let (_dir, index) = bio_index();
        let results = query(&index, "chr5", 5999000, 6007000, 0);
        let gap = results.iter().find(|r| r.event_type == EventType::Insertion).unwrap();
        assert_eq!(gap.asm_start, 996000);
        assert_eq!(gap.asm_end, 1000000);
    }

    #[test]
    fn test_translocation_with_cigar() {
        let (_dir, index) = bio_index();
        let results = query(&index, "chr6", 6999000, 7012000, 0);
        let types: Vec<_> = results.iter().map(|r| &r.event_type).collect();
        assert_eq!(
            types,
            vec![&EventType::Alignment, &EventType::Translocation, &EventType::Alignment]
        );
    }

    #[test]
    fn test_translocation_different_contigs() {
        let (_dir, index) = bio_index();
        let results = query(&index, "chr6", 6999000, 7012000, 0);
        let alns: Vec<_> = results.iter().filter(|r| r.event_type == EventType::Alignment).collect();
        assert_eq!(alns[0].asm_chrom, "ctg6a");
        assert_eq!(alns[1].asm_chrom, "ctg6b");
    }

    // -----------------------------------------------------------------------
    // Edge cases: single-base, zero-width, boundaries, mapq filtering
    // -----------------------------------------------------------------------

    #[test]
    fn test_single_base_query() {
        let (_dir, index) = bio_index();
        let results = query(&index, "chr1", 1002000, 1002001, 0);
        assert_eq!(results.len(), 1);
        assert_eq!(results[0].ref_start, 1002000);
        assert_eq!(results[0].ref_end, 1002001);
        assert_eq!(results[0].asm_start, 502000);
        assert_eq!(results[0].asm_end, 502001);
    }

    #[test]
    fn test_zero_width_query() {
        let (_dir, index) = bio_index();
        let results = query(&index, "chr1", 1002000, 1002000, 0);
        assert_eq!(results.len(), 0);
    }

    #[test]
    fn test_query_at_block_start() {
        let (_dir, index) = bio_index();
        let results = query(&index, "chr1", 1000000, 1000100, 0);
        assert_eq!(results.len(), 1);
        assert_eq!(results[0].ref_start, 1000000);
        assert_eq!(results[0].asm_start, 500000);
    }

    #[test]
    fn test_query_at_block_end() {
        let (_dir, index) = bio_index();
        let results = query(&index, "chr1", 1004900, 1005000, 0);
        assert_eq!(results.len(), 1);
        assert_eq!(results[0].ref_end, 1005000);
        assert_eq!(results[0].asm_end, 505000);
    }

    #[test]
    fn test_query_just_past_block_end() {
        let (_dir, index) = bio_index();
        let results = query(&index, "chr1", 1005000, 1005001, 0);
        let alns: Vec<_> = results.iter().filter(|r| r.event_type == EventType::Alignment).collect();
        assert_eq!(alns.len(), 0);
    }

    #[test]
    fn test_no_overlap() {
        let (_dir, index) = bio_index();
        let results = query(&index, "chr1", 50000000, 50001000, 0);
        assert_eq!(results.len(), 0);
    }

    #[test]
    fn test_unknown_chrom() {
        let (_dir, index) = bio_index();
        let results = query(&index, "chrX", 0, 1000000, 0);
        assert_eq!(results.len(), 0);
    }

    #[test]
    fn test_query_exactly_matching_block() {
        let (_dir, index) = bio_index();
        let results = query(&index, "chr1", 1000000, 1005000, 0);
        assert_eq!(results.len(), 1);
        assert_eq!(results[0].ref_start, 1000000);
        assert_eq!(results[0].ref_end, 1005000);
        assert_eq!(results[0].asm_start, 500000);
        assert_eq!(results[0].asm_end, 505000);
    }

    #[test]
    fn test_wide_query_spans_deletion_and_insertion() {
        let (_dir, index) = bio_index();
        let results = query(&index, "chr1", 999000, 2011000, 0);
        let types: Vec<_> = results.iter().map(|r| &r.event_type).collect();
        assert_eq!(types.iter().filter(|t| ***t == EventType::Alignment).count(), 4);
        assert!(types.contains(&&EventType::Deletion));
        assert!(types.contains(&&EventType::Insertion));
    }

    // -----------------------------------------------------------------------
    // Mapq filtering
    // -----------------------------------------------------------------------

    fn mapq_index() -> (tempfile::TempDir, MappingIndex) {
        let dir = tempfile::tempdir().unwrap();
        let data = make_index_json(&[
            ("chr1", 100000, 105000, "ctg1", 10000, 15000, "+", 60, None),
            ("chr1", 100000, 105000, "ctg1", 10000, 15000, "+", 0, None),
            ("chr1", 110000, 115000, "ctg1", 20000, 25000, "+", 20, None),
            ("chr1", 120000, 125000, "ctg1", 30000, 35000, "+", 60, None),
        ]);
        let path = write_gz_json(&dir, &data);
        let index = load_index(&path).unwrap();
        (dir, index)
    }

    #[test]
    fn test_mapq_default_no_filtering() {
        let (_dir, index) = mapq_index();
        let results = query(&index, "chr1", 99000, 126000, 0);
        let alns: Vec<_> = results.iter().filter(|r| r.event_type == EventType::Alignment).collect();
        assert_eq!(alns.len(), 4);
    }

    #[test]
    fn test_mapq_1_excludes_zero() {
        let (_dir, index) = mapq_index();
        let results = query(&index, "chr1", 99000, 126000, 1);
        let alns: Vec<_> = results.iter().filter(|r| r.event_type == EventType::Alignment).collect();
        assert_eq!(alns.len(), 3);
        for a in &alns {
            assert!(a.mapq >= 1);
        }
    }

    #[test]
    fn test_mapq_30_excludes_low_quality() {
        let (_dir, index) = mapq_index();
        let results = query(&index, "chr1", 99000, 126000, 30);
        let alns: Vec<_> = results.iter().filter(|r| r.event_type == EventType::Alignment).collect();
        assert_eq!(alns.len(), 2);
        for a in &alns {
            assert!(a.mapq >= 30);
        }
    }

    #[test]
    fn test_mapq_60_only_primary() {
        let (_dir, index) = mapq_index();
        let results = query(&index, "chr1", 99000, 126000, 60);
        let alns: Vec<_> = results.iter().filter(|r| r.event_type == EventType::Alignment).collect();
        assert_eq!(alns.len(), 2);
        for a in &alns {
            assert_eq!(a.mapq, 60);
        }
    }

    #[test]
    fn test_mapq_too_high_returns_empty() {
        let (_dir, index) = mapq_index();
        let results = query(&index, "chr1", 99000, 126000, 61);
        assert_eq!(results.len(), 0);
    }

    #[test]
    fn test_mapq_filtering_with_gaps() {
        let (_dir, index) = mapq_index();
        let results = query(&index, "chr1", 99000, 126000, 30);
        let gaps: Vec<_> = results.iter().filter(|r| r.event_type != EventType::Alignment).collect();
        assert_eq!(gaps.len(), 1);
        assert_eq!(gaps[0].ref_start, 105000);
        assert_eq!(gaps[0].ref_end, 120000);
    }

    // -----------------------------------------------------------------------
    // CIGAR projection edge cases (via query)
    // -----------------------------------------------------------------------

    fn cigar_edge_index() -> (tempfile::TempDir, MappingIndex) {
        let dir = tempfile::tempdir().unwrap();
        let data = make_index_json(&[
            // Block 1: large deletion inside single block (1000=5000D1000=)
            ("chr7", 100000, 107000, "ctg7", 200000, 202000, "+", 60, Some("1000=5000D1000=")),
            // Block 2: multiple small indels (100=3I50=2D100=5I50=100=)
            ("chr7", 200000, 200402, "ctg7", 300000, 300408, "+", 60, Some("100=3I50=2D100=5I50=100=")),
            // Block 3: extended CIGAR (500=10X490=5D500=10I500=)
            ("chr7", 300000, 302005, "ctg7", 400000, 402010, "+", 60, Some("500=10X490=5D500=10I500=")),
            // Block 4: minus strand with insertion (2000=500I3000=)
            ("chr8", 400000, 405000, "ctg8", 500000, 505500, "-", 60, Some("2000=500I3000=")),
        ]);
        let path = write_gz_json(&dir, &data);
        let index = load_index(&path).unwrap();
        (dir, index)
    }

    #[test]
    fn test_large_del_query_before() {
        let (_dir, index) = cigar_edge_index();
        let results = query(&index, "chr7", 100000, 100500, 0);
        assert_eq!(results.len(), 1);
        assert_eq!(results[0].asm_start, 200000);
        assert_eq!(results[0].asm_end, 200500);
    }

    #[test]
    fn test_large_del_query_spanning() {
        let (_dir, index) = cigar_edge_index();
        let results = query(&index, "chr7", 100500, 106500, 0);
        assert_eq!(results.len(), 1);
        assert_eq!(results[0].asm_start, 200500);
        assert_eq!(results[0].asm_end, 201500);
    }

    #[test]
    fn test_large_del_query_exactly_in_deletion() {
        let (_dir, index) = cigar_edge_index();
        let results = query(&index, "chr7", 101000, 106000, 0);
        assert_eq!(results.len(), 1);
        assert_eq!(results[0].asm_start, 201000);
        assert_eq!(results[0].asm_end, 201000);
    }

    #[test]
    fn test_large_del_query_after() {
        let (_dir, index) = cigar_edge_index();
        let results = query(&index, "chr7", 106500, 107000, 0);
        assert_eq!(results.len(), 1);
        assert_eq!(results[0].asm_start, 201500);
        assert_eq!(results[0].asm_end, 202000);
    }

    #[test]
    fn test_multi_indel_full_block() {
        let (_dir, index) = cigar_edge_index();
        let results = query(&index, "chr7", 200000, 200402, 0);
        assert_eq!(results.len(), 1);
        assert_eq!(results[0].asm_start, 300000);
        assert_eq!(results[0].asm_end, 300408);
    }

    #[test]
    fn test_multi_indel_query_before_first_insertion() {
        let (_dir, index) = cigar_edge_index();
        let results = query(&index, "chr7", 200000, 200050, 0);
        assert_eq!(results.len(), 1);
        assert_eq!(results[0].asm_start, 300000);
        assert_eq!(results[0].asm_end, 300050);
    }

    #[test]
    fn test_multi_indel_query_spanning_first_insertion() {
        let (_dir, index) = cigar_edge_index();
        let results = query(&index, "chr7", 200090, 200110, 0);
        assert_eq!(results.len(), 1);
        assert_eq!(results[0].asm_start, 300090);
        assert_eq!(results[0].asm_end, 300113);
    }

    #[test]
    fn test_multi_indel_query_spanning_deletion() {
        let (_dir, index) = cigar_edge_index();
        let results = query(&index, "chr7", 200145, 200160, 0);
        assert_eq!(results.len(), 1);
        assert_eq!(results[0].asm_start, 300148);
        assert_eq!(results[0].asm_end, 300161);
    }

    #[test]
    fn test_eqx_cigar_mismatch_region() {
        let (_dir, index) = cigar_edge_index();
        let results = query(&index, "chr7", 300495, 300515, 0);
        assert_eq!(results.len(), 1);
        assert_eq!(results[0].asm_start, 400495);
        assert_eq!(results[0].asm_end, 400515);
    }

    #[test]
    fn test_eqx_cigar_spanning_del_and_ins() {
        let (_dir, index) = cigar_edge_index();
        let results = query(&index, "chr7", 300990, 301510, 0);
        assert_eq!(results.len(), 1);
        assert_eq!(results[0].asm_start, 400990);
        assert_eq!(results[0].asm_end, 401515);
    }

    #[test]
    fn test_minus_strand_insertion_in_cigar_full() {
        let (_dir, index) = cigar_edge_index();
        let results = query(&index, "chr8", 400000, 405000, 0);
        assert_eq!(results.len(), 1);
        assert_eq!(results[0].asm_start, 500000);
        assert_eq!(results[0].asm_end, 505500);
    }

    #[test]
    fn test_minus_strand_insertion_in_cigar_partial_before() {
        let (_dir, index) = cigar_edge_index();
        let results = query(&index, "chr8", 400000, 401000, 0);
        assert_eq!(results.len(), 1);
        assert_eq!(results[0].asm_start, 504500);
        assert_eq!(results[0].asm_end, 505500);
    }

    #[test]
    fn test_minus_strand_insertion_in_cigar_straddling() {
        let (_dir, index) = cigar_edge_index();
        let results = query(&index, "chr8", 401500, 402500, 0);
        assert_eq!(results.len(), 1);
        assert_eq!(results[0].asm_start, 502500);
        assert_eq!(results[0].asm_end, 504000);
    }

    #[test]
    fn test_minus_strand_insertion_in_cigar_after() {
        let (_dir, index) = cigar_edge_index();
        let results = query(&index, "chr8", 403000, 404000, 0);
        assert_eq!(results.len(), 1);
        assert_eq!(results[0].asm_start, 501000);
        assert_eq!(results[0].asm_end, 502000);
    }

    // -----------------------------------------------------------------------
    // project_cigar unit-level edge cases
    // -----------------------------------------------------------------------

    #[test]
    fn test_query_exactly_at_cigar_end() {
        let ops = parse_cigar("100M50D100M");
        let (s, e) = project_cigar(&ops, 0, 250, "+", 0, 200, None);
        assert_eq!((s, e), (0, 200));
    }

    #[test]
    fn test_query_beyond_cigar_end() {
        let ops = parse_cigar("100M");
        let (s, e) = project_cigar(&ops, 0, 200, "+", 0, 100, None);
        assert_eq!((s, e), (0, 100));
    }

    #[test]
    fn test_single_base_at_deletion_start() {
        let ops = parse_cigar("50M10D50M");
        let (s, e) = project_cigar(&ops, 50, 51, "+", 0, 100, None);
        assert_eq!(s, 50);
        assert_eq!(e, 50);
    }

    #[test]
    fn test_single_base_at_deletion_end() {
        let ops = parse_cigar("50M10D50M");
        let (s, e) = project_cigar(&ops, 59, 60, "+", 0, 100, None);
        assert_eq!(s, 50);
        assert_eq!(e, 50);
    }

    #[test]
    fn test_single_base_right_after_deletion() {
        let ops = parse_cigar("50M10D50M");
        let (s, e) = project_cigar(&ops, 60, 61, "+", 0, 100, None);
        assert_eq!(s, 50);
        assert_eq!(e, 51);
    }

    #[test]
    fn test_insertion_at_cigar_start() {
        let ops = parse_cigar("10I90M");
        let (s, e) = project_cigar(&ops, 0, 10, "+", 0, 100, None);
        assert_eq!(s, 10);
        assert_eq!(e, 20);
    }

    #[test]
    fn test_consecutive_deletions() {
        let ops = parse_cigar("50M5D10D50M");
        let (s, e) = project_cigar(&ops, 45, 70, "+", 0, 100, None);
        assert_eq!(s, 45);
        assert_eq!(e, 55);
    }

    #[test]
    fn test_consecutive_insertions() {
        let ops = parse_cigar("50M10I5I50M");
        let (s, e) = project_cigar(&ops, 45, 55, "+", 0, 115, None);
        assert_eq!(s, 45);
        assert_eq!(e, 70);
    }

    #[test]
    fn test_minus_strand_with_insertion_at_cigar_start() {
        let ops = parse_cigar("10I100M");
        let (s, e) = project_cigar(&ops, 0, 100, "-", 200, 310, None);
        assert_eq!((s, e), (200, 300));
    }

    // -----------------------------------------------------------------------
    // Soft clip / hard clip: S is asm-consuming, H/P consume nothing
    // -----------------------------------------------------------------------

    #[test]
    fn test_soft_clip_asm_consuming() {
        // 10S90M: S advances asm by 10, ref consumed = 90, asm consumed = 100
        let ops = parse_cigar("10S90M");
        let idx = build_cigar_index(&ops);
        assert_eq!(idx.ref_cumul, vec![0, 0, 90]);
        assert_eq!(idx.asm_cumul, vec![0, 10, 100]);
    }

    #[test]
    fn test_hard_clip_consumes_nothing() {
        let ops = parse_cigar("5H100M3H");
        let idx = build_cigar_index(&ops);
        assert_eq!(idx.ref_cumul, vec![0, 0, 100, 100]);
        assert_eq!(idx.asm_cumul, vec![0, 0, 100, 100]);
    }

    #[test]
    fn test_padding_consumes_nothing() {
        let ops = parse_cigar("100M5P50M");
        let idx = build_cigar_index(&ops);
        assert_eq!(idx.ref_cumul, vec![0, 100, 100, 150]);
        assert_eq!(idx.asm_cumul, vec![0, 100, 100, 150]);
    }

    // -----------------------------------------------------------------------
    // project_cigar with/without index produce identical results
    // -----------------------------------------------------------------------

    fn assert_same_result(cigar: &str, ref_start: u64, ref_end: u64, strand: &str, asm_s: u64, asm_e: u64) {
        let ops = parse_cigar(cigar);
        let idx = build_cigar_index(&ops);
        let r_linear = project_cigar(&ops, ref_start, ref_end, strand, asm_s, asm_e, None);
        let r_indexed = project_cigar(&ops, ref_start, ref_end, strand, asm_s, asm_e, Some(&idx));
        assert_eq!(r_linear, r_indexed, "CIGAR={cigar} ref=[{ref_start},{ref_end}) strand={strand}");
    }

    #[test]
    fn test_index_vs_linear_simple_match() {
        assert_same_result("500M", 100, 200, "+", 0, 1000);
    }

    #[test]
    fn test_index_vs_linear_with_deletion() {
        assert_same_result("100M50D200M", 90, 160, "+", 0, 1000);
    }

    #[test]
    fn test_index_vs_linear_with_insertion() {
        assert_same_result("100M10I200M", 90, 200, "+", 0, 1000);
    }

    #[test]
    fn test_index_vs_linear_query_inside_deletion() {
        assert_same_result("100M50D200M", 110, 140, "+", 0, 1000);
    }

    #[test]
    fn test_index_vs_linear_query_spans_deletion() {
        assert_same_result("100M50D200M", 50, 200, "+", 0, 1000);
    }

    #[test]
    fn test_index_vs_linear_minus_strand() {
        assert_same_result("100M5D100M", 50, 150, "-", 0, 200);
    }

    #[test]
    fn test_index_vs_linear_complex_various_offsets() {
        let cigar = "100M5D200M3I50M10D100M";
        for (s, e) in [(0, 50), (95, 110), (200, 350), (300, 465)] {
            assert_same_result(cigar, s, e, "+", 0, 500);
        }
    }

    #[test]
    fn test_index_vs_linear_long_cigar_deep_offset() {
        let cigar: String = (0..250).map(|_| "1000M5D").collect::<Vec<_>>().join("");
        let ops = parse_cigar(&cigar);
        let idx = build_cigar_index(&ops);
        let r_linear = project_cigar(&ops, 200000, 200500, "+", 0, 300000, None);
        let r_indexed = project_cigar(&ops, 200000, 200500, "+", 0, 300000, Some(&idx));
        assert_eq!(r_linear, r_indexed);
    }

    // -----------------------------------------------------------------------
    // No-CIGAR fallback (linear interpolation)
    // -----------------------------------------------------------------------

    #[test]
    fn test_no_cigar_plus_strand_linear() {
        let dir = tempfile::tempdir().unwrap();
        let data = make_index_json(&[
            ("chr1", 10000, 13000, "ctg1", 1000, 4000, "+", 60, None),
            ("chr1", 20000, 23000, "ctg1", 5000, 8000, "-", 50, None),
        ]);
        let path = write_gz_json(&dir, &data);
        let index = load_index(&path).unwrap();

        let results = query(&index, "chr1", 11000, 12000, 0);
        assert_eq!(results.len(), 1);
        assert_eq!(results[0].asm_start, 2000);
        assert_eq!(results[0].asm_end, 3000);
    }

    #[test]
    fn test_no_cigar_minus_strand_linear() {
        let dir = tempfile::tempdir().unwrap();
        let data = make_index_json(&[
            ("chr1", 10000, 13000, "ctg1", 1000, 4000, "+", 60, None),
            ("chr1", 20000, 23000, "ctg1", 5000, 8000, "-", 50, None),
        ]);
        let path = write_gz_json(&dir, &data);
        let index = load_index(&path).unwrap();

        let results = query(&index, "chr1", 21000, 22000, 0);
        assert_eq!(results.len(), 1);
        assert_eq!(results[0].asm_start, 6000);
        assert_eq!(results[0].asm_end, 7000);
    }

    // -----------------------------------------------------------------------
    // Phantom gap regression
    // -----------------------------------------------------------------------

    fn phantom_gap_index() -> (tempfile::TempDir, MappingIndex) {
        let dir = tempfile::tempdir().unwrap();
        let data = make_index_json(&[
            ("chr11", 100, 1000, "ctg11", 5000, 5900, "+", 60, None),
            ("chr11", 200, 300, "ctg11", 8000, 8100, "+", 5, None),
            ("chr11", 400, 500, "ctg11", 5300, 5400, "+", 60, None),
            ("chr11", 1200, 1300, "ctg11", 5900, 5950, "+", 60, None),
        ]);
        let path = write_gz_json(&dir, &data);
        let index = load_index(&path).unwrap();
        (dir, index)
    }

    #[test]
    fn test_no_phantom_gap_from_nested_supplementary() {
        let (_dir, index) = phantom_gap_index();
        let results = query(&index, "chr11", 0, 600, 0);
        let gaps: Vec<_> = results.iter().filter(|r| r.event_type != EventType::Alignment).collect();
        assert_eq!(gaps.len(), 0, "Unexpected gap events: {gaps:?}");
    }

    #[test]
    fn test_genuine_gap_still_detected() {
        let (_dir, index) = phantom_gap_index();
        let results = query(&index, "chr11", 0, 1500, 0);
        let gaps: Vec<_> = results.iter().filter(|r| r.event_type != EventType::Alignment).collect();
        assert_eq!(gaps.len(), 1, "Expected 1 gap event, got: {gaps:?}");
        assert_eq!(gaps[0].ref_start, 1000);
        assert_eq!(gaps[0].ref_end, 1200);
    }

    #[test]
    fn test_genuine_gap_classified_as_deletion() {
        let (_dir, index) = phantom_gap_index();
        let results = query(&index, "chr11", 0, 1500, 0);
        let gap = results.iter().find(|r| r.event_type != EventType::Alignment).unwrap();
        assert_eq!(gap.event_type, EventType::Deletion);
        assert_eq!(gap.ref_gap_size, Some(200));
    }

    #[test]
    fn test_all_alignment_blocks_returned() {
        let (_dir, index) = phantom_gap_index();
        let results = query(&index, "chr11", 0, 1500, 0);
        let alns: Vec<_> = results.iter().filter(|r| r.event_type == EventType::Alignment).collect();
        assert_eq!(alns.len(), 4);
    }

    #[test]
    fn test_no_phantom_gap_with_mapq_filter() {
        let (_dir, index) = phantom_gap_index();
        let results = query(&index, "chr11", 0, 600, 10);
        let gaps: Vec<_> = results.iter().filter(|r| r.event_type != EventType::Alignment).collect();
        assert_eq!(gaps.len(), 0, "Unexpected gap events: {gaps:?}");
    }

    // -----------------------------------------------------------------------
    // Interleaving contigs
    // -----------------------------------------------------------------------

    fn interleaving_index() -> (tempfile::TempDir, MappingIndex) {
        let dir = tempfile::tempdir().unwrap();
        let data = make_index_json(&[
            ("chr12", 100, 500, "ctgX", 1000, 1400, "+", 60, None),
            ("chr12", 300, 800, "ctgY", 2000, 2500, "+", 60, None),
            ("chr12", 1000, 1200, "ctgX", 1400, 1600, "+", 60, None),
        ]);
        let path = write_gz_json(&dir, &data);
        let index = load_index(&path).unwrap();
        (dir, index)
    }

    #[test]
    fn test_interleaving_all_blocks_returned() {
        let (_dir, index) = interleaving_index();
        let results = query(&index, "chr12", 0, 1500, 0);
        let alns: Vec<_> = results.iter().filter(|r| r.event_type == EventType::Alignment).collect();
        assert_eq!(alns.len(), 3);
    }

    #[test]
    fn test_interleaving_gap_uses_high_water_mark() {
        let (_dir, index) = interleaving_index();
        let results = query(&index, "chr12", 0, 1500, 0);
        let gaps: Vec<_> = results.iter().filter(|r| r.event_type != EventType::Alignment).collect();
        assert_eq!(gaps.len(), 1);
        assert_eq!(gaps[0].ref_start, 800);
        assert_eq!(gaps[0].ref_end, 1000);
    }

    #[test]
    fn test_interleaving_gap_classified_as_translocation() {
        let (_dir, index) = interleaving_index();
        let results = query(&index, "chr12", 0, 1500, 0);
        let gaps: Vec<_> = results.iter().filter(|r| r.event_type != EventType::Alignment).collect();
        assert_eq!(gaps[0].event_type, EventType::Translocation);
    }

    #[test]
    fn test_no_phantom_gap_in_overlap_zone() {
        let (_dir, index) = interleaving_index();
        let results = query(&index, "chr12", 300, 500, 0);
        let gaps: Vec<_> = results.iter().filter(|r| r.event_type != EventType::Alignment).collect();
        assert_eq!(gaps.len(), 0);
    }

    // -----------------------------------------------------------------------
    // CIGAR insertion at exact query boundary
    // -----------------------------------------------------------------------

    fn insertion_boundary_index() -> (tempfile::TempDir, MappingIndex) {
        let dir = tempfile::tempdir().unwrap();
        let data = make_index_json(&[
            ("chr13", 100, 115, "ctg_ins", 200, 265, "+", 60, Some("10=50I5=")),
        ]);
        let path = write_gz_json(&dir, &data);
        let index = load_index(&path).unwrap();
        (dir, index)
    }

    #[test]
    fn test_straddle_insertion_includes_inserted_bases() {
        let (_dir, index) = insertion_boundary_index();
        let results = query(&index, "chr13", 105, 113, 0);
        assert_eq!(results.len(), 1);
        assert_eq!(results[0].asm_start, 205);
        assert_eq!(results[0].asm_end, 263);
    }

    #[test]
    fn test_query_starting_at_insertion_skips_insert() {
        let (_dir, index) = insertion_boundary_index();
        let results = query(&index, "chr13", 110, 115, 0);
        assert_eq!(results.len(), 1);
        assert_eq!(results[0].asm_start, 260);
        assert_eq!(results[0].asm_end, 265);
    }

    #[test]
    fn test_query_ending_at_insertion_excludes_insert() {
        let (_dir, index) = insertion_boundary_index();
        let results = query(&index, "chr13", 100, 110, 0);
        assert_eq!(results.len(), 1);
        assert_eq!(results[0].asm_start, 200);
        assert_eq!(results[0].asm_end, 210);
    }

    #[test]
    fn test_one_base_before_insertion_captures_it() {
        let (_dir, index) = insertion_boundary_index();
        let results = query(&index, "chr13", 109, 115, 0);
        assert_eq!(results.len(), 1);
        assert_eq!(results[0].asm_start, 209);
        assert_eq!(results[0].asm_end, 265);
    }

    // -----------------------------------------------------------------------
    // Tandem duplication complex (adjacent ref, overlapping asm)
    // -----------------------------------------------------------------------

    #[test]
    fn test_tandem_dup_complex_classification() {
        let dir = tempfile::tempdir().unwrap();
        let data = make_index_json(&[
            ("chr14", 5000, 6000, "ctgT", 1000, 2000, "+", 60, None),
            ("chr14", 6000, 7000, "ctgT", 1500, 2500, "+", 60, None),
        ]);
        let path = write_gz_json(&dir, &data);
        let index = load_index(&path).unwrap();

        let results = query(&index, "chr14", 4500, 7500, 0);
        let gaps: Vec<_> = results.iter().filter(|r| r.event_type != EventType::Alignment).collect();
        assert_eq!(gaps.len(), 1);
        assert_eq!(gaps[0].event_type, EventType::Complex);
        assert_eq!(gaps[0].asm_gap_size, Some(-500));
        assert_eq!(gaps[0].ref_gap_size, Some(0));
    }

    // -----------------------------------------------------------------------
    // Tandem dup: identical ref coordinates, no gap emitted
    // -----------------------------------------------------------------------

    #[test]
    fn test_tandem_dup_identical_ref_no_gap() {
        let dir = tempfile::tempdir().unwrap();
        let data = make_index_json(&[
            ("chr15", 8000, 9000, "ctgD", 3000, 4000, "+", 60, None),
            ("chr15", 8000, 9000, "ctgD", 4000, 5000, "+", 60, None),
        ]);
        let path = write_gz_json(&dir, &data);
        let index = load_index(&path).unwrap();

        let results = query(&index, "chr15", 7500, 9500, 0);
        let gaps: Vec<_> = results.iter().filter(|r| r.event_type != EventType::Alignment).collect();
        assert_eq!(gaps.len(), 0);
        let alns: Vec<_> = results.iter().filter(|r| r.event_type == EventType::Alignment).collect();
        assert_eq!(alns.len(), 2);
        let mut asm_starts: Vec<u64> = alns.iter().map(|r| r.asm_start).collect();
        asm_starts.sort();
        assert_eq!(asm_starts, vec![3000, 4000]);
    }

    // -----------------------------------------------------------------------
    // Overlapping / supplementary alignment blocks
    // -----------------------------------------------------------------------

    fn overlap_index() -> (tempfile::TempDir, MappingIndex) {
        let dir = tempfile::tempdir().unwrap();
        let data = make_index_json(&[
            ("chr9", 1000, 3000, "ctg9", 10000, 12000, "+", 60, None),
            ("chr9", 2500, 5000, "ctg9", 12000, 14500, "+", 5, None),
            ("chr10", 100000, 102000, "ctg10a", 200000, 202000, "+", 60, None),
            ("chr10", 100000, 102000, "ctg10b", 300000, 302000, "+", 60, None),
        ]);
        let path = write_gz_json(&dir, &data);
        let index = load_index(&path).unwrap();
        (dir, index)
    }

    #[test]
    fn test_overlapping_blocks_both_returned() {
        let (_dir, index) = overlap_index();
        let results = query(&index, "chr9", 1000, 5000, 0);
        let alns: Vec<_> = results.iter().filter(|r| r.event_type == EventType::Alignment).collect();
        assert_eq!(alns.len(), 2);
    }

    #[test]
    fn test_overlapping_blocks_ref_overlap_region() {
        let (_dir, index) = overlap_index();
        let results = query(&index, "chr9", 2600, 2900, 0);
        let alns: Vec<_> = results.iter().filter(|r| r.event_type == EventType::Alignment).collect();
        assert_eq!(alns.len(), 2);
        for a in &alns {
            assert_eq!(a.ref_start, 2600);
            assert_eq!(a.ref_end, 2900);
        }
    }

    #[test]
    fn test_overlapping_blocks_different_asm_coords() {
        let (_dir, index) = overlap_index();
        let results = query(&index, "chr9", 2600, 2900, 0);
        let alns: Vec<_> = results.iter().filter(|r| r.event_type == EventType::Alignment).collect();
        let asm_starts: Vec<u64> = alns.iter().map(|r| r.asm_start).collect();
        assert_eq!(asm_starts.len(), 2);
        assert_ne!(asm_starts[0], asm_starts[1]);
    }

    #[test]
    fn test_overlapping_filtered_by_mapq() {
        let (_dir, index) = overlap_index();
        let results = query(&index, "chr9", 2600, 2900, 10);
        let alns: Vec<_> = results.iter().filter(|r| r.event_type == EventType::Alignment).collect();
        assert_eq!(alns.len(), 1);
        assert_eq!(alns[0].mapq, 60);
    }

    #[test]
    fn test_segdup_two_contigs_same_ref() {
        let (_dir, index) = overlap_index();
        let results = query(&index, "chr10", 100000, 102000, 0);
        let alns: Vec<_> = results.iter().filter(|r| r.event_type == EventType::Alignment).collect();
        assert_eq!(alns.len(), 2);
        let mut contigs: Vec<_> = alns.iter().map(|r| r.asm_chrom.as_str()).collect();
        contigs.sort();
        assert_eq!(contigs, vec!["ctg10a", "ctg10b"]);
    }

    #[test]
    fn test_segdup_different_asm_ranges() {
        let (_dir, index) = overlap_index();
        let results = query(&index, "chr10", 100500, 101500, 0);
        let alns: Vec<_> = results.iter().filter(|r| r.event_type == EventType::Alignment).collect();
        assert_eq!(alns.len(), 2);
        let mut asm_ranges: Vec<_> = alns.iter().map(|r| (r.asm_start, r.asm_end)).collect();
        asm_ranges.sort();
        assert_eq!(asm_ranges, vec![(200500, 201500), (300500, 301500)]);
    }
}
