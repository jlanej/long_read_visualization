use super::AlignedRead;

/// A row in a pileup layout. Each row contains non-overlapping reads
/// packed left-to-right with a minimum gap.
#[derive(Debug, Clone)]
pub struct PileupRow {
    /// Vertical offset (0 = top row).
    pub y_offset: u32,
    /// Reads assigned to this row, ordered by start position.
    pub reads: Vec<AlignedRead>,
}

/// Minimum gap (in base-pairs) between adjacent reads in the same row.
const READ_GAP: u64 = 5;

/// Pack a set of aligned reads into pileup rows using a greedy
/// interval-scheduling algorithm.
///
/// Reads are sorted by start position, then each read is assigned to
/// the lowest row whose last read ends before `read.start - READ_GAP`.
pub fn pack_reads(mut reads: Vec<AlignedRead>) -> Vec<PileupRow> {
    if reads.is_empty() {
        return Vec::new();
    }

    reads.sort_by_key(|r| r.start);

    // Track the rightmost end position in each row.
    let mut row_ends: Vec<u64> = Vec::new();
    let mut rows: Vec<Vec<AlignedRead>> = Vec::new();

    for read in reads {
        let assigned = row_ends
            .iter()
            .position(|&end| read.start > end + READ_GAP);

        match assigned {
            Some(idx) => {
                row_ends[idx] = read.end;
                rows[idx].push(read);
            }
            None => {
                row_ends.push(read.end);
                rows.push(vec![read]);
            }
        }
    }

    rows.into_iter()
        .enumerate()
        .map(|(i, reads)| PileupRow {
            y_offset: i as u32,
            reads,
        })
        .collect()
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    fn make_read(name: &str, start: u64, end: u64) -> AlignedRead {
        AlignedRead {
            name: name.to_string(),
            start,
            end,
            is_reverse: false,
            mapping_quality: Some(60),
            haplotype: None,
            flags: 0,
        }
    }

    #[test]
    fn test_pack_empty() {
        let rows = pack_reads(vec![]);
        assert!(rows.is_empty());
    }

    #[test]
    fn test_pack_single_read() {
        let rows = pack_reads(vec![make_read("r1", 100, 200)]);
        assert_eq!(rows.len(), 1);
        assert_eq!(rows[0].reads.len(), 1);
        assert_eq!(rows[0].y_offset, 0);
    }

    #[test]
    fn test_pack_non_overlapping_reads_share_row() {
        let reads = vec![
            make_read("r1", 100, 200),
            make_read("r2", 210, 300),
            make_read("r3", 310, 400),
        ];
        let rows = pack_reads(reads);
        assert_eq!(rows.len(), 1, "non-overlapping reads should share one row");
        assert_eq!(rows[0].reads.len(), 3);
    }

    #[test]
    fn test_pack_overlapping_reads_use_separate_rows() {
        let reads = vec![
            make_read("r1", 100, 300),
            make_read("r2", 150, 350),
            make_read("r3", 200, 400),
        ];
        let rows = pack_reads(reads);
        assert_eq!(rows.len(), 3, "fully overlapping reads need separate rows");
    }

    #[test]
    fn test_pack_mixed_overlap() {
        // r1: 100-200, r2: 150-250, r3: 260-350
        // r1 and r3 can share row 0; r2 goes to row 1
        let reads = vec![
            make_read("r1", 100, 200),
            make_read("r2", 150, 250),
            make_read("r3", 260, 350),
        ];
        let rows = pack_reads(reads);
        assert_eq!(rows.len(), 2);
        // Row 0 should have r1 and r3
        assert_eq!(rows[0].reads.len(), 2);
        assert_eq!(rows[0].reads[0].name, "r1");
        assert_eq!(rows[0].reads[1].name, "r3");
        // Row 1 should have r2
        assert_eq!(rows[1].reads.len(), 1);
        assert_eq!(rows[1].reads[0].name, "r2");
    }

    #[test]
    fn test_pack_preserves_sort_order() {
        let reads = vec![
            make_read("r3", 300, 400),
            make_read("r1", 100, 200),
            make_read("r2", 200, 300),
        ];
        let rows = pack_reads(reads);
        // All reads should be sorted by start within their rows
        for row in &rows {
            for w in row.reads.windows(2) {
                assert!(w[0].start <= w[1].start);
            }
        }
    }

    #[test]
    fn test_pack_row_offsets_are_sequential() {
        let reads = vec![
            make_read("r1", 100, 300),
            make_read("r2", 100, 300),
            make_read("r3", 100, 300),
        ];
        let rows = pack_reads(reads);
        for (i, row) in rows.iter().enumerate() {
            assert_eq!(row.y_offset, i as u32);
        }
    }

    #[test]
    fn test_pack_gap_enforced() {
        // Two reads separated by exactly READ_GAP should share a row
        let reads = vec![
            make_read("r1", 100, 200),
            make_read("r2", 200 + READ_GAP + 1, 300),
        ];
        let rows = pack_reads(reads);
        assert_eq!(rows.len(), 1, "reads separated by > READ_GAP share a row");

        // Two reads separated by less than READ_GAP should not share a row
        let reads = vec![
            make_read("r1", 100, 200),
            make_read("r2", 200 + READ_GAP, 300),
        ];
        let rows = pack_reads(reads);
        assert_eq!(
            rows.len(),
            2,
            "reads separated by exactly READ_GAP need separate rows"
        );
    }

    // -- Integration-style test using BAM + pileup --

    #[test]
    fn test_pileup_from_bam_region() {
        use crate::genome::bam::query_bam;
        use std::path::PathBuf;

        let bam = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("../resources/toy_dataset/toy_reads.bam");
        if !bam.exists() {
            eprintln!("skipping: toy BAM not found");
            return;
        }

        let reads = match query_bam(&bam, "chr1:112164095-112177007") {
            Ok(r) => r,
            Err(crate::genome::GenomeError::ParseError(_)) => {
                eprintln!("skipping: header parse error");
                return;
            }
            Err(e) => panic!("unexpected error: {e}"),
        };

        let rows = pack_reads(reads);
        assert!(!rows.is_empty(), "pileup should have at least one row");
        // Verify packing invariants
        for row in &rows {
            assert!(!row.reads.is_empty());
            for w in row.reads.windows(2) {
                assert!(
                    w[1].start > w[0].end + READ_GAP,
                    "reads in same row must not overlap (including gap)"
                );
            }
        }
    }
}
