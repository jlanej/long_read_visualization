//! Pileup layout: pack aligned reads into rows for rendering.
//!
//! Uses a greedy interval-scheduling algorithm to minimize the number of rows
//! while placing reads on the first row where they fit without overlap.

use crate::genome::bam::AlignedRead;
use crate::pileup::{PileupData, ReadRow};

/// Minimum gap (bp) between reads on the same row.
const MIN_GAP: u64 = 5;

/// Pack reads into rows using a greedy first-fit algorithm.
///
/// Reads are sorted by start position, then each read is placed in the first
/// row where it doesn't overlap with the last placed read. This produces a
/// compact vertical layout.
pub fn layout_reads(
    mut reads: Vec<AlignedRead>,
    region_start: u64,
    region_end: u64,
    sort_by_hp: bool,
) -> PileupData {
    let total_reads = reads.len();

    // Sort by HP tag first (if sorting by haplotype), then by start position
    if sort_by_hp {
        reads.sort_by(|a, b| a.hp_tag.cmp(&b.hp_tag).then(a.start.cmp(&b.start)));
    } else {
        reads.sort_by_key(|r| r.start);
    }

    let mut rows: Vec<ReadRow> = Vec::new();
    // Track the end position of the last read in each row
    let mut row_ends: Vec<u64> = Vec::new();

    for read in reads {
        // Find the first row where this read fits
        let mut placed = false;
        for (i, row_end) in row_ends.iter_mut().enumerate() {
            if read.start >= *row_end + MIN_GAP {
                rows[i].reads.push(read.clone());
                *row_end = read.end;
                placed = true;
                break;
            }
        }

        if !placed {
            // Create a new row
            let end = read.end;
            let mut row = ReadRow::new();
            row.reads.push(read);
            rows.push(row);
            row_ends.push(end);
        }
    }

    PileupData {
        rows,
        total_reads,
        region_start,
        region_end,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use noodles::sam::alignment::record::cigar::op::Kind as CigarOpKind;

    fn make_read(name: &str, start: u64, end: u64, hp: u8) -> AlignedRead {
        AlignedRead {
            name: name.to_string(),
            chrom: "chr1".to_string(),
            start,
            end,
            hp_tag: hp,
            mapq: 60,
            is_reverse: false,
            cigar_ops: vec![((end - start) as u32, CigarOpKind::Match)],
            flags: 0,
        }
    }

    #[test]
    fn test_layout_empty() {
        let pileup = layout_reads(vec![], 0, 1000, false);
        assert_eq!(pileup.rows.len(), 0);
        assert_eq!(pileup.total_reads, 0);
    }

    #[test]
    fn test_layout_single_read() {
        let reads = vec![make_read("r1", 100, 200, 0)];
        let pileup = layout_reads(reads, 0, 1000, false);
        assert_eq!(pileup.rows.len(), 1);
        assert_eq!(pileup.total_reads, 1);
        assert_eq!(pileup.rows[0].reads.len(), 1);
    }

    #[test]
    fn test_layout_non_overlapping() {
        // Two non-overlapping reads should fit on one row
        let reads = vec![
            make_read("r1", 100, 200, 0),
            make_read("r2", 300, 400, 0),
        ];
        let pileup = layout_reads(reads, 0, 1000, false);
        assert_eq!(pileup.rows.len(), 1);
        assert_eq!(pileup.rows[0].reads.len(), 2);
    }

    #[test]
    fn test_layout_overlapping() {
        // Two overlapping reads need two rows
        let reads = vec![
            make_read("r1", 100, 300, 0),
            make_read("r2", 200, 400, 0),
        ];
        let pileup = layout_reads(reads, 0, 1000, false);
        assert_eq!(pileup.rows.len(), 2);
    }

    #[test]
    fn test_layout_many_overlapping() {
        // Five fully overlapping reads need five rows
        let reads: Vec<_> = (0..5)
            .map(|i| make_read(&format!("r{i}"), 100, 200, 0))
            .collect();
        let pileup = layout_reads(reads, 0, 1000, false);
        assert_eq!(pileup.rows.len(), 5);
        assert_eq!(pileup.total_reads, 5);
    }

    #[test]
    fn test_layout_hp_sorted() {
        // Reads should group by HP tag when sort_by_hp is true
        let reads = vec![
            make_read("r1", 100, 200, 2),
            make_read("r2", 100, 200, 1),
            make_read("r3", 100, 200, 0),
        ];
        let pileup = layout_reads(reads, 0, 1000, true);
        // Each read overlaps, so they all get different rows
        assert_eq!(pileup.rows.len(), 3);
        // First row should have HP=0 (sorted first)
        assert_eq!(pileup.rows[0].reads[0].hp_tag, 0);
        assert_eq!(pileup.rows[1].reads[0].hp_tag, 1);
        assert_eq!(pileup.rows[2].reads[0].hp_tag, 2);
    }

    #[test]
    fn test_layout_staircase() {
        // Staircase pattern: each read starts slightly after the previous
        let reads: Vec<_> = (0..10)
            .map(|i| make_read(&format!("r{i}"), i * 50, i * 50 + 100, 0))
            .collect();
        let pileup = layout_reads(reads, 0, 1000, false);
        // Should pack efficiently into just a few rows
        assert!(pileup.rows.len() <= 3, "Staircase should pack into ≤3 rows");
        assert_eq!(pileup.total_reads, 10);
    }
}
