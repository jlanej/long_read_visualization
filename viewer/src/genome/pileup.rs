use std::cmp::Reverse;
use std::collections::BinaryHeap;

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

/// Default indel threshold for suppression (base pairs).
///
/// Set to 50 bp for long-read data where small indels are common sequencing
/// noise and structural variants of interest are typically ≥50 bp.
pub const DEFAULT_INDEL_THRESHOLD: u32 = 50;

/// Display configuration for a pileup panel.
#[derive(Debug, Clone)]
pub struct PileupDisplayConfig {
    /// If true, use squished (compact, single-pixel-height) rows.
    /// If false, use expanded (taller, multi-pixel-height) rows.
    pub squished: bool,
    /// If true, suppress display of indels ≤ `indel_threshold` bp.
    pub hide_small_indels: bool,
    /// Indel size threshold (bp) below which indels are hidden.
    pub indel_threshold: u32,
    /// If true, group reads by HP tag (hap1 → hap2 → unphased) before packing.
    pub sort_by_haplotype: bool,
    /// If true, render soft-clipped bases as semi-transparent extensions.
    pub show_soft_clips: bool,
    /// If true, render base-level mismatches with standard nucleotide coloring.
    pub show_mismatches: bool,
    /// If true, show coverage histogram above reads.
    pub show_coverage: bool,
}

impl Default for PileupDisplayConfig {
    fn default() -> Self {
        Self {
            squished: true,
            hide_small_indels: true,
            indel_threshold: DEFAULT_INDEL_THRESHOLD,
            sort_by_haplotype: false,
            show_soft_clips: false,
            show_mismatches: false,
            show_coverage: false,
        }
    }
}

impl PileupDisplayConfig {
    /// Return a copy with `squished` overridden.
    pub fn with_squished(&self, squished: bool) -> Self {
        let mut cfg = self.clone();
        cfg.squished = squished;
        cfg
    }
}

/// Row heights (in pixels) for the two display modes.
pub const SQUISHED_ROW_HEIGHT: f32 = 3.0;
pub const EXPANDED_ROW_HEIGHT: f32 = 10.0;
/// Gap between rows in pixels.
pub const ROW_SPACING: f32 = 1.0;

/// Map an HP (haplotype) tag value to a display color.
///
/// - HP 1 → green
/// - HP 2 → orange
/// - None / other → grey (ambiguous)
#[inline]
pub fn hp_color(haplotype: Option<u8>) -> [u8; 3] {
    match haplotype {
        Some(1) => [100, 200, 100], // green
        Some(2) => [230, 160, 60],  // orange
        _ => [160, 160, 160],       // grey (ambiguous / unphased)
    }
}

/// Return the visible indels for a read given the current display config.
///
/// If `config.hide_small_indels` is true, indels with length ≤ threshold
/// are filtered out. Otherwise all indels are returned.
pub fn visible_indels(read: &AlignedRead, config: &PileupDisplayConfig) -> Vec<super::Indel> {
    if config.hide_small_indels {
        read.indels
            .iter()
            .filter(|indel| indel.length > config.indel_threshold)
            .cloned()
            .collect()
    } else {
        read.indels.clone()
    }
}

/// Pack a set of aligned reads into pileup rows using a greedy
/// interval-scheduling algorithm with a min-heap for O(n log n) performance.
///
/// Reads are sorted by start position, then each read is assigned to
/// the row whose last read ended earliest (if the gap constraint is met).
/// If no existing row has room, a new row is created.
///
/// Complexity: O(n log n) sort + O(n log R) packing where R = row count.
pub fn pack_reads(mut reads: Vec<AlignedRead>) -> Vec<PileupRow> {
    if reads.is_empty() {
        return Vec::new();
    }

    reads.sort_by_key(|r| r.start);

    // Min-heap of (end_position, row_index) — earliest-ending row on top.
    let mut heap: BinaryHeap<Reverse<(u64, usize)>> = BinaryHeap::new();
    let mut rows: Vec<Vec<AlignedRead>> = Vec::new();

    for read in reads {
        if let Some(&Reverse((end, idx))) = heap.peek()
            && read.start > end + READ_GAP
        {
            // Reuse the row with the earliest end.
            heap.pop();
            heap.push(Reverse((read.end, idx)));
            rows[idx].push(read);
            continue;
        }
        // No row available — create a new one.
        let idx = rows.len();
        heap.push(Reverse((read.end, idx)));
        rows.push(vec![read]);
    }

    rows.into_iter()
        .enumerate()
        .map(|(i, reads)| PileupRow {
            y_offset: i as u32,
            reads,
        })
        .collect()
}

/// Pack reads grouped by haplotype: HP=1 first, then HP=2, then unphased.
///
/// Within each haplotype group, reads are packed using the standard greedy
/// interval-scheduling algorithm.  A separator row is inserted between groups
/// when the group is non-empty.
pub fn pack_reads_by_haplotype(reads: Vec<AlignedRead>) -> Vec<PileupRow> {
    if reads.is_empty() {
        return Vec::new();
    }

    let mut hap1 = Vec::new();
    let mut hap2 = Vec::new();
    let mut unphased = Vec::new();

    for read in reads {
        match read.haplotype {
            Some(1) => hap1.push(read),
            Some(2) => hap2.push(read),
            _ => unphased.push(read),
        }
    }

    let mut all_rows: Vec<PileupRow> = Vec::new();
    let mut y = 0u32;

    for group in [hap1, hap2, unphased] {
        if group.is_empty() {
            continue;
        }
        if y > 0 {
            // Visual separator: skip a y-offset value so the renderer leaves
            // a blank row gap between haplotype groups.
            y += 1;
        }
        let packed = pack_reads(group);
        for row in packed {
            all_rows.push(PileupRow {
                y_offset: y,
                reads: row.reads,
            });
            y += 1;
        }
    }

    all_rows
}

/// Map a nucleotide base to a display color following standard genomics conventions.
///
/// - A → green
/// - C → blue
/// - G → orange/yellow
/// - T → red
/// - Other → grey
#[inline]
pub fn nucleotide_color(base: u8) -> [u8; 3] {
    match base.to_ascii_uppercase() {
        b'A' => [0, 180, 0],    // green
        b'C' => [0, 0, 200],    // blue
        b'G' => [209, 159, 0],  // orange/yellow
        b'T' => [200, 0, 0],    // red
        _ => [128, 128, 128],   // grey (N or unknown)
    }
}

/// Compute the pixel rectangles for all reads in the pileup.
///
/// Returns `(x, y, width, height, color_rgb)` tuples suitable for
/// rendering with any immediate-mode graphics API (e.g. egui painter).
///
/// `view_start` / `view_end`: the genomic coordinate range visible in
/// the panel. `panel_width`: the pixel width of the drawing area.
pub fn layout_read_rects(
    rows: &[PileupRow],
    config: &PileupDisplayConfig,
    view_start: u64,
    view_end: u64,
    panel_width: f32,
) -> Vec<ReadRect> {
    let bp_span = (view_end.saturating_sub(view_start)).max(1) as f32;
    let bp_per_px = panel_width / bp_span;
    let row_h = if config.squished {
        SQUISHED_ROW_HEIGHT
    } else {
        EXPANDED_ROW_HEIGHT
    };

    let mut rects = Vec::new();
    for row in rows {
        let y = row.y_offset as f32 * (row_h + ROW_SPACING);
        for read in &row.reads {
            // Clamp to visible window
            let r_start = read.start.max(view_start);
            let r_end = read.end.min(view_end);
            if r_start > r_end {
                continue;
            }
            let x = (r_start - view_start) as f32 * bp_per_px;
            let w = ((r_end - r_start + 1) as f32 * bp_per_px).max(1.0);
            let color = hp_color(read.haplotype);
            rects.push(ReadRect {
                x,
                y,
                width: w,
                height: row_h,
                color,
                read_name: read.name.clone(),
                tooltip: Some(ReadTooltipInfo::from_read(read)),
            });

            // Draw soft-clip overlays (semi-transparent extensions beyond alignment)
            if config.show_soft_clips {
                for clip in &read.soft_clips {
                    let clip_color = [180, 180, 220]; // light blue-grey
                    if clip.is_leading {
                        // Leading clip: extends to the left of the alignment start
                        let clip_end = clip.ref_pos;
                        let clip_start = clip_end.saturating_sub(clip.length as u64);
                        let cs = clip_start.max(view_start);
                        let ce = clip_end.min(view_end);
                        if cs < ce {
                            let cx = (cs - view_start) as f32 * bp_per_px;
                            let cw = ((ce - cs) as f32 * bp_per_px).max(1.0);
                            rects.push(ReadRect {
                                x: cx,
                                y,
                                width: cw,
                                height: row_h,
                                color: clip_color,
                                read_name: read.name.clone(),
                                tooltip: None,
                            });
                        }
                    } else {
                        // Trailing clip: extends to the right of the alignment end
                        let clip_start = clip.ref_pos;
                        let clip_end = clip_start + clip.length as u64;
                        let cs = clip_start.max(view_start);
                        let ce = clip_end.min(view_end);
                        if cs < ce {
                            let cx = (cs - view_start) as f32 * bp_per_px;
                            let cw = ((ce - cs) as f32 * bp_per_px).max(1.0);
                            rects.push(ReadRect {
                                x: cx,
                                y,
                                width: cw,
                                height: row_h,
                                color: clip_color,
                                read_name: read.name.clone(),
                                tooltip: None,
                            });
                        }
                    }
                }
            }

            // Draw mismatch markers (nucleotide-colored overlays)
            if config.show_mismatches {
                for mm in &read.mismatches {
                    if mm.ref_pos < view_start || mm.ref_pos > view_end {
                        continue;
                    }
                    let mx = (mm.ref_pos - view_start) as f32 * bp_per_px;
                    let mw = bp_per_px.max(1.0);
                    let mc = nucleotide_color(mm.read_base);
                    rects.push(ReadRect {
                        x: mx,
                        y,
                        width: mw,
                        height: row_h,
                        color: mc,
                        read_name: read.name.clone(),
                        tooltip: None,
                    });
                }
            }

            // Draw indel markers
            let vis_indels = visible_indels(read, config);
            for indel in &vis_indels {
                if indel.ref_pos < view_start || indel.ref_pos > view_end {
                    continue;
                }
                let ix = (indel.ref_pos - view_start) as f32 * bp_per_px;
                let iw = match indel.kind {
                    super::IndelKind::Insertion => 2.0_f32.max(bp_per_px),
                    super::IndelKind::Deletion => (indel.length as f32 * bp_per_px).max(1.0),
                };
                let ic = match indel.kind {
                    super::IndelKind::Insertion => [180, 80, 220], // purple
                    super::IndelKind::Deletion => [30, 30, 30],    // dark
                };
                rects.push(ReadRect {
                    x: ix,
                    y,
                    width: iw,
                    height: row_h,
                    color: ic,
                    read_name: read.name.clone(),
                    tooltip: None,
                });
            }
        }
    }
    rects
}

/// A positioned rectangle for a read or indel marker in the pileup display.
#[derive(Debug, Clone)]
pub struct ReadRect {
    pub x: f32,
    pub y: f32,
    pub width: f32,
    pub height: f32,
    /// RGB color.
    pub color: [u8; 3],
    /// Source read name (for tooltips / identification).
    #[allow(dead_code)]
    pub read_name: String,
    /// Optional rich tooltip metadata (populated for read body rects).
    pub tooltip: Option<ReadTooltipInfo>,
}

/// Rich metadata for hover tooltips on reads.
#[derive(Debug, Clone)]
pub struct ReadTooltipInfo {
    /// Read name (QNAME).
    pub name: String,
    /// Mapping quality.
    pub mapq: Option<u8>,
    /// Haplotype assignment (HP tag).
    pub haplotype: Option<u8>,
    /// True if reverse strand.
    pub is_reverse: bool,
    /// Alignment start (1-based).
    pub start: u64,
    /// Alignment end (1-based, inclusive).
    pub end: u64,
    /// Number of indels in the read.
    pub indel_count: usize,
    /// Alignment length in bases.
    pub alignment_length: u64,
}

impl ReadTooltipInfo {
    /// Build tooltip info from an `AlignedRead`.
    pub fn from_read(read: &super::AlignedRead) -> Self {
        Self {
            name: read.name.clone(),
            mapq: read.mapping_quality,
            haplotype: read.haplotype,
            is_reverse: read.is_reverse,
            start: read.start,
            end: read.end,
            indel_count: read.indels.len(),
            alignment_length: read.end.saturating_sub(read.start) + 1,
        }
    }

    /// Format as multi-line tooltip text.
    pub fn tooltip_text(&self) -> String {
        let strand = if self.is_reverse { "reverse (−)" } else { "forward (+)" };
        let hp = match self.haplotype {
            Some(1) => "HP:1 (hap1)".to_string(),
            Some(2) => "HP:2 (hap2)".to_string(),
            Some(v) => format!("HP:{v}"),
            None => "unphased".to_string(),
        };
        let mapq = self.mapq.map_or("N/A".to_string(), |q| q.to_string());
        format!(
            "{}\n  {}:{}-{} ({} bp)\n  MAPQ: {} | {} | {}\n  Indels: {}",
            self.name,
            strand,
            self.start,
            self.end,
            self.alignment_length,
            mapq,
            hp,
            if self.is_reverse { "rev" } else { "fwd" },
            self.indel_count,
        )
    }
}

/// A single bin in a coverage histogram.
#[derive(Debug, Clone, PartialEq)]
pub struct CoverageBin {
    /// Start position of this bin (genomic coordinate).
    pub start: u64,
    /// End position of this bin (genomic coordinate, exclusive).
    pub end: u64,
    /// Total coverage depth in this bin.
    pub depth: u32,
    /// HP=1 coverage depth.
    pub hap1_depth: u32,
    /// HP=2 coverage depth.
    pub hap2_depth: u32,
}

/// Compute binned coverage from pileup rows.
///
/// `view_start` / `view_end`: the visible genomic coordinate range.
/// `num_bins`: number of bins to divide the view into.
///
/// Returns a vector of `CoverageBin` entries with total and per-HP depth.
pub fn compute_coverage(
    rows: &[PileupRow],
    view_start: u64,
    view_end: u64,
    num_bins: usize,
) -> Vec<CoverageBin> {
    let span = view_end.saturating_sub(view_start);
    if span == 0 || num_bins == 0 {
        return Vec::new();
    }

    let bin_size = (span as f64 / num_bins as f64).ceil() as u64;
    let bin_size = bin_size.max(1);

    let mut bins: Vec<CoverageBin> = (0..num_bins)
        .map(|i| {
            let bin_start = view_start + i as u64 * bin_size;
            let bin_end = (bin_start + bin_size).min(view_end);
            CoverageBin {
                start: bin_start,
                end: bin_end,
                depth: 0,
                hap1_depth: 0,
                hap2_depth: 0,
            }
        })
        .collect();

    for row in rows {
        for read in &row.reads {
            // Determine which bins this read overlaps
            if read.end < view_start || read.start > view_end {
                continue;
            }
            let r_start = read.start.max(view_start);
            let r_end = read.end.min(view_end);

            let first_bin = ((r_start - view_start) / bin_size) as usize;
            let last_bin = (((r_end - view_start).saturating_sub(1)) / bin_size) as usize;

            for bin_idx in first_bin..=last_bin.min(bins.len() - 1) {
                bins[bin_idx].depth += 1;
                match read.haplotype {
                    Some(1) => bins[bin_idx].hap1_depth += 1,
                    Some(2) => bins[bin_idx].hap2_depth += 1,
                    _ => {}
                }
            }
        }
    }

    bins
}

/// Default height for coverage track in pixels.
pub const COVERAGE_TRACK_HEIGHT: f32 = 40.0;

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::genome::{Indel, IndelKind};

    fn make_read(name: &str, start: u64, end: u64) -> AlignedRead {
        AlignedRead {
            name: name.to_string(),
            start,
            end,
            is_reverse: false,
            mapping_quality: Some(60),
            haplotype: None,
            flags: 0,
            indels: Vec::new(),
            mismatches: Vec::new(),
            soft_clips: Vec::new(),
        }
    }

    fn make_read_hp(name: &str, start: u64, end: u64, hp: Option<u8>) -> AlignedRead {
        AlignedRead {
            name: name.to_string(),
            start,
            end,
            is_reverse: false,
            mapping_quality: Some(60),
            haplotype: hp,
            flags: 0,
            indels: Vec::new(),
            mismatches: Vec::new(),
            soft_clips: Vec::new(),
        }
    }

    // -----------------------------------------------------------------------
    // pack_reads tests (existing)
    // -----------------------------------------------------------------------

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

    // -----------------------------------------------------------------------
    // Stress tests for row-packer (10k+ reads)
    // -----------------------------------------------------------------------

    #[test]
    fn test_pack_stress_10k_reads_non_overlapping() {
        // 10,000 non-overlapping reads → should all fit in 1 row
        let reads: Vec<AlignedRead> = (0..10_000)
            .map(|i| {
                let start = i * 1000;
                let end = start + 500;
                make_read(&format!("r{i}"), start, end)
            })
            .collect();
        let rows = pack_reads(reads);
        assert_eq!(
            rows.len(),
            1,
            "10k non-overlapping reads should share 1 row"
        );
        assert_eq!(rows[0].reads.len(), 10_000);
    }

    #[test]
    fn test_pack_stress_10k_reads_fully_overlapping() {
        // 10,000 reads all at the same position → 10,000 rows
        let reads: Vec<AlignedRead> = (0..10_000)
            .map(|i| make_read(&format!("r{i}"), 1000, 2000))
            .collect();
        let rows = pack_reads(reads);
        assert_eq!(rows.len(), 10_000);
        for (i, row) in rows.iter().enumerate() {
            assert_eq!(row.y_offset, i as u32);
            assert_eq!(row.reads.len(), 1);
        }
    }

    #[test]
    fn test_pack_stress_10k_mixed_coverage() {
        // Simulate realistic coverage: reads of ~1000bp across a 50kb region
        // with ~30x coverage, generating many reads.
        let region_size = 50_000u64;
        let read_len = 1_000u64;
        let n_reads = 10_000usize;
        let reads: Vec<AlignedRead> = (0..n_reads)
            .map(|i| {
                let start = (i as u64 * 5) % region_size;
                let end = start + read_len;
                make_read(&format!("r{i}"), start, end)
            })
            .collect();
        let rows = pack_reads(reads);

        // Verify packing invariants
        let total_reads: usize = rows.iter().map(|r| r.reads.len()).sum();
        assert_eq!(total_reads, n_reads, "all reads must be assigned");
        for row in &rows {
            for w in row.reads.windows(2) {
                assert!(
                    w[1].start > w[0].end + READ_GAP,
                    "reads in same row must not overlap (including gap)"
                );
            }
        }
    }

    #[test]
    fn test_pack_stress_50k_reads() {
        // 50,000 reads: ensures the packer scales well
        let reads: Vec<AlignedRead> = (0..50_000)
            .map(|i| {
                let start = (i as u64 * 3) % 100_000;
                let end = start + 500;
                make_read(&format!("r{i}"), start, end)
            })
            .collect();
        let rows = pack_reads(reads);
        let total_reads: usize = rows.iter().map(|r| r.reads.len()).sum();
        assert_eq!(total_reads, 50_000);
        for row in &rows {
            for w in row.reads.windows(2) {
                assert!(w[1].start > w[0].end + READ_GAP);
            }
        }
    }

    // -----------------------------------------------------------------------
    // HP coloring tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_hp_color_hap1_green() {
        let c = hp_color(Some(1));
        // Green channel should be dominant
        assert!(c[1] > c[0], "HP1 green channel should exceed red");
        assert!(c[1] > c[2], "HP1 green channel should exceed blue");
    }

    #[test]
    fn test_hp_color_hap2_orange() {
        let c = hp_color(Some(2));
        // Red channel should be high, blue low → orange
        assert!(c[0] > c[2], "HP2 red should exceed blue");
    }

    #[test]
    fn test_hp_color_ambiguous_grey() {
        let c = hp_color(None);
        // All channels should be equal → grey
        assert_eq!(c[0], c[1]);
        assert_eq!(c[1], c[2]);
    }

    #[test]
    fn test_hp_color_unknown_value_grey() {
        let c = hp_color(Some(3));
        assert_eq!(c, hp_color(None), "unknown HP values should be grey");
    }

    #[test]
    fn test_hp_colors_distinct() {
        let c1 = hp_color(Some(1));
        let c2 = hp_color(Some(2));
        let cg = hp_color(None);
        assert_ne!(c1, c2);
        assert_ne!(c1, cg);
        assert_ne!(c2, cg);
    }

    // -----------------------------------------------------------------------
    // Indel filtering tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_visible_indels_all_shown_when_disabled() {
        let mut read = make_read("r1", 100, 500);
        read.indels = vec![
            Indel {
                ref_pos: 150,
                length: 1,
                kind: IndelKind::Insertion,
            },
            Indel {
                ref_pos: 200,
                length: 2,
                kind: IndelKind::Deletion,
            },
            Indel {
                ref_pos: 300,
                length: 10,
                kind: IndelKind::Insertion,
            },
        ];
        let config = PileupDisplayConfig {
            hide_small_indels: false,
            ..Default::default()
        };
        let vis = visible_indels(&read, &config);
        assert_eq!(vis.len(), 3, "all indels visible when filtering disabled");
    }

    #[test]
    fn test_visible_indels_small_hidden() {
        let mut read = make_read("r1", 100, 500);
        read.indels = vec![
            Indel {
                ref_pos: 150,
                length: 1,
                kind: IndelKind::Insertion,
            },
            Indel {
                ref_pos: 200,
                length: 3,
                kind: IndelKind::Deletion,
            },
            Indel {
                ref_pos: 300,
                length: 4,
                kind: IndelKind::Insertion,
            },
            Indel {
                ref_pos: 400,
                length: 10,
                kind: IndelKind::Deletion,
            },
        ];
        let config = PileupDisplayConfig {
            hide_small_indels: true,
            indel_threshold: 3,
            ..Default::default()
        };
        let vis = visible_indels(&read, &config);
        assert_eq!(vis.len(), 2, "only indels > 3bp should be visible");
        assert_eq!(vis[0].length, 4);
        assert_eq!(vis[1].length, 10);
    }

    #[test]
    fn test_visible_indels_threshold_boundary() {
        let mut read = make_read("r1", 100, 500);
        read.indels = vec![
            Indel {
                ref_pos: 150,
                length: 3,
                kind: IndelKind::Insertion,
            },
            Indel {
                ref_pos: 200,
                length: 4,
                kind: IndelKind::Insertion,
            },
        ];
        let config = PileupDisplayConfig {
            hide_small_indels: true,
            indel_threshold: 3,
            ..Default::default()
        };
        let vis = visible_indels(&read, &config);
        assert_eq!(vis.len(), 1, "length==threshold should be hidden (≤)");
        assert_eq!(vis[0].length, 4);
    }

    #[test]
    fn test_visible_indels_empty_read() {
        let read = make_read("r1", 100, 500);
        let config = PileupDisplayConfig::default();
        let vis = visible_indels(&read, &config);
        assert!(vis.is_empty());
    }

    // -----------------------------------------------------------------------
    // Display config tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_display_config_defaults() {
        let cfg = PileupDisplayConfig::default();
        assert!(cfg.squished, "default display should be squished");
        assert!(cfg.hide_small_indels, "small indels hidden by default");
        assert_eq!(cfg.indel_threshold, DEFAULT_INDEL_THRESHOLD);
    }

    // -----------------------------------------------------------------------
    // layout_read_rects tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_layout_empty() {
        let rows: Vec<PileupRow> = vec![];
        let config = PileupDisplayConfig::default();
        let rects = layout_read_rects(&rows, &config, 100, 200, 800.0);
        assert!(rects.is_empty());
    }

    #[test]
    fn test_layout_single_read() {
        let rows = pack_reads(vec![make_read_hp("r1", 100, 200, Some(1))]);
        let config = PileupDisplayConfig::default();
        let rects = layout_read_rects(&rows, &config, 100, 200, 800.0);
        assert!(!rects.is_empty());
        let r = &rects[0];
        assert!(r.x >= 0.0);
        assert!(r.width > 0.0);
        assert_eq!(r.color, hp_color(Some(1)));
    }

    #[test]
    fn test_layout_hp_colors_correct() {
        let reads = vec![
            make_read_hp("hap1", 100, 200, Some(1)),
            make_read_hp("hap2", 300, 400, Some(2)),
            make_read_hp("none", 500, 600, None),
        ];
        let rows = pack_reads(reads);
        let config = PileupDisplayConfig::default();
        let rects = layout_read_rects(&rows, &config, 100, 600, 1000.0);
        assert_eq!(rects.len(), 3);
        assert_eq!(rects[0].color, hp_color(Some(1)));
        assert_eq!(rects[1].color, hp_color(Some(2)));
        assert_eq!(rects[2].color, hp_color(None));
    }

    #[test]
    fn test_layout_squished_vs_expanded_height() {
        let rows = pack_reads(vec![make_read("r1", 100, 200)]);
        let squished_cfg = PileupDisplayConfig {
            squished: true,
            ..Default::default()
        };
        let expanded_cfg = PileupDisplayConfig {
            squished: false,
            ..Default::default()
        };
        let sq_rects = layout_read_rects(&rows, &squished_cfg, 100, 200, 800.0);
        let ex_rects = layout_read_rects(&rows, &expanded_cfg, 100, 200, 800.0);
        assert!(
            sq_rects[0].height < ex_rects[0].height,
            "squished should be shorter"
        );
        assert_eq!(sq_rects[0].height, SQUISHED_ROW_HEIGHT);
        assert_eq!(ex_rects[0].height, EXPANDED_ROW_HEIGHT);
    }

    #[test]
    fn test_layout_clips_to_view() {
        // Read extends well beyond view
        let rows = pack_reads(vec![make_read("r1", 50, 500)]);
        let config = PileupDisplayConfig::default();
        let rects = layout_read_rects(&rows, &config, 100, 200, 800.0);
        assert_eq!(rects.len(), 1);
        let r = &rects[0];
        // x should start at 0 (clamped to view_start)
        assert!((r.x - 0.0).abs() < 0.01, "should start at panel left edge");
    }

    #[test]
    fn test_layout_read_outside_view() {
        let rows = pack_reads(vec![make_read("r1", 500, 600)]);
        let config = PileupDisplayConfig::default();
        let rects = layout_read_rects(&rows, &config, 100, 200, 800.0);
        assert!(rects.is_empty(), "read outside view should be skipped");
    }

    #[test]
    fn test_layout_indel_markers_shown() {
        let mut read = make_read_hp("r1", 100, 500, Some(1));
        read.indels = vec![Indel {
            ref_pos: 200,
            length: 10,
            kind: IndelKind::Insertion,
        }];
        let rows = pack_reads(vec![read]);
        let config = PileupDisplayConfig {
            hide_small_indels: true,
            indel_threshold: 3,
            ..Default::default()
        };
        let rects = layout_read_rects(&rows, &config, 100, 500, 800.0);
        // Should have 1 read rect + 1 indel marker
        assert_eq!(rects.len(), 2);
    }

    #[test]
    fn test_layout_small_indels_hidden() {
        let mut read = make_read_hp("r1", 100, 500, Some(1));
        read.indels = vec![Indel {
            ref_pos: 200,
            length: 2,
            kind: IndelKind::Insertion,
        }];
        let rows = pack_reads(vec![read]);
        let config = PileupDisplayConfig {
            hide_small_indels: true,
            indel_threshold: 3,
            ..Default::default()
        };
        let rects = layout_read_rects(&rows, &config, 100, 500, 800.0);
        // Should have only the read rect, no indel marker
        assert_eq!(rects.len(), 1);
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

        let reads = match query_bam(&bam, "chr1:112064095-112277008:1-212914") {
            Ok(r) => r,
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

    // -----------------------------------------------------------------------
    // Integration: BAM reads + HP → correct colors in layout
    // -----------------------------------------------------------------------

    #[test]
    fn test_integration_hp_coloring_layout() {
        use crate::genome::bam::query_bam;
        use std::path::PathBuf;

        let bam = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("../resources/toy_dataset/toy_reads.bam");
        if !bam.exists() {
            eprintln!("skipping: toy BAM not found");
            return;
        }

        let reads = match query_bam(&bam, "chr1:112064095-112277008:1-212914") {
            Ok(r) => r,
            Err(e) => panic!("unexpected error: {e}"),
        };

        if reads.is_empty() {
            eprintln!("skipping: no reads in region");
            return;
        }

        // Build expected color map from reads
        let expected_colors: Vec<(String, [u8; 3])> = reads
            .iter()
            .map(|r| (r.name.clone(), hp_color(r.haplotype)))
            .collect();

        let rows = pack_reads(reads);
        let config = PileupDisplayConfig::default();
        let rects = layout_read_rects(&rows, &config, 1, 212914, 1000.0);

        // Every read should appear in the layout with the correct color
        for (name, expected_color) in &expected_colors {
            let matching: Vec<_> = rects
                .iter()
                .filter(|r| &r.read_name == name && r.color == *expected_color)
                .collect();
            assert!(
                !matching.is_empty(),
                "read {name} should have color {expected_color:?} in layout"
            );
        }
    }

    // -----------------------------------------------------------------------
    // Long-read indel threshold (DEFAULT_INDEL_THRESHOLD = 50)
    // -----------------------------------------------------------------------

    #[test]
    fn test_default_indel_threshold_is_50() {
        assert_eq!(
            DEFAULT_INDEL_THRESHOLD, 50,
            "default indel threshold should be 50 bp for long-read data"
        );
    }

    #[test]
    fn test_long_read_indel_threshold_hides_small_indels() {
        let mut read = make_read("r1", 100, 1000);
        read.indels = vec![
            Indel {
                ref_pos: 200,
                length: 5,
                kind: IndelKind::Insertion,
            },
            Indel {
                ref_pos: 400,
                length: 49,
                kind: IndelKind::Deletion,
            },
            Indel {
                ref_pos: 600,
                length: 50,
                kind: IndelKind::Insertion,
            },
            Indel {
                ref_pos: 800,
                length: 100,
                kind: IndelKind::Deletion,
            },
        ];
        let config = PileupDisplayConfig::default(); // threshold = 50
        let vis = visible_indels(&read, &config);
        // Only indels > 50 bp should be visible (length 100)
        assert_eq!(vis.len(), 1, "only indels > 50bp visible with default");
        assert_eq!(vis[0].length, 100);
    }

    // -----------------------------------------------------------------------
    // Haplotype-based read sorting
    // -----------------------------------------------------------------------

    #[test]
    fn test_pack_reads_by_haplotype_empty() {
        let rows = pack_reads_by_haplotype(vec![]);
        assert!(rows.is_empty());
    }

    #[test]
    fn test_pack_reads_by_haplotype_groups_correctly() {
        let reads = vec![
            make_read_hp("u1", 100, 200, None),
            make_read_hp("h2a", 100, 200, Some(2)),
            make_read_hp("h1a", 100, 200, Some(1)),
            make_read_hp("h1b", 300, 400, Some(1)),
            make_read_hp("h2b", 300, 400, Some(2)),
            make_read_hp("u2", 300, 400, None),
        ];
        let rows = pack_reads_by_haplotype(reads);

        // Collect all read names in order
        let all_names: Vec<String> = rows
            .iter()
            .flat_map(|r| r.reads.iter().map(|rd| rd.name.clone()))
            .collect();

        // HP1 reads should come first, then HP2, then unphased
        let h1_pos = all_names.iter().position(|n| n == "h1a").unwrap();
        let h2_pos = all_names.iter().position(|n| n == "h2a").unwrap();
        let u_pos = all_names.iter().position(|n| n == "u1").unwrap();
        assert!(
            h1_pos < h2_pos,
            "HP1 reads should come before HP2 reads"
        );
        assert!(
            h2_pos < u_pos,
            "HP2 reads should come before unphased reads"
        );
    }

    #[test]
    fn test_pack_reads_by_haplotype_separator_rows() {
        // Two groups of non-overlapping reads → each fits in 1 row.
        // With separator, expect: row0 (hp1), gap, row2 (hp2) = y_offsets 0, 2
        let reads = vec![
            make_read_hp("h1", 100, 200, Some(1)),
            make_read_hp("h2", 100, 200, Some(2)),
        ];
        let rows = pack_reads_by_haplotype(reads);
        assert_eq!(rows.len(), 2);
        assert_eq!(rows[0].y_offset, 0);
        assert_eq!(rows[1].y_offset, 2, "separator gap between groups");
    }

    #[test]
    fn test_pack_reads_by_haplotype_single_group() {
        // Only HP1 reads → no separator needed
        let reads = vec![
            make_read_hp("h1a", 100, 200, Some(1)),
            make_read_hp("h1b", 300, 400, Some(1)),
        ];
        let rows = pack_reads_by_haplotype(reads);
        assert_eq!(rows.len(), 1);
        assert_eq!(rows[0].y_offset, 0);
        assert_eq!(rows[0].reads.len(), 2);
    }

    #[test]
    fn test_pack_reads_by_haplotype_preserves_all_reads() {
        let reads = vec![
            make_read_hp("u1", 100, 200, None),
            make_read_hp("h2", 100, 200, Some(2)),
            make_read_hp("h1", 100, 200, Some(1)),
            make_read_hp("u2", 300, 400, Some(3)), // Unknown HP → unphased
        ];
        let total: usize = pack_reads_by_haplotype(reads).iter().map(|r| r.reads.len()).sum();
        assert_eq!(total, 4, "all reads must be preserved");
    }

    #[test]
    fn test_pack_reads_by_haplotype_performance_10k() {
        // Verify O(n log n) scalability with 10k mixed-HP reads
        let reads: Vec<AlignedRead> = (0..10_000)
            .map(|i| {
                let hp = match i % 3 {
                    0 => Some(1),
                    1 => Some(2),
                    _ => None,
                };
                let start = (i as u64 * 5) % 50_000;
                make_read_hp(&format!("r{i}"), start, start + 1000, hp)
            })
            .collect();
        let rows = pack_reads_by_haplotype(reads);
        let total: usize = rows.iter().map(|r| r.reads.len()).sum();
        assert_eq!(total, 10_000);
    }

    // -----------------------------------------------------------------------
    // Nucleotide color tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_nucleotide_colors_distinct() {
        let a = nucleotide_color(b'A');
        let c = nucleotide_color(b'C');
        let g = nucleotide_color(b'G');
        let t = nucleotide_color(b'T');
        assert_ne!(a, c);
        assert_ne!(a, g);
        assert_ne!(a, t);
        assert_ne!(c, g);
        assert_ne!(c, t);
        assert_ne!(g, t);
    }

    #[test]
    fn test_nucleotide_color_case_insensitive() {
        assert_eq!(nucleotide_color(b'a'), nucleotide_color(b'A'));
        assert_eq!(nucleotide_color(b't'), nucleotide_color(b'T'));
    }

    #[test]
    fn test_nucleotide_color_unknown_is_grey() {
        let n = nucleotide_color(b'N');
        assert_eq!(n[0], n[1]);
        assert_eq!(n[1], n[2]);
    }

    // -----------------------------------------------------------------------
    // Soft-clip rendering tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_layout_soft_clips_hidden_by_default() {
        use crate::genome::SoftClip;

        let mut read = make_read("r1", 100, 500);
        read.soft_clips = vec![SoftClip {
            ref_pos: 100,
            length: 50,
            is_leading: true,
        }];
        let rows = pack_reads(vec![read]);
        let config = PileupDisplayConfig::default(); // show_soft_clips = false
        let rects = layout_read_rects(&rows, &config, 50, 600, 1000.0);
        // Should have only 1 rect (the read), no soft-clip overlay
        assert_eq!(rects.len(), 1);
    }

    #[test]
    fn test_layout_soft_clips_shown_when_enabled() {
        use crate::genome::SoftClip;

        let mut read = make_read("r1", 100, 500);
        read.soft_clips = vec![
            SoftClip {
                ref_pos: 100,
                length: 50,
                is_leading: true,
            },
            SoftClip {
                ref_pos: 501,
                length: 30,
                is_leading: false,
            },
        ];
        let rows = pack_reads(vec![read]);
        let config = PileupDisplayConfig {
            show_soft_clips: true,
            ..Default::default()
        };
        let rects = layout_read_rects(&rows, &config, 50, 600, 1000.0);
        // Should have 3 rects: 1 read + 2 soft-clip overlays
        assert_eq!(rects.len(), 3);
        // Soft clips should have the clip color [180, 180, 220]
        assert_eq!(rects[1].color, [180, 180, 220]);
        assert_eq!(rects[2].color, [180, 180, 220]);
    }

    // -----------------------------------------------------------------------
    // Mismatch rendering tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_layout_mismatches_hidden_by_default() {
        use crate::genome::Mismatch;

        let mut read = make_read("r1", 100, 500);
        read.mismatches = vec![Mismatch {
            ref_pos: 200,
            read_base: b'A',
        }];
        let rows = pack_reads(vec![read]);
        let config = PileupDisplayConfig::default(); // show_mismatches = false
        let rects = layout_read_rects(&rows, &config, 100, 500, 800.0);
        assert_eq!(rects.len(), 1, "mismatches should be hidden by default");
    }

    #[test]
    fn test_layout_mismatches_shown_when_enabled() {
        use crate::genome::Mismatch;

        let mut read = make_read("r1", 100, 500);
        read.mismatches = vec![
            Mismatch {
                ref_pos: 200,
                read_base: b'A',
            },
            Mismatch {
                ref_pos: 300,
                read_base: b'T',
            },
        ];
        let rows = pack_reads(vec![read]);
        let config = PileupDisplayConfig {
            show_mismatches: true,
            ..Default::default()
        };
        let rects = layout_read_rects(&rows, &config, 100, 500, 800.0);
        // 1 read rect + 2 mismatch markers
        assert_eq!(rects.len(), 3);
        // Mismatch colors should correspond to nucleotide coloring
        assert_eq!(rects[1].color, nucleotide_color(b'A'));
        assert_eq!(rects[2].color, nucleotide_color(b'T'));
    }

    #[test]
    fn test_layout_mismatches_outside_view_skipped() {
        use crate::genome::Mismatch;

        let mut read = make_read("r1", 100, 500);
        read.mismatches = vec![
            Mismatch {
                ref_pos: 50,
                read_base: b'C',
            }, // Before view
            Mismatch {
                ref_pos: 600,
                read_base: b'G',
            }, // After view
            Mismatch {
                ref_pos: 300,
                read_base: b'T',
            }, // In view
        ];
        let rows = pack_reads(vec![read]);
        let config = PileupDisplayConfig {
            show_mismatches: true,
            ..Default::default()
        };
        let rects = layout_read_rects(&rows, &config, 100, 500, 800.0);
        // 1 read rect + 1 visible mismatch (only ref_pos=300 is in view)
        assert_eq!(rects.len(), 2);
    }

    // -----------------------------------------------------------------------
    // Display config new fields tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_display_config_new_defaults() {
        let cfg = PileupDisplayConfig::default();
        assert!(!cfg.sort_by_haplotype, "HP sort should be off by default");
        assert!(!cfg.show_soft_clips, "soft clips should be hidden by default");
        assert!(!cfg.show_mismatches, "mismatches should be hidden by default");
    }

    // -----------------------------------------------------------------------
    // Min-heap packer correctness under various conditions
    // -----------------------------------------------------------------------

    #[test]
    fn test_pack_heap_optimal_row_count() {
        // 3 reads all at same position → needs 3 rows (optimal)
        let reads = vec![
            make_read("a", 100, 200),
            make_read("b", 100, 200),
            make_read("c", 100, 200),
        ];
        let rows = pack_reads(reads);
        assert_eq!(rows.len(), 3, "heap packer should use optimal row count");
    }

    #[test]
    fn test_pack_heap_reuses_earliest_ending_row() {
        // r1: 100-200, r2: 100-500, r3: 210-300
        // r3 should reuse r1's row (ends at 200) not r2's (ends at 500)
        let reads = vec![
            make_read("r1", 100, 200),
            make_read("r2", 100, 500),
            make_read("r3", 210, 300),
        ];
        let rows = pack_reads(reads);
        assert_eq!(rows.len(), 2, "should only need 2 rows");
        // r3 should be in same row as r1 (the earlier-ending row)
        let row_with_r1 = rows.iter().find(|r| r.reads.iter().any(|rd| rd.name == "r1")).unwrap();
        assert!(
            row_with_r1.reads.iter().any(|rd| rd.name == "r3"),
            "r3 should share row with r1 (earlier end)"
        );
    }

    #[test]
    fn test_pack_heap_stress_50k_reads() {
        let reads: Vec<AlignedRead> = (0..50_000)
            .map(|i| {
                let start = (i as u64 * 3) % 100_000;
                let end = start + 500;
                make_read(&format!("r{i}"), start, end)
            })
            .collect();
        let rows = pack_reads(reads);
        let total: usize = rows.iter().map(|r| r.reads.len()).sum();
        assert_eq!(total, 50_000);
        // Verify packing invariant
        for row in &rows {
            for w in row.reads.windows(2) {
                assert!(w[1].start > w[0].end + READ_GAP);
            }
        }
    }

    // -----------------------------------------------------------------------
    // Coverage computation tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_compute_coverage_empty() {
        let bins = compute_coverage(&[], 100, 200, 10);
        assert!(bins.is_empty() || bins.iter().all(|b| b.depth == 0));
    }

    #[test]
    fn test_compute_coverage_zero_span() {
        let rows = pack_reads(vec![make_read("r1", 100, 200)]);
        let bins = compute_coverage(&rows, 100, 100, 10);
        assert!(bins.is_empty());
    }

    #[test]
    fn test_compute_coverage_single_read() {
        let rows = pack_reads(vec![make_read("r1", 100, 200)]);
        let bins = compute_coverage(&rows, 100, 200, 10);
        assert_eq!(bins.len(), 10);
        // All bins should have depth 1 since the read spans the entire view
        for bin in &bins {
            assert_eq!(bin.depth, 1, "bin {:?} should have depth 1", bin);
        }
    }

    #[test]
    fn test_compute_coverage_hp_stratified() {
        let reads = vec![
            make_read_hp("h1a", 100, 200, Some(1)),
            make_read_hp("h2a", 100, 200, Some(2)),
            make_read_hp("u1", 100, 200, None),
        ];
        let rows = pack_reads(reads);
        let bins = compute_coverage(&rows, 100, 200, 5);
        for bin in &bins {
            assert_eq!(bin.depth, 3);
            assert_eq!(bin.hap1_depth, 1);
            assert_eq!(bin.hap2_depth, 1);
        }
    }

    #[test]
    fn test_compute_coverage_partial_overlap() {
        // Read covers only left half of view
        let rows = pack_reads(vec![make_read("r1", 100, 150)]);
        let bins = compute_coverage(&rows, 100, 200, 10);
        // Bins in first half should have depth 1, second half depth 0
        let nonzero: Vec<_> = bins.iter().filter(|b| b.depth > 0).collect();
        let zero: Vec<_> = bins.iter().filter(|b| b.depth == 0).collect();
        assert!(!nonzero.is_empty());
        assert!(!zero.is_empty());
    }

    #[test]
    fn test_compute_coverage_read_outside_view() {
        let rows = pack_reads(vec![make_read("r1", 500, 600)]);
        let bins = compute_coverage(&rows, 100, 200, 10);
        for bin in &bins {
            assert_eq!(bin.depth, 0);
        }
    }

    // -----------------------------------------------------------------------
    // Tooltip info tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_tooltip_info_from_read() {
        let mut read = make_read_hp("test_read", 100, 500, Some(1));
        read.is_reverse = true;
        read.mapping_quality = Some(42);
        read.indels = vec![
            Indel { ref_pos: 200, length: 10, kind: IndelKind::Insertion },
            Indel { ref_pos: 300, length: 20, kind: IndelKind::Deletion },
        ];
        let info = ReadTooltipInfo::from_read(&read);
        assert_eq!(info.name, "test_read");
        assert_eq!(info.mapq, Some(42));
        assert_eq!(info.haplotype, Some(1));
        assert!(info.is_reverse);
        assert_eq!(info.start, 100);
        assert_eq!(info.end, 500);
        assert_eq!(info.indel_count, 2);
        assert_eq!(info.alignment_length, 401);
    }

    #[test]
    fn test_tooltip_text_content() {
        let read = make_read_hp("my_read", 1000, 2000, Some(2));
        let info = ReadTooltipInfo::from_read(&read);
        let text = info.tooltip_text();
        assert!(text.contains("my_read"));
        assert!(text.contains("HP:2"));
        assert!(text.contains("forward"));
        assert!(text.contains("1000"));
        assert!(text.contains("2000"));
    }

    #[test]
    fn test_read_rect_has_tooltip_on_body() {
        let rows = pack_reads(vec![make_read("r1", 100, 200)]);
        let config = PileupDisplayConfig::default();
        let rects = layout_read_rects(&rows, &config, 100, 200, 800.0);
        assert_eq!(rects.len(), 1);
        assert!(rects[0].tooltip.is_some(), "body rect should have tooltip");
    }

    // -----------------------------------------------------------------------
    // Display config coverage field test
    // -----------------------------------------------------------------------

    #[test]
    fn test_display_config_coverage_default_off() {
        let cfg = PileupDisplayConfig::default();
        assert!(!cfg.show_coverage, "coverage should be off by default");
    }
}
