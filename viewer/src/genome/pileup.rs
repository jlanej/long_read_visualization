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
pub const DEFAULT_INDEL_THRESHOLD: u32 = 3;

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
}

impl Default for PileupDisplayConfig {
    fn default() -> Self {
        Self {
            squished: true,
            hide_small_indels: true,
            indel_threshold: DEFAULT_INDEL_THRESHOLD,
        }
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
            });

            // Draw indel markers
            let vis_indels = visible_indels(read, config);
            for indel in &vis_indels {
                if indel.ref_pos < view_start || indel.ref_pos > view_end {
                    continue;
                }
                let ix = (indel.ref_pos - view_start) as f32 * bp_per_px;
                let iw = match indel.kind {
                    super::IndelKind::Insertion => 2.0_f32.max(bp_per_px),
                    super::IndelKind::Deletion => {
                        (indel.length as f32 * bp_per_px).max(1.0)
                    }
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
    pub read_name: String,
}

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
        assert_eq!(rows.len(), 1, "10k non-overlapping reads should share 1 row");
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
            Indel { ref_pos: 150, length: 1, kind: IndelKind::Insertion },
            Indel { ref_pos: 200, length: 2, kind: IndelKind::Deletion },
            Indel { ref_pos: 300, length: 10, kind: IndelKind::Insertion },
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
            Indel { ref_pos: 150, length: 1, kind: IndelKind::Insertion },
            Indel { ref_pos: 200, length: 3, kind: IndelKind::Deletion },
            Indel { ref_pos: 300, length: 4, kind: IndelKind::Insertion },
            Indel { ref_pos: 400, length: 10, kind: IndelKind::Deletion },
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
            Indel { ref_pos: 150, length: 3, kind: IndelKind::Insertion },
            Indel { ref_pos: 200, length: 4, kind: IndelKind::Insertion },
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
        let squished_cfg = PileupDisplayConfig { squished: true, ..Default::default() };
        let expanded_cfg = PileupDisplayConfig { squished: false, ..Default::default() };
        let sq_rects = layout_read_rects(&rows, &squished_cfg, 100, 200, 800.0);
        let ex_rects = layout_read_rects(&rows, &expanded_cfg, 100, 200, 800.0);
        assert!(sq_rects[0].height < ex_rects[0].height, "squished should be shorter");
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
        read.indels = vec![
            Indel { ref_pos: 200, length: 10, kind: IndelKind::Insertion },
        ];
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
        read.indels = vec![
            Indel { ref_pos: 200, length: 2, kind: IndelKind::Insertion },
        ];
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

        let reads = match query_bam(&bam, "chr1:112164095-112177007") {
            Ok(r) => r,
            Err(crate::genome::GenomeError::ParseError(_)) => {
                eprintln!("skipping: header parse error");
                return;
            }
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
        let rects = layout_read_rects(&rows, &config, 112164095, 112177007, 1000.0);

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
}
