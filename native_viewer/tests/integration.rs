//! Integration tests for multi-panel synchronization and viewer state.

use lrv_native_viewer::app::ViewerState;
use lrv_native_viewer::config::AppConfig;
use lrv_native_viewer::coordinate::CoordinateTranslator;
use lrv_native_viewer::filters::DisplayFilters;
use lrv_native_viewer::genome::bam::AlignedRead;
use lrv_native_viewer::genome::region::GenomicRegion;
use lrv_native_viewer::pileup::layout::layout_reads;
use noodles::sam::alignment::record::cigar::op::Kind as CigarOpKind;
use std::path::Path;

/// Helper to create a minimal ViewerState for testing.
fn make_viewer_state() -> ViewerState {
    let config = AppConfig::from_cli(None, None, None);
    ViewerState::new(config, None, 3, true)
}

/// Helper to make a test read.
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

// ── Navigation Tests ──────────────────────────────────────────────────────

#[test]
fn test_region_navigation_forward() {
    let manifest_path =
        Path::new(env!("CARGO_MANIFEST_DIR")).join("../resources/toy_dataset/toy_manifest.json");
    if !manifest_path.exists() {
        eprintln!("Skipping: toy_manifest.json not found");
        return;
    }

    let config = AppConfig::from_cli(None, None, Some(manifest_path.to_str().unwrap()));
    let mut state = ViewerState::new(config, None, 3, true);

    assert!(state.config.region_count() > 0, "Should have regions");
    assert!(state.current_region_idx.is_none());

    // Navigate to first region
    state.next_region();
    assert_eq!(state.current_region_idx, Some(0));
    assert!(state.ref_region.is_some());

    let first_region = state.ref_region.clone().unwrap();

    // Navigate forward
    state.next_region();
    assert_eq!(state.current_region_idx, Some(1));

    let second_region = state.ref_region.clone().unwrap();
    assert_ne!(
        format!("{first_region}"),
        format!("{second_region}"),
        "Regions should differ"
    );
}

#[test]
fn test_region_navigation_backward() {
    let manifest_path =
        Path::new(env!("CARGO_MANIFEST_DIR")).join("../resources/toy_dataset/toy_manifest.json");
    if !manifest_path.exists() {
        eprintln!("Skipping: toy_manifest.json not found");
        return;
    }

    let config = AppConfig::from_cli(None, None, Some(manifest_path.to_str().unwrap()));
    let mut state = ViewerState::new(config, None, 3, true);

    // Go to first region
    state.next_region();
    assert_eq!(state.current_region_idx, Some(0));

    // Go backward should wrap to last
    state.prev_region();
    assert_eq!(
        state.current_region_idx,
        Some(state.config.region_count() - 1)
    );
}

#[test]
fn test_region_navigation_wraps() {
    let manifest_path =
        Path::new(env!("CARGO_MANIFEST_DIR")).join("../resources/toy_dataset/toy_manifest.json");
    if !manifest_path.exists() {
        eprintln!("Skipping: toy_manifest.json not found");
        return;
    }

    let config = AppConfig::from_cli(None, None, Some(manifest_path.to_str().unwrap()));
    let mut state = ViewerState::new(config, None, 3, true);
    let count = state.config.region_count();

    // Navigate through all regions
    for i in 0..count {
        state.next_region();
        assert_eq!(state.current_region_idx, Some(i));
    }

    // Next should wrap to 0
    state.next_region();
    assert_eq!(state.current_region_idx, Some(0));
}

// ── Pan/Zoom Tests ────────────────────────────────────────────────────────

#[test]
fn test_pan_updates_region() {
    let mut state = make_viewer_state();
    state.ref_region = Some(GenomicRegion::new("chr1", 1000, 2000));

    state.pan(0.5);
    let region = state.ref_region.as_ref().unwrap();
    assert_eq!(region.start, 1500);
    assert_eq!(region.end, 2500);
    assert_eq!(region.span(), 1000); // span preserved
}

#[test]
fn test_zoom_preserves_center() {
    let mut state = make_viewer_state();
    state.ref_region = Some(GenomicRegion::new("chr1", 1000, 2000));
    let original_center = state.ref_region.as_ref().unwrap().center();

    state.zoom(2.0); // zoom out
    let region = state.ref_region.as_ref().unwrap();
    assert_eq!(region.center(), original_center);
    assert_eq!(region.span(), 2000); // doubled span
}

#[test]
fn test_pan_clears_pileup_data() {
    let mut state = make_viewer_state();
    state.ref_region = Some(GenomicRegion::new("chr1", 1000, 2000));

    // Simulate having pileup data
    state.ref_pileup = Some(lrv_native_viewer::pileup::PileupData {
        rows: vec![],
        total_reads: 0,
        region_start: 1000,
        region_end: 2000,
    });

    state.pan(0.1);
    assert!(
        state.ref_pileup.is_none(),
        "Pan should clear pileup to trigger reload"
    );
}

// ── Display Filter Tests ──────────────────────────────────────────────────

#[test]
fn test_display_filter_toggle() {
    let mut state = make_viewer_state();

    assert!(state.filters.squished);
    state.filters.squished = false;
    assert!(!state.filters.squished);

    assert!(state.filters.hide_small_indels);
    state.filters.hide_small_indels = false;
    assert!(!state.filters.hide_small_indels);
}

#[test]
fn test_indel_threshold_filtering() {
    let filters = DisplayFilters::new(3, true);

    // Indels at or below threshold hidden
    assert!(!filters.show_indel(1));
    assert!(!filters.show_indel(2));
    assert!(!filters.show_indel(3));

    // Indels above threshold shown
    assert!(filters.show_indel(4));
    assert!(filters.show_indel(50));
    assert!(filters.show_indel(100));
}

// ── Pileup Layout Integration Tests ──────────────────────────────────────

#[test]
fn test_pileup_hp_grouping() {
    let reads = vec![
        make_read("r1", 100, 300, 1),
        make_read("r2", 100, 300, 2),
        make_read("r3", 100, 300, 1),
        make_read("r4", 100, 300, 0),
        make_read("r5", 100, 300, 2),
    ];

    let pileup = layout_reads(reads, 0, 1000, true);
    assert_eq!(pileup.total_reads, 5);

    // With HP sorting, reads should be grouped
    // All reads overlap so each gets its own row
    assert_eq!(pileup.rows.len(), 5);

    // Check ordering: HP0, HP1, HP1, HP2, HP2
    let hp_order: Vec<u8> = pileup.rows.iter().map(|r| r.reads[0].hp_tag).collect();
    assert_eq!(hp_order, vec![0, 1, 1, 2, 2]);
}

#[test]
fn test_pileup_large_dataset() {
    // Simulate a realistic pileup with many reads
    let reads: Vec<_> = (0..1000)
        .map(|i| {
            let start = (i * 100) % 50_000;
            let end = start + 5000 + (i % 3) * 1000;
            let hp = (i % 3) as u8;
            make_read(&format!("r{i}"), start, end, hp)
        })
        .collect();

    let pileup = layout_reads(reads, 0, 100_000, true);
    assert_eq!(pileup.total_reads, 1000);
    assert!(
        pileup.rows.len() < 1000,
        "Layout should pack reads efficiently"
    );
}

// ── Coordinate Translation Tests ──────────────────────────────────────────

#[test]
fn test_translator_empty() {
    let translator = CoordinateTranslator::new();
    let result = translator.translate(
        "sample1",
        &GenomicRegion::new("chr1", 1000, 2000),
        0,
    );
    assert!(result.hap1.is_empty());
    assert!(result.hap2.is_empty());
}

#[test]
fn test_translator_has_index() {
    let translator = CoordinateTranslator::new();
    assert!(!translator.has_index("sample1", "hap1"));
    assert!(!translator.has_index("sample1", "hap2"));
}

// ── Config Loading Tests ──────────────────────────────────────────────────

#[test]
fn test_manifest_loading_toy_dataset() {
    let manifest_path =
        Path::new(env!("CARGO_MANIFEST_DIR")).join("../resources/toy_dataset/toy_manifest.json");
    if !manifest_path.exists() {
        eprintln!("Skipping: toy_manifest.json not found");
        return;
    }

    let config = AppConfig::from_cli(None, None, Some(manifest_path.to_str().unwrap()));
    assert_eq!(
        config.regions.len(),
        10,
        "Toy manifest should have 10 regions"
    );

    // Verify first region
    let first = &config.regions[0];
    assert!(first.chrom.starts_with("chr"), "Region should have chr prefix");
    assert!(first.span() > 0, "Region should have non-zero span");
}

#[test]
fn test_config_from_tsv() {
    let dir = tempfile::tempdir().unwrap();
    let tsv_path = dir.path().join("config.tsv");
    std::fs::write(
        &tsv_path,
        "# sample_id\toutput_dir\treference\n\
         sample1\t/tmp/output\t/tmp/ref.fa.gz\n",
    )
    .unwrap();

    let config = AppConfig::from_tsv(tsv_path.to_str().unwrap()).unwrap();
    assert_eq!(config.samples.len(), 1);
    assert_eq!(config.samples[0].sample_id, "sample1");
}

// ── State Sync Tests ──────────────────────────────────────────────────────

#[test]
fn test_navigate_to_region_clears_all_pileups() {
    let manifest_path =
        Path::new(env!("CARGO_MANIFEST_DIR")).join("../resources/toy_dataset/toy_manifest.json");
    if !manifest_path.exists() {
        eprintln!("Skipping: toy_manifest.json not found");
        return;
    }

    let config = AppConfig::from_cli(None, None, Some(manifest_path.to_str().unwrap()));
    let mut state = ViewerState::new(config, None, 3, true);

    // Set fake pileup data for all three panels
    let fake_pileup = lrv_native_viewer::pileup::PileupData {
        rows: vec![],
        total_reads: 42,
        region_start: 0,
        region_end: 1000,
    };
    state.ref_pileup = Some(fake_pileup.clone());
    state.hap1_pileup = Some(fake_pileup.clone());
    state.hap2_pileup = Some(fake_pileup);

    // Navigate to a region
    state.navigate_to_region(0);

    // All pileups should be cleared (triggers reload)
    assert!(state.ref_pileup.is_none());
    assert!(state.hap1_pileup.is_none());
    assert!(state.hap2_pileup.is_none());
}

#[test]
fn test_zoom_clears_all_pileups() {
    let mut state = make_viewer_state();
    state.ref_region = Some(GenomicRegion::new("chr1", 1000, 2000));

    let fake_pileup = lrv_native_viewer::pileup::PileupData {
        rows: vec![],
        total_reads: 0,
        region_start: 1000,
        region_end: 2000,
    };
    state.ref_pileup = Some(fake_pileup.clone());
    state.hap1_pileup = Some(fake_pileup.clone());
    state.hap2_pileup = Some(fake_pileup);

    state.zoom(2.0);

    assert!(state.ref_pileup.is_none());
    assert!(state.hap1_pileup.is_none());
    assert!(state.hap2_pileup.is_none());
}
