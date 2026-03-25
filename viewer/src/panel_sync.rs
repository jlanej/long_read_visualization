//! Synchronized per-panel view state management.
//!
//! Manages independent view ranges (pan/zoom) for the three panels
//! (Reference, Haplotype 1, Haplotype 2) and provides synchronized
//! updates: when navigation occurs in one panel the others can be
//! updated via coordinate translation.

use crate::region::GenomicRegion;

// ---------------------------------------------------------------------------
// Types
// ---------------------------------------------------------------------------

/// Identifies one of the three genome browser panels.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum PanelId {
    Reference,
    Haplotype1,
    Haplotype2,
}

impl PanelId {
    /// All panels in display order.
    pub const ALL: [PanelId; 3] = [PanelId::Reference, PanelId::Haplotype1, PanelId::Haplotype2];
}

/// View state for a single panel.
#[derive(Debug, Clone, PartialEq)]
pub struct PanelView {
    /// Current genomic region displayed in this panel.
    pub region: Option<GenomicRegion>,
    /// View start (may be panned/zoomed relative to the region).
    pub view_start: u64,
    /// View end (may be panned/zoomed relative to the region).
    pub view_end: u64,
}

impl Default for PanelView {
    fn default() -> Self {
        Self {
            region: None,
            view_start: 0,
            view_end: 1000,
        }
    }
}

impl PanelView {
    /// Current view span in base pairs.
    pub fn span(&self) -> u64 {
        self.view_end.saturating_sub(self.view_start)
    }

    /// Midpoint of the current view.
    pub fn center(&self) -> u64 {
        self.view_start + self.span() / 2
    }
}

/// Manages synchronized view state for all three panels.
#[derive(Debug, Clone)]
pub struct PanelSyncManager {
    panels: [PanelView; 3],
    /// When `true`, navigation in any panel propagates to the others.
    pub sync_enabled: bool,
}

impl Default for PanelSyncManager {
    fn default() -> Self {
        Self {
            panels: [
                PanelView::default(),
                PanelView::default(),
                PanelView::default(),
            ],
            sync_enabled: true,
        }
    }
}

impl PanelSyncManager {
    pub fn new() -> Self {
        Self::default()
    }

    // -- Accessors ----------------------------------------------------------

    fn idx(panel: PanelId) -> usize {
        match panel {
            PanelId::Reference => 0,
            PanelId::Haplotype1 => 1,
            PanelId::Haplotype2 => 2,
        }
    }

    /// Get the view state for a panel.
    pub fn view(&self, panel: PanelId) -> &PanelView {
        &self.panels[Self::idx(panel)]
    }

    /// Get a mutable reference to the view state for a panel.
    #[allow(dead_code)] // Public API used in tests; available for direct panel manipulation.
    pub fn view_mut(&mut self, panel: PanelId) -> &mut PanelView {
        &mut self.panels[Self::idx(panel)]
    }

    // -- Region loading (from navigator) ------------------------------------

    /// Set all panels to their respective regions from the manifest entry.
    /// This is the typical path when stepping through regions with the
    /// region navigator.
    pub fn set_regions(
        &mut self,
        ref_region: &GenomicRegion,
        hap1_region: Option<&GenomicRegion>,
        hap2_region: Option<&GenomicRegion>,
    ) {
        // Reference panel
        let pv = &mut self.panels[0];
        pv.view_start = ref_region.start;
        pv.view_end = ref_region.end;
        pv.region = Some(ref_region.clone());

        // Hap1 panel
        let pv = &mut self.panels[1];
        if let Some(r) = hap1_region {
            pv.view_start = r.start;
            pv.view_end = r.end;
            pv.region = Some(r.clone());
        } else {
            // Mirror the reference region when no hap region is available
            pv.view_start = ref_region.start;
            pv.view_end = ref_region.end;
            pv.region = None;
        }

        // Hap2 panel
        let pv = &mut self.panels[2];
        if let Some(r) = hap2_region {
            pv.view_start = r.start;
            pv.view_end = r.end;
            pv.region = Some(r.clone());
        } else {
            pv.view_start = ref_region.start;
            pv.view_end = ref_region.end;
            pv.region = None;
        }
    }

    // -- Pan / Zoom ---------------------------------------------------------

    /// Zoom into the given panel, keeping the view centered.
    /// `factor` > 1.0 zooms in (narrower view), < 1.0 zooms out.
    /// The minimum span is 10 bp.
    pub fn zoom(&mut self, panel: PanelId, factor: f64) {
        let idx = Self::idx(panel);
        let pv = &self.panels[idx];
        let center = pv.center();
        let span = pv.span();
        let new_span = ((span as f64 / factor).round() as u64).max(10);
        let half = new_span / 2;
        let new_start = center.saturating_sub(half);
        let new_end = new_start + new_span;
        self.panels[idx].view_start = new_start;
        self.panels[idx].view_end = new_end;

        if self.sync_enabled {
            self.propagate_zoom(panel, factor);
        }
    }

    /// Pan the given panel by `delta` base pairs (positive = right).
    pub fn pan(&mut self, panel: PanelId, delta: i64) {
        let idx = Self::idx(panel);
        let pv = &mut self.panels[idx];
        if delta >= 0 {
            pv.view_start = pv.view_start.saturating_add(delta as u64);
            pv.view_end = pv.view_end.saturating_add(delta as u64);
        } else {
            let abs = (-delta) as u64;
            pv.view_start = pv.view_start.saturating_sub(abs);
            pv.view_end = pv.view_end.saturating_sub(abs);
        }

        if self.sync_enabled {
            self.propagate_pan(panel, delta);
        }
    }

    /// Apply matching zoom to all *other* panels.
    fn propagate_zoom(&mut self, source: PanelId, factor: f64) {
        for &target in &PanelId::ALL {
            if target == source {
                continue;
            }
            let idx = Self::idx(target);
            let pv = &self.panels[idx];
            let center = pv.center();
            let span = pv.span();
            let new_span = ((span as f64 / factor).round() as u64).max(10);
            let half = new_span / 2;
            let new_start = center.saturating_sub(half);
            let new_end = new_start + new_span;
            self.panels[idx].view_start = new_start;
            self.panels[idx].view_end = new_end;
        }
    }

    /// Apply matching pan delta to all *other* panels.
    fn propagate_pan(&mut self, source: PanelId, delta: i64) {
        for &target in &PanelId::ALL {
            if target == source {
                continue;
            }
            let idx = Self::idx(target);
            let pv = &mut self.panels[idx];
            if delta >= 0 {
                pv.view_start = pv.view_start.saturating_add(delta as u64);
                pv.view_end = pv.view_end.saturating_add(delta as u64);
            } else {
                let abs = (-delta) as u64;
                pv.view_start = pv.view_start.saturating_sub(abs);
                pv.view_end = pv.view_end.saturating_sub(abs);
            }
        }
    }

    /// Check whether all panels have matching view spans (useful in tests).
    pub fn all_spans_match(&self) -> bool {
        let s0 = self.panels[0].span();
        self.panels[1].span() == s0 && self.panels[2].span() == s0
    }

    /// Return the current highlight region for a panel (its view range as a
    /// `GenomicRegion`).  Returns `None` when the panel has no region set.
    pub fn highlight_region(&self, panel: PanelId) -> Option<GenomicRegion> {
        let pv = self.view(panel);
        pv.region.as_ref().map(|r| GenomicRegion {
            chrom: r.chrom.clone(),
            start: pv.view_start,
            end: pv.view_end,
        })
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    fn make_region(chrom: &str, start: u64, end: u64) -> GenomicRegion {
        GenomicRegion {
            chrom: chrom.to_string(),
            start,
            end,
        }
    }

    // -- Basic construction / defaults --------------------------------------

    #[test]
    fn test_default_sync_enabled() {
        let mgr = PanelSyncManager::new();
        assert!(mgr.sync_enabled);
    }

    #[test]
    fn test_default_views() {
        let mgr = PanelSyncManager::new();
        for &p in &PanelId::ALL {
            assert_eq!(mgr.view(p).view_start, 0);
            assert_eq!(mgr.view(p).view_end, 1000);
            assert!(mgr.view(p).region.is_none());
        }
    }

    // -- set_regions --------------------------------------------------------

    #[test]
    fn test_set_regions_all_present() {
        let mut mgr = PanelSyncManager::new();
        let ref_r = make_region("chr1", 1000, 2000);
        let h1 = make_region("ctg1", 5000, 6000);
        let h2 = make_region("ctg2", 7000, 8000);
        mgr.set_regions(&ref_r, Some(&h1), Some(&h2));

        assert_eq!(mgr.view(PanelId::Reference).view_start, 1000);
        assert_eq!(mgr.view(PanelId::Reference).view_end, 2000);
        assert_eq!(mgr.view(PanelId::Haplotype1).view_start, 5000);
        assert_eq!(mgr.view(PanelId::Haplotype1).view_end, 6000);
        assert_eq!(mgr.view(PanelId::Haplotype2).view_start, 7000);
        assert_eq!(mgr.view(PanelId::Haplotype2).view_end, 8000);
    }

    #[test]
    fn test_set_regions_missing_hap() {
        let mut mgr = PanelSyncManager::new();
        let ref_r = make_region("chr1", 1000, 2000);
        mgr.set_regions(&ref_r, None, None);

        // Hap panels fall back to reference coordinates
        assert_eq!(mgr.view(PanelId::Haplotype1).view_start, 1000);
        assert_eq!(mgr.view(PanelId::Haplotype2).view_start, 1000);
        assert!(mgr.view(PanelId::Haplotype1).region.is_none());
        assert!(mgr.view(PanelId::Haplotype2).region.is_none());
    }

    // -- Navigation: all panels match expected region after region step ------

    #[test]
    fn test_region_step_all_panels_match() {
        let mut mgr = PanelSyncManager::new();
        let ref1 = make_region("chr1", 1000, 2000);
        let h1_1 = make_region("ctg1", 5000, 6000);
        let h2_1 = make_region("ctg2", 7000, 8000);
        mgr.set_regions(&ref1, Some(&h1_1), Some(&h2_1));

        // All panels should show their respective regions
        assert_eq!(mgr.view(PanelId::Reference).span(), 1000);
        assert_eq!(mgr.view(PanelId::Haplotype1).span(), 1000);
        assert_eq!(mgr.view(PanelId::Haplotype2).span(), 1000);

        // Step to region 2
        let ref2 = make_region("chr2", 3000, 5000);
        let h1_2 = make_region("ctg3", 10000, 12000);
        let h2_2 = make_region("ctg4", 15000, 17000);
        mgr.set_regions(&ref2, Some(&h1_2), Some(&h2_2));

        assert_eq!(mgr.view(PanelId::Reference).view_start, 3000);
        assert_eq!(mgr.view(PanelId::Reference).view_end, 5000);
        assert_eq!(mgr.view(PanelId::Haplotype1).view_start, 10000);
        assert_eq!(mgr.view(PanelId::Haplotype1).view_end, 12000);
        assert_eq!(mgr.view(PanelId::Haplotype2).view_start, 15000);
        assert_eq!(mgr.view(PanelId::Haplotype2).view_end, 17000);
    }

    // -- Pan ----------------------------------------------------------------

    #[test]
    fn test_pan_right_synced() {
        let mut mgr = PanelSyncManager::new();
        let ref_r = make_region("chr1", 1000, 2000);
        let h1 = make_region("ctg1", 5000, 6000);
        let h2 = make_region("ctg2", 7000, 8000);
        mgr.set_regions(&ref_r, Some(&h1), Some(&h2));

        mgr.pan(PanelId::Reference, 500);

        // Source panel
        assert_eq!(mgr.view(PanelId::Reference).view_start, 1500);
        assert_eq!(mgr.view(PanelId::Reference).view_end, 2500);
        // Synced panels also moved
        assert_eq!(mgr.view(PanelId::Haplotype1).view_start, 5500);
        assert_eq!(mgr.view(PanelId::Haplotype1).view_end, 6500);
        assert_eq!(mgr.view(PanelId::Haplotype2).view_start, 7500);
        assert_eq!(mgr.view(PanelId::Haplotype2).view_end, 8500);
    }

    #[test]
    fn test_pan_left_synced() {
        let mut mgr = PanelSyncManager::new();
        let ref_r = make_region("chr1", 1000, 2000);
        let h1 = make_region("ctg1", 5000, 6000);
        mgr.set_regions(&ref_r, Some(&h1), None);

        mgr.pan(PanelId::Haplotype1, -200);

        assert_eq!(mgr.view(PanelId::Haplotype1).view_start, 4800);
        assert_eq!(mgr.view(PanelId::Haplotype1).view_end, 5800);
        // Reference panel also moved
        assert_eq!(mgr.view(PanelId::Reference).view_start, 800);
        assert_eq!(mgr.view(PanelId::Reference).view_end, 1800);
    }

    #[test]
    fn test_pan_no_sync() {
        let mut mgr = PanelSyncManager::new();
        mgr.sync_enabled = false;
        let ref_r = make_region("chr1", 1000, 2000);
        let h1 = make_region("ctg1", 5000, 6000);
        mgr.set_regions(&ref_r, Some(&h1), None);

        mgr.pan(PanelId::Reference, 500);

        // Source moved
        assert_eq!(mgr.view(PanelId::Reference).view_start, 1500);
        // Others did NOT move
        assert_eq!(mgr.view(PanelId::Haplotype1).view_start, 5000);
    }

    #[test]
    fn test_pan_saturates_at_zero() {
        let mut mgr = PanelSyncManager::new();
        let ref_r = make_region("chr1", 100, 200);
        mgr.set_regions(&ref_r, None, None);

        mgr.pan(PanelId::Reference, -500);

        assert_eq!(mgr.view(PanelId::Reference).view_start, 0);
    }

    // -- Zoom ---------------------------------------------------------------

    #[test]
    fn test_zoom_in_synced() {
        let mut mgr = PanelSyncManager::new();
        let ref_r = make_region("chr1", 1000, 3000);
        let h1 = make_region("ctg1", 5000, 7000);
        let h2 = make_region("ctg2", 9000, 11000);
        mgr.set_regions(&ref_r, Some(&h1), Some(&h2));

        // Zoom in 2× on reference panel
        mgr.zoom(PanelId::Reference, 2.0);

        let ref_span = mgr.view(PanelId::Reference).span();
        assert_eq!(ref_span, 1000); // 2000 / 2 = 1000
        // All panels should have the same span
        assert!(mgr.all_spans_match());
    }

    #[test]
    fn test_zoom_out_synced() {
        let mut mgr = PanelSyncManager::new();
        let ref_r = make_region("chr1", 1000, 2000);
        let h1 = make_region("ctg1", 5000, 6000);
        mgr.set_regions(&ref_r, Some(&h1), None);

        // Zoom out 2×
        mgr.zoom(PanelId::Reference, 0.5);

        let ref_span = mgr.view(PanelId::Reference).span();
        assert_eq!(ref_span, 2000); // 1000 / 0.5 = 2000
        assert!(mgr.all_spans_match());
    }

    #[test]
    fn test_zoom_minimum_span() {
        let mut mgr = PanelSyncManager::new();
        let ref_r = make_region("chr1", 1000, 1020);
        mgr.set_regions(&ref_r, None, None);

        // Aggressive zoom that would go below 10bp
        mgr.zoom(PanelId::Reference, 100.0);

        assert!(mgr.view(PanelId::Reference).span() >= 10);
    }

    #[test]
    fn test_zoom_no_sync() {
        let mut mgr = PanelSyncManager::new();
        mgr.sync_enabled = false;
        let ref_r = make_region("chr1", 1000, 3000);
        let h1 = make_region("ctg1", 5000, 7000);
        mgr.set_regions(&ref_r, Some(&h1), None);

        mgr.zoom(PanelId::Reference, 2.0);

        assert_eq!(mgr.view(PanelId::Reference).span(), 1000);
        // Hap1 unchanged
        assert_eq!(mgr.view(PanelId::Haplotype1).span(), 2000);
    }

    // -- Cross-navigation: changing region in one panel updates others ------

    #[test]
    fn test_cross_nav_pan_from_hap2_updates_all() {
        let mut mgr = PanelSyncManager::new();
        let ref_r = make_region("chr1", 1000, 2000);
        let h1 = make_region("ctg1", 5000, 6000);
        let h2 = make_region("ctg2", 9000, 10000);
        mgr.set_regions(&ref_r, Some(&h1), Some(&h2));

        // Pan from Hap2
        mgr.pan(PanelId::Haplotype2, 300);

        // Verify all panels moved by the same delta
        assert_eq!(mgr.view(PanelId::Reference).view_start, 1300);
        assert_eq!(mgr.view(PanelId::Reference).view_end, 2300);
        assert_eq!(mgr.view(PanelId::Haplotype1).view_start, 5300);
        assert_eq!(mgr.view(PanelId::Haplotype1).view_end, 6300);
        assert_eq!(mgr.view(PanelId::Haplotype2).view_start, 9300);
        assert_eq!(mgr.view(PanelId::Haplotype2).view_end, 10300);
    }

    #[test]
    fn test_cross_nav_zoom_from_hap1_updates_all() {
        let mut mgr = PanelSyncManager::new();
        let ref_r = make_region("chr1", 1000, 5000);
        let h1 = make_region("ctg1", 10000, 14000);
        let h2 = make_region("ctg2", 20000, 24000);
        mgr.set_regions(&ref_r, Some(&h1), Some(&h2));

        // All start at span 4000
        assert!(mgr.all_spans_match());

        // Zoom in 2× from hap1
        mgr.zoom(PanelId::Haplotype1, 2.0);

        // All panels should now have span 2000
        assert_eq!(mgr.view(PanelId::Reference).span(), 2000);
        assert_eq!(mgr.view(PanelId::Haplotype1).span(), 2000);
        assert_eq!(mgr.view(PanelId::Haplotype2).span(), 2000);
        assert!(mgr.all_spans_match());
    }

    // -- Highlight region ---------------------------------------------------

    #[test]
    fn test_highlight_region_with_data() {
        let mut mgr = PanelSyncManager::new();
        let ref_r = make_region("chr1", 1000, 2000);
        mgr.set_regions(&ref_r, None, None);

        let h = mgr.highlight_region(PanelId::Reference).unwrap();
        assert_eq!(h.chrom, "chr1");
        assert_eq!(h.start, 1000);
        assert_eq!(h.end, 2000);
    }

    #[test]
    fn test_highlight_region_after_pan() {
        let mut mgr = PanelSyncManager::new();
        let ref_r = make_region("chr1", 1000, 2000);
        mgr.set_regions(&ref_r, None, None);

        mgr.pan(PanelId::Reference, 500);

        let h = mgr.highlight_region(PanelId::Reference).unwrap();
        assert_eq!(h.start, 1500);
        assert_eq!(h.end, 2500);
    }

    #[test]
    fn test_highlight_region_none_when_no_region() {
        let mgr = PanelSyncManager::new();
        // Haplotype1 has no region set in default state
        assert!(mgr.highlight_region(PanelId::Haplotype1).is_none());
    }

    // -- Integration: full navigation sequence ------------------------------

    #[test]
    fn test_full_navigation_sequence() {
        let mut mgr = PanelSyncManager::new();

        // Step 1: Load initial region
        let ref1 = make_region("chr1", 10000, 20000);
        let h1_1 = make_region("ctg1", 50000, 60000);
        let h2_1 = make_region("ctg2", 80000, 90000);
        mgr.set_regions(&ref1, Some(&h1_1), Some(&h2_1));
        assert!(mgr.all_spans_match());

        // Step 2: Zoom in
        mgr.zoom(PanelId::Reference, 2.0);
        assert_eq!(mgr.view(PanelId::Reference).span(), 5000);
        assert!(mgr.all_spans_match());

        // Step 3: Pan right
        mgr.pan(PanelId::Reference, 1000);
        let ref_v = mgr.view(PanelId::Reference);
        let expected_center = 15000 + 1000; // original center was 15000
        assert_eq!(ref_v.center(), expected_center);
        assert!(mgr.all_spans_match());

        // Step 4: Step to new region (resets all views)
        let ref2 = make_region("chr2", 30000, 40000);
        let h1_2 = make_region("ctg3", 70000, 80000);
        let h2_2 = make_region("ctg4", 90000, 100000);
        mgr.set_regions(&ref2, Some(&h1_2), Some(&h2_2));

        assert_eq!(mgr.view(PanelId::Reference).view_start, 30000);
        assert_eq!(mgr.view(PanelId::Haplotype1).view_start, 70000);
        assert_eq!(mgr.view(PanelId::Haplotype2).view_start, 90000);
        assert!(mgr.all_spans_match());
    }

    // -- Edge cases ---------------------------------------------------------

    #[test]
    fn test_pan_view_helper_methods() {
        let pv = PanelView {
            region: Some(make_region("chr1", 1000, 3000)),
            view_start: 1000,
            view_end: 3000,
        };
        assert_eq!(pv.span(), 2000);
        assert_eq!(pv.center(), 2000);
    }

    #[test]
    fn test_repeated_zoom_in_out_returns_close() {
        let mut mgr = PanelSyncManager::new();
        let ref_r = make_region("chr1", 10000, 20000);
        mgr.set_regions(&ref_r, None, None);

        let orig_start = mgr.view(PanelId::Reference).view_start;
        let orig_end = mgr.view(PanelId::Reference).view_end;

        // Zoom in then out
        mgr.zoom(PanelId::Reference, 2.0);
        mgr.zoom(PanelId::Reference, 0.5);

        // Should be approximately the same (rounding may differ by 1)
        let delta_start =
            (mgr.view(PanelId::Reference).view_start as i64 - orig_start as i64).unsigned_abs();
        let delta_end =
            (mgr.view(PanelId::Reference).view_end as i64 - orig_end as i64).unsigned_abs();
        assert!(delta_start <= 1, "start drifted by {delta_start}");
        assert!(delta_end <= 1, "end drifted by {delta_end}");
    }
}
