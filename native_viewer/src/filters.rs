/// Display filter settings controlling read rendering.
#[derive(Debug, Clone)]
pub struct DisplayFilters {
    /// Hide indels ≤ this size (bp). Default: 3.
    pub indel_threshold: u32,
    /// Whether small-indel hiding is enabled.
    pub hide_small_indels: bool,
    /// Squished (compact) display mode.
    pub squished: bool,
    /// Show soft-clipped bases.
    pub show_soft_clips: bool,
    /// Show base mismatches (SNVs).
    pub show_mismatches: bool,
}

impl DisplayFilters {
    pub fn new(indel_threshold: u32, squished: bool) -> Self {
        Self {
            indel_threshold,
            hide_small_indels: true,
            squished,
            show_soft_clips: false,
            show_mismatches: true,
        }
    }

    /// Whether a given indel of `size` bp should be displayed.
    pub fn show_indel(&self, size: u32) -> bool {
        if !self.hide_small_indels {
            return true;
        }
        size > self.indel_threshold
    }

    /// Row height in pixels based on display mode.
    pub fn row_height(&self) -> f32 {
        if self.squished { 4.0 } else { 12.0 }
    }

    /// Gap between rows in pixels.
    pub fn row_gap(&self) -> f32 {
        if self.squished { 1.0 } else { 2.0 }
    }
}

impl Default for DisplayFilters {
    fn default() -> Self {
        Self::new(3, true)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_default_filters() {
        let f = DisplayFilters::default();
        assert_eq!(f.indel_threshold, 3);
        assert!(f.hide_small_indels);
        assert!(f.squished);
        assert!(!f.show_soft_clips);
        assert!(f.show_mismatches);
    }

    #[test]
    fn test_show_indel() {
        let f = DisplayFilters::default();
        // Indels ≤ 3bp hidden
        assert!(!f.show_indel(1));
        assert!(!f.show_indel(2));
        assert!(!f.show_indel(3));
        // Indels > 3bp shown
        assert!(f.show_indel(4));
        assert!(f.show_indel(100));
    }

    #[test]
    fn test_show_indel_disabled() {
        let mut f = DisplayFilters::default();
        f.hide_small_indels = false;
        assert!(f.show_indel(1));
        assert!(f.show_indel(3));
    }

    #[test]
    fn test_row_height() {
        let mut f = DisplayFilters::default();
        assert_eq!(f.row_height(), 4.0); // squished
        f.squished = false;
        assert_eq!(f.row_height(), 12.0); // expanded
    }
}
