//! Haplotype-phase coloring scheme.
//!
//! Matches the IGV.js default palette for HP tags:
//! - HP=1 (Haplotype 1): Blue tones
//! - HP=2 (Haplotype 2): Red/pink tones
//! - Unphased (HP=0): Grey

use egui::Color32;

/// Default colors for haplotype phase tags.
pub const HP1_COLOR: Color32 = Color32::from_rgb(70, 130, 180);   // Steel blue
pub const HP2_COLOR: Color32 = Color32::from_rgb(205, 92, 92);    // Indian red
pub const UNPHASED_COLOR: Color32 = Color32::from_rgb(160, 160, 160); // Grey

/// Lighter variants for squished mode.
pub const HP1_LIGHT: Color32 = Color32::from_rgb(135, 175, 210);
pub const HP2_LIGHT: Color32 = Color32::from_rgb(225, 150, 150);
pub const UNPHASED_LIGHT: Color32 = Color32::from_rgb(200, 200, 200);

/// Get the color for a read based on its HP tag and display mode.
pub fn read_color(hp_tag: u8, squished: bool) -> Color32 {
    if squished {
        match hp_tag {
            1 => HP1_LIGHT,
            2 => HP2_LIGHT,
            _ => UNPHASED_LIGHT,
        }
    } else {
        match hp_tag {
            1 => HP1_COLOR,
            2 => HP2_COLOR,
            _ => UNPHASED_COLOR,
        }
    }
}

/// Color for deletion markers in reads.
pub const DELETION_COLOR: Color32 = Color32::from_rgb(50, 50, 50);

/// Color for insertion markers in reads.
pub const INSERTION_COLOR: Color32 = Color32::from_rgb(128, 0, 128);

/// Color for soft-clipped segments.
pub const SOFTCLIP_COLOR: Color32 = Color32::from_rgb(200, 200, 50);

/// Color for mismatches.
pub const MISMATCH_COLOR: Color32 = Color32::from_rgb(255, 0, 0);

/// Panel header background.
pub const PANEL_HEADER_BG: Color32 = Color32::from_rgb(40, 44, 52);

/// Panel background.
pub const PANEL_BG: Color32 = Color32::from_rgb(255, 255, 255);

/// Reference sequence track color.
pub const REF_TRACK_COLOR: Color32 = Color32::from_rgb(200, 200, 200);

/// Base colors for the reference sequence.
pub fn base_color(base: u8) -> Color32 {
    match base {
        b'A' | b'a' => Color32::from_rgb(0, 150, 0),    // Green
        b'C' | b'c' => Color32::from_rgb(0, 0, 200),    // Blue
        b'G' | b'g' => Color32::from_rgb(200, 160, 0),  // Gold
        b'T' | b't' => Color32::from_rgb(200, 0, 0),    // Red
        b'N' | b'n' => Color32::from_rgb(128, 128, 128), // Grey
        _ => Color32::from_rgb(200, 200, 200),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_color_hp_tags() {
        assert_eq!(read_color(1, false), HP1_COLOR);
        assert_eq!(read_color(2, false), HP2_COLOR);
        assert_eq!(read_color(0, false), UNPHASED_COLOR);
    }

    #[test]
    fn test_read_color_squished() {
        assert_eq!(read_color(1, true), HP1_LIGHT);
        assert_eq!(read_color(2, true), HP2_LIGHT);
        assert_eq!(read_color(0, true), UNPHASED_LIGHT);
    }

    #[test]
    fn test_base_colors_differ() {
        // Each base should have a distinct color
        let a = base_color(b'A');
        let c = base_color(b'C');
        let g = base_color(b'G');
        let t = base_color(b'T');
        assert_ne!(a, c);
        assert_ne!(a, g);
        assert_ne!(a, t);
        assert_ne!(c, g);
        assert_ne!(c, t);
        assert_ne!(g, t);
    }
}
