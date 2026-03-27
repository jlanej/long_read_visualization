//! Genomic coordinate ruler with adaptive tick marks.
//!
//! Renders a coordinate axis at the top of each panel showing chromosome
//! position with automatically scaled tick intervals (1bp → Mbp).

/// A single tick mark on the coordinate ruler.
#[derive(Debug, Clone, PartialEq)]
pub struct Tick {
    /// Genomic coordinate of this tick.
    pub position: u64,
    /// Pixel X offset from the left edge of the panel.
    pub x: f32,
    /// Whether this is a major tick (with a label) or a minor tick.
    pub is_major: bool,
    /// Human-readable label for major ticks (e.g., "1,234,567").
    pub label: Option<String>,
}

/// Height of the ruler track in pixels.
pub const RULER_HEIGHT: f32 = 20.0;

/// Desired minimum spacing between major tick labels in pixels.
const MIN_MAJOR_TICK_PX: f32 = 80.0;

/// Compute a "nice" tick interval for the given view span and panel width.
///
/// Returns `(major_interval, minor_interval)` in base pairs.  The major
/// interval is chosen so that labels are at least `MIN_MAJOR_TICK_PX` apart;
/// the minor interval divides the major interval into sub-ticks.
pub fn tick_intervals(view_span: u64, panel_width: f32) -> (u64, u64) {
    if view_span == 0 || panel_width <= 0.0 {
        return (1, 1);
    }

    // How many bp one pixel covers
    let bp_per_px = view_span as f64 / panel_width as f64;
    // Minimum interval in bp that would give enough spacing
    let min_interval_bp = bp_per_px * MIN_MAJOR_TICK_PX as f64;

    // Round up to a "nice" number: 1, 2, 5, 10, 20, 50, 100, …
    let magnitude = 10.0_f64.powf(min_interval_bp.log10().floor());
    let residual = min_interval_bp / magnitude;

    let major = if residual <= 1.0 {
        magnitude as u64
    } else if residual <= 2.0 {
        (2.0 * magnitude) as u64
    } else if residual <= 5.0 {
        (5.0 * magnitude) as u64
    } else {
        (10.0 * magnitude) as u64
    }
    .max(1);

    // Minor ticks: split major interval into 5 parts (or 2 if too small)
    let minor = if major >= 10 { major / 5 } else { major / 2 }.max(1);

    (major, minor)
}

/// Generate tick marks for the current view.
///
/// `view_start` / `view_end`: visible genomic coordinate range.
/// `panel_width`: pixel width of the drawing area.
pub fn compute_ticks(view_start: u64, view_end: u64, panel_width: f32) -> Vec<Tick> {
    let span = view_end.saturating_sub(view_start);
    if span == 0 || panel_width <= 0.0 {
        return Vec::new();
    }

    let (major, minor) = tick_intervals(span, panel_width);
    let bp_per_px = panel_width / span as f32;

    let mut ticks = Vec::new();

    // Start from the first minor tick at or after view_start
    let first_minor = if view_start.is_multiple_of(minor) {
        view_start
    } else {
        view_start + (minor - view_start % minor)
    };

    let mut pos = first_minor;
    while pos <= view_end {
        let x = (pos - view_start) as f32 * bp_per_px;
        let is_major = pos % major == 0;
        let label = if is_major {
            Some(format_position(pos))
        } else {
            None
        };
        ticks.push(Tick {
            position: pos,
            x,
            is_major,
            label,
        });
        pos = pos.saturating_add(minor);
        if pos <= first_minor {
            break; // overflow guard
        }
    }

    ticks
}

/// Format a genomic position for display with SI suffixes.
fn format_position(pos: u64) -> String {
    if pos >= 1_000_000 && pos.is_multiple_of(1_000_000) {
        format!("{}M", pos / 1_000_000)
    } else if pos >= 1_000_000 {
        format!("{:.1}M", pos as f64 / 1_000_000.0)
    } else if pos >= 1_000 && pos.is_multiple_of(1_000) {
        format!("{}k", pos / 1_000)
    } else if pos >= 10_000 {
        format!("{:.1}k", pos as f64 / 1_000.0)
    } else {
        format_with_commas(pos)
    }
}

/// Format a number with comma separators.
fn format_with_commas(n: u64) -> String {
    let s = n.to_string();
    let mut result = String::with_capacity(s.len() + s.len() / 3);
    for (i, c) in s.chars().enumerate() {
        if i > 0 && (s.len() - i).is_multiple_of(3) {
            result.push(',');
        }
        result.push(c);
    }
    result
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_tick_intervals_basic() {
        // 1000 bp in 800 px → ~1.25 bp/px → min_interval ≈ 100 bp
        let (major, minor) = tick_intervals(1000, 800.0);
        assert!(major >= 100, "major={major}");
        assert!(minor >= 1, "minor={minor}");
        assert!(major >= minor);
    }

    #[test]
    fn test_tick_intervals_wide_view() {
        // 1M bp in 1000 px → 1000 bp/px → min_interval ≈ 80000 bp
        let (major, _minor) = tick_intervals(1_000_000, 1000.0);
        assert!(major >= 50_000, "major={major}");
    }

    #[test]
    fn test_tick_intervals_narrow_view() {
        // 100 bp in 1000 px → 0.1 bp/px → min_interval ≈ 8 bp
        let (major, minor) = tick_intervals(100, 1000.0);
        assert!(major <= 20, "major={major}");
        assert!(minor >= 1);
    }

    #[test]
    fn test_tick_intervals_zero_span() {
        let (major, minor) = tick_intervals(0, 1000.0);
        assert_eq!(major, 1);
        assert_eq!(minor, 1);
    }

    #[test]
    fn test_tick_intervals_zero_width() {
        let (major, minor) = tick_intervals(1000, 0.0);
        assert_eq!(major, 1);
        assert_eq!(minor, 1);
    }

    #[test]
    fn test_compute_ticks_basic() {
        let ticks = compute_ticks(1000, 2000, 1000.0);
        assert!(!ticks.is_empty());
        // All ticks should be within view range
        for t in &ticks {
            assert!(t.position >= 1000 && t.position <= 2000);
            assert!(t.x >= 0.0 && t.x <= 1000.0);
        }
        // Should have at least one major tick
        assert!(ticks.iter().any(|t| t.is_major));
    }

    #[test]
    fn test_compute_ticks_empty_span() {
        let ticks = compute_ticks(1000, 1000, 800.0);
        assert!(ticks.is_empty());
    }

    #[test]
    fn test_compute_ticks_major_have_labels() {
        let ticks = compute_ticks(0, 10_000, 1000.0);
        for t in &ticks {
            if t.is_major {
                assert!(
                    t.label.is_some(),
                    "Major tick at {} has no label",
                    t.position
                );
            } else {
                assert!(t.label.is_none());
            }
        }
    }

    #[test]
    fn test_format_position_small() {
        assert_eq!(format_position(0), "0");
        assert_eq!(format_position(42), "42");
        assert_eq!(format_position(999), "999");
    }

    #[test]
    fn test_format_position_thousands() {
        assert_eq!(format_position(1000), "1k");
        assert_eq!(format_position(5000), "5k");
        assert_eq!(format_position(1_234), "1,234");
        assert_eq!(format_position(10_500), "10.5k");
    }

    #[test]
    fn test_format_position_millions() {
        assert_eq!(format_position(1_000_000), "1M");
        assert_eq!(format_position(2_500_000), "2.5M");
    }

    #[test]
    fn test_format_with_commas() {
        assert_eq!(format_with_commas(0), "0");
        assert_eq!(format_with_commas(999), "999");
        assert_eq!(format_with_commas(1000), "1,000");
        assert_eq!(format_with_commas(1_234_567), "1,234,567");
    }

    #[test]
    fn test_ticks_are_sorted_by_position() {
        let ticks = compute_ticks(5000, 15000, 800.0);
        for w in ticks.windows(2) {
            assert!(w[0].position <= w[1].position);
        }
    }

    #[test]
    fn test_ticks_x_increases_monotonically() {
        let ticks = compute_ticks(10000, 20000, 600.0);
        for w in ticks.windows(2) {
            assert!(w[0].x <= w[1].x);
        }
    }
}
