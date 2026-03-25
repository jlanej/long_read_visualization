use serde::{Deserialize, Serialize};
use std::fmt;

/// A genomic interval: chromosome name + 1-based start/end coordinates.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct GenomicRegion {
    pub chrom: String,
    pub start: u64,
    pub end: u64,
}

impl GenomicRegion {
    pub fn new(chrom: impl Into<String>, start: u64, end: u64) -> Self {
        Self {
            chrom: chrom.into(),
            start,
            end,
        }
    }

    /// Parse a region string like "chr1:1000-2000".
    pub fn parse(s: &str) -> Option<Self> {
        let (chrom, coords) = s.rsplit_once(':')?;
        let (start_s, end_s) = coords.split_once('-')?;
        let start: u64 = start_s.trim().replace(',', "").parse().ok()?;
        let end: u64 = end_s.trim().replace(',', "").parse().ok()?;
        if start > end {
            return None;
        }
        Some(Self {
            chrom: chrom.to_string(),
            start,
            end,
        })
    }

    /// Span (length) of this region in base pairs.
    pub fn span(&self) -> u64 {
        self.end.saturating_sub(self.start)
    }

    /// Center position of this region.
    pub fn center(&self) -> u64 {
        (self.start + self.end) / 2
    }

    /// Pan the region by a fraction of its span.
    /// Positive fraction = move right, negative = move left.
    pub fn pan(&mut self, fraction: f64) {
        let shift = (self.span() as f64 * fraction) as i64;
        let new_start = (self.start as i64 + shift).max(0) as u64;
        let span = self.span();
        self.start = new_start;
        self.end = new_start + span;
    }

    /// Zoom by a factor. >1 = zoom out, <1 = zoom in.
    pub fn zoom(&mut self, factor: f64) {
        let center = self.center();
        let half_span = ((self.span() as f64 * factor) / 2.0).max(1.0) as u64;
        self.start = center.saturating_sub(half_span);
        self.end = center + half_span;
    }

    /// Check if this region overlaps another region on the same chromosome.
    pub fn overlaps(&self, other: &GenomicRegion) -> bool {
        self.chrom == other.chrom && self.start < other.end && other.start < self.end
    }

    /// Check if a position falls within this region.
    pub fn contains(&self, chrom: &str, pos: u64) -> bool {
        self.chrom == chrom && pos >= self.start && pos < self.end
    }
}

impl fmt::Display for GenomicRegion {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}:{}-{}", self.chrom, self.start, self.end)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_region() {
        let r = GenomicRegion::parse("chr1:1000-2000").unwrap();
        assert_eq!(r.chrom, "chr1");
        assert_eq!(r.start, 1000);
        assert_eq!(r.end, 2000);
    }

    #[test]
    fn test_parse_region_with_commas() {
        let r = GenomicRegion::parse("chr1:1,000-2,000").unwrap();
        assert_eq!(r.start, 1000);
        assert_eq!(r.end, 2000);
    }

    #[test]
    fn test_parse_invalid() {
        assert!(GenomicRegion::parse("").is_none());
        assert!(GenomicRegion::parse("chr1").is_none());
        assert!(GenomicRegion::parse("chr1:2000-1000").is_none()); // start > end
    }

    #[test]
    fn test_span() {
        let r = GenomicRegion::new("chr1", 1000, 2000);
        assert_eq!(r.span(), 1000);
    }

    #[test]
    fn test_pan() {
        let mut r = GenomicRegion::new("chr1", 1000, 2000);
        r.pan(0.5); // move right by half the span
        assert_eq!(r.start, 1500);
        assert_eq!(r.end, 2500);
    }

    #[test]
    fn test_pan_left_clamp() {
        let mut r = GenomicRegion::new("chr1", 100, 200);
        r.pan(-2.0); // try to move left past zero
        assert_eq!(r.start, 0);
        assert_eq!(r.end, 100);
    }

    #[test]
    fn test_zoom_in() {
        let mut r = GenomicRegion::new("chr1", 1000, 2000);
        r.zoom(0.5);
        assert_eq!(r.start, 1250);
        assert_eq!(r.end, 1750);
    }

    #[test]
    fn test_zoom_out() {
        let mut r = GenomicRegion::new("chr1", 1000, 2000);
        r.zoom(2.0);
        assert_eq!(r.start, 500);
        assert_eq!(r.end, 2500);
    }

    #[test]
    fn test_overlaps() {
        let r1 = GenomicRegion::new("chr1", 1000, 2000);
        let r2 = GenomicRegion::new("chr1", 1500, 2500);
        let r3 = GenomicRegion::new("chr1", 3000, 4000);
        let r4 = GenomicRegion::new("chr2", 1000, 2000);

        assert!(r1.overlaps(&r2));
        assert!(!r1.overlaps(&r3));
        assert!(!r1.overlaps(&r4)); // different chrom
    }

    #[test]
    fn test_contains() {
        let r = GenomicRegion::new("chr1", 1000, 2000);
        assert!(r.contains("chr1", 1000));
        assert!(r.contains("chr1", 1500));
        assert!(!r.contains("chr1", 2000)); // half-open [start, end)
        assert!(!r.contains("chr2", 1500));
    }

    #[test]
    fn test_display() {
        let r = GenomicRegion::new("chr1", 1000, 2000);
        assert_eq!(format!("{r}"), "chr1:1000-2000");
    }
}
