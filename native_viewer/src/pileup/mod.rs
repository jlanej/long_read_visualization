pub mod layout;

use crate::genome::bam::AlignedRead;

/// Pileup data for a single panel: pre-computed read layout.
#[derive(Debug, Clone)]
pub struct PileupData {
    /// Reads packed into rows for rendering.
    pub rows: Vec<ReadRow>,
    /// Total number of reads before layout.
    pub total_reads: usize,
    /// Visible region start (0-based).
    pub region_start: u64,
    /// Visible region end (0-based, exclusive).
    pub region_end: u64,
}

/// A row of non-overlapping reads in the pileup display.
#[derive(Debug, Clone)]
pub struct ReadRow {
    pub reads: Vec<AlignedRead>,
}

impl ReadRow {
    pub fn new() -> Self {
        Self { reads: Vec::new() }
    }
}

impl Default for ReadRow {
    fn default() -> Self {
        Self::new()
    }
}
