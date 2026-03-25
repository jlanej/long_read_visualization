//! Indexed FASTA reader using noodles for reference sequence retrieval.

use std::fs::File;
use std::io::BufReader;
use std::path::Path;

use anyhow::{Context, Result};
use noodles::core::Region;
use noodles::fasta;
use noodles::fasta::fai;

use crate::genome::region::GenomicRegion;

/// Result of a FASTA query: the reference sequence for a region.
#[derive(Debug, Clone)]
pub struct ReferenceSequence {
    /// The chromosome/contig name.
    pub chrom: String,
    /// 0-based start position.
    pub start: u64,
    /// The nucleotide sequence (uppercase).
    pub sequence: Vec<u8>,
}

impl ReferenceSequence {
    /// Get the base at a given 0-based position relative to the region start.
    pub fn base_at(&self, offset: usize) -> Option<u8> {
        self.sequence.get(offset).copied()
    }

    /// Get the base at an absolute genomic position.
    pub fn base_at_pos(&self, pos: u64) -> Option<u8> {
        if pos < self.start {
            return None;
        }
        let offset = (pos - self.start) as usize;
        self.base_at(offset)
    }

    /// Length of the retrieved sequence.
    pub fn len(&self) -> usize {
        self.sequence.len()
    }

    pub fn is_empty(&self) -> bool {
        self.sequence.is_empty()
    }
}

/// Query a FASTA file for a specific region.
///
/// Supports both plain and bgzip-compressed FASTA files.
/// Requires a .fai index file alongside the FASTA.
pub fn query_fasta(
    fasta_path: &Path,
    region: &GenomicRegion,
) -> Result<ReferenceSequence> {
    // Load the FASTA index
    let fai_path = {
        let mut p = fasta_path.as_os_str().to_owned();
        p.push(".fai");
        std::path::PathBuf::from(p)
    };

    let index = fai::fs::read(&fai_path)
        .with_context(|| format!("Failed to read FASTA index: {}", fai_path.display()))?;

    // Open the FASTA file
    let file = File::open(fasta_path)
        .with_context(|| format!("Failed to open FASTA: {}", fasta_path.display()))?;

    let mut reader = fasta::io::IndexedReader::new(BufReader::new(file), index);

    // Build region string (1-based for noodles)
    let region_str = format!("{}:{}-{}", region.chrom, region.start + 1, region.end);
    let query_region: Region = region_str
        .parse()
        .with_context(|| format!("Invalid FASTA region: {region_str}"))?;

    // Query the sequence
    let record = reader.query(&query_region)?;
    let seq_bytes: Vec<u8> = record
        .sequence()
        .as_ref()
        .iter()
        .map(|&b| b.to_ascii_uppercase())
        .collect();

    Ok(ReferenceSequence {
        chrom: region.chrom.clone(),
        start: region.start,
        sequence: seq_bytes,
    })
}

/// List all sequences (chromosomes) in a FASTA index.
pub fn list_sequences(fasta_path: &Path) -> Result<Vec<(String, u64)>> {
    let fai_path = {
        let mut p = fasta_path.as_os_str().to_owned();
        p.push(".fai");
        std::path::PathBuf::from(p)
    };

    let index = fai::fs::read(&fai_path)
        .with_context(|| format!("Failed to read FASTA index: {}", fai_path.display()))?;

    let seqs: Vec<(String, u64)> = index
        .as_ref()
        .iter()
        .map(|record: &fai::Record| {
            (
                String::from_utf8_lossy(record.name()).to_string(),
                record.length(),
            )
        })
        .collect();

    Ok(seqs)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_reference_sequence_base_at() {
        let seq = ReferenceSequence {
            chrom: "chr1".to_string(),
            start: 100,
            sequence: b"ACGTACGT".to_vec(),
        };
        assert_eq!(seq.base_at(0), Some(b'A'));
        assert_eq!(seq.base_at(3), Some(b'T'));
        assert_eq!(seq.base_at(8), None);
        assert_eq!(seq.base_at_pos(100), Some(b'A'));
        assert_eq!(seq.base_at_pos(103), Some(b'T'));
        assert_eq!(seq.base_at_pos(99), None);
        assert_eq!(seq.len(), 8);
    }
}
