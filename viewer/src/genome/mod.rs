pub mod bam;
pub mod fasta;
pub mod pileup;

use std::fmt;
use std::io;
use std::path::PathBuf;

/// Error type for genome file operations.
#[derive(Debug)]
pub enum GenomeError {
    /// I/O error (file not found, permission denied, etc.)
    Io(io::Error),
    /// Index file not found at expected path.
    IndexNotFound(PathBuf),
    /// Invalid or unparseable region string.
    InvalidRegion(String),
    /// Parse error in the genomic file format.
    ParseError(String),
}

impl fmt::Display for GenomeError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            GenomeError::Io(e) => write!(f, "I/O error: {e}"),
            GenomeError::IndexNotFound(p) => write!(f, "Index not found: {}", p.display()),
            GenomeError::InvalidRegion(r) => write!(f, "Invalid region: {r}"),
            GenomeError::ParseError(msg) => write!(f, "Parse error: {msg}"),
        }
    }
}

impl std::error::Error for GenomeError {}

impl From<io::Error> for GenomeError {
    fn from(e: io::Error) -> Self {
        GenomeError::Io(e)
    }
}

/// A single aligned read extracted from a BAM or CRAM file.
#[derive(Debug, Clone)]
pub struct AlignedRead {
    /// Read name (QNAME).
    pub name: String,
    /// 1-based alignment start position on the reference.
    pub start: u64,
    /// 1-based alignment end position on the reference (inclusive).
    pub end: u64,
    /// Whether the read maps to the reverse strand.
    pub is_reverse: bool,
    /// Mapping quality (0–255).
    pub mapping_quality: Option<u8>,
    /// Haplotype assignment from the HP tag (typically 1 or 2).
    pub haplotype: Option<u8>,
    /// Raw SAM flag bits.
    pub flags: u16,
}

/// A reference or assembly sequence extracted from a FASTA file.
#[derive(Debug, Clone)]
pub struct FastaSequence {
    /// Sequence name (e.g., "chr1").
    pub name: String,
    /// 1-based start of the returned sequence.
    pub start: u64,
    /// 1-based end of the returned sequence (inclusive).
    pub end: u64,
    /// The nucleotide sequence as an uppercase ASCII string.
    pub sequence: String,
}
