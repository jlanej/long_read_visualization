pub mod bam;
pub mod coordinate_mapper;
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

/// Type of insertion or deletion in a read alignment.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IndelKind {
    /// Insertion into the read (bases present in read but not reference).
    Insertion,
    /// Deletion from the reference (bases present in reference but not read).
    Deletion,
}

/// A single insertion or deletion event within an aligned read.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Indel {
    /// Reference position (1-based) where the indel occurs.
    pub ref_pos: u64,
    /// Length of the indel in base pairs.
    pub length: u32,
    /// Whether this is an insertion or deletion.
    pub kind: IndelKind,
}

/// A single base-level mismatch between a read and the reference.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Mismatch {
    /// Reference position (1-based) of the mismatch.
    pub ref_pos: u64,
    /// The read base at this position (uppercase ASCII).
    pub read_base: u8,
}

/// A soft-clipped region at the start or end of a read.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SoftClip {
    /// Reference position (1-based) where the clip occurs.
    /// For leading clips this equals alignment start; for trailing clips, alignment end + 1.
    pub ref_pos: u64,
    /// Number of soft-clipped bases.
    pub length: u32,
    /// Whether this is a leading (true) or trailing (false) soft clip.
    pub is_leading: bool,
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
    /// Indels extracted from the CIGAR string.
    pub indels: Vec<Indel>,
    /// Base-level mismatches (from `X` CIGAR ops with `--eqx` encoding).
    pub mismatches: Vec<Mismatch>,
    /// Soft-clipped regions at the read boundaries.
    pub soft_clips: Vec<SoftClip>,
}

/// A reference or assembly sequence extracted from a FASTA file.
#[derive(Debug, Clone)]
pub struct FastaSequence {
    /// Sequence name (e.g., "chr1").
    pub name: String,
    /// 1-based start of the returned sequence.
    #[allow(dead_code)]
    pub start: u64,
    /// 1-based end of the returned sequence (inclusive).
    #[allow(dead_code)]
    pub end: u64,
    /// The nucleotide sequence as an uppercase ASCII string.
    pub sequence: String,
}
