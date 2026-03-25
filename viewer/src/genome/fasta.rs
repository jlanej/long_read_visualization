use std::io;
use std::path::Path;

use noodles::core::{Position, Region};
use noodles::core::region::Interval;
use noodles::fasta as noodles_fasta;

use super::{FastaSequence, GenomeError};

/// Query a FASTA file using a region string like `"chr1:100-200"`.
///
/// Works when sequence names do **not** contain colons. For FASTA files
/// whose sequence names include colons (e.g., the toy dataset), use
/// [`query_fasta_region`] instead.
pub fn query_fasta(fasta_path: &Path, region_str: &str) -> Result<FastaSequence, GenomeError> {
    let region: Region = region_str
        .parse()
        .map_err(|e| GenomeError::InvalidRegion(format!("{region_str}: {e}")))?;

    query_fasta_with_region(fasta_path, &region)
}

/// Query a FASTA file by exact sequence name and 1-based coordinates.
///
/// This avoids the ambiguity of noodles' region string parser when
/// sequence names contain colons (common in extracted sub-FASTA files).
pub fn query_fasta_region(
    fasta_path: &Path,
    name: &str,
    start: u64,
    end: u64,
) -> Result<FastaSequence, GenomeError> {
    let s = Position::try_from(start as usize)
        .map_err(|e| GenomeError::InvalidRegion(format!("start {start}: {e}")))?;
    let e = Position::try_from(end as usize)
        .map_err(|e| GenomeError::InvalidRegion(format!("end {end}: {e}")))?;

    let region = Region::new(name, Interval::from(s..=e));
    query_fasta_with_region(fasta_path, &region)
}

/// List sequences available in a FASTA index.
///
/// Returns `(name, length)` pairs read from the `.fai` sidecar file.
pub fn list_sequences(fasta_path: &Path) -> Result<Vec<(String, u64)>, GenomeError> {
    let fai_path = fai_path_for(fasta_path);
    let content = std::fs::read_to_string(&fai_path).map_err(|e| match e.kind() {
        io::ErrorKind::NotFound => GenomeError::IndexNotFound(fai_path),
        _ => GenomeError::Io(e),
    })?;

    let mut seqs = Vec::new();
    for line in content.lines() {
        let mut cols = line.split('\t');
        if let (Some(name), Some(len_str)) = (cols.next(), cols.next()) {
            let len = len_str
                .parse::<u64>()
                .map_err(|e| GenomeError::ParseError(format!("fai length: {e}")))?;
            seqs.push((name.to_string(), len));
        }
    }
    Ok(seqs)
}

// ---------------------------------------------------------------------------
// Internal helpers
// ---------------------------------------------------------------------------

fn query_fasta_with_region(
    fasta_path: &Path,
    region: &Region,
) -> Result<FastaSequence, GenomeError> {
    if !fasta_path.exists() {
        return Err(GenomeError::Io(io::Error::new(
            io::ErrorKind::NotFound,
            format!("FASTA file not found: {}", fasta_path.display()),
        )));
    }

    let mut reader = noodles_fasta::io::indexed_reader::Builder::default()
        .build_from_path(fasta_path)
        .map_err(|e| match e.kind() {
            io::ErrorKind::NotFound => {
                GenomeError::IndexNotFound(fai_path_for(fasta_path))
            }
            _ => GenomeError::Io(e),
        })?;

    let record = reader
        .query(region)
        .map_err(|e| GenomeError::ParseError(format!("FASTA query: {e}")))?;

    let name = String::from_utf8_lossy(record.name()).into_owned();
    let seq_bytes = record.sequence().as_ref();
    let sequence: String = seq_bytes
        .iter()
        .map(|&b| (b as char).to_ascii_uppercase())
        .collect();

    let (start, end) = region_bounds(region);

    Ok(FastaSequence {
        name,
        start,
        end,
        sequence,
    })
}

/// Expected `.fai` path for a given FASTA file.
fn fai_path_for(fasta_path: &Path) -> std::path::PathBuf {
    let mut p = fasta_path.as_os_str().to_owned();
    p.push(".fai");
    p.into()
}

/// Extract 1-based start/end from a [`Region`].
fn region_bounds(region: &Region) -> (u64, u64) {
    let interval = region.interval();

    let start = interval
        .start()
        .map(|p| usize::from(p) as u64)
        .unwrap_or(1);

    let end = interval
        .end()
        .map(|p| usize::from(p) as u64)
        .unwrap_or(start);

    (start, end)
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    fn toy_reference() -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("../resources/toy_dataset/toy_reference.fa.gz")
    }

    fn toy_hap1() -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("../resources/toy_dataset/toy_hap1.fa.gz")
    }

    // -- query_fasta_region tests (exact sequence name) --

    #[test]
    fn test_query_fasta_region_returns_sequence() {
        let fa = toy_reference();
        if !fa.exists() {
            eprintln!("skipping: toy reference not found");
            return;
        }
        // The toy reference has sequence named "chr1:112064095-112277008".
        // Query local positions 1-100 within that sequence.
        let result = query_fasta_region(
            &fa,
            "chr1:112064095-112277008",
            1,
            100,
        );
        match result {
            Ok(seq) => {
                assert_eq!(seq.sequence.len(), 100);
                assert!(
                    seq.sequence.chars().all(|c| "ACGTNRYSWKMBDHV".contains(c)),
                    "unexpected bases in sequence: {}",
                    &seq.sequence[..20.min(seq.sequence.len())]
                );
            }
            Err(e) => panic!("unexpected error: {e}"),
        }
    }

    #[test]
    fn test_query_fasta_region_hap1() {
        let fa = toy_hap1();
        if !fa.exists() {
            return;
        }
        // First hap1 sequence
        let result = query_fasta_region(
            &fa,
            "NA21110#1#CM089663.1:111967298-112167366",
            1,
            100,
        );
        match result {
            Ok(seq) => {
                assert_eq!(seq.sequence.len(), 100);
            }
            Err(e) => panic!("unexpected error: {e}"),
        }
    }

    #[test]
    fn test_query_fasta_region_length_matches() {
        let fa = toy_reference();
        if !fa.exists() {
            return;
        }
        let result = query_fasta_region(
            &fa,
            "chr1:112064095-112277008",
            500,
            599,
        );
        match result {
            Ok(seq) => {
                assert_eq!(
                    seq.sequence.len(),
                    100,
                    "expected 100 bases for a 100bp window"
                );
            }
            Err(e) => panic!("unexpected error: {e}"),
        }
    }

    // -- list_sequences --

    #[test]
    fn test_list_sequences() {
        let fa = toy_reference();
        if !fa.exists() {
            return;
        }
        let seqs = list_sequences(&fa).expect("should read index");
        assert_eq!(seqs.len(), 10, "toy reference has 10 sequences");
        // First entry should be the chr1 region
        assert!(
            seqs[0].0.starts_with("chr1:"),
            "first seq should be chr1 region, got {}",
            seqs[0].0
        );
        assert!(seqs[0].1 > 0, "sequence length should be positive");
    }

    // -- Error handling --

    #[test]
    fn test_query_fasta_missing_file() {
        let result = query_fasta(Path::new("/nonexistent/ref.fa"), "chr1:1-100");
        assert!(result.is_err());
    }

    #[test]
    fn test_query_fasta_invalid_region() {
        let fa = toy_reference();
        if !fa.exists() {
            return;
        }
        // A malformed interval that the region parser should reject
        let result = query_fasta(&fa, "chr1:end-start");
        assert!(
            result.is_err(),
            "malformed interval should produce an error"
        );
    }

    #[test]
    fn test_query_fasta_missing_index() {
        let tmp = tempfile::tempdir().unwrap();
        let fa_copy = tmp.path().join("no_index.fa.gz");
        std::fs::copy(toy_reference(), &fa_copy).ok();
        if !fa_copy.exists() {
            return;
        }
        let result = query_fasta(&fa_copy, "chr1:1-100");
        assert!(result.is_err(), "should fail without .fai");
    }

    #[test]
    fn test_query_fasta_region_missing_file() {
        let result = query_fasta_region(Path::new("/nonexistent/ref.fa"), "chr1", 1, 100);
        assert!(result.is_err());
    }
}
