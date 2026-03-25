use std::io;
use std::io::Read;
use std::path::Path;

use noodles::bam;
use noodles::core::Region;
use noodles::cram;
use noodles::fasta;
use noodles::fasta::repository::adapters::IndexedReader;
use noodles::sam;
use noodles::sam::alignment::record::cigar::Cigar as _;
use noodles::sam::alignment::record::data::field::tag::Tag;

use super::{AlignedRead, GenomeError, Indel, IndelKind};

// ---------------------------------------------------------------------------
// Lenient SAM header parsing helpers
// ---------------------------------------------------------------------------

/// Remove duplicate `@PG` lines from SAM header text, keeping only the first
/// occurrence of each program ID.  This works around real-world BAM/CRAM files
/// produced by samtools pipelines that repeat `@PG` entries.
fn deduplicate_pg_ids(header_text: &str) -> String {
    use std::collections::HashSet;

    let mut seen: HashSet<String> = HashSet::new();
    let mut out = String::with_capacity(header_text.len());

    for line in header_text.lines() {
        if line.starts_with("@PG") {
            let id = line
                .split('\t')
                .find_map(|f| f.strip_prefix("ID:"))
                .unwrap_or("");
            if !seen.insert(id.to_string()) {
                continue; // duplicate – skip
            }
        }
        out.push_str(line);
        out.push('\n');
    }

    out
}

/// Read the raw SAM header text from a BAM file by decoding its BGZF blocks.
///
/// BAM layout: `magic(4)  header_len(u32-LE)  header_text(header_len bytes) …`
fn read_bam_raw_header(bam_path: &Path) -> Result<String, GenomeError> {
    let file = std::fs::File::open(bam_path).map_err(GenomeError::Io)?;
    let mut gz = flate2::read::MultiGzDecoder::new(file);

    let mut magic = [0u8; 4];
    gz.read_exact(&mut magic).map_err(GenomeError::Io)?;

    let mut len_buf = [0u8; 4];
    gz.read_exact(&mut len_buf).map_err(GenomeError::Io)?;
    let header_len = u32::from_le_bytes(len_buf) as usize;

    let mut header_bytes = vec![0u8; header_len];
    gz.read_exact(&mut header_bytes).map_err(GenomeError::Io)?;

    String::from_utf8(header_bytes)
        .map(|s| s.trim_end_matches('\0').to_string())
        .map_err(|e| GenomeError::ParseError(format!("BAM header not UTF-8: {e}")))
}

/// Try to parse a SAM header leniently: if the raw text has duplicate `@PG`
/// IDs, remove duplicates and re-parse.
fn parse_header_lenient(raw: &str) -> Result<sam::Header, GenomeError> {
    let sanitized = deduplicate_pg_ids(raw);
    sanitized
        .parse::<sam::Header>()
        .map_err(|e| GenomeError::ParseError(format!("SAM header (lenient): {e}")))
}

// ---------------------------------------------------------------------------
// BAM
// ---------------------------------------------------------------------------

/// Query a BAM file for aligned reads overlapping `region_str`.
///
/// The BAM must have an associated `.bai` index discoverable by noodles
/// (i.e. `<bam_path>.bai` or replacing the `.bam` extension with `.bai`).
pub fn query_bam(bam_path: &Path, region_str: &str) -> Result<Vec<AlignedRead>, GenomeError> {
    if !bam_path.exists() {
        return Err(GenomeError::Io(io::Error::new(
            io::ErrorKind::NotFound,
            format!("BAM file not found: {}", bam_path.display()),
        )));
    }

    let region: Region = region_str
        .parse()
        .map_err(|e| GenomeError::InvalidRegion(format!("{region_str}: {e}")))?;

    let mut reader = bam::io::indexed_reader::Builder::default()
        .build_from_path(bam_path)
        .map_err(|e| match e.kind() {
            io::ErrorKind::NotFound => {
                let bai = bam_path.with_extension("bam.bai");
                GenomeError::IndexNotFound(bai)
            }
            _ => GenomeError::Io(e),
        })?;

    let header = match reader.read_header() {
        Ok(h) => h,
        Err(e) if e.to_string().contains("duplicate program") => {
            // Noodles rejects headers with duplicate @PG IDs, which are
            // common in samtools-produced BAM files.  Fall back to reading
            // the raw SAM header text, deduplicating @PG entries, and
            // re-parsing.  A fresh reader is needed because the failed
            // read_header() leaves the internal cursor in an undefined
            // position; indexed queries use absolute virtual-file-offset
            // seeks, so skipping read_header() on the new reader is safe.
            let raw = read_bam_raw_header(bam_path)?;
            let header = parse_header_lenient(&raw)?;
            reader = bam::io::indexed_reader::Builder::default()
                .build_from_path(bam_path)
                .map_err(GenomeError::Io)?;
            header
        }
        Err(e) => return Err(GenomeError::ParseError(format!("BAM header: {e}"))),
    };

    let mut reads = Vec::new();

    for result in reader.query(&header, &region)?.records() {
        match result {
            Ok(record) => {
                if let Some(r) = extract_bam_read(&record)? {
                    reads.push(r);
                }
            }
            Err(_) => continue, // skip unparseable records
        }
    }

    Ok(reads)
}

/// Extract an [`AlignedRead`] from a BAM record (returns `None` for unmapped).
fn extract_bam_read(record: &bam::Record) -> Result<Option<AlignedRead>, GenomeError> {
    let flags = record.flags();

    if flags.is_unmapped() {
        return Ok(None);
    }

    let start = match record.alignment_start() {
        Some(Ok(pos)) => usize::from(pos) as u64,
        _ => return Ok(None),
    };

    let cigar = record.cigar();
    let span = cigar
        .alignment_span()
        .map_err(|e| GenomeError::ParseError(format!("CIGAR: {e}")))?;
    let end = if span > 0 {
        start + span as u64 - 1
    } else {
        start
    };

    let indels = extract_indels_from_bam_cigar(&cigar, start)?;

    let name = record.name().map(|n| format!("{n}")).unwrap_or_default();

    let mq = record.mapping_quality().map(u8::from);
    let hp = hp_from_bam(record);

    Ok(Some(AlignedRead {
        name,
        start,
        end,
        is_reverse: flags.is_reverse_complemented(),
        mapping_quality: mq,
        haplotype: hp,
        flags: u16::from(flags),
        indels,
    }))
}

/// Walk a BAM CIGAR string and collect insertion/deletion events with their
/// reference positions and lengths.
fn extract_indels_from_bam_cigar(
    cigar: &noodles::bam::record::Cigar<'_>,
    alignment_start: u64,
) -> Result<Vec<Indel>, GenomeError> {
    use noodles::sam::alignment::record::cigar::op::Kind;

    let mut indels = Vec::new();
    let mut ref_pos = alignment_start;

    for result in cigar.iter() {
        let op = result.map_err(|e| GenomeError::ParseError(format!("CIGAR op: {e}")))?;
        let len = op.len();
        match op.kind() {
            Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                ref_pos += len as u64;
            }
            Kind::Insertion => {
                indels.push(Indel {
                    ref_pos,
                    length: len as u32,
                    kind: IndelKind::Insertion,
                });
                // Insertions don't consume reference bases
            }
            Kind::Deletion => {
                indels.push(Indel {
                    ref_pos,
                    length: len as u32,
                    kind: IndelKind::Deletion,
                });
                ref_pos += len as u64;
            }
            Kind::SoftClip | Kind::HardClip | Kind::Pad => {
                // Soft/hard clips and pads don't consume reference positions
            }
            Kind::Skip => {
                ref_pos += len as u64;
            }
        }
    }

    Ok(indels)
}

/// Read the HP (haplotype) auxiliary tag from a BAM record.
fn hp_from_bam(record: &bam::Record) -> Option<u8> {
    let tag = Tag::new(b'H', b'P');
    record.data().get(&tag).and_then(|res| {
        res.ok().and_then(|v| {
            use noodles::sam::alignment::record::data::field::Value;
            match v {
                Value::UInt8(n) => Some(n),
                Value::Int8(n) => Some(n as u8),
                Value::UInt16(n) => Some(n as u8),
                Value::Int16(n) => Some(n as u8),
                Value::UInt32(n) => Some(n as u8),
                Value::Int32(n) => Some(n as u8),
                _ => None,
            }
        })
    })
}

// ---------------------------------------------------------------------------
// CRAM
// ---------------------------------------------------------------------------

/// Read the raw SAM header text from a CRAM file.
///
/// The header container typically stores the SAM text uncompressed.  We scan
/// the first portion of the file for the `@HD`/`@SQ` marker that begins the
/// SAM header block.
fn read_cram_raw_header(cram_path: &Path) -> Result<String, GenomeError> {
    let mut file = std::fs::File::open(cram_path).map_err(GenomeError::Io)?;

    // Read the first 2 MiB – more than enough for any SAM header.
    let mut buf = vec![0u8; 2 * 1024 * 1024];
    let n = file.read(&mut buf).map_err(GenomeError::Io)?;
    buf.truncate(n);

    // Skip the 26-byte CRAM file definition and look for the SAM header.
    let search_start = 26.min(n);
    let header_start = buf[search_start..]
        .windows(3)
        .position(|w| w == b"@HD" || w == b"@SQ" || w == b"@RG" || w == b"@PG")
        .map(|p| search_start + p)
        .ok_or_else(|| GenomeError::ParseError("SAM header not found in CRAM".into()))?;

    // The SAM text ends at a NUL byte or a control byte below TAB.
    let header_end = buf[header_start..]
        .iter()
        .position(|&b| b == 0 || b < b'\t')
        .map(|p| header_start + p)
        .unwrap_or(n);

    String::from_utf8(buf[header_start..header_end].to_vec())
        .map_err(|e| GenomeError::ParseError(format!("CRAM header not UTF-8: {e}")))
}

/// Query a CRAM file for aligned reads overlapping `region_str`.
///
/// If `reference_path` is provided it is used for CRAM decoding; otherwise
/// an empty repository is used (works when the CRAM embeds reference blocks).
pub fn query_cram(
    cram_path: &Path,
    reference_path: Option<&Path>,
    region_str: &str,
) -> Result<Vec<AlignedRead>, GenomeError> {
    if !cram_path.exists() {
        return Err(GenomeError::Io(io::Error::new(
            io::ErrorKind::NotFound,
            format!("CRAM file not found: {}", cram_path.display()),
        )));
    }

    let region: Region = region_str
        .parse()
        .map_err(|e| GenomeError::InvalidRegion(format!("{region_str}: {e}")))?;

    let build_repo = || -> io::Result<fasta::Repository> {
        Ok(reference_path
            .map(|p| fasta::io::indexed_reader::Builder::default().build_from_path(p))
            .transpose()?
            .map(IndexedReader::new)
            .map(fasta::Repository::new)
            .unwrap_or_default())
    };

    let repo = build_repo()?;

    let mut reader = cram::io::indexed_reader::Builder::default()
        .set_reference_sequence_repository(repo)
        .build_from_path(cram_path)
        .map_err(|e| match e.kind() {
            io::ErrorKind::NotFound => {
                let crai = cram_path.with_extension("cram.crai");
                GenomeError::IndexNotFound(crai)
            }
            _ => GenomeError::Io(e),
        })?;

    let header = match reader.read_header() {
        Ok(h) => h,
        Err(e) if e.to_string().contains("duplicate program") => {
            let raw = read_cram_raw_header(cram_path)?;
            let header = parse_header_lenient(&raw)?;
            let repo = build_repo()?;
            reader = cram::io::indexed_reader::Builder::default()
                .set_reference_sequence_repository(repo)
                .build_from_path(cram_path)
                .map_err(GenomeError::Io)?;
            header
        }
        Err(e) => return Err(GenomeError::ParseError(format!("CRAM header: {e}"))),
    };

    let mut reads = Vec::new();

    for result in reader.query(&header, &region)? {
        match result {
            Ok(record) => {
                if let Some(r) = extract_cram_read(&record) {
                    reads.push(r);
                }
            }
            Err(_) => continue,
        }
    }

    Ok(reads)
}

/// Extract an [`AlignedRead`] from a CRAM `RecordBuf`.
fn extract_cram_read(record: &sam::alignment::RecordBuf) -> Option<AlignedRead> {
    let flags = record.flags();

    if flags.is_unmapped() {
        return None;
    }

    let start = record.alignment_start().map(|p| usize::from(p) as u64)?;

    let end = record
        .alignment_end()
        .map(|p| usize::from(p) as u64)
        .unwrap_or(start);

    let name = record.name().map(|n| n.to_string()).unwrap_or_default();

    let mq = record.mapping_quality().map(u8::from);
    let hp = hp_from_record_buf(record);

    let indels = extract_indels_from_cigar_buf(record.cigar(), start);

    Some(AlignedRead {
        name,
        start,
        end,
        is_reverse: flags.is_reverse_complemented(),
        mapping_quality: mq,
        haplotype: hp,
        flags: u16::from(flags),
        indels,
    })
}

/// Walk a CIGAR from a RecordBuf and collect insertion/deletion events.
fn extract_indels_from_cigar_buf(
    cigar: &noodles::sam::alignment::record_buf::Cigar,
    alignment_start: u64,
) -> Vec<Indel> {
    use noodles::sam::alignment::record::cigar::op::Kind;

    let mut indels = Vec::new();
    let mut ref_pos = alignment_start;

    for op in cigar.as_ref() {
        let len = op.len();
        match op.kind() {
            Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                ref_pos += len as u64;
            }
            Kind::Insertion => {
                indels.push(Indel {
                    ref_pos,
                    length: len as u32,
                    kind: IndelKind::Insertion,
                });
            }
            Kind::Deletion => {
                indels.push(Indel {
                    ref_pos,
                    length: len as u32,
                    kind: IndelKind::Deletion,
                });
                ref_pos += len as u64;
            }
            Kind::SoftClip | Kind::HardClip | Kind::Pad => {}
            Kind::Skip => {
                ref_pos += len as u64;
            }
        }
    }

    indels
}

/// Read the HP tag from a fully-parsed `RecordBuf`.
fn hp_from_record_buf(record: &sam::alignment::RecordBuf) -> Option<u8> {
    use noodles::sam::alignment::record_buf::data::field::Value;
    let tag = Tag::new(b'H', b'P');
    record.data().get(&tag).and_then(|v| match v {
        Value::UInt8(n) => Some(*n),
        Value::Int8(n) => Some(*n as u8),
        Value::UInt16(n) => Some(*n as u8),
        Value::Int16(n) => Some(*n as u8),
        Value::UInt32(n) => Some(*n as u8),
        Value::Int32(n) => Some(*n as u8),
        _ => None,
    })
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    fn toy_bam() -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../resources/toy_dataset/toy_reads.bam")
    }

    // -- BAM unit tests --

    /// The toy BAM uses region-sliced contig names (e.g.
    /// "chr1:112064095-112277008").  Querying by simple "chr1" fails because
    /// that contig doesn't exist; tests that need reads must use the full
    /// sliced name with local coordinates (1-based within the slice).
    const TOY_CONTIG: &str = "chr1:112064095-112277008";

    #[test]
    fn test_query_bam_returns_reads() {
        let bam = toy_bam();
        if !bam.exists() {
            eprintln!("skipping: toy BAM not found");
            return;
        }
        // Query the first sliced contig with local coordinates.
        let result = query_bam(&bam, &format!("{TOY_CONTIG}:1-212914"));
        match result {
            Ok(reads) => {
                assert!(!reads.is_empty(), "expected reads in toy contig");
                for r in &reads {
                    assert!(r.start >= 1, "start should be 1-based");
                    assert!(r.end >= r.start, "end >= start");
                }
            }
            Err(e) => panic!("unexpected error: {e}"),
        }
    }

    #[test]
    fn test_query_bam_region_correctness() {
        let bam = toy_bam();
        if !bam.exists() {
            return;
        }
        // Query a small window inside the first sliced contig.
        let result = query_bam(&bam, &format!("{TOY_CONTIG}:1-100000"));
        match result {
            Ok(reads) => {
                for r in &reads {
                    assert!(
                        r.end >= 1 && r.start <= 100000,
                        "read {}:{}-{} does not overlap query",
                        r.name,
                        r.start,
                        r.end
                    );
                }
            }
            Err(e) => panic!("unexpected error: {e}"),
        }
    }

    #[test]
    fn test_query_bam_reads_have_valid_fields() {
        let bam = toy_bam();
        if !bam.exists() {
            return;
        }
        let result = query_bam(&bam, &format!("{TOY_CONTIG}:1-212914"));
        match result {
            Ok(reads) => {
                for r in &reads {
                    assert!(!r.name.is_empty(), "read name should not be empty");
                    // mapping quality can be None (255 in BAM → None)
                    if let Some(hp) = r.haplotype {
                        assert!(hp <= 2, "HP tag should be 0, 1 or 2, got {hp}");
                    }
                }
            }
            Err(e) => panic!("unexpected error: {e}"),
        }
    }

    #[test]
    fn test_query_bam_empty_region() {
        let bam = toy_bam();
        if !bam.exists() {
            return;
        }
        // Query a region with no reads (chrUn doesn't exist in the BAM)
        let result = query_bam(&bam, "chrUn:1-100");
        if let Ok(reads) = result {
            assert!(reads.is_empty(), "no reads expected in chrUn");
        }
        // Either empty or error is acceptable for non-existent contig
    }

    // -- Error handling tests --

    #[test]
    fn test_query_bam_missing_file() {
        let result = query_bam(Path::new("/nonexistent/file.bam"), "chr1:1-100");
        assert!(result.is_err());
        let msg = format!("{}", result.unwrap_err());
        assert!(
            msg.contains("not found") || msg.contains("No such file"),
            "expected not-found error, got: {msg}"
        );
    }

    #[test]
    fn test_query_bam_invalid_region() {
        let bam = toy_bam();
        if !bam.exists() {
            return;
        }
        // A malformed interval that the region parser should reject
        let result = query_bam(&bam, "chr1:end-start");
        assert!(
            result.is_err(),
            "malformed interval should produce an error"
        );
    }

    #[test]
    fn test_query_bam_missing_index() {
        // Create a temp copy of the BAM without an index
        let tmp = tempfile::tempdir().unwrap();
        let bam_copy = tmp.path().join("no_index.bam");
        std::fs::copy(toy_bam(), &bam_copy).ok(); // may fail if toy absent
        if !bam_copy.exists() {
            return;
        }
        let result = query_bam(&bam_copy, "chr1:1-100");
        assert!(result.is_err(), "should fail without index");
    }

    // -- CRAM error tests (no CRAM in toy dataset, test error paths) --

    #[test]
    fn test_query_cram_missing_file() {
        let result = query_cram(Path::new("/nonexistent/file.cram"), None, "chr1:1-100");
        assert!(result.is_err());
    }

    // -- Lenient header parsing tests --

    #[test]
    fn test_deduplicate_pg_ids_removes_duplicates() {
        let header = "\
@HD\tVN:1.6\tSO:coordinate\n\
@SQ\tSN:chr1\tLN:100000\n\
@PG\tID:minimap2\tPN:minimap2\n\
@PG\tID:samtools\tPN:samtools\tPP:minimap2\n\
@PG\tID:samtools.1\tPN:samtools\tPP:samtools\n\
@PG\tID:samtools\tPN:samtools\tPP:samtools.1\n\
@PG\tID:samtools\tPN:samtools\tPP:samtools.1\n";

        let result = deduplicate_pg_ids(header);
        let pg_lines: Vec<&str> = result.lines().filter(|l| l.starts_with("@PG")).collect();
        // Only 3 unique IDs: minimap2, samtools, samtools.1
        assert_eq!(pg_lines.len(), 3);
        assert!(pg_lines[0].contains("ID:minimap2"));
        assert!(pg_lines[1].contains("ID:samtools\t"));
        assert!(pg_lines[2].contains("ID:samtools.1"));
    }

    #[test]
    fn test_deduplicate_pg_ids_no_duplicates() {
        let header = "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:100\n@PG\tID:bwa\tPN:bwa\n";
        let result = deduplicate_pg_ids(header);
        assert_eq!(result.lines().count(), header.lines().count());
    }

    #[test]
    fn test_parse_header_lenient_succeeds() {
        let header = "\
@HD\tVN:1.6\tSO:coordinate\n\
@SQ\tSN:chr1\tLN:100000\n\
@PG\tID:samtools\tPN:samtools\n\
@PG\tID:samtools\tPN:samtools\n";

        // Strict parsing would fail on duplicate PG ID
        assert!(header.parse::<sam::Header>().is_err());
        // Lenient parsing should succeed
        let h = parse_header_lenient(header).unwrap();
        assert!(h.reference_sequences().contains_key(&b"chr1"[..]));
    }

    #[test]
    fn test_read_bam_raw_header_toy_dataset() {
        let bam = toy_bam();
        if !bam.exists() {
            return;
        }
        let raw = read_bam_raw_header(&bam).unwrap();
        // Must contain @HD and at least one @SQ line
        assert!(raw.contains("@HD"));
        assert!(raw.contains("@SQ"));
    }

    #[test]
    fn test_toy_bam_lenient_parse() {
        let bam = toy_bam();
        if !bam.exists() {
            return;
        }
        let raw = read_bam_raw_header(&bam).unwrap();
        let header = parse_header_lenient(&raw).unwrap();
        assert!(
            !header.reference_sequences().is_empty(),
            "header should contain reference sequences"
        );
    }
}
