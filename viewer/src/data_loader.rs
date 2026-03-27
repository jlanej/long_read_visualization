use std::path::PathBuf;
use std::sync::Arc;
use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};
use std::sync::mpsc;
use std::thread;

use crate::genome;
use crate::genome::FastaSequence;
use crate::genome::pileup::{self, PileupRow};
use crate::region::GenomicRegion;
use crate::region_cache::{CachedAssemblyTrack, CachedPanelData, CachedRegion};

// ---------------------------------------------------------------------------
// Load request & result types
// ---------------------------------------------------------------------------

/// Snapshot of all paths needed for one load operation.
///
/// This is `Send + 'static` so it can be shipped to a background thread.
#[derive(Debug, Clone)]
pub struct LoadPaths {
    pub reads_bam: Option<PathBuf>,
    pub reference_fasta: Option<PathBuf>,
    pub hap1_fasta: Option<PathBuf>,
    pub hap2_fasta: Option<PathBuf>,
    pub reads_to_hap1_bam: Option<PathBuf>,
    pub reads_to_hap2_bam: Option<PathBuf>,
    pub cram_ref: Option<PathBuf>,
    // Cross-alignment BAMs
    pub hap1_to_ref_bam: Option<PathBuf>,
    pub hap2_to_ref_bam: Option<PathBuf>,
    pub ref_to_hap1_bam: Option<PathBuf>,
    pub ref_to_hap2_bam: Option<PathBuf>,
    pub hap1_to_hap2_bam: Option<PathBuf>,
    pub hap2_to_hap1_bam: Option<PathBuf>,
}

/// A request to load data for one region.
#[derive(Debug, Clone)]
pub struct LoadRequest {
    /// Unique monotonic ID — used to detect stale results.
    pub id: u64,
    /// Reference region to load.
    pub ref_region: GenomicRegion,
    /// Haplotype 1 region (if available).
    pub hap1_region: Option<GenomicRegion>,
    /// Haplotype 2 region (if available).
    pub hap2_region: Option<GenomicRegion>,
    /// Whether to sort reads by haplotype.
    pub sort_by_haplotype: bool,
    /// File paths snapshot.
    pub paths: LoadPaths,
    /// Cache key for storing the result.
    pub cache_key: String,
}

/// Result of a background load operation.
#[derive(Debug)]
#[allow(dead_code)]
pub struct LoadResult {
    /// Request ID — used to match results with requests.
    pub id: u64,
    /// Cache key for the loaded region.
    pub cache_key: String,
    /// The loaded data.
    pub data: CachedRegion,
}

// ---------------------------------------------------------------------------
// Background data loader
// ---------------------------------------------------------------------------

/// Manages a background thread for non-blocking BAM/CRAM/FASTA loading.
///
/// The UI thread sends [`LoadRequest`]s via `request_tx` and polls
/// [`LoadResult`]s from `result_rx` every frame.  A cancellation token
/// (`cancel_flag`) tells the background thread to abandon stale work.
pub struct DataLoader {
    request_tx: mpsc::Sender<LoadRequest>,
    result_rx: mpsc::Receiver<LoadResult>,
    cancel_flag: Arc<AtomicBool>,
    current_id: Arc<AtomicU64>,
    /// egui context for requesting repaints when data arrives.
    ctx: egui::Context,
}

impl DataLoader {
    /// Spawn a new background loader thread tied to the given egui context.
    pub fn new(ctx: egui::Context) -> Self {
        let (request_tx, request_rx) = mpsc::channel::<LoadRequest>();
        let (result_tx, result_rx) = mpsc::channel::<LoadResult>();
        let cancel_flag = Arc::new(AtomicBool::new(false));
        let current_id = Arc::new(AtomicU64::new(0));

        let cancel = cancel_flag.clone();
        let id_tracker = current_id.clone();
        let repaint_ctx = ctx.clone();

        thread::spawn(move || {
            loader_thread(request_rx, result_tx, cancel, id_tracker, repaint_ctx);
        });

        Self {
            request_tx,
            result_rx,
            cancel_flag,
            current_id,
            ctx,
        }
    }

    /// Submit a new load request.  Any in-flight request with an older ID is
    /// implicitly cancelled.
    pub fn submit(&self, request: LoadRequest) {
        // Set the latest expected ID so the thread can skip stale work.
        self.current_id.store(request.id, Ordering::SeqCst);
        // Clear cancellation flag for the new request.
        self.cancel_flag.store(false, Ordering::SeqCst);
        let _ = self.request_tx.send(request);
        self.ctx.request_repaint();
    }

    /// Cancel any currently in-flight request.
    pub fn cancel(&self) {
        self.cancel_flag.store(true, Ordering::SeqCst);
    }

    /// Poll for a completed result.  Non-blocking — returns `None` if
    /// nothing is ready yet.
    pub fn try_recv(&self) -> Option<LoadResult> {
        self.result_rx.try_recv().ok()
    }

    /// Generate the next monotonically increasing request ID.
    pub fn next_id(&self) -> u64 {
        self.current_id.fetch_add(1, Ordering::SeqCst) + 1
    }
}

impl std::fmt::Debug for DataLoader {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("DataLoader")
            .field("current_id", &self.current_id.load(Ordering::Relaxed))
            .finish()
    }
}

// ---------------------------------------------------------------------------
// Background thread logic
// ---------------------------------------------------------------------------

use eframe::egui;

fn loader_thread(
    rx: mpsc::Receiver<LoadRequest>,
    tx: mpsc::Sender<LoadResult>,
    cancel: Arc<AtomicBool>,
    current_id: Arc<AtomicU64>,
    ctx: egui::Context,
) {
    while let Ok(req) = rx.recv() {
        // Skip if a newer request has already been submitted.
        if req.id < current_id.load(Ordering::SeqCst) {
            continue;
        }
        if cancel.load(Ordering::SeqCst) {
            continue;
        }

        let result = execute_load(&req, &cancel);

        // Only send if still the current request.
        if !cancel.load(Ordering::SeqCst) && req.id >= current_id.load(Ordering::SeqCst) {
            let _ = tx.send(LoadResult {
                id: req.id,
                cache_key: req.cache_key.clone(),
                data: result,
            });
            ctx.request_repaint();
        }
    }
}

/// Execute the actual data loading (BAM queries, FASTA loads, packing).
fn execute_load(req: &LoadRequest, cancel: &Arc<AtomicBool>) -> CachedRegion {
    let mut ref_data = CachedPanelData::default();
    let mut hap1_data = CachedPanelData::default();
    let mut hap2_data = CachedPanelData::default();
    let mut status_message = String::new();

    // -- Reference panel --
    if !cancel.load(Ordering::SeqCst)
        && let Some(reads_path) = &req.paths.reads_bam
    {
        let region_str = req.ref_region.to_string();
        let cram_ref = req
            .paths
            .cram_ref
            .as_deref()
            .or(req.paths.reference_fasta.as_deref());
        let result = if reads_path.extension().is_some_and(|e| e == "cram") {
            genome::bam::query_cram(reads_path, cram_ref, &region_str)
        } else {
            genome::bam::query_bam(reads_path, &region_str)
        };
        match result {
            Ok(reads) => {
                let n_reads = reads.len();
                let n_reverse = reads.iter().filter(|r| r.is_reverse).count();
                let avg_mapq = {
                    let mapq_iter: Vec<u64> = reads
                        .iter()
                        .filter_map(|r| r.mapping_quality)
                        .map(u64::from)
                        .collect();
                    let n = mapq_iter.len() as u64;
                    mapq_iter.into_iter().sum::<u64>().checked_div(n)
                };
                let n_flagged = reads.iter().filter(|r| r.flags & 0x100 != 0).count();

                ref_data.rows = pack_for_config(reads, req.sort_by_haplotype);

                status_message = format!(
                    "{n_reads} reads ({n_reverse} rev, {n_flagged} secondary, avg MAPQ {})",
                    avg_mapq.map_or("N/A".to_string(), |q| q.to_string()),
                );
            }
            Err(e) => {
                status_message = format!("Reads error: {e}");
            }
        }
    }
    if !cancel.load(Ordering::SeqCst)
        && let Some(fasta_path) = &req.paths.reference_fasta
    {
        ref_data.sequence = load_fasta_for_region(fasta_path, &req.ref_region);
    }
    // Assembly tracks for reference panel
    if !cancel.load(Ordering::SeqCst) {
        let region_str = req.ref_region.to_string();
        ref_data.assembly_tracks = load_assembly_tracks_cached(
            &[
                (
                    "Hap1 → Ref",
                    [60, 160, 80],
                    req.paths.hap1_to_ref_bam.as_deref(),
                ),
                (
                    "Hap2 → Ref",
                    [180, 100, 60],
                    req.paths.hap2_to_ref_bam.as_deref(),
                ),
            ],
            &region_str,
        );
    }

    // -- Haplotype 1 panel --
    if let Some(region) = &req.hap1_region {
        if !cancel.load(Ordering::SeqCst) {
            if let Some(bam_path) = &req.paths.reads_to_hap1_bam {
                let region_str = region.to_string();
                if let Ok(reads) = genome::bam::query_bam(bam_path, &region_str) {
                    hap1_data.rows = pack_for_config(reads, req.sort_by_haplotype);
                }
            }
            if let Some(fasta_path) = &req.paths.hap1_fasta {
                hap1_data.sequence = load_fasta_for_region(fasta_path, region);
            }
        }
        if !cancel.load(Ordering::SeqCst) {
            let region_str = region.to_string();
            hap1_data.assembly_tracks = load_assembly_tracks_cached(
                &[
                    (
                        "Ref → Hap1",
                        [70, 130, 180],
                        req.paths.ref_to_hap1_bam.as_deref(),
                    ),
                    (
                        "Hap2 → Hap1",
                        [180, 100, 60],
                        req.paths.hap2_to_hap1_bam.as_deref(),
                    ),
                ],
                &region_str,
            );
        }
    }

    // -- Haplotype 2 panel --
    if let Some(region) = &req.hap2_region {
        if !cancel.load(Ordering::SeqCst) {
            if let Some(bam_path) = &req.paths.reads_to_hap2_bam {
                let region_str = region.to_string();
                if let Ok(reads) = genome::bam::query_bam(bam_path, &region_str) {
                    hap2_data.rows = pack_for_config(reads, req.sort_by_haplotype);
                }
            }
            if let Some(fasta_path) = &req.paths.hap2_fasta {
                hap2_data.sequence = load_fasta_for_region(fasta_path, region);
            }
        }
        if !cancel.load(Ordering::SeqCst) {
            let region_str = region.to_string();
            hap2_data.assembly_tracks = load_assembly_tracks_cached(
                &[
                    (
                        "Ref → Hap2",
                        [70, 130, 180],
                        req.paths.ref_to_hap2_bam.as_deref(),
                    ),
                    (
                        "Hap1 → Hap2",
                        [60, 160, 80],
                        req.paths.hap1_to_hap2_bam.as_deref(),
                    ),
                ],
                &region_str,
            );
        }
    }

    CachedRegion {
        ref_data,
        hap1_data,
        hap2_data,
        status_message,
    }
}

// ---------------------------------------------------------------------------
// Helpers (duplicated from app.rs to avoid circular dependency — these are
// pure functions that only depend on genome::* types)
// ---------------------------------------------------------------------------

fn pack_for_config(reads: Vec<genome::AlignedRead>, sort_by_haplotype: bool) -> Vec<PileupRow> {
    if sort_by_haplotype {
        pileup::pack_reads_by_haplotype(reads)
    } else {
        pileup::pack_reads(reads)
    }
}

fn load_fasta_for_region(
    fasta_path: &std::path::Path,
    region: &GenomicRegion,
) -> Option<FastaSequence> {
    if region.end < region.start {
        return None;
    }
    let full_name = region.to_string();
    if let Ok(seqs) = genome::fasta::list_sequences(fasta_path) {
        for (name, len) in &seqs {
            if *name == full_name {
                let end = (*len).min(region.end - region.start + 1);
                return genome::fasta::query_fasta_region(fasta_path, name, 1, end).ok();
            }
            if *name == region.chrom {
                return genome::fasta::query_fasta_region(
                    fasta_path,
                    name,
                    region.start,
                    region.end,
                )
                .ok();
            }
            if name.starts_with(&region.chrom)
                && name[region.chrom.len()..].starts_with(|c: char| !c.is_ascii_alphanumeric())
            {
                let end = (*len).min(region.end - region.start + 1);
                return genome::fasta::query_fasta_region(fasta_path, name, 1, end).ok();
            }
        }
    }
    genome::fasta::query_fasta(fasta_path, &region.to_string()).ok()
}

fn load_assembly_tracks_cached(
    track_specs: &[(&str, [u8; 3], Option<&std::path::Path>)],
    region_str: &str,
) -> Vec<CachedAssemblyTrack> {
    let mut tracks = Vec::new();
    for &(label, color, bam_path) in track_specs {
        let Some(path) = bam_path else { continue };
        if let Ok(reads) = genome::bam::query_bam(path, region_str)
            && !reads.is_empty()
        {
            let rows = pileup::pack_reads(reads);
            tracks.push(CachedAssemblyTrack {
                label: label.to_string(),
                color,
                rows,
            });
        }
    }
    tracks
}

/// Build a cache key from the region strings and sort mode.
pub fn make_cache_key(
    ref_region: &GenomicRegion,
    hap1_region: Option<&GenomicRegion>,
    hap2_region: Option<&GenomicRegion>,
    sort_by_haplotype: bool,
) -> String {
    let h1 = hap1_region
        .map(|r| r.to_string())
        .unwrap_or_else(|| "-".to_string());
    let h2 = hap2_region
        .map(|r| r.to_string())
        .unwrap_or_else(|| "-".to_string());
    let sort_tag = if sort_by_haplotype { "hp" } else { "pos" };
    format!("{ref_region}|{h1}|{h2}|{sort_tag}")
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_make_cache_key_with_all_regions() {
        let ref_r = GenomicRegion {
            chrom: "chr1".to_string(),
            start: 100,
            end: 200,
        };
        let h1 = GenomicRegion {
            chrom: "ctg1".to_string(),
            start: 500,
            end: 600,
        };
        let h2 = GenomicRegion {
            chrom: "ctg2".to_string(),
            start: 800,
            end: 900,
        };
        let key = make_cache_key(&ref_r, Some(&h1), Some(&h2), false);
        assert_eq!(key, "chr1:100-200|ctg1:500-600|ctg2:800-900|pos");
    }

    #[test]
    fn test_make_cache_key_no_hap_regions() {
        let ref_r = GenomicRegion {
            chrom: "chr1".to_string(),
            start: 100,
            end: 200,
        };
        let key = make_cache_key(&ref_r, None, None, false);
        assert_eq!(key, "chr1:100-200|-|-|pos");
    }

    #[test]
    fn test_make_cache_key_haplotype_sorted() {
        let ref_r = GenomicRegion {
            chrom: "chr1".to_string(),
            start: 100,
            end: 200,
        };
        let key = make_cache_key(&ref_r, None, None, true);
        assert_eq!(key, "chr1:100-200|-|-|hp");
    }

    #[test]
    fn test_load_paths_clone() {
        let paths = LoadPaths {
            reads_bam: Some(PathBuf::from("/reads.bam")),
            reference_fasta: None,
            hap1_fasta: None,
            hap2_fasta: None,
            reads_to_hap1_bam: None,
            reads_to_hap2_bam: None,
            cram_ref: None,
            hap1_to_ref_bam: None,
            hap2_to_ref_bam: None,
            ref_to_hap1_bam: None,
            ref_to_hap2_bam: None,
            hap1_to_hap2_bam: None,
            hap2_to_hap1_bam: None,
        };
        let cloned = paths.clone();
        assert_eq!(cloned.reads_bam, Some(PathBuf::from("/reads.bam")));
    }

    #[test]
    fn test_pack_for_config_position() {
        let reads = vec![
            genome::AlignedRead {
                name: "r1".to_string(),
                start: 100,
                end: 200,
                is_reverse: false,
                mapping_quality: Some(60),
                haplotype: None,
                flags: 0,
                indels: vec![],
                mismatches: vec![],
                soft_clips: vec![],
            },
            genome::AlignedRead {
                name: "r2".to_string(),
                start: 150,
                end: 250,
                is_reverse: false,
                mapping_quality: Some(60),
                haplotype: None,
                flags: 0,
                indels: vec![],
                mismatches: vec![],
                soft_clips: vec![],
            },
        ];
        let rows = pack_for_config(reads, false);
        assert!(!rows.is_empty());
    }

    #[test]
    fn test_pack_for_config_haplotype() {
        let reads = vec![
            genome::AlignedRead {
                name: "r1".to_string(),
                start: 100,
                end: 200,
                is_reverse: false,
                mapping_quality: Some(60),
                haplotype: Some(1),
                flags: 0,
                indels: vec![],
                mismatches: vec![],
                soft_clips: vec![],
            },
            genome::AlignedRead {
                name: "r2".to_string(),
                start: 150,
                end: 250,
                is_reverse: false,
                mapping_quality: Some(60),
                haplotype: Some(2),
                flags: 0,
                indels: vec![],
                mismatches: vec![],
                soft_clips: vec![],
            },
        ];
        let rows = pack_for_config(reads, true);
        assert!(!rows.is_empty());
    }

    // -- Cancellation / stale-request logic ---------------------------------

    #[test]
    fn test_cancel_flag_skips_work() {
        let cancel = Arc::new(AtomicBool::new(true));
        let req = LoadRequest {
            id: 1,
            ref_region: GenomicRegion {
                chrom: "chr1".to_string(),
                start: 100,
                end: 200,
            },
            hap1_region: None,
            hap2_region: None,
            sort_by_haplotype: false,
            paths: LoadPaths {
                reads_bam: None,
                reference_fasta: None,
                hap1_fasta: None,
                hap2_fasta: None,
                reads_to_hap1_bam: None,
                reads_to_hap2_bam: None,
                cram_ref: None,
                hap1_to_ref_bam: None,
                hap2_to_ref_bam: None,
                ref_to_hap1_bam: None,
                ref_to_hap2_bam: None,
                hap1_to_hap2_bam: None,
                hap2_to_hap1_bam: None,
            },
            cache_key: "test".to_string(),
        };
        // Even with cancellation set, execute_load should not panic — it
        // just skips the actual file I/O and returns empty data.
        let result = execute_load(&req, &cancel);
        assert!(result.ref_data.rows.is_empty());
        assert!(result.status_message.is_empty());
    }

    #[test]
    fn test_execute_load_no_files() {
        let cancel = Arc::new(AtomicBool::new(false));
        let req = LoadRequest {
            id: 1,
            ref_region: GenomicRegion {
                chrom: "chr1".to_string(),
                start: 100,
                end: 200,
            },
            hap1_region: None,
            hap2_region: None,
            sort_by_haplotype: false,
            paths: LoadPaths {
                reads_bam: None,
                reference_fasta: None,
                hap1_fasta: None,
                hap2_fasta: None,
                reads_to_hap1_bam: None,
                reads_to_hap2_bam: None,
                cram_ref: None,
                hap1_to_ref_bam: None,
                hap2_to_ref_bam: None,
                ref_to_hap1_bam: None,
                ref_to_hap2_bam: None,
                hap1_to_hap2_bam: None,
                hap2_to_hap1_bam: None,
            },
            cache_key: "test".to_string(),
        };
        let result = execute_load(&req, &cancel);
        // No files → empty data, no error
        assert!(result.ref_data.rows.is_empty());
        assert!(result.hap1_data.rows.is_empty());
        assert!(result.hap2_data.rows.is_empty());
    }

    #[test]
    fn test_load_assembly_tracks_cached_no_bams() {
        let tracks = load_assembly_tracks_cached(
            &[
                ("Hap1 → Ref", [60, 160, 80], None),
                ("Hap2 → Ref", [180, 100, 60], None),
            ],
            "chr1:100-200",
        );
        assert!(tracks.is_empty());
    }
}
