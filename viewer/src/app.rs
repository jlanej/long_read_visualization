use std::collections::HashSet;
use std::path::PathBuf;
use std::time::Instant;

use eframe::egui;

use crate::data_loader::{self, DataLoader, LoadPaths, LoadRequest};
use crate::dot_plot::{self, DEFAULT_K, DotPlotComparison, DotPlotResult};
use crate::genome::pileup::{self, PileupDisplayConfig, PileupRow, ReadRect};
use crate::genome::{self, FastaSequence};
use crate::panel_sync::{PanelId, PanelSyncManager};
use crate::region::RegionNavigator;
use crate::region_cache::{
    CachedAssemblyTrack, CachedPanelData, CachedRegion, DEFAULT_CACHE_CAPACITY, RegionCache,
};
use crate::ruler;

// ---------------------------------------------------------------------------
// Data file paths
// ---------------------------------------------------------------------------

/// Paths to data files for the viewer.
#[derive(Debug, Clone, Default)]
pub struct DataPaths {
    /// Reference genome FASTA (indexed, .fa.gz with .fai/.gzi).
    pub reference_fasta: Option<PathBuf>,
    /// Haplotype 1 assembly FASTA.
    pub hap1_fasta: Option<PathBuf>,
    /// Haplotype 2 assembly FASTA.
    pub hap2_fasta: Option<PathBuf>,
    /// Reads BAM/CRAM file (indexed).
    pub reads_bam: Option<PathBuf>,
    /// Optional coordinate mapping index (gzipped JSON) — legacy single-index.
    pub coordinate_index: Option<PathBuf>,
    /// Optional path to the manifest/regions JSON loaded from config.
    pub regions: Option<PathBuf>,
    /// Hap1-to-ref coordinate mapping index (gzipped JSON).
    pub hap1_coord_index: Option<PathBuf>,
    /// Hap2-to-ref coordinate mapping index (gzipped JSON).
    pub hap2_coord_index: Option<PathBuf>,
    /// Reads aligned to haplotype 1 assembly (BAM, indexed).
    pub reads_to_hap1_bam: Option<PathBuf>,
    /// Reads aligned to haplotype 2 assembly (BAM, indexed).
    pub reads_to_hap2_bam: Option<PathBuf>,
    /// Reference FASTA used for CRAM encoding (when different from panel ref).
    pub cram_ref: Option<PathBuf>,
    // -- Cross-alignment (assembly) BAMs --
    /// Hap1 assembly aligned to reference (BAM, indexed).
    pub hap1_to_ref_bam: Option<PathBuf>,
    /// Hap2 assembly aligned to reference (BAM, indexed).
    pub hap2_to_ref_bam: Option<PathBuf>,
    /// Reference aligned to hap1 assembly (BAM, indexed).
    pub ref_to_hap1_bam: Option<PathBuf>,
    /// Reference aligned to hap2 assembly (BAM, indexed).
    pub ref_to_hap2_bam: Option<PathBuf>,
    /// Hap1 aligned to hap2 (BAM, indexed) — cross-haplotype track in hap2 panel.
    pub hap1_to_hap2_bam: Option<PathBuf>,
    /// Hap2 aligned to hap1 (BAM, indexed) — cross-haplotype track in hap1 panel.
    pub hap2_to_hap1_bam: Option<PathBuf>,
}

impl DataPaths {
    /// Load data paths from a server config TSV file (the format produced by
    /// `generate_server_config.py` and consumed by the Python `app.py`).
    ///
    /// Expected TSV columns: `sample_id`, `output_dir`, `reference`,
    /// `hap1_assembly`, `hap2_assembly`, `reads_bam` (optional),
    /// `regions` (optional), `cram_ref` (optional).
    ///
    /// If the TSV contains multiple samples, only the first is loaded.
    pub fn from_tsv(tsv_path: &std::path::Path) -> Result<Self, String> {
        let content = std::fs::read_to_string(tsv_path)
            .map_err(|e| format!("Failed to read config TSV: {e}"))?;

        let mut header: Option<Vec<String>> = None;
        let mut first_row: Option<std::collections::HashMap<String, String>> = None;

        for line in content.lines() {
            let line = line.trim();
            if line.is_empty() {
                continue;
            }
            if line.starts_with('#') {
                if header.is_none() {
                    // Header line starting with '#'
                    header = Some(
                        line.trim_start_matches('#')
                            .trim()
                            .split('\t')
                            .map(|s| s.trim().to_string())
                            .collect(),
                    );
                }
                continue;
            }
            if header.is_none() {
                // First non-comment line is the header
                header = Some(line.split('\t').map(|s| s.trim().to_string()).collect());
                continue;
            }
            // Data row
            if first_row.is_none() {
                let cols: Vec<&str> = line.split('\t').collect();
                let h = header.as_ref().unwrap();
                let mut row = std::collections::HashMap::new();
                for (i, col_name) in h.iter().enumerate() {
                    let val = cols.get(i).unwrap_or(&"").trim().to_string();
                    row.insert(col_name.clone(), val);
                }
                first_row = Some(row);
                break; // Only load first sample
            }
        }

        let row = first_row.ok_or_else(|| "No data rows found in TSV".to_string())?;

        let non_empty = |key: &str| -> Option<PathBuf> {
            row.get(key).filter(|v| !v.is_empty()).map(PathBuf::from)
        };

        // Discover preprocessing output files from output_dir + sample_id.
        let output_dir = non_empty("output_dir");
        let sample_id = non_empty("sample_id");

        let discover = |suffix: &str| -> Option<PathBuf> {
            let dir = output_dir.as_ref()?;
            let sid = sample_id.as_ref()?.to_string_lossy().to_string();
            let candidate = dir.join(format!("{sid}{suffix}"));
            if candidate.is_file() {
                Some(candidate)
            } else {
                None
            }
        };

        Ok(Self {
            reference_fasta: non_empty("reference"),
            hap1_fasta: non_empty("hap1_assembly"),
            hap2_fasta: non_empty("hap2_assembly"),
            reads_bam: non_empty("reads_bam"),
            coordinate_index: None,
            regions: non_empty("regions"),
            hap1_coord_index: discover("_hap1_to_ref.mapping.json.gz"),
            hap2_coord_index: discover("_hap2_to_ref.mapping.json.gz"),
            reads_to_hap1_bam: discover("_reads_to_hap1.bam"),
            reads_to_hap2_bam: discover("_reads_to_hap2.bam"),
            cram_ref: non_empty("cram_ref"),
            // Cross-alignment BAMs
            hap1_to_ref_bam: discover("_hap1_to_ref.bam"),
            hap2_to_ref_bam: discover("_hap2_to_ref.bam"),
            ref_to_hap1_bam: discover("_ref_to_hap1.bam"),
            ref_to_hap2_bam: discover("_ref_to_hap2.bam"),
            hap1_to_hap2_bam: discover("_hap1_to_hap2.bam"),
            hap2_to_hap1_bam: discover("_hap2_to_hap1.bam"),
        })
    }
}

/// Parse ALL samples from a TSV config file, returning (sample_id, DataPaths)
/// for each data row.
pub fn parse_all_samples(tsv_path: &std::path::Path) -> Result<Vec<(String, DataPaths)>, String> {
    let content =
        std::fs::read_to_string(tsv_path).map_err(|e| format!("Failed to read config TSV: {e}"))?;

    let mut header: Option<Vec<String>> = None;
    let mut samples = Vec::new();

    for line in content.lines() {
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        if line.starts_with('#') {
            if header.is_none() {
                header = Some(
                    line.trim_start_matches('#')
                        .trim()
                        .split('\t')
                        .map(|s| s.trim().to_string())
                        .collect(),
                );
            }
            continue;
        }
        if header.is_none() {
            header = Some(line.split('\t').map(|s| s.trim().to_string()).collect());
            continue;
        }
        // Data row
        let cols: Vec<&str> = line.split('\t').collect();
        let h = header.as_ref().unwrap();
        let mut row = std::collections::HashMap::new();
        for (i, col_name) in h.iter().enumerate() {
            let val = cols.get(i).unwrap_or(&"").trim().to_string();
            row.insert(col_name.clone(), val);
        }

        let non_empty = |key: &str| -> Option<PathBuf> {
            row.get(key).filter(|v| !v.is_empty()).map(PathBuf::from)
        };

        let sample_id = non_empty("sample_id")
            .map(|p| p.to_string_lossy().to_string())
            .unwrap_or_default();

        let output_dir = non_empty("output_dir");
        let sid_for_discover = non_empty("sample_id");

        let discover = |suffix: &str| -> Option<PathBuf> {
            let dir = output_dir.as_ref()?;
            let sid = sid_for_discover.as_ref()?.to_string_lossy().to_string();
            let candidate = dir.join(format!("{sid}{suffix}"));
            if candidate.is_file() {
                Some(candidate)
            } else {
                None
            }
        };

        let dp = DataPaths {
            reference_fasta: non_empty("reference"),
            hap1_fasta: non_empty("hap1_assembly"),
            hap2_fasta: non_empty("hap2_assembly"),
            reads_bam: non_empty("reads_bam"),
            coordinate_index: None,
            regions: non_empty("regions"),
            hap1_coord_index: discover("_hap1_to_ref.mapping.json.gz"),
            hap2_coord_index: discover("_hap2_to_ref.mapping.json.gz"),
            reads_to_hap1_bam: discover("_reads_to_hap1.bam"),
            reads_to_hap2_bam: discover("_reads_to_hap2.bam"),
            cram_ref: non_empty("cram_ref"),
            hap1_to_ref_bam: discover("_hap1_to_ref.bam"),
            hap2_to_ref_bam: discover("_hap2_to_ref.bam"),
            ref_to_hap1_bam: discover("_ref_to_hap1.bam"),
            ref_to_hap2_bam: discover("_ref_to_hap2.bam"),
            hap1_to_hap2_bam: discover("_hap1_to_hap2.bam"),
            hap2_to_hap1_bam: discover("_hap2_to_hap1.bam"),
        };
        samples.push((sample_id, dp));
    }

    if samples.is_empty() {
        return Err("No data rows found in TSV".to_string());
    }
    Ok(samples)
}

// ---------------------------------------------------------------------------
// Per-panel loaded data
// ---------------------------------------------------------------------------

/// Data loaded for a single panel (reads and/or sequence).
#[derive(Debug, Clone, Default)]
struct PanelData {
    /// Packed pileup rows from BAM query.
    rows: Vec<PileupRow>,
    /// FASTA sequence for this region.
    sequence: Option<FastaSequence>,
    /// Assembly cross-alignment tracks (each entry: label, color, packed rows).
    assembly_tracks: Vec<AssemblyTrack>,
}

/// A single assembly cross-alignment track for display below reads.
#[derive(Debug, Clone)]
struct AssemblyTrack {
    /// Display label (e.g. "Hap1 → Ref").
    label: String,
    /// RGB color for this track.
    color: [u8; 3],
    /// Packed pileup rows from the assembly BAM query.
    rows: Vec<PileupRow>,
}

impl AssemblyTrack {
    fn from_cached(t: &CachedAssemblyTrack) -> Self {
        Self {
            label: t.label.clone(),
            color: t.color,
            rows: t.rows.clone(),
        }
    }
}

impl PanelData {
    fn from_cached(c: &CachedPanelData) -> Self {
        Self {
            rows: c.rows.clone(),
            sequence: c.sequence.clone(),
            assembly_tracks: c
                .assembly_tracks
                .iter()
                .map(AssemblyTrack::from_cached)
                .collect(),
        }
    }
}

// ---------------------------------------------------------------------------
// Dot plot state
// ---------------------------------------------------------------------------

/// State for the dot plot comparison view.
#[derive(Debug, Clone)]
struct DotPlotState {
    /// Which pair of sequences to compare.
    comparison: DotPlotComparison,
    /// Computed result (if sequences are available).
    result: Option<DotPlotResult>,
    /// Whether the dot plot panel is visible.
    show: bool,
    /// K-mer size for dot plot.
    k: usize,
}

impl Default for DotPlotState {
    fn default() -> Self {
        Self {
            comparison: DotPlotComparison::RefVsHap1,
            result: None,
            show: false,
            k: DEFAULT_K,
        }
    }
}

// ---------------------------------------------------------------------------
// Panel identifiers
// ---------------------------------------------------------------------------

/// Panel identifiers for the 3-panel layout.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Panel {
    Reference,
    Haplotype1,
    Haplotype2,
}

impl Panel {
    pub fn label(self) -> &'static str {
        match self {
            Panel::Reference => "Reference",
            Panel::Haplotype1 => "Haplotype 1",
            Panel::Haplotype2 => "Haplotype 2",
        }
    }

    pub fn color(self) -> egui::Color32 {
        match self {
            Panel::Reference => egui::Color32::from_rgb(70, 130, 180),
            Panel::Haplotype1 => egui::Color32::from_rgb(60, 160, 80),
            Panel::Haplotype2 => egui::Color32::from_rgb(180, 100, 60),
        }
    }

    /// Convert to the panel_sync PanelId.
    #[allow(dead_code)]
    pub fn sync_id(self) -> PanelId {
        match self {
            Panel::Reference => PanelId::Reference,
            Panel::Haplotype1 => PanelId::Haplotype1,
            Panel::Haplotype2 => PanelId::Haplotype2,
        }
    }
}

// ---------------------------------------------------------------------------
// Compare Reads result
// ---------------------------------------------------------------------------

/// Statistics from comparing read names across the three panels.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct CompareReadsResult {
    pub ref_only: usize,
    pub hap1_only: usize,
    pub hap2_only: usize,
    pub ref_and_hap1: usize,
    pub ref_and_hap2: usize,
    pub hap1_and_hap2: usize,
    pub all_three: usize,
    pub total_unique: usize,
}

// ---------------------------------------------------------------------------
// SV annotation overlay
// ---------------------------------------------------------------------------

/// An SV gap event annotation for inter-panel display.
#[derive(Debug, Clone)]
pub struct SvAnnotation {
    pub event_type: genome::coordinate_mapper::EventType,
    pub ref_start: u64,
    pub ref_end: u64,
    #[allow(dead_code)]
    pub asm_start: u64,
    #[allow(dead_code)]
    pub asm_end: u64,
    #[allow(dead_code)]
    pub gap_size: Option<u64>,
    pub label: String,
}

impl SvAnnotation {
    /// Return the display color for this SV event type.
    pub fn color(&self) -> [u8; 3] {
        sv_event_color(&self.event_type)
    }
}

/// Map an SV event type to a display color.
pub fn sv_event_color(event_type: &genome::coordinate_mapper::EventType) -> [u8; 3] {
    use genome::coordinate_mapper::EventType;
    match event_type {
        EventType::Deletion => [200, 60, 60],       // red
        EventType::Insertion => [60, 100, 200],     // blue
        EventType::Inversion => [150, 60, 200],     // purple
        EventType::Translocation => [200, 200, 60], // yellow
        EventType::Complex => [140, 140, 140],      // grey
        EventType::Alignment => [100, 100, 100],    // dark grey (shouldn't appear)
    }
}

// ---------------------------------------------------------------------------
// Main application
// ---------------------------------------------------------------------------

/// Main application state.
pub struct ViewerApp {
    pub navigator: RegionNavigator,
    pub status_message: String,
    /// Display configuration (shared across panels).
    pub display_config: PileupDisplayConfig,
    /// Synchronized per-panel view state (pan/zoom).
    pub sync_manager: PanelSyncManager,
    /// Paths to data files.
    pub data_paths: DataPaths,
    /// Per-panel loaded data.
    ref_data: PanelData,
    hap1_data: PanelData,
    hap2_data: PanelData,
    /// Dot plot comparison state.
    dot_plot: DotPlotState,
    /// Legacy single coordinate mapping index (backward compat).
    coord_index: Option<genome::coordinate_mapper::MappingIndex>,
    /// Hap1-to-ref coordinate mapping index.
    hap1_coord_index: Option<genome::coordinate_mapper::MappingIndex>,
    /// Hap2-to-ref coordinate mapping index.
    hap2_coord_index: Option<genome::coordinate_mapper::MappingIndex>,
    /// Background data loader (None in tests, Some when GUI is active).
    loader: Option<DataLoader>,
    /// Whether a background load is currently in flight.
    is_loading: bool,
    /// LRU cache for previously loaded region data.
    region_cache: RegionCache,
    /// Timestamp of last pan/zoom input, used for debounce.
    last_pan_zoom_time: Option<Instant>,
    /// Whether a debounced reload is pending.
    pending_debounce_reload: bool,
    // -- Phase 4 fields --
    /// All available samples from TSV config (sample_id → DataPaths).
    samples: Vec<(String, DataPaths)>,
    /// Index of the currently selected sample in `samples`.
    selected_sample_idx: usize,
    /// Whether to show the compare-reads modal.
    show_compare_reads: bool,
    /// Cached compare-reads result.
    compare_result: Option<CompareReadsResult>,
    /// Whether to show SV overlay markers between panels.
    show_sv_overlay: bool,
    /// SV gap annotations for the current region.
    sv_annotations: Vec<SvAnnotation>,
    /// Whether a PNG screenshot export is pending.
    pending_screenshot: bool,
}

impl Default for ViewerApp {
    fn default() -> Self {
        Self {
            navigator: RegionNavigator::new(),
            status_message:
                "No regions loaded. Use File > Load Regions to open a manifest JSON or VCF."
                    .to_string(),
            display_config: PileupDisplayConfig::default(),
            sync_manager: PanelSyncManager::new(),
            data_paths: DataPaths::default(),
            ref_data: PanelData::default(),
            hap1_data: PanelData::default(),
            hap2_data: PanelData::default(),
            dot_plot: DotPlotState::default(),
            coord_index: None,
            hap1_coord_index: None,
            hap2_coord_index: None,
            loader: None,
            is_loading: false,
            region_cache: RegionCache::new(DEFAULT_CACHE_CAPACITY),
            last_pan_zoom_time: None,
            pending_debounce_reload: false,
            samples: Vec::new(),
            selected_sample_idx: 0,
            show_compare_reads: false,
            compare_result: None,
            show_sv_overlay: false,
            sv_annotations: Vec::new(),
            pending_screenshot: false,
        }
    }
}

/// Parameters for rendering a single panel.
struct PanelRenderParams<'a> {
    panel: Panel,
    region_text: &'a str,
    data: &'a PanelData,
    config: &'a PileupDisplayConfig,
    view_start: u64,
    view_end: u64,
    coord_info: Option<&'a str>,
    /// Optional metadata to show in the header (genotype, size).
    metadata: Option<&'a str>,
}

/// Result of mouse interaction within a panel.
#[derive(Debug, Default)]
struct PanelInteraction {
    /// Pan delta in base pairs (positive = right).
    pan_delta: i64,
    /// Zoom factor (>1 = zoom in) and cursor fraction for cursor-centric zoom.
    zoom: Option<(f64, f64)>,
}

impl ViewerApp {
    pub fn new(
        cc: &eframe::CreationContext<'_>,
        manifest_path: Option<&std::path::Path>,
        data_paths: DataPaths,
        samples: Vec<(String, DataPaths)>,
    ) -> Self {
        configure_fonts(&cc.egui_ctx);
        let mut app = Self {
            data_paths,
            loader: Some(DataLoader::new(cc.egui_ctx.clone())),
            samples,
            ..Self::default()
        };

        // Load coordinate mapping index if provided (legacy single index).
        if let Some(idx_path) = &app.data_paths.coordinate_index {
            match genome::coordinate_mapper::load_index(&idx_path.to_string_lossy()) {
                Ok(index) => {
                    app.coord_index = Some(index);
                }
                Err(e) => {
                    app.status_message = format!("Warning: failed to load coord index: {e}");
                }
            }
        }

        // Load hap1 coordinate mapping index.
        if let Some(idx_path) = &app.data_paths.hap1_coord_index {
            match genome::coordinate_mapper::load_index(&idx_path.to_string_lossy()) {
                Ok(index) => {
                    app.hap1_coord_index = Some(index);
                }
                Err(e) => {
                    eprintln!("Warning: failed to load hap1 coord index: {e}");
                }
            }
        }

        // Load hap2 coordinate mapping index.
        if let Some(idx_path) = &app.data_paths.hap2_coord_index {
            match genome::coordinate_mapper::load_index(&idx_path.to_string_lossy()) {
                Ok(index) => {
                    app.hap2_coord_index = Some(index);
                }
                Err(e) => {
                    eprintln!("Warning: failed to load hap2 coord index: {e}");
                }
            }
        }

        // Load manifest/regions if provided.
        if let Some(path) = manifest_path {
            match app.navigator.load_regions(path) {
                Ok(()) => {
                    app.status_message = format!(
                        "Loaded {} regions from {}",
                        app.navigator.len(),
                        path.display()
                    );
                    // Load data for the first region.
                    app.load_region_data();
                }
                Err(e) => {
                    app.status_message = format!("Error loading manifest: {e}");
                }
            }
        }
        app
    }

    // -- Data loading -------------------------------------------------------

    /// Pack reads using the current display config (position-based or
    /// haplotype-grouped).
    fn pack_reads_for_config(&self, reads: Vec<genome::AlignedRead>) -> Vec<pileup::PileupRow> {
        if self.display_config.sort_by_haplotype {
            pileup::pack_reads_by_haplotype(reads)
        } else {
            pileup::pack_reads(reads)
        }
    }

    /// Load BAM reads and FASTA sequences for the current region.
    ///
    /// When a background loader is available, this first checks the LRU cache
    /// for an instant hit, then dispatches the work to the background thread.
    /// When no loader is present (unit tests), it falls back to synchronous
    /// loading so existing tests continue to pass without a GUI context.
    fn load_region_data(&mut self) {
        let entry = match self.navigator.current() {
            Some(e) => e.clone(),
            None => return,
        };

        // -- Compute haplotype regions via coordinate mapper when missing --
        let hap1_region = entry.hap1_region.clone().or_else(|| {
            let index = self.hap1_coord_index.as_ref()?;
            collapse_asm_regions(&genome::coordinate_mapper::query(
                index,
                &entry.ref_region.chrom,
                entry.ref_region.start,
                entry.ref_region.end,
                0,
            ))
        });
        let hap2_region = entry.hap2_region.clone().or_else(|| {
            let index = self.hap2_coord_index.as_ref()?;
            collapse_asm_regions(&genome::coordinate_mapper::query(
                index,
                &entry.ref_region.chrom,
                entry.ref_region.start,
                entry.ref_region.end,
                0,
            ))
        });

        // Sync panel views to the new region.
        self.sync_manager.set_regions(
            &entry.ref_region,
            hap1_region.as_ref(),
            hap2_region.as_ref(),
        );

        let cache_key = data_loader::make_cache_key(
            &entry.ref_region,
            hap1_region.as_ref(),
            hap2_region.as_ref(),
            self.display_config.sort_by_haplotype,
        );

        // -- Check LRU cache first --
        if let Some(cached) = self.region_cache.get(&cache_key).cloned() {
            self.apply_cached_region(&cached);
            return;
        }

        // -- Dispatch to background loader if available --
        if let Some(loader) = &self.loader {
            let id = loader.next_id();
            let paths = LoadPaths {
                reads_bam: self.data_paths.reads_bam.clone(),
                reference_fasta: self.data_paths.reference_fasta.clone(),
                hap1_fasta: self.data_paths.hap1_fasta.clone(),
                hap2_fasta: self.data_paths.hap2_fasta.clone(),
                reads_to_hap1_bam: self.data_paths.reads_to_hap1_bam.clone(),
                reads_to_hap2_bam: self.data_paths.reads_to_hap2_bam.clone(),
                cram_ref: self.data_paths.cram_ref.clone(),
                hap1_to_ref_bam: self.data_paths.hap1_to_ref_bam.clone(),
                hap2_to_ref_bam: self.data_paths.hap2_to_ref_bam.clone(),
                ref_to_hap1_bam: self.data_paths.ref_to_hap1_bam.clone(),
                ref_to_hap2_bam: self.data_paths.ref_to_hap2_bam.clone(),
                hap1_to_hap2_bam: self.data_paths.hap1_to_hap2_bam.clone(),
                hap2_to_hap1_bam: self.data_paths.hap2_to_hap1_bam.clone(),
            };
            let request = LoadRequest {
                id,
                ref_region: entry.ref_region.clone(),
                hap1_region,
                hap2_region,
                sort_by_haplotype: self.display_config.sort_by_haplotype,
                paths,
                cache_key,
            };
            loader.submit(request);
            self.is_loading = true;
            self.status_message = "Loading…".to_string();
            return;
        }

        // -- Fallback: synchronous load (unit tests without GUI context) --
        self.load_region_data_sync(
            &entry.ref_region,
            hap1_region.as_ref(),
            hap2_region.as_ref(),
        );
    }

    /// Apply a cached region result to the panel data.
    fn apply_cached_region(&mut self, cached: &CachedRegion) {
        self.ref_data = PanelData::from_cached(&cached.ref_data);
        self.hap1_data = PanelData::from_cached(&cached.hap1_data);
        self.hap2_data = PanelData::from_cached(&cached.hap2_data);
        self.status_message = cached.status_message.clone();
        if self.dot_plot.show {
            self.recompute_dot_plot();
        }
        self.sv_annotations = self.collect_sv_annotations();
    }

    /// Poll the background loader for completed results.  Called each frame.
    fn poll_load_results(&mut self) {
        let result = match &self.loader {
            Some(loader) => loader.try_recv(),
            None => None,
        };
        if let Some(result) = result {
            self.is_loading = false;
            // Store in cache before applying
            self.region_cache
                .insert(result.cache_key.clone(), result.data.clone());
            self.apply_cached_region(&result.data);
        }
    }

    /// Synchronous data loading — used only in tests (no GUI context).
    fn load_region_data_sync(
        &mut self,
        ref_region: &crate::region::GenomicRegion,
        hap1_region: Option<&crate::region::GenomicRegion>,
        hap2_region: Option<&crate::region::GenomicRegion>,
    ) {
        // -- Reference panel: BAM/CRAM reads + FASTA sequence --
        self.ref_data = PanelData::default();
        if let Some(reads_path) = &self.data_paths.reads_bam {
            let region_str = ref_region.to_string();
            let cram_ref = self
                .data_paths
                .cram_ref
                .as_deref()
                .or(self.data_paths.reference_fasta.as_deref());
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

                    self.ref_data.rows = self.pack_reads_for_config(reads);

                    self.status_message = format!(
                        "{n_reads} reads ({n_reverse} rev, {n_flagged} secondary, avg MAPQ {})",
                        avg_mapq.map_or("N/A".to_string(), |q| q.to_string()),
                    );
                }
                Err(e) => {
                    self.status_message = format!("Reads error: {e}");
                }
            }
        }
        if let Some(fasta_path) = &self.data_paths.reference_fasta {
            self.ref_data.sequence = load_fasta_for_region(fasta_path, ref_region);
        }
        // Assembly cross-alignment tracks for reference panel
        {
            let region_str = ref_region.to_string();
            self.ref_data.assembly_tracks = load_assembly_tracks(
                &[
                    (
                        "Hap1 → Ref",
                        [60, 160, 80],
                        self.data_paths.hap1_to_ref_bam.as_deref(),
                    ),
                    (
                        "Hap2 → Ref",
                        [180, 100, 60],
                        self.data_paths.hap2_to_ref_bam.as_deref(),
                    ),
                ],
                &region_str,
            );
        }

        // -- Haplotype 1 panel: BAM reads + FASTA sequence --
        self.hap1_data = PanelData::default();
        if let Some(region) = hap1_region {
            if let Some(bam_path) = &self.data_paths.reads_to_hap1_bam {
                let region_str = region.to_string();
                if let Ok(reads) = genome::bam::query_bam(bam_path, &region_str) {
                    self.hap1_data.rows = self.pack_reads_for_config(reads);
                }
            }
            if let Some(fasta_path) = &self.data_paths.hap1_fasta {
                self.hap1_data.sequence = load_fasta_for_region(fasta_path, region);
            }
            // Assembly cross-alignment tracks for hap1 panel
            let region_str = region.to_string();
            self.hap1_data.assembly_tracks = load_assembly_tracks(
                &[
                    (
                        "Ref → Hap1",
                        [70, 130, 180],
                        self.data_paths.ref_to_hap1_bam.as_deref(),
                    ),
                    (
                        "Hap2 → Hap1",
                        [180, 100, 60],
                        self.data_paths.hap2_to_hap1_bam.as_deref(),
                    ),
                ],
                &region_str,
            );
        }

        // -- Haplotype 2 panel: BAM reads + FASTA sequence --
        self.hap2_data = PanelData::default();
        if let Some(region) = hap2_region {
            if let Some(bam_path) = &self.data_paths.reads_to_hap2_bam {
                let region_str = region.to_string();
                if let Ok(reads) = genome::bam::query_bam(bam_path, &region_str) {
                    self.hap2_data.rows = self.pack_reads_for_config(reads);
                }
            }
            if let Some(fasta_path) = &self.data_paths.hap2_fasta {
                self.hap2_data.sequence = load_fasta_for_region(fasta_path, region);
            }
            // Assembly cross-alignment tracks for hap2 panel
            let region_str = region.to_string();
            self.hap2_data.assembly_tracks = load_assembly_tracks(
                &[
                    (
                        "Ref → Hap2",
                        [70, 130, 180],
                        self.data_paths.ref_to_hap2_bam.as_deref(),
                    ),
                    (
                        "Hap1 → Hap2",
                        [60, 160, 80],
                        self.data_paths.hap1_to_hap2_bam.as_deref(),
                    ),
                ],
                &region_str,
            );
        }

        // -- Dot plot: recompute if visible --
        if self.dot_plot.show {
            self.recompute_dot_plot();
        }

        // -- Collect SV annotations for overlay --
        self.sv_annotations = self.collect_sv_annotations();
    }

    /// Recompute the dot plot from loaded FASTA sequences.
    fn recompute_dot_plot(&mut self) {
        let (seq1, seq2) = match self.dot_plot.comparison {
            DotPlotComparison::RefVsHap1 => (
                self.ref_data.sequence.as_ref(),
                self.hap1_data.sequence.as_ref(),
            ),
            DotPlotComparison::RefVsHap2 => (
                self.ref_data.sequence.as_ref(),
                self.hap2_data.sequence.as_ref(),
            ),
            DotPlotComparison::Hap1VsHap2 => (
                self.hap1_data.sequence.as_ref(),
                self.hap2_data.sequence.as_ref(),
            ),
        };

        self.dot_plot.result = match (seq1, seq2) {
            (Some(s1), Some(s2)) => Some(dot_plot::compute_dotplot(
                &s1.sequence,
                &s2.sequence,
                self.dot_plot.k,
            )),
            _ => None,
        };
    }

    /// Compare read names across the three panels and compute set statistics.
    fn compare_reads(&self) -> CompareReadsResult {
        let collect_names = |data: &PanelData| -> HashSet<String> {
            data.rows
                .iter()
                .flat_map(|row| row.reads.iter().map(|r| r.name.clone()))
                .collect()
        };
        let ref_names = collect_names(&self.ref_data);
        let h1_names = collect_names(&self.hap1_data);
        let h2_names = collect_names(&self.hap2_data);

        let all_names: HashSet<&String> = ref_names
            .iter()
            .chain(h1_names.iter())
            .chain(h2_names.iter())
            .collect();

        let mut result = CompareReadsResult::default();
        for name in &all_names {
            let in_ref = ref_names.contains(*name);
            let in_h1 = h1_names.contains(*name);
            let in_h2 = h2_names.contains(*name);
            match (in_ref, in_h1, in_h2) {
                (true, true, true) => result.all_three += 1,
                (true, true, false) => result.ref_and_hap1 += 1,
                (true, false, true) => result.ref_and_hap2 += 1,
                (false, true, true) => result.hap1_and_hap2 += 1,
                (true, false, false) => result.ref_only += 1,
                (false, true, false) => result.hap1_only += 1,
                (false, false, true) => result.hap2_only += 1,
                (false, false, false) => {} // impossible
            }
        }
        result.total_unique = all_names.len();
        result
    }

    /// Collect SV annotations from coordinate mapper indices for the current region.
    fn collect_sv_annotations(&self) -> Vec<SvAnnotation> {
        let entry = match self.navigator.current() {
            Some(e) => e,
            None => return Vec::new(),
        };

        let mut annotations = Vec::new();
        let indices: Vec<(&str, &genome::coordinate_mapper::MappingIndex)> = [
            ("hap1", self.hap1_coord_index.as_ref()),
            ("hap2", self.hap2_coord_index.as_ref()),
            ("", self.coord_index.as_ref()),
        ]
        .into_iter()
        .filter_map(|(label, idx)| idx.map(|i| (label, i)))
        .collect();

        for (label_prefix, index) in indices {
            let results = genome::coordinate_mapper::query(
                index,
                &entry.ref_region.chrom,
                entry.ref_region.start,
                entry.ref_region.end,
                0,
            );
            for r in &results {
                if r.event_type == genome::coordinate_mapper::EventType::Alignment {
                    continue;
                }
                let gap_size = r.ref_gap_size;
                let label = if label_prefix.is_empty() {
                    format!("{:?}", r.event_type)
                } else {
                    format!("{label_prefix}: {:?}", r.event_type)
                };
                annotations.push(SvAnnotation {
                    event_type: r.event_type.clone(),
                    ref_start: r.ref_start,
                    ref_end: r.ref_end,
                    asm_start: r.asm_start,
                    asm_end: r.asm_end,
                    gap_size,
                    label,
                });
            }
        }
        annotations
    }

    /// Paint SV overlay markers between panels.
    fn paint_sv_overlay(&self, ui: &mut egui::Ui) {
        if self.sv_annotations.is_empty() {
            return;
        }
        let panel_width = ui.available_width();
        let ref_view = self.sync_manager.view(PanelId::Reference);
        let view_span = ref_view.view_end.saturating_sub(ref_view.view_start).max(1) as f32;
        let bar_height = 16.0_f32;

        let (response, painter) =
            ui.allocate_painter(egui::vec2(panel_width, bar_height), egui::Sense::hover());
        let origin = response.rect.left_top();

        painter.rect_filled(
            response.rect,
            0.0,
            egui::Color32::from_rgba_premultiplied(30, 30, 40, 200),
        );

        for ann in &self.sv_annotations {
            let x_start = (ann.ref_start.saturating_sub(ref_view.view_start)) as f32 / view_span
                * panel_width;
            let x_end =
                (ann.ref_end.saturating_sub(ref_view.view_start)) as f32 / view_span * panel_width;
            let w = (x_end - x_start).max(3.0);
            let c = ann.color();
            let rect = egui::Rect::from_min_size(
                origin + egui::vec2(x_start, 1.0),
                egui::vec2(w, bar_height - 2.0),
            );
            painter.rect_filled(rect, 2.0, egui::Color32::from_rgb(c[0], c[1], c[2]));
            // Label
            if w > 20.0 {
                painter.text(
                    rect.center(),
                    egui::Align2::CENTER_CENTER,
                    &ann.label,
                    egui::FontId::monospace(8.0),
                    egui::Color32::WHITE,
                );
            }
        }
    }

    // -- Toolbar ------------------------------------------------------------

    /// Render the top toolbar with navigation controls.
    fn show_toolbar(&mut self, ui: &mut egui::Ui) {
        ui.horizontal(|ui| {
            // File menu
            ui.menu_button("File", |ui| {
                if ui.button("Load Regions…").clicked() {
                    if let Some(path) = rfd::FileDialog::new()
                        .add_filter("Regions", &["json", "vcf", "gz"])
                        .pick_file()
                    {
                        match self.navigator.load_regions(&path) {
                            Ok(()) => {
                                self.status_message = format!(
                                    "Loaded {} regions from {}",
                                    self.navigator.len(),
                                    path.display()
                                );
                                self.load_region_data();
                            }
                            Err(e) => {
                                self.status_message = format!("Error: {e}");
                            }
                        }
                    }
                    ui.close_menu();
                }
                if ui.button("Load BAM…").clicked() {
                    if let Some(path) = rfd::FileDialog::new()
                        .add_filter("BAM", &["bam"])
                        .add_filter("CRAM", &["cram"])
                        .pick_file()
                    {
                        self.data_paths.reads_bam = Some(path);
                        self.region_cache.clear();
                        self.load_region_data();
                    }
                    ui.close_menu();
                }
                if ui.button("Load Reference FASTA…").clicked() {
                    if let Some(path) = rfd::FileDialog::new()
                        .add_filter("FASTA", &["fa", "fa.gz", "fasta", "fasta.gz"])
                        .pick_file()
                    {
                        self.data_paths.reference_fasta = Some(path);
                        self.region_cache.clear();
                        self.load_region_data();
                    }
                    ui.close_menu();
                }
                if ui.button("Load Hap1 FASTA…").clicked() {
                    if let Some(path) = rfd::FileDialog::new()
                        .add_filter("FASTA", &["fa", "fa.gz", "fasta", "fasta.gz"])
                        .pick_file()
                    {
                        self.data_paths.hap1_fasta = Some(path);
                        self.region_cache.clear();
                        self.load_region_data();
                    }
                    ui.close_menu();
                }
                if ui.button("Load Hap2 FASTA…").clicked() {
                    if let Some(path) = rfd::FileDialog::new()
                        .add_filter("FASTA", &["fa", "fa.gz", "fasta", "fasta.gz"])
                        .pick_file()
                    {
                        self.data_paths.hap2_fasta = Some(path);
                        self.region_cache.clear();
                        self.load_region_data();
                    }
                    ui.close_menu();
                }
                ui.separator();
                if ui.button("Export PNG…").clicked() {
                    self.pending_screenshot = true;
                    ui.close_menu();
                }
                ui.separator();
                if ui.button("Quit").clicked() {
                    ui.ctx().send_viewport_cmd(egui::ViewportCommand::Close);
                }
            });

            ui.separator();

            // Sample selector (multi-sample support)
            if self.samples.len() > 1 {
                let current_sample = self
                    .samples
                    .get(self.selected_sample_idx)
                    .map(|(id, _)| id.clone())
                    .unwrap_or_default();
                let prev_idx = self.selected_sample_idx;
                egui::ComboBox::from_id_salt("sample_selector")
                    .selected_text(format!("Sample: {current_sample}"))
                    .width(150.0)
                    .show_ui(ui, |ui| {
                        for (i, (id, _)) in self.samples.iter().enumerate() {
                            if ui
                                .selectable_label(i == self.selected_sample_idx, id)
                                .clicked()
                            {
                                self.selected_sample_idx = i;
                            }
                        }
                    });
                if self.selected_sample_idx != prev_idx {
                    // Switch sample: update data_paths, reload coord indices, clear cache
                    if let Some((_, dp)) = self.samples.get(self.selected_sample_idx) {
                        self.data_paths = dp.clone();
                        self.coord_index = None;
                        self.hap1_coord_index =
                            try_load_coord_index(self.data_paths.hap1_coord_index.as_deref());
                        self.hap2_coord_index =
                            try_load_coord_index(self.data_paths.hap2_coord_index.as_deref());
                        self.region_cache.clear();
                        self.load_region_data();
                    }
                }
                ui.separator();
            }

            // Navigation controls
            let has_regions = !self.navigator.is_empty();

            if ui
                .add_enabled(has_regions, egui::Button::new("◀ Prev"))
                .on_hover_text("Previous region (← or [)")
                .clicked()
            {
                self.navigator.prev();
                self.load_region_data();
                self.update_status_for_current_region();
            }

            // Region selector dropdown
            if has_regions {
                let current_label = self
                    .navigator
                    .current()
                    .map(|r| r.label.clone())
                    .unwrap_or_default();
                let labels = self.navigator.labels();
                egui::ComboBox::from_id_salt("region_selector")
                    .selected_text(&current_label)
                    .width(250.0)
                    .show_ui(ui, |ui| {
                        for (i, label) in labels.iter().enumerate() {
                            if ui
                                .selectable_label(i == self.navigator.current_index(), label)
                                .clicked()
                            {
                                self.navigator.go_to(i);
                                self.load_region_data();
                                self.update_status_for_current_region();
                            }
                        }
                    });
            } else {
                ui.label("No regions");
            }

            if ui
                .add_enabled(has_regions, egui::Button::new("Next ▶"))
                .on_hover_text("Next region (→ or ])")
                .clicked()
            {
                self.navigator.next();
                self.load_region_data();
                self.update_status_for_current_region();
            }

            ui.label(self.navigator.counter_text());

            ui.separator();

            // Display toggle: squished / expanded
            let squish_label = if self.display_config.squished {
                "Squished"
            } else {
                "Expanded"
            };
            if ui
                .button(squish_label)
                .on_hover_text("Toggle squished/expanded read display")
                .clicked()
            {
                self.display_config.squished = !self.display_config.squished;
            }

            // Display toggle: hide small indels
            let indel_label = if self.display_config.hide_small_indels {
                format!("Indels ≤{}bp: hidden", self.display_config.indel_threshold)
            } else {
                "Indels: shown".to_string()
            };
            if ui
                .button(&indel_label)
                .on_hover_text("Toggle small indel display")
                .clicked()
            {
                self.display_config.hide_small_indels = !self.display_config.hide_small_indels;
            }

            // Display toggle: show mismatches
            let mm_label = if self.display_config.show_mismatches {
                "SNVs: on"
            } else {
                "SNVs: off"
            };
            if ui
                .button(mm_label)
                .on_hover_text("Toggle base-level mismatch/SNV display")
                .clicked()
            {
                self.display_config.show_mismatches = !self.display_config.show_mismatches;
            }

            // Display toggle: show soft clips
            let sc_label = if self.display_config.show_soft_clips {
                "Clips: on"
            } else {
                "Clips: off"
            };
            if ui
                .button(sc_label)
                .on_hover_text("Toggle soft-clip overlay display")
                .clicked()
            {
                self.display_config.show_soft_clips = !self.display_config.show_soft_clips;
            }

            // Display toggle: sort by haplotype
            let hp_label = if self.display_config.sort_by_haplotype {
                "HP Sort: on"
            } else {
                "HP Sort: off"
            };
            if ui
                .button(hp_label)
                .on_hover_text("Toggle haplotype-based read grouping")
                .clicked()
            {
                self.display_config.sort_by_haplotype = !self.display_config.sort_by_haplotype;
                // Re-pack reads with the new sort mode
                self.load_region_data();
            }

            // Display toggle: show coverage histogram
            let cov_label = if self.display_config.show_coverage {
                "Coverage: on"
            } else {
                "Coverage: off"
            };
            if ui
                .button(cov_label)
                .on_hover_text("Toggle coverage depth histogram")
                .clicked()
            {
                self.display_config.show_coverage = !self.display_config.show_coverage;
            }

            // Display toggle: colorblind palette
            let cb_label = if self.display_config.use_colorblind_palette {
                "🎨 Colorblind"
            } else {
                "🎨 Default"
            };
            if ui
                .button(cb_label)
                .on_hover_text("Toggle colorblind-safe palette")
                .clicked()
            {
                self.display_config.use_colorblind_palette =
                    !self.display_config.use_colorblind_palette;
            }

            // Display toggle: SV overlay
            let sv_label = if self.show_sv_overlay {
                "SV Overlay: on"
            } else {
                "SV Overlay: off"
            };
            if ui
                .button(sv_label)
                .on_hover_text("Toggle SV gap event overlay between panels")
                .clicked()
            {
                self.show_sv_overlay = !self.show_sv_overlay;
            }

            ui.separator();

            // Sync toggle
            let sync_label = if self.sync_manager.sync_enabled {
                "🔗 Sync"
            } else {
                "🔗̸ Unsync"
            };
            if ui
                .button(sync_label)
                .on_hover_text("Toggle cross-panel synchronization (S)")
                .clicked()
            {
                self.sync_manager.sync_enabled = !self.sync_manager.sync_enabled;
            }

            // Zoom controls
            if ui.button("🔍+").on_hover_text("Zoom in (+)").clicked() {
                self.sync_manager.zoom(PanelId::Reference, 2.0);
                self.schedule_debounced_reload();
            }
            if ui.button("🔍−").on_hover_text("Zoom out (-)").clicked() {
                self.sync_manager.zoom(PanelId::Reference, 0.5);
                self.schedule_debounced_reload();
            }

            // Pan controls
            let pan_delta = self.sync_manager.view(PanelId::Reference).span() as i64 / 4;
            if ui.button("◀◀").on_hover_text("Pan left").clicked() {
                self.sync_manager.pan(PanelId::Reference, -pan_delta);
                self.schedule_debounced_reload();
            }
            if ui.button("▶▶").on_hover_text("Pan right").clicked() {
                self.sync_manager.pan(PanelId::Reference, pan_delta);
                self.schedule_debounced_reload();
            }

            ui.separator();

            // Dot plot toggle
            let dp_label = if self.dot_plot.show {
                "Dot Plot ▼"
            } else {
                "Dot Plot ▶"
            };
            if ui
                .button(dp_label)
                .on_hover_text("Toggle dot plot panel")
                .clicked()
            {
                self.dot_plot.show = !self.dot_plot.show;
                if self.dot_plot.show {
                    self.recompute_dot_plot();
                }
            }

            // Compare Reads button
            if ui
                .button("Compare Reads")
                .on_hover_text("Compare read names across all three panels")
                .clicked()
            {
                self.compare_result = Some(self.compare_reads());
                self.show_compare_reads = true;
            }
        });
    }

    // -- Panel rendering ----------------------------------------------------

    /// Render a single panel with pileup reads drawn on a canvas.
    /// Returns mouse interaction events (pan/zoom) for the app to process.
    fn show_panel(ui: &mut egui::Ui, params: &PanelRenderParams<'_>) -> PanelInteraction {
        let mut interaction = PanelInteraction::default();
        let header_color = params.panel.color();

        // Panel header
        let header_rect = ui.allocate_space(egui::vec2(ui.available_width(), 24.0)).1;
        ui.painter().rect_filled(header_rect, 0.0, header_color);

        let mut header_text = format!("{}  |  {}", params.panel.label(), params.region_text);
        if let Some(meta) = params.metadata
            && !meta.is_empty()
        {
            header_text.push_str(&format!("  {meta}"));
        }
        if let Some(info) = params.coord_info {
            header_text.push_str(&format!("  [{info}]"));
        }
        ui.painter().text(
            header_rect.left_center() + egui::vec2(8.0, 0.0),
            egui::Align2::LEFT_CENTER,
            header_text,
            egui::FontId::proportional(13.0),
            egui::Color32::WHITE,
        );

        // Coordinate ruler
        let panel_width = ui.available_width();
        Self::paint_ruler(ui, params.view_start, params.view_end, panel_width);

        // Panel body: dark background with pileup rendering
        let body = egui::Frame::new()
            .fill(egui::Color32::from_gray(32))
            .inner_margin(egui::Margin::same(4));

        body.show(ui, |ui| {
            ui.set_min_height(80.0);

            let has_reads = !params.data.rows.is_empty();
            let has_assembly = !params.data.assembly_tracks.is_empty();
            let has_seq = params.data.sequence.is_some();

            if !has_reads && !has_assembly && !has_seq {
                // Placeholder when no data is loaded
                ui.centered_and_justified(|ui| {
                    ui.label(
                        egui::RichText::new(format!("{} – no data loaded", params.panel.label()))
                            .color(egui::Color32::from_gray(120))
                            .italics(),
                    );
                });
            } else {
                // Render base-level sequence at high zoom
                if let Some(seq) = &params.data.sequence {
                    let view_span = params.view_end.saturating_sub(params.view_start);
                    let bp_per_px = if view_span > 0 {
                        ui.available_width() / view_span as f32
                    } else {
                        0.0
                    };

                    if bp_per_px >= 5.0 {
                        // High zoom: render individual nucleotide letters
                        Self::paint_base_sequence(ui, seq, params.view_start, params.view_end);
                    } else {
                        // Low zoom: show text summary
                        let seq_info = format!(
                            "{}: {} ({} bp)",
                            seq.name,
                            if seq.sequence.len() > 40 {
                                format!("{}…", &seq.sequence[..40])
                            } else {
                                seq.sequence.clone()
                            },
                            seq.sequence.len()
                        );
                        ui.label(
                            egui::RichText::new(seq_info)
                                .small()
                                .color(egui::Color32::from_gray(170))
                                .monospace(),
                        );
                    }
                }

                // Coverage histogram
                if params.config.show_coverage && has_reads {
                    let cov_panel_width = ui.available_width();
                    let num_bins = (cov_panel_width as usize).clamp(10, 500);
                    let bins = pileup::compute_coverage(
                        &params.data.rows,
                        params.view_start,
                        params.view_end,
                        num_bins,
                    );
                    Self::paint_coverage(ui, &bins, cov_panel_width);
                }

                // Render pileup reads
                if has_reads {
                    let pw = ui.available_width();
                    let rects = pileup::layout_read_rects(
                        &params.data.rows,
                        params.config,
                        params.view_start,
                        params.view_end,
                        pw,
                    );
                    Self::paint_pileup_with_tooltips(ui, &rects);
                }

                // Render assembly cross-alignment tracks below reads
                for track in &params.data.assembly_tracks {
                    if track.rows.is_empty() {
                        continue;
                    }
                    // Small label for the assembly track
                    ui.label(egui::RichText::new(&track.label).small().color(
                        egui::Color32::from_rgb(track.color[0], track.color[1], track.color[2]),
                    ));
                    // Assembly tracks always use squished display
                    let asm_config = params.config.with_squished(true);
                    let asm_pw = ui.available_width();
                    let rects = pileup::layout_read_rects(
                        &track.rows,
                        &asm_config,
                        params.view_start,
                        params.view_end,
                        asm_pw,
                    );
                    // Override colors with the track's distinct color
                    let colored_rects: Vec<ReadRect> = rects
                        .into_iter()
                        .map(|mut r| {
                            r.color = track.color;
                            r
                        })
                        .collect();
                    Self::paint_pileup_with_tooltips(ui, &colored_rects);
                }
            }

            // Mouse interaction: capture drag (pan) and scroll (zoom) on the panel body
            let panel_rect = ui.min_rect();
            let panel_resp = ui.interact(
                panel_rect,
                ui.id().with("mouse"),
                egui::Sense::click_and_drag(),
            );

            // Drag → pan
            if panel_resp.dragged() {
                let drag_dx = panel_resp.drag_delta().x;
                let body_width = panel_rect.width();
                if body_width > 0.0 {
                    let view_span = params.view_end.saturating_sub(params.view_start);
                    let bp_per_px = view_span as f64 / body_width as f64;
                    // Drag right → view moves left (pan_delta < 0)
                    interaction.pan_delta = -(drag_dx as f64 * bp_per_px) as i64;
                }
            }

            // Scroll → zoom (cursor-centric)
            let scroll_delta = ui.input(|i| i.smooth_scroll_delta.y);
            if scroll_delta.abs() > 0.1
                && let Some(hover_pos) = ui.input(|i| i.pointer.hover_pos())
            {
                let frac =
                    ((hover_pos.x - panel_rect.left()) / panel_rect.width()).clamp(0.0, 1.0) as f64;
                // Scroll up → zoom in, scroll down → zoom out
                let factor = if scroll_delta > 0.0 { 1.2 } else { 1.0 / 1.2 };
                interaction.zoom = Some((factor, frac));
            }
        });
        interaction
    }

    /// Paint pre-computed read rectangles onto the UI with hover tooltips.
    fn paint_pileup_with_tooltips(ui: &mut egui::Ui, rects: &[ReadRect]) {
        if rects.is_empty() {
            return;
        }

        // Compute the total height needed
        let max_y = rects.iter().map(|r| r.y + r.height).fold(0.0_f32, f32::max);
        let total_height = max_y + 4.0; // small padding

        let (response, painter) = ui.allocate_painter(
            egui::vec2(ui.available_width(), total_height),
            egui::Sense::hover(),
        );

        let origin = response.rect.left_top();

        for rect in rects {
            let min = origin + egui::vec2(rect.x, rect.y);
            let max = min + egui::vec2(rect.width, rect.height);
            let color = egui::Color32::from_rgb(rect.color[0], rect.color[1], rect.color[2]);
            painter.rect_filled(egui::Rect::from_min_max(min, max), 0.0, color);
        }

        // Hover tooltip: find the read under cursor
        if let Some(hover_pos) = response.hover_pos() {
            for rect in rects {
                if let Some(ref tip) = rect.tooltip {
                    let min = origin + egui::vec2(rect.x, rect.y);
                    let max = min + egui::vec2(rect.width, rect.height);
                    let r = egui::Rect::from_min_max(min, max);
                    if r.contains(hover_pos) {
                        egui::show_tooltip_at_pointer(
                            ui.ctx(),
                            ui.layer_id(),
                            ui.id().with("read_tooltip"),
                            |ui| {
                                ui.label(
                                    egui::RichText::new(tip.tooltip_text())
                                        .monospace()
                                        .size(11.0),
                                );
                            },
                        );
                        break; // show only one tooltip
                    }
                }
            }
        }
    }

    /// Paint the genomic coordinate ruler.
    fn paint_ruler(ui: &mut egui::Ui, view_start: u64, view_end: u64, panel_width: f32) {
        let ticks = ruler::compute_ticks(view_start, view_end, panel_width);
        if ticks.is_empty() {
            return;
        }

        let (response, painter) = ui.allocate_painter(
            egui::vec2(panel_width, ruler::RULER_HEIGHT),
            egui::Sense::hover(),
        );
        let origin = response.rect.left_top();
        let bottom = response.rect.bottom();

        // Background
        painter.rect_filled(response.rect, 0.0, egui::Color32::from_gray(40));

        for tick in &ticks {
            let x = origin.x + tick.x;
            let tick_top = if tick.is_major {
                bottom - 10.0
            } else {
                bottom - 5.0
            };
            let tick_color = if tick.is_major {
                egui::Color32::from_gray(200)
            } else {
                egui::Color32::from_gray(100)
            };
            painter.line_segment(
                [egui::pos2(x, tick_top), egui::pos2(x, bottom)],
                egui::Stroke::new(1.0, tick_color),
            );
            if let Some(ref label) = tick.label {
                painter.text(
                    egui::pos2(x + 2.0, origin.y + 2.0),
                    egui::Align2::LEFT_TOP,
                    label,
                    egui::FontId::monospace(9.0),
                    egui::Color32::from_gray(200),
                );
            }
        }
    }

    /// Paint the coverage depth histogram.
    fn paint_coverage(ui: &mut egui::Ui, bins: &[pileup::CoverageBin], panel_width: f32) {
        if bins.is_empty() {
            return;
        }

        let track_height = pileup::COVERAGE_TRACK_HEIGHT;
        let max_depth = bins.iter().map(|b| b.depth).max().unwrap_or(1).max(1);
        let bin_width = panel_width / bins.len() as f32;

        let (response, painter) =
            ui.allocate_painter(egui::vec2(panel_width, track_height), egui::Sense::hover());
        let origin = response.rect.left_top();
        let bottom = response.rect.bottom();

        // Background
        painter.rect_filled(
            response.rect,
            0.0,
            egui::Color32::from_rgba_premultiplied(20, 20, 30, 200),
        );

        for (i, bin) in bins.iter().enumerate() {
            let x = origin.x + i as f32 * bin_width;
            let total_h = (bin.depth as f32 / max_depth as f32) * track_height;

            // HP-stratified: hap1 (green) at bottom, hap2 (orange) on top, grey for rest
            let h1_h = (bin.hap1_depth as f32 / max_depth as f32) * track_height;
            let h2_h = (bin.hap2_depth as f32 / max_depth as f32) * track_height;
            let other_h = total_h - h1_h - h2_h;

            // Draw HP1 (green) at bottom
            if h1_h > 0.1 {
                let r = egui::Rect::from_min_max(
                    egui::pos2(x, bottom - h1_h),
                    egui::pos2(x + bin_width, bottom),
                );
                painter.rect_filled(
                    r,
                    0.0,
                    egui::Color32::from_rgba_premultiplied(80, 180, 80, 180),
                );
            }
            // Draw HP2 (orange) above HP1
            if h2_h > 0.1 {
                let r = egui::Rect::from_min_max(
                    egui::pos2(x, bottom - h1_h - h2_h),
                    egui::pos2(x + bin_width, bottom - h1_h),
                );
                painter.rect_filled(
                    r,
                    0.0,
                    egui::Color32::from_rgba_premultiplied(200, 140, 50, 180),
                );
            }
            // Draw unphased (grey) on top
            if other_h > 0.1 {
                let r = egui::Rect::from_min_max(
                    egui::pos2(x, bottom - total_h),
                    egui::pos2(x + bin_width, bottom - h1_h - h2_h),
                );
                painter.rect_filled(
                    r,
                    0.0,
                    egui::Color32::from_rgba_premultiplied(140, 140, 140, 160),
                );
            }
        }

        // Depth label
        painter.text(
            origin + egui::vec2(2.0, 1.0),
            egui::Align2::LEFT_TOP,
            format!("max: {max_depth}×"),
            egui::FontId::monospace(9.0),
            egui::Color32::from_gray(200),
        );
    }

    /// Paint per-base nucleotide letters for a reference sequence at high zoom.
    fn paint_base_sequence(ui: &mut egui::Ui, seq: &FastaSequence, view_start: u64, view_end: u64) {
        let panel_width = ui.available_width();
        let view_span = view_end.saturating_sub(view_start);
        if view_span == 0 {
            return;
        }
        let bp_per_px = panel_width / view_span as f32;
        let base_height = 14.0_f32;

        let (response, painter) = ui.allocate_painter(
            egui::vec2(panel_width, base_height + 2.0),
            egui::Sense::hover(),
        );
        let origin = response.rect.left_top();

        // Background
        painter.rect_filled(response.rect, 0.0, egui::Color32::from_gray(24));

        // Determine which bases from the sequence fall in the view
        let seq_start = seq.start;
        let seq_bytes = seq.sequence.as_bytes();

        for pos in view_start..view_end {
            if pos < seq_start {
                continue;
            }
            let idx = (pos - seq_start) as usize;
            if idx >= seq_bytes.len() {
                break;
            }
            let base = seq_bytes[idx];
            let x = (pos - view_start) as f32 * bp_per_px;

            let color = pileup::nucleotide_color(base);
            let base_char = (base as char).to_ascii_uppercase();

            if bp_per_px >= 8.0 {
                // Draw letter
                painter.text(
                    origin + egui::vec2(x + bp_per_px * 0.5, base_height * 0.5 + 1.0),
                    egui::Align2::CENTER_CENTER,
                    base_char.to_string(),
                    egui::FontId::monospace(11.0),
                    egui::Color32::from_rgb(color[0], color[1], color[2]),
                );
            } else {
                // Draw colored bar
                let r = egui::Rect::from_min_size(
                    origin + egui::vec2(x, 1.0),
                    egui::vec2(bp_per_px.max(1.0), base_height),
                );
                painter.rect_filled(
                    r,
                    0.0,
                    egui::Color32::from_rgb(color[0], color[1], color[2]),
                );
            }
        }
    }

    /// Render the dot plot panel.
    fn show_dot_plot(&mut self, ui: &mut egui::Ui) {
        ui.horizontal(|ui| {
            ui.label("Comparison:");
            for cmp in DotPlotComparison::ALL {
                if ui
                    .selectable_label(self.dot_plot.comparison == cmp, cmp.label())
                    .clicked()
                {
                    self.dot_plot.comparison = cmp;
                    self.recompute_dot_plot();
                }
            }
        });

        match &self.dot_plot.result {
            Some(result) => {
                let total = result.total_matches();
                ui.label(format!(
                    "{} matches (fwd: {}, rev: {}, pal: {}){}",
                    total,
                    result.forward.len(),
                    result.reverse.len(),
                    result.palindrome.len(),
                    if result.truncated { " [TRUNCATED]" } else { "" }
                ));

                // Draw dot plot canvas
                if total > 0 {
                    let size = 200.0_f32;
                    let (response, painter) =
                        ui.allocate_painter(egui::vec2(size, size), egui::Sense::hover());
                    let rect = response.rect;

                    // Background
                    painter.rect_filled(rect, 0.0, egui::Color32::from_gray(20));

                    let x_scale = size / result.seq1_len.max(1) as f32;
                    let y_scale = size / result.seq2_len.max(1) as f32;

                    // Draw forward matches (blue)
                    for m in &result.forward {
                        let p = rect.left_top()
                            + egui::vec2(m.x as f32 * x_scale, m.y as f32 * y_scale);
                        painter.rect_filled(
                            egui::Rect::from_min_size(p, egui::vec2(1.0, 1.0)),
                            0.0,
                            egui::Color32::from_rgb(80, 140, 255),
                        );
                    }

                    // Draw reverse matches (red)
                    for m in &result.reverse {
                        let p = rect.left_top()
                            + egui::vec2(m.x as f32 * x_scale, m.y as f32 * y_scale);
                        painter.rect_filled(
                            egui::Rect::from_min_size(p, egui::vec2(1.0, 1.0)),
                            0.0,
                            egui::Color32::from_rgb(255, 80, 80),
                        );
                    }

                    // Draw palindromic matches (green)
                    for m in &result.palindrome {
                        let p = rect.left_top()
                            + egui::vec2(m.x as f32 * x_scale, m.y as f32 * y_scale);
                        painter.rect_filled(
                            egui::Rect::from_min_size(p, egui::vec2(1.0, 1.0)),
                            0.0,
                            egui::Color32::from_rgb(80, 200, 80),
                        );
                    }
                }
            }
            None => {
                ui.label(
                    egui::RichText::new("Load FASTA files to compute dot plot")
                        .color(egui::Color32::from_gray(120))
                        .italics(),
                );
            }
        }
    }

    // -- Status / navigation helpers ----------------------------------------

    /// Update status message to reflect the current region.
    fn update_status_for_current_region(&mut self) {
        if let Some(entry) = self.navigator.current() {
            let mut msg = format!(
                "Region {}: {} → ref {}",
                self.navigator.counter_text(),
                entry.label,
                entry.ref_region
            );

            // Append variant metadata.
            if !entry.genotype.is_empty() {
                msg.push_str(&format!(" | GT: {}", entry.genotype));
            }
            if entry.sv_size > 0 {
                msg.push_str(&format!(" | size: {} bp", entry.sv_size));
            }

            // Show coordinate mapper info if available (legacy single index
            // or hap-specific indices).
            let indices: Vec<(&str, &genome::coordinate_mapper::MappingIndex)> = [
                ("", self.coord_index.as_ref()),
                ("hap1", self.hap1_coord_index.as_ref()),
                ("hap2", self.hap2_coord_index.as_ref()),
            ]
            .into_iter()
            .filter_map(|(label, idx)| idx.map(|i| (label, i)))
            .collect();

            for (label, index) in &indices {
                let results = genome::coordinate_mapper::query(
                    index,
                    &entry.ref_region.chrom,
                    entry.ref_region.start,
                    entry.ref_region.end,
                    0,
                );
                if !results.is_empty() {
                    let prefix = if label.is_empty() {
                        "mapped".to_string()
                    } else {
                        label.to_string()
                    };
                    let events: Vec<String> = results
                        .iter()
                        .filter(|r| r.event_type == genome::coordinate_mapper::EventType::Alignment)
                        .map(|r| format!("{}:{}-{}", r.asm_chrom, r.asm_start, r.asm_end))
                        .collect();
                    if !events.is_empty() {
                        msg.push_str(&format!(" | {prefix}: {}", events.join(", ")));
                    }
                }
            }

            // Append sync/view state info.
            if self.sync_manager.all_spans_match() {
                msg.push_str(" | spans synced");
            }
            if let Some(hr) = self.sync_manager.highlight_region(PanelId::Reference) {
                msg.push_str(&format!(" | view: {hr}"));
            }

            self.status_message = msg;
        }
    }

    /// Schedule a debounced data reload after pan/zoom.
    ///
    /// Each call resets the debounce timer.  The actual reload fires in
    /// `update()` once `PAN_ZOOM_DEBOUNCE_MS` has elapsed without another
    /// pan/zoom event.
    fn schedule_debounced_reload(&mut self) {
        self.last_pan_zoom_time = Some(Instant::now());
        self.pending_debounce_reload = true;
        // Cancel any in-flight background load for stale view.
        if let Some(loader) = &self.loader {
            loader.cancel();
        }
    }

    /// Handle keyboard shortcuts for region navigation and zoom.
    fn handle_keyboard(&mut self, ctx: &egui::Context) {
        if ctx.wants_keyboard_input() {
            return; // Don't capture keys when a text field is focused
        }

        let prev = ctx.input(|i| {
            i.key_pressed(egui::Key::ArrowLeft) || i.key_pressed(egui::Key::OpenBracket)
        });
        let next = ctx.input(|i| {
            i.key_pressed(egui::Key::ArrowRight) || i.key_pressed(egui::Key::CloseBracket)
        });
        let zoom_in =
            ctx.input(|i| i.key_pressed(egui::Key::Plus) || i.key_pressed(egui::Key::Equals));
        let zoom_out = ctx.input(|i| i.key_pressed(egui::Key::Minus));
        let toggle_sync = ctx.input(|i| i.key_pressed(egui::Key::S));

        if prev && !self.navigator.is_empty() {
            self.navigator.prev();
            self.load_region_data();
            self.update_status_for_current_region();
        }
        if next && !self.navigator.is_empty() {
            self.navigator.next();
            self.load_region_data();
            self.update_status_for_current_region();
        }
        if zoom_in {
            self.sync_manager.zoom(PanelId::Reference, 2.0);
            self.schedule_debounced_reload();
        }
        if zoom_out {
            self.sync_manager.zoom(PanelId::Reference, 0.5);
            self.schedule_debounced_reload();
        }
        if toggle_sync {
            self.sync_manager.sync_enabled = !self.sync_manager.sync_enabled;
        }
    }
}

/// Debounce interval for pan/zoom before triggering data reload (ms).
const PAN_ZOOM_DEBOUNCE_MS: u128 = 200;

impl eframe::App for ViewerApp {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        // Poll for background load results.
        self.poll_load_results();

        // Handle debounced pan/zoom reload.
        if self.pending_debounce_reload
            && self
                .last_pan_zoom_time
                .is_some_and(|t| t.elapsed().as_millis() >= PAN_ZOOM_DEBOUNCE_MS)
        {
            self.pending_debounce_reload = false;
            self.load_region_data();
        }

        // Request continuous repaints while loading or debouncing.
        if self.is_loading || self.pending_debounce_reload {
            ctx.request_repaint();
        }

        self.handle_keyboard(ctx);

        // Handle pending screenshot export.
        if self.pending_screenshot {
            self.pending_screenshot = false;
            ctx.send_viewport_cmd(egui::ViewportCommand::Screenshot(egui::UserData::default()));
        }
        // Check for screenshot result in events.
        let screenshot_image = ctx.input(|i| {
            for event in &i.raw.events {
                if let egui::Event::Screenshot { image, .. } = event {
                    return Some(image.clone());
                }
            }
            None
        });
        if let Some(image) = screenshot_image
            && let Some(path) = rfd::FileDialog::new()
                .add_filter("PPM Image", &["ppm"])
                .set_file_name("screenshot.ppm")
                .save_file()
        {
            let w = image.width() as u32;
            let h = image.height() as u32;
            let mut data = Vec::new();
            data.extend_from_slice(format!("P6\n{w} {h}\n255\n").as_bytes());
            for pixel in image.pixels.iter() {
                data.push(pixel.r());
                data.push(pixel.g());
                data.push(pixel.b());
            }
            if let Err(e) = std::fs::write(&path, &data) {
                self.status_message = format!("Export failed: {e}");
            } else {
                self.status_message = format!("Exported screenshot to {}", path.display());
            }
        }

        // Compare Reads modal window
        if self.show_compare_reads {
            let mut open = self.show_compare_reads;
            egui::Window::new("Compare Reads")
                .open(&mut open)
                .resizable(false)
                .show(ctx, |ui| {
                    if let Some(ref result) = self.compare_result {
                        ui.label(format!("Total unique reads: {}", result.total_unique));
                        ui.separator();
                        ui.label(format!("Ref only: {}", result.ref_only));
                        ui.label(format!("Hap1 only: {}", result.hap1_only));
                        ui.label(format!("Hap2 only: {}", result.hap2_only));
                        ui.separator();
                        ui.label(format!("Ref ∩ Hap1: {}", result.ref_and_hap1));
                        ui.label(format!("Ref ∩ Hap2: {}", result.ref_and_hap2));
                        ui.label(format!("Hap1 ∩ Hap2: {}", result.hap1_and_hap2));
                        ui.separator();
                        ui.label(format!("All three: {}", result.all_three));
                    } else {
                        ui.label("No data");
                    }
                });
            self.show_compare_reads = open;
        }

        // Top toolbar
        egui::TopBottomPanel::top("toolbar").show(ctx, |ui| {
            self.show_toolbar(ui);
        });

        // Bottom status bar
        egui::TopBottomPanel::bottom("status_bar").show(ctx, |ui| {
            ui.horizontal(|ui| {
                if self.is_loading {
                    ui.spinner();
                    ui.label(
                        egui::RichText::new("Loading…")
                            .small()
                            .color(egui::Color32::from_gray(180)),
                    );
                } else {
                    ui.label(
                        egui::RichText::new(&self.status_message)
                            .small()
                            .color(egui::Color32::from_gray(180)),
                    );
                }
            });
        });

        // Clone display_config for immutable borrow inside closure
        let config = self.display_config.clone();

        // Build metadata string for Reference panel header
        let metadata_text: Option<String> = self.navigator.current().and_then(|entry| {
            let mut parts = Vec::new();
            if !entry.genotype.is_empty() {
                parts.push(format!("GT:{}", entry.genotype));
            }
            if entry.sv_size > 0 {
                parts.push(format!("{} bp", entry.sv_size));
            }
            if parts.is_empty() {
                None
            } else {
                Some(parts.join(" | "))
            }
        });

        // Central area: 3 panels + optional SV overlay + optional dot plot
        let mut interactions: Vec<(PanelId, PanelInteraction)> = Vec::new();
        egui::CentralPanel::default().show(ctx, |ui| {
            let current = self.navigator.current();
            let ref_text = current
                .map(|e| e.ref_region.to_string())
                .unwrap_or_else(|| "—".to_string());

            let h1_view = self.sync_manager.view(PanelId::Haplotype1);
            let h2_view = self.sync_manager.view(PanelId::Haplotype2);
            let hap1_text = if h1_view.view_start > 0 || h1_view.view_end > 0 {
                self.sync_manager
                    .highlight_region(PanelId::Haplotype1)
                    .map(|r| r.to_string())
                    .unwrap_or_else(|| "—".to_string())
            } else {
                "—".to_string()
            };
            let hap2_text = if h2_view.view_start > 0 || h2_view.view_end > 0 {
                self.sync_manager
                    .highlight_region(PanelId::Haplotype2)
                    .map(|r| r.to_string())
                    .unwrap_or_else(|| "—".to_string())
            } else {
                "—".to_string()
            };

            let ref_view = self.sync_manager.view(PanelId::Reference);

            let (ref_start, ref_end) = (ref_view.view_start, ref_view.view_end);
            let (h1_start, h1_end) = (h1_view.view_start, h1_view.view_end);
            let (h2_start, h2_end) = (h2_view.view_start, h2_view.view_end);

            let ref_coord_info = self.coord_index.as_ref().and_then(|index| {
                let entry = self.navigator.current()?;
                let results = genome::coordinate_mapper::query(
                    index,
                    &entry.ref_region.chrom,
                    entry.ref_region.start,
                    entry.ref_region.end,
                    0,
                );
                if results.is_empty() {
                    None
                } else {
                    Some(format!("{} blocks", results.len()))
                }
            });

            let hap1_coord_info: Option<String> = if !self.hap1_data.rows.is_empty() {
                let n: usize = self.hap1_data.rows.iter().map(|r| r.reads.len()).sum();
                Some(format!("{n} reads"))
            } else {
                None
            };
            let hap2_coord_info: Option<String> = if !self.hap2_data.rows.is_empty() {
                let n: usize = self.hap2_data.rows.iter().map(|r| r.reads.len()).sum();
                Some(format!("{n} reads"))
            } else {
                None
            };

            let sv_overlay_height = if self.show_sv_overlay && !self.sv_annotations.is_empty() {
                20.0
            } else {
                0.0
            };
            let dot_plot_height = if self.dot_plot.show { 280.0 } else { 0.0 };
            // *2.0 because the overlay appears between Ref↔Hap1 and Hap1↔Hap2
            let available = ui.available_height() - dot_plot_height - sv_overlay_height * 2.0;
            let panel_height = (available - 16.0) / 3.0;

            ui.allocate_ui(egui::vec2(ui.available_width(), panel_height), |ui| {
                let inter = Self::show_panel(
                    ui,
                    &PanelRenderParams {
                        panel: Panel::Reference,
                        region_text: &ref_text,
                        data: &self.ref_data,
                        config: &config,
                        view_start: ref_start,
                        view_end: ref_end,
                        coord_info: ref_coord_info.as_deref(),
                        metadata: metadata_text.as_deref(),
                    },
                );
                interactions.push((PanelId::Reference, inter));
            });

            // SV overlay between Reference and Hap1
            if self.show_sv_overlay {
                self.paint_sv_overlay(ui);
            }

            ui.add_space(4.0);

            ui.allocate_ui(egui::vec2(ui.available_width(), panel_height), |ui| {
                let inter = Self::show_panel(
                    ui,
                    &PanelRenderParams {
                        panel: Panel::Haplotype1,
                        region_text: &hap1_text,
                        data: &self.hap1_data,
                        config: &config,
                        view_start: h1_start,
                        view_end: h1_end,
                        coord_info: hap1_coord_info.as_deref(),
                        metadata: None,
                    },
                );
                interactions.push((PanelId::Haplotype1, inter));
            });

            // SV overlay between Hap1 and Hap2
            if self.show_sv_overlay {
                self.paint_sv_overlay(ui);
            }

            ui.add_space(4.0);

            ui.allocate_ui(egui::vec2(ui.available_width(), panel_height), |ui| {
                let inter = Self::show_panel(
                    ui,
                    &PanelRenderParams {
                        panel: Panel::Haplotype2,
                        region_text: &hap2_text,
                        data: &self.hap2_data,
                        config: &config,
                        view_start: h2_start,
                        view_end: h2_end,
                        coord_info: hap2_coord_info.as_deref(),
                        metadata: None,
                    },
                );
                interactions.push((PanelId::Haplotype2, inter));
            });

            // Dot plot panel
            if self.dot_plot.show {
                ui.add_space(4.0);
                ui.group(|ui| {
                    ui.set_min_height(250.0);
                    self.show_dot_plot(ui);
                });
            }
        });

        // Process mouse interactions from panels
        for (panel_id, inter) in interactions {
            let mut needs_reload = false;
            if inter.pan_delta != 0 {
                self.sync_manager.pan(panel_id, inter.pan_delta);
                needs_reload = true;
            }
            if let Some((factor, cursor_frac)) = inter.zoom {
                self.sync_manager.zoom_at(panel_id, factor, cursor_frac);
                needs_reload = true;
            }
            if needs_reload {
                self.schedule_debounced_reload();
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/// Merge coordinate-mapper query results into a single assembly region.
///
/// Port of the Python `_collapse_asm_regions()` function.  Only alignment
/// events (not SV gap events) are considered.  The intervals are merged
/// per contig and the first merged region is returned as a [`GenomicRegion`].
fn collapse_asm_regions(
    results: &[genome::coordinate_mapper::QueryResult],
) -> Option<crate::region::GenomicRegion> {
    use std::collections::BTreeMap;

    let mut by_contig: BTreeMap<&str, Vec<(u64, u64)>> = BTreeMap::new();
    for r in results {
        if r.event_type != genome::coordinate_mapper::EventType::Alignment {
            continue;
        }
        let (s, e) = if r.asm_start <= r.asm_end {
            (r.asm_start, r.asm_end)
        } else {
            (r.asm_end, r.asm_start)
        };
        by_contig.entry(&r.asm_chrom).or_default().push((s, e));
    }

    // Return the first merged region (covering the most contiguous aligned span).
    for (contig, mut intervals) in by_contig {
        if intervals.is_empty() {
            continue;
        }
        intervals.sort();
        let mut merged: Vec<(u64, u64)> = vec![intervals[0]];
        for &(s, e) in &intervals[1..] {
            // Safety: merged always has at least one element (initialised above).
            let last = merged.last_mut().unwrap();
            if s <= last.1 {
                last.1 = last.1.max(e);
            } else {
                merged.push((s, e));
            }
        }
        // Use the bounding interval across all merged segments on this contig.
        // Safety: merged was initialised with intervals[0] so it is non-empty.
        let start = merged.first().unwrap().0;
        let end = merged.last().unwrap().1;
        // Ensure start is at least 1 for GenomicRegion (1-based).
        let start = start.max(1);
        if end >= start {
            return Some(crate::region::GenomicRegion {
                chrom: contig.to_string(),
                start,
                end,
            });
        }
    }
    None
}

/// Load a FASTA sequence for a genomic region.
///
/// Uses `query_fasta_region()` with the region's chrom as the exact sequence
/// name and local coordinates, to handle FASTAs with colons in sequence names.
fn load_fasta_for_region(
    fasta_path: &std::path::Path,
    region: &crate::region::GenomicRegion,
) -> Option<FastaSequence> {
    if region.end < region.start {
        return None;
    }

    // First try: use the full "chrom:start-end" as a sequence name (for
    // extracted sub-FASTA files where the name includes coordinates).
    let full_name = region.to_string();
    if let Ok(seqs) = genome::fasta::list_sequences(fasta_path) {
        // Check if any sequence name matches the region's chrom or full string.
        for (name, len) in &seqs {
            if *name == full_name {
                // The full region string IS the sequence name; query local coords.
                let end = (*len).min(region.end - region.start + 1);
                return genome::fasta::query_fasta_region(fasta_path, name, 1, end).ok();
            }
            if *name == region.chrom {
                // Standard case: chrom matches exactly, use region coordinates.
                return genome::fasta::query_fasta_region(
                    fasta_path,
                    name,
                    region.start,
                    region.end,
                )
                .ok();
            }
            // Handle assembly contig names (e.g., "NA21110#1#CM089663.1:111967298-112167366")
            // where the region chrom is the full contig name.  Use exact prefix
            // match with a delimiter boundary to avoid false positives like
            // "chr1" matching "chr10".
            if name.starts_with(&region.chrom)
                && name[region.chrom.len()..].starts_with(|c: char| !c.is_ascii_alphanumeric())
            {
                let end = (*len).min(region.end - region.start + 1);
                return genome::fasta::query_fasta_region(fasta_path, name, 1, end).ok();
            }
        }
    }
    // Fallback: try direct query (may fail for colon-containing names).
    genome::fasta::query_fasta(fasta_path, &region.to_string()).ok()
}

/// Load assembly cross-alignment tracks from a list of BAM file descriptors.
///
/// Each entry is `(label, color, optional_bam_path)`. BAM files that are
/// `None` or fail to query are silently skipped (graceful degradation).
/// Assembly tracks always use squished (compact) display.
fn load_assembly_tracks(
    track_specs: &[(&str, [u8; 3], Option<&std::path::Path>)],
    region_str: &str,
) -> Vec<AssemblyTrack> {
    let mut tracks = Vec::new();
    for &(label, color, bam_path) in track_specs {
        let Some(path) = bam_path else { continue };
        if let Ok(reads) = genome::bam::query_bam(path, region_str)
            && !reads.is_empty()
        {
            let rows = pileup::pack_reads(reads);
            tracks.push(AssemblyTrack {
                label: label.to_string(),
                color,
                rows,
            });
        }
    }
    tracks
}

/// Attempt to load a coordinate mapping index from an optional path.
fn try_load_coord_index(
    path: Option<&std::path::Path>,
) -> Option<genome::coordinate_mapper::MappingIndex> {
    let p = path?;
    genome::coordinate_mapper::load_index(&p.to_string_lossy()).ok()
}

/// Configure default fonts/styles.
fn configure_fonts(ctx: &egui::Context) {
    let mut style = (*ctx.style()).clone();
    style
        .text_styles
        .insert(egui::TextStyle::Body, egui::FontId::proportional(14.0));
    style
        .text_styles
        .insert(egui::TextStyle::Button, egui::FontId::proportional(13.0));
    ctx.set_style(style);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_panel_labels() {
        assert_eq!(Panel::Reference.label(), "Reference");
        assert_eq!(Panel::Haplotype1.label(), "Haplotype 1");
        assert_eq!(Panel::Haplotype2.label(), "Haplotype 2");
    }

    #[test]
    fn test_panel_colors_distinct() {
        let c1 = Panel::Reference.color();
        let c2 = Panel::Haplotype1.color();
        let c3 = Panel::Haplotype2.color();
        assert_ne!(c1, c2);
        assert_ne!(c2, c3);
        assert_ne!(c1, c3);
    }

    #[test]
    fn test_panel_sync_ids() {
        assert_eq!(Panel::Reference.sync_id(), PanelId::Reference);
        assert_eq!(Panel::Haplotype1.sync_id(), PanelId::Haplotype1);
        assert_eq!(Panel::Haplotype2.sync_id(), PanelId::Haplotype2);
    }

    #[test]
    fn test_viewer_app_default() {
        let app = ViewerApp::default();
        assert!(app.navigator.is_empty());
        assert!(app.status_message.contains("No regions loaded"));
        assert!(app.display_config.squished);
        assert!(app.display_config.hide_small_indels);
        assert!(app.sync_manager.sync_enabled);
    }

    #[test]
    fn test_update_status_for_current_region() {
        let mut app = ViewerApp::default();
        let json = r#"{
          "variants": [{
            "chrom": "chr1", "pos": 100, "size": 50,
            "ref_region": "chr1:100-150"
          }]
        }"#;
        app.navigator.load_manifest_str(json).unwrap();
        app.load_region_data();
        app.update_status_for_current_region();
        assert!(app.status_message.contains("1/1"));
        assert!(app.status_message.contains("chr1:100-150"));
        // Verify sync manager was updated
        assert_eq!(app.sync_manager.view(PanelId::Reference).view_start, 100);
        assert_eq!(app.sync_manager.view(PanelId::Reference).view_end, 150);
    }

    #[test]
    fn test_display_config_toggle() {
        let mut app = ViewerApp::default();
        assert!(app.display_config.squished);
        app.display_config.squished = false;
        assert!(!app.display_config.squished);

        assert!(app.display_config.hide_small_indels);
        app.display_config.hide_small_indels = false;
        assert!(!app.display_config.hide_small_indels);
    }

    // -- Integration: navigation updates sync_manager -----------------------

    #[test]
    fn test_navigation_updates_sync_manager() {
        let mut app = ViewerApp::default();
        let json = r#"{
          "variants": [
            {
              "chrom": "chr1", "pos": 1000, "size": 500,
              "ref_region": "chr1:1000-1500",
              "hap1_regions": ["ctg1:5000-5500"],
              "hap2_regions": ["ctg2:8000-8500"]
            },
            {
              "chrom": "chr2", "pos": 2000, "size": 300,
              "ref_region": "chr2:2000-2300",
              "hap1_regions": ["ctg3:6000-6300"],
              "hap2_regions": ["ctg4:9000-9300"]
            }
          ]
        }"#;
        app.navigator.load_manifest_str(json).unwrap();
        app.load_region_data();
        app.update_status_for_current_region();

        // Verify region 1
        assert_eq!(app.sync_manager.view(PanelId::Reference).view_start, 1000);
        assert_eq!(app.sync_manager.view(PanelId::Haplotype1).view_start, 5000);
        assert_eq!(app.sync_manager.view(PanelId::Haplotype2).view_start, 8000);

        // Navigate to next region
        app.navigator.next();
        app.load_region_data();
        app.update_status_for_current_region();

        // Verify region 2 — all panels updated
        assert_eq!(app.sync_manager.view(PanelId::Reference).view_start, 2000);
        assert_eq!(app.sync_manager.view(PanelId::Reference).view_end, 2300);
        assert_eq!(app.sync_manager.view(PanelId::Haplotype1).view_start, 6000);
        assert_eq!(app.sync_manager.view(PanelId::Haplotype1).view_end, 6300);
        assert_eq!(app.sync_manager.view(PanelId::Haplotype2).view_start, 9000);
        assert_eq!(app.sync_manager.view(PanelId::Haplotype2).view_end, 9300);
    }

    #[test]
    fn test_zoom_updates_all_panels_synced() {
        let mut app = ViewerApp::default();
        let json = r#"{
          "variants": [{
            "chrom": "chr1", "pos": 1000, "size": 1000,
            "ref_region": "chr1:1000-2000",
            "hap1_regions": ["ctg1:5000-6000"],
            "hap2_regions": ["ctg2:8000-9000"]
          }]
        }"#;
        app.navigator.load_manifest_str(json).unwrap();
        app.load_region_data();
        app.update_status_for_current_region();

        // Zoom in on reference panel
        app.sync_manager.zoom(PanelId::Reference, 2.0);

        // All panels should have matching spans
        assert!(app.sync_manager.all_spans_match());
        assert_eq!(app.sync_manager.view(PanelId::Reference).span(), 500);
    }

    #[test]
    fn test_pan_updates_all_panels_synced() {
        let mut app = ViewerApp::default();
        let json = r#"{
          "variants": [{
            "chrom": "chr1", "pos": 1000, "size": 1000,
            "ref_region": "chr1:1000-2000",
            "hap1_regions": ["ctg1:5000-6000"],
            "hap2_regions": ["ctg2:8000-9000"]
          }]
        }"#;
        app.navigator.load_manifest_str(json).unwrap();
        app.load_region_data();
        app.update_status_for_current_region();

        // Pan right from reference
        app.sync_manager.pan(PanelId::Reference, 200);

        // All panels moved
        assert_eq!(app.sync_manager.view(PanelId::Reference).view_start, 1200);
        assert_eq!(app.sync_manager.view(PanelId::Haplotype1).view_start, 5200);
        assert_eq!(app.sync_manager.view(PanelId::Haplotype2).view_start, 8200);
    }

    #[test]
    fn test_data_paths_default() {
        let dp = DataPaths::default();
        assert!(dp.reference_fasta.is_none());
        assert!(dp.hap1_fasta.is_none());
        assert!(dp.hap2_fasta.is_none());
        assert!(dp.reads_bam.is_none());
        assert!(dp.coordinate_index.is_none());
    }

    #[test]
    fn test_dot_plot_state_default() {
        let app = ViewerApp::default();
        assert!(!app.dot_plot.show);
        assert_eq!(app.dot_plot.comparison, DotPlotComparison::RefVsHap1);
        assert!(app.dot_plot.result.is_none());
    }

    #[test]
    fn test_recompute_dot_plot_no_sequences() {
        let mut app = ViewerApp::default();
        app.recompute_dot_plot();
        assert!(app.dot_plot.result.is_none());
    }

    #[test]
    fn test_recompute_dot_plot_with_sequences() {
        let mut app = ViewerApp::default();
        app.ref_data.sequence = Some(FastaSequence {
            name: "test".to_string(),
            start: 1,
            end: 100,
            sequence: "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT".to_string(),
        });
        app.hap1_data.sequence = Some(FastaSequence {
            name: "test2".to_string(),
            start: 1,
            end: 100,
            sequence: "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT".to_string(),
        });
        app.dot_plot.k = 5;
        app.recompute_dot_plot();
        assert!(app.dot_plot.result.is_some());
        let result = app.dot_plot.result.as_ref().unwrap();
        assert!(result.total_matches() > 0);
    }

    #[test]
    fn test_load_region_data_without_files() {
        let mut app = ViewerApp::default();
        let json = r#"{
          "variants": [{
            "chrom": "chr1", "pos": 1000, "size": 500,
            "ref_region": "chr1:1000-1500"
          }]
        }"#;
        app.navigator.load_manifest_str(json).unwrap();
        app.load_region_data();
        // No files loaded, panels should be empty
        assert!(app.ref_data.rows.is_empty());
        assert!(app.ref_data.sequence.is_none());
    }

    #[test]
    fn test_data_paths_from_tsv() {
        let tmp = tempfile::tempdir().unwrap();
        let tsv_path = tmp.path().join("config.tsv");
        std::fs::write(
            &tsv_path,
            "sample_id\treference\thap1_assembly\thap2_assembly\treads_bam\tregions\n\
             sample1\t/path/to/ref.fa.gz\t/path/to/hap1.fa.gz\t/path/to/hap2.fa.gz\t/path/to/reads.bam\t/path/to/manifest.json\n",
        )
        .unwrap();

        let dp = DataPaths::from_tsv(&tsv_path).unwrap();
        assert_eq!(
            dp.reference_fasta.as_deref(),
            Some(std::path::Path::new("/path/to/ref.fa.gz"))
        );
        assert_eq!(
            dp.hap1_fasta.as_deref(),
            Some(std::path::Path::new("/path/to/hap1.fa.gz"))
        );
        assert_eq!(
            dp.hap2_fasta.as_deref(),
            Some(std::path::Path::new("/path/to/hap2.fa.gz"))
        );
        assert_eq!(
            dp.reads_bam.as_deref(),
            Some(std::path::Path::new("/path/to/reads.bam"))
        );
        assert_eq!(
            dp.regions.as_deref(),
            Some(std::path::Path::new("/path/to/manifest.json"))
        );
    }

    #[test]
    fn test_data_paths_from_tsv_with_comment_header() {
        let tmp = tempfile::tempdir().unwrap();
        let tsv_path = tmp.path().join("config.tsv");
        std::fs::write(
            &tsv_path,
            "#sample_id\treference\thap1_assembly\thap2_assembly\treads_bam\n\
             sample1\t/ref.fa.gz\t/hap1.fa.gz\t/hap2.fa.gz\t/reads.bam\n",
        )
        .unwrap();

        let dp = DataPaths::from_tsv(&tsv_path).unwrap();
        assert_eq!(
            dp.reference_fasta.as_deref(),
            Some(std::path::Path::new("/ref.fa.gz"))
        );
        assert_eq!(
            dp.hap1_fasta.as_deref(),
            Some(std::path::Path::new("/hap1.fa.gz"))
        );
    }

    #[test]
    fn test_data_paths_from_tsv_empty_optional_fields() {
        let tmp = tempfile::tempdir().unwrap();
        let tsv_path = tmp.path().join("config.tsv");
        std::fs::write(
            &tsv_path,
            "sample_id\treference\thap1_assembly\thap2_assembly\treads_bam\tregions\n\
             sample1\t/ref.fa.gz\t/hap1.fa.gz\t/hap2.fa.gz\t\t\n",
        )
        .unwrap();

        let dp = DataPaths::from_tsv(&tsv_path).unwrap();
        assert_eq!(
            dp.reference_fasta.as_deref(),
            Some(std::path::Path::new("/ref.fa.gz"))
        );
        assert!(dp.reads_bam.is_none());
        assert!(dp.regions.is_none());
    }

    #[test]
    fn test_data_paths_from_tsv_missing_file() {
        let result = DataPaths::from_tsv(std::path::Path::new("/nonexistent/config.tsv"));
        assert!(result.is_err());
    }

    #[test]
    fn test_data_paths_default_new_fields() {
        let dp = DataPaths::default();
        assert!(dp.hap1_coord_index.is_none());
        assert!(dp.hap2_coord_index.is_none());
        assert!(dp.reads_to_hap1_bam.is_none());
        assert!(dp.reads_to_hap2_bam.is_none());
        assert!(dp.cram_ref.is_none());
    }

    #[test]
    fn test_data_paths_from_tsv_discovers_output_files() {
        let tmp = tempfile::tempdir().unwrap();
        let output_dir = tmp.path().join("output").join("NA21110");
        std::fs::create_dir_all(&output_dir).unwrap();

        // Create the files that from_tsv should discover
        for suffix in &[
            "_hap1_to_ref.mapping.json.gz",
            "_hap2_to_ref.mapping.json.gz",
            "_reads_to_hap1.bam",
            "_reads_to_hap2.bam",
        ] {
            let path = output_dir.join(format!("NA21110{suffix}"));
            std::fs::write(&path, b"dummy").unwrap();
        }

        let tsv_path = tmp.path().join("config.tsv");
        std::fs::write(
            &tsv_path,
            format!(
                "#sample_id\toutput_dir\treference\thap1_assembly\thap2_assembly\treads_bam\tregions\tcram_ref\n\
                 NA21110\t{}\t/ref.fa.gz\t/hap1.fa.gz\t/hap2.fa.gz\t/reads.cram\t/regions.vcf.gz\t/cram_ref.fa.gz\n",
                output_dir.display()
            ),
        )
        .unwrap();

        let dp = DataPaths::from_tsv(&tsv_path).unwrap();
        assert!(dp.hap1_coord_index.is_some());
        assert!(dp.hap2_coord_index.is_some());
        assert!(dp.reads_to_hap1_bam.is_some());
        assert!(dp.reads_to_hap2_bam.is_some());
        assert_eq!(
            dp.cram_ref.as_deref(),
            Some(std::path::Path::new("/cram_ref.fa.gz"))
        );
        // Verify discovered file names contain expected suffix
        assert!(
            dp.hap1_coord_index
                .as_ref()
                .unwrap()
                .to_string_lossy()
                .ends_with("_hap1_to_ref.mapping.json.gz")
        );
        assert!(
            dp.reads_to_hap2_bam
                .as_ref()
                .unwrap()
                .to_string_lossy()
                .ends_with("_reads_to_hap2.bam")
        );
    }

    #[test]
    fn test_data_paths_from_tsv_no_output_dir_no_discovery() {
        let tmp = tempfile::tempdir().unwrap();
        let tsv_path = tmp.path().join("config.tsv");
        std::fs::write(
            &tsv_path,
            "sample_id\treference\thap1_assembly\thap2_assembly\treads_bam\n\
             sample1\t/ref.fa.gz\t/hap1.fa.gz\t/hap2.fa.gz\t/reads.bam\n",
        )
        .unwrap();

        let dp = DataPaths::from_tsv(&tsv_path).unwrap();
        // Without output_dir column, discovery should yield None
        assert!(dp.hap1_coord_index.is_none());
        assert!(dp.hap2_coord_index.is_none());
        assert!(dp.reads_to_hap1_bam.is_none());
        assert!(dp.reads_to_hap2_bam.is_none());
    }

    #[test]
    fn test_collapse_asm_regions_empty() {
        assert!(collapse_asm_regions(&[]).is_none());
    }

    #[test]
    fn test_collapse_asm_regions_single_alignment() {
        use crate::genome::coordinate_mapper::{EventType, QueryResult};
        let results = vec![QueryResult {
            ref_chrom: "chr1".to_string(),
            ref_start: 1000,
            ref_end: 2000,
            asm_chrom: "ctg1".to_string(),
            asm_start: 5000,
            asm_end: 6000,
            strand: "+".to_string(),
            mapq: 60,
            event_type: EventType::Alignment,
            ref_gap_size: None,
            asm_gap_size: None,
        }];
        let region = collapse_asm_regions(&results).unwrap();
        assert_eq!(region.chrom, "ctg1");
        assert_eq!(region.start, 5000);
        assert_eq!(region.end, 6000);
    }

    #[test]
    fn test_collapse_asm_regions_skips_sv_events() {
        use crate::genome::coordinate_mapper::{EventType, QueryResult};
        let results = vec![
            QueryResult {
                ref_chrom: "chr1".to_string(),
                ref_start: 1000,
                ref_end: 2000,
                asm_chrom: "ctg1".to_string(),
                asm_start: 5000,
                asm_end: 6000,
                strand: "+".to_string(),
                mapq: 60,
                event_type: EventType::Alignment,
                ref_gap_size: None,
                asm_gap_size: None,
            },
            QueryResult {
                ref_chrom: "chr1".to_string(),
                ref_start: 2000,
                ref_end: 2500,
                asm_chrom: "ctg1".to_string(),
                asm_start: 6000,
                asm_end: 6000,
                strand: "+".to_string(),
                mapq: 60,
                event_type: EventType::Deletion,
                ref_gap_size: Some(500),
                asm_gap_size: Some(0),
            },
            QueryResult {
                ref_chrom: "chr1".to_string(),
                ref_start: 2500,
                ref_end: 3500,
                asm_chrom: "ctg1".to_string(),
                asm_start: 6000,
                asm_end: 7000,
                strand: "+".to_string(),
                mapq: 60,
                event_type: EventType::Alignment,
                ref_gap_size: None,
                asm_gap_size: None,
            },
        ];
        // Should merge two alignment blocks, skip the deletion
        let region = collapse_asm_regions(&results).unwrap();
        assert_eq!(region.chrom, "ctg1");
        assert_eq!(region.start, 5000);
        assert_eq!(region.end, 7000);
    }

    #[test]
    fn test_collapse_asm_regions_only_sv_events() {
        use crate::genome::coordinate_mapper::{EventType, QueryResult};
        let results = vec![QueryResult {
            ref_chrom: "chr1".to_string(),
            ref_start: 1000,
            ref_end: 2000,
            asm_chrom: "ctg1".to_string(),
            asm_start: 5000,
            asm_end: 5000,
            strand: "+".to_string(),
            mapq: 60,
            event_type: EventType::Insertion,
            ref_gap_size: Some(1000),
            asm_gap_size: Some(5000),
        }];
        // No alignment events → None
        assert!(collapse_asm_regions(&results).is_none());
    }

    // -----------------------------------------------------------------------
    // DataPaths cross-alignment fields
    // -----------------------------------------------------------------------

    #[test]
    fn test_data_paths_default_has_no_cross_alignment_bams() {
        let dp = DataPaths::default();
        assert!(dp.hap1_to_ref_bam.is_none());
        assert!(dp.hap2_to_ref_bam.is_none());
        assert!(dp.ref_to_hap1_bam.is_none());
        assert!(dp.ref_to_hap2_bam.is_none());
        assert!(dp.hap1_to_hap2_bam.is_none());
        assert!(dp.hap2_to_hap1_bam.is_none());
    }

    #[test]
    fn test_data_paths_cross_alignment_fields_settable() {
        let dp = DataPaths {
            hap1_to_ref_bam: Some(PathBuf::from("/tmp/h1_to_ref.bam")),
            ref_to_hap1_bam: Some(PathBuf::from("/tmp/ref_to_h1.bam")),
            hap2_to_hap1_bam: Some(PathBuf::from("/tmp/h2_to_h1.bam")),
            ..Default::default()
        };
        assert!(dp.hap1_to_ref_bam.is_some());
        assert!(dp.ref_to_hap1_bam.is_some());
        assert!(dp.hap2_to_hap1_bam.is_some());
    }

    // -----------------------------------------------------------------------
    // PanelData assembly tracks
    // -----------------------------------------------------------------------

    #[test]
    fn test_panel_data_default_has_empty_assembly_tracks() {
        let pd = PanelData::default();
        assert!(pd.assembly_tracks.is_empty());
    }

    // -----------------------------------------------------------------------
    // Assembly track loading with missing files (graceful degradation)
    // -----------------------------------------------------------------------

    #[test]
    fn test_load_assembly_tracks_with_no_bams() {
        let tracks = load_assembly_tracks(
            &[
                ("Hap1 → Ref", [60, 160, 80], None),
                ("Hap2 → Ref", [180, 100, 60], None),
            ],
            "chr1:100-200",
        );
        assert!(
            tracks.is_empty(),
            "should gracefully return empty when no BAMs"
        );
    }

    #[test]
    fn test_load_assembly_tracks_with_nonexistent_bam() {
        let tracks = load_assembly_tracks(
            &[(
                "Hap1 → Ref",
                [60, 160, 80],
                Some(std::path::Path::new("/nonexistent/file.bam")),
            )],
            "chr1:100-200",
        );
        assert!(
            tracks.is_empty(),
            "should gracefully skip nonexistent BAM files"
        );
    }

    // -----------------------------------------------------------------------
    // New toolbar toggle defaults
    // -----------------------------------------------------------------------

    #[test]
    fn test_viewer_app_default_new_toggles() {
        let app = ViewerApp::default();
        assert!(
            !app.display_config.show_mismatches,
            "mismatches off by default"
        );
        assert!(
            !app.display_config.show_soft_clips,
            "soft clips off by default"
        );
        assert!(
            !app.display_config.sort_by_haplotype,
            "HP sort off by default"
        );
        assert_eq!(
            app.display_config.indel_threshold,
            pileup::DEFAULT_INDEL_THRESHOLD
        );
        assert_eq!(app.display_config.indel_threshold, 50);
    }

    // -----------------------------------------------------------------------
    // Phase 4: Colorblind palette defaults
    // -----------------------------------------------------------------------

    #[test]
    fn test_viewer_app_default_colorblind_off() {
        let app = ViewerApp::default();
        assert!(
            !app.display_config.use_colorblind_palette,
            "colorblind palette off by default"
        );
    }

    // -----------------------------------------------------------------------
    // Phase 4: Multi-sample support
    // -----------------------------------------------------------------------

    #[test]
    fn test_parse_all_samples_single() {
        let tmp = tempfile::tempdir().unwrap();
        let tsv_path = tmp.path().join("config.tsv");
        std::fs::write(
            &tsv_path,
            "sample_id\treference\thap1_assembly\thap2_assembly\n\
             sample1\t/ref.fa.gz\t/hap1.fa.gz\t/hap2.fa.gz\n",
        )
        .unwrap();

        let samples = parse_all_samples(&tsv_path).unwrap();
        assert_eq!(samples.len(), 1);
        assert_eq!(samples[0].0, "sample1");
        assert_eq!(
            samples[0].1.reference_fasta.as_deref(),
            Some(std::path::Path::new("/ref.fa.gz"))
        );
    }

    #[test]
    fn test_parse_all_samples_multiple() {
        let tmp = tempfile::tempdir().unwrap();
        let tsv_path = tmp.path().join("config.tsv");
        std::fs::write(
            &tsv_path,
            "#sample_id\treference\thap1_assembly\thap2_assembly\n\
             sample1\t/ref1.fa.gz\t/hap1a.fa.gz\t/hap2a.fa.gz\n\
             sample2\t/ref2.fa.gz\t/hap1b.fa.gz\t/hap2b.fa.gz\n\
             sample3\t/ref3.fa.gz\t/hap1c.fa.gz\t/hap2c.fa.gz\n",
        )
        .unwrap();

        let samples = parse_all_samples(&tsv_path).unwrap();
        assert_eq!(samples.len(), 3);
        assert_eq!(samples[0].0, "sample1");
        assert_eq!(samples[1].0, "sample2");
        assert_eq!(samples[2].0, "sample3");
        assert_eq!(
            samples[1].1.reference_fasta.as_deref(),
            Some(std::path::Path::new("/ref2.fa.gz"))
        );
    }

    #[test]
    fn test_parse_all_samples_empty_file() {
        let tmp = tempfile::tempdir().unwrap();
        let tsv_path = tmp.path().join("config.tsv");
        std::fs::write(&tsv_path, "sample_id\treference\n").unwrap();

        let result = parse_all_samples(&tsv_path);
        assert!(result.is_err());
    }

    #[test]
    fn test_viewer_app_default_samples_empty() {
        let app = ViewerApp::default();
        assert!(app.samples.is_empty());
        assert_eq!(app.selected_sample_idx, 0);
    }

    // -----------------------------------------------------------------------
    // Phase 4: Compare Reads
    // -----------------------------------------------------------------------

    #[test]
    fn test_compare_reads_empty_panels() {
        let app = ViewerApp::default();
        let result = app.compare_reads();
        assert_eq!(result.total_unique, 0);
        assert_eq!(result.ref_only, 0);
        assert_eq!(result.hap1_only, 0);
        assert_eq!(result.hap2_only, 0);
        assert_eq!(result.all_three, 0);
    }

    #[test]
    fn test_compare_reads_with_data() {
        use crate::genome::AlignedRead;

        let make_read = |name: &str| AlignedRead {
            name: name.to_string(),
            start: 100,
            end: 200,
            is_reverse: false,
            mapping_quality: Some(60),
            haplotype: None,
            flags: 0,
            indels: vec![],
            soft_clips: vec![],
            mismatches: vec![],
        };

        let mut app = ViewerApp::default();
        // Ref panel: reads A, B, C
        app.ref_data.rows = vec![PileupRow {
            y_offset: 0,
            reads: vec![make_read("A"), make_read("B"), make_read("C")],
        }];
        // Hap1 panel: reads B, D
        app.hap1_data.rows = vec![PileupRow {
            y_offset: 0,
            reads: vec![make_read("B"), make_read("D")],
        }];
        // Hap2 panel: reads C, D, E
        app.hap2_data.rows = vec![PileupRow {
            y_offset: 0,
            reads: vec![make_read("C"), make_read("D"), make_read("E")],
        }];

        let result = app.compare_reads();
        assert_eq!(result.total_unique, 5); // A, B, C, D, E
        assert_eq!(result.ref_only, 1); // A
        assert_eq!(result.hap1_only, 0); // (B shared with ref, D shared with hap2)
        assert_eq!(result.hap2_only, 1); // E
        assert_eq!(result.ref_and_hap1, 1); // B
        assert_eq!(result.ref_and_hap2, 1); // C
        assert_eq!(result.hap1_and_hap2, 1); // D
        assert_eq!(result.all_three, 0);
    }

    #[test]
    fn test_compare_reads_all_shared() {
        use crate::genome::AlignedRead;

        let make_read = |name: &str| AlignedRead {
            name: name.to_string(),
            start: 100,
            end: 200,
            is_reverse: false,
            mapping_quality: Some(60),
            haplotype: None,
            flags: 0,
            indels: vec![],
            soft_clips: vec![],
            mismatches: vec![],
        };

        let mut app = ViewerApp::default();
        let reads = vec![make_read("X"), make_read("Y")];
        app.ref_data.rows = vec![PileupRow {
            y_offset: 0,
            reads: reads.clone(),
        }];
        app.hap1_data.rows = vec![PileupRow {
            y_offset: 0,
            reads: reads.clone(),
        }];
        app.hap2_data.rows = vec![PileupRow { y_offset: 0, reads }];

        let result = app.compare_reads();
        assert_eq!(result.total_unique, 2);
        assert_eq!(result.all_three, 2);
        assert_eq!(result.ref_only, 0);
        assert_eq!(result.hap1_only, 0);
        assert_eq!(result.hap2_only, 0);
    }

    #[test]
    fn test_compare_reads_result_default() {
        let result = CompareReadsResult::default();
        assert_eq!(result.total_unique, 0);
        assert_eq!(result.ref_only, 0);
    }

    #[test]
    fn test_viewer_app_default_compare_reads_off() {
        let app = ViewerApp::default();
        assert!(!app.show_compare_reads);
        assert!(app.compare_result.is_none());
    }

    // -----------------------------------------------------------------------
    // Phase 4: SV Annotation overlay
    // -----------------------------------------------------------------------

    #[test]
    fn test_sv_event_color_deletion() {
        use crate::genome::coordinate_mapper::EventType;
        assert_eq!(sv_event_color(&EventType::Deletion), [200, 60, 60]);
    }

    #[test]
    fn test_sv_event_color_insertion() {
        use crate::genome::coordinate_mapper::EventType;
        assert_eq!(sv_event_color(&EventType::Insertion), [60, 100, 200]);
    }

    #[test]
    fn test_sv_event_color_inversion() {
        use crate::genome::coordinate_mapper::EventType;
        assert_eq!(sv_event_color(&EventType::Inversion), [150, 60, 200]);
    }

    #[test]
    fn test_sv_event_color_translocation() {
        use crate::genome::coordinate_mapper::EventType;
        assert_eq!(sv_event_color(&EventType::Translocation), [200, 200, 60]);
    }

    #[test]
    fn test_sv_event_color_complex() {
        use crate::genome::coordinate_mapper::EventType;
        assert_eq!(sv_event_color(&EventType::Complex), [140, 140, 140]);
    }

    #[test]
    fn test_sv_annotation_color_method() {
        use crate::genome::coordinate_mapper::EventType;
        let ann = SvAnnotation {
            event_type: EventType::Deletion,
            ref_start: 100,
            ref_end: 200,
            asm_start: 100,
            asm_end: 100,
            gap_size: Some(100),
            label: "test".to_string(),
        };
        assert_eq!(ann.color(), [200, 60, 60]);
    }

    #[test]
    fn test_viewer_app_default_sv_overlay_off() {
        let app = ViewerApp::default();
        assert!(!app.show_sv_overlay);
        assert!(app.sv_annotations.is_empty());
    }

    #[test]
    fn test_collect_sv_annotations_no_index() {
        let mut app = ViewerApp::default();
        let json = r#"{
          "variants": [{
            "chrom": "chr1", "pos": 100, "size": 50,
            "ref_region": "chr1:100-150"
          }]
        }"#;
        app.navigator.load_manifest_str(json).unwrap();
        let annotations = app.collect_sv_annotations();
        assert!(annotations.is_empty());
    }

    // -----------------------------------------------------------------------
    // Phase 4: Export/Screenshot defaults
    // -----------------------------------------------------------------------

    #[test]
    fn test_viewer_app_default_pending_screenshot() {
        let app = ViewerApp::default();
        assert!(!app.pending_screenshot);
    }

    // -----------------------------------------------------------------------
    // Phase 4: Variant metadata in status message
    // -----------------------------------------------------------------------

    #[test]
    fn test_status_message_includes_metadata() {
        let mut app = ViewerApp::default();
        let json = r#"{
          "variants": [{
            "chrom": "chr1", "pos": 100, "size": 5000,
            "genotype": "0|1",
            "ref_region": "chr1:100-5100"
          }]
        }"#;
        app.navigator.load_manifest_str(json).unwrap();
        app.load_region_data();
        app.update_status_for_current_region();
        assert!(
            app.status_message.contains("GT: 0|1"),
            "status should contain genotype, got: {}",
            app.status_message
        );
        assert!(
            app.status_message.contains("size: 5000 bp"),
            "status should contain size, got: {}",
            app.status_message
        );
    }
}
