use std::path::PathBuf;

use eframe::egui;

use crate::dot_plot::{self, DEFAULT_K, DotPlotComparison, DotPlotResult};
use crate::genome::pileup::{self, PileupDisplayConfig, PileupRow, ReadRect};
use crate::genome::{self, FastaSequence};
use crate::panel_sync::{PanelId, PanelSyncManager};
use crate::region::RegionNavigator;

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
    /// Optional coordinate mapping index (gzipped JSON).
    pub coordinate_index: Option<PathBuf>,
    /// Optional path to the manifest/regions JSON loaded from config.
    pub regions: Option<PathBuf>,
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

        Ok(Self {
            reference_fasta: non_empty("reference"),
            hap1_fasta: non_empty("hap1_assembly"),
            hap2_fasta: non_empty("hap2_assembly"),
            reads_bam: non_empty("reads_bam"),
            coordinate_index: None, // Not in TSV format; discovered via output_dir
            regions: non_empty("regions"),
        })
    }
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
    /// Optional coordinate mapping index.
    coord_index: Option<genome::coordinate_mapper::MappingIndex>,
}

impl Default for ViewerApp {
    fn default() -> Self {
        Self {
            navigator: RegionNavigator::new(),
            status_message:
                "No regions loaded. Use File > Load Manifest to open a region manifest JSON."
                    .to_string(),
            display_config: PileupDisplayConfig::default(),
            sync_manager: PanelSyncManager::new(),
            data_paths: DataPaths::default(),
            ref_data: PanelData::default(),
            hap1_data: PanelData::default(),
            hap2_data: PanelData::default(),
            dot_plot: DotPlotState::default(),
            coord_index: None,
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
}

impl ViewerApp {
    pub fn new(
        cc: &eframe::CreationContext<'_>,
        manifest_path: Option<&std::path::Path>,
        data_paths: DataPaths,
    ) -> Self {
        configure_fonts(&cc.egui_ctx);
        let mut app = Self {
            data_paths,
            ..Self::default()
        };

        // Load coordinate mapping index if provided.
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

        // Load manifest if provided.
        if let Some(path) = manifest_path {
            match app.navigator.load_manifest(path) {
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

    /// Load BAM reads and FASTA sequences for the current region.
    fn load_region_data(&mut self) {
        let entry = match self.navigator.current() {
            Some(e) => e.clone(),
            None => return,
        };

        // Sync panel views to the new region.
        self.sync_manager.set_regions(
            &entry.ref_region,
            entry.hap1_region.as_ref(),
            entry.hap2_region.as_ref(),
        );

        // -- Reference panel: BAM/CRAM reads + FASTA sequence --
        self.ref_data = PanelData::default();
        if let Some(reads_path) = &self.data_paths.reads_bam {
            let region_str = entry.ref_region.to_string();
            let result = if reads_path.extension().is_some_and(|e| e == "cram") {
                genome::bam::query_cram(
                    reads_path,
                    self.data_paths.reference_fasta.as_deref(),
                    &region_str,
                )
            } else {
                genome::bam::query_bam(reads_path, &region_str)
            };
            match result {
                Ok(reads) => {
                    // Read field summary for status: count reads, note reverse-strand and MAPQ.
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

                    self.ref_data.rows = pileup::pack_reads(reads);

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
            self.ref_data.sequence = load_fasta_for_region(fasta_path, &entry.ref_region);
        }

        // -- Haplotype 1 panel: FASTA sequence --
        self.hap1_data = PanelData::default();
        if let Some(fasta_path) = &self.data_paths.hap1_fasta
            && let Some(region) = &entry.hap1_region
        {
            self.hap1_data.sequence = load_fasta_for_region(fasta_path, region);
        }

        // -- Haplotype 2 panel: FASTA sequence --
        self.hap2_data = PanelData::default();
        if let Some(fasta_path) = &self.data_paths.hap2_fasta
            && let Some(region) = &entry.hap2_region
        {
            self.hap2_data.sequence = load_fasta_for_region(fasta_path, region);
        }

        // -- Dot plot: recompute if visible --
        if self.dot_plot.show {
            self.recompute_dot_plot();
        }
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

    // -- Toolbar ------------------------------------------------------------

    /// Render the top toolbar with navigation controls.
    fn show_toolbar(&mut self, ui: &mut egui::Ui) {
        ui.horizontal(|ui| {
            // File menu
            ui.menu_button("File", |ui| {
                if ui.button("Load Manifest…").clicked() {
                    if let Some(path) = rfd::FileDialog::new()
                        .add_filter("JSON Manifest", &["json"])
                        .pick_file()
                    {
                        match self.navigator.load_manifest(&path) {
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
                        self.load_region_data();
                    }
                    ui.close_menu();
                }
                ui.separator();
                if ui.button("Quit").clicked() {
                    ui.ctx().send_viewport_cmd(egui::ViewportCommand::Close);
                }
            });

            ui.separator();

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
            }
            if ui.button("🔍−").on_hover_text("Zoom out (-)").clicked() {
                self.sync_manager.zoom(PanelId::Reference, 0.5);
            }

            // Pan controls
            let pan_delta = self.sync_manager.view(PanelId::Reference).span() as i64 / 4;
            if ui.button("◀◀").on_hover_text("Pan left").clicked() {
                self.sync_manager.pan(PanelId::Reference, -pan_delta);
            }
            if ui.button("▶▶").on_hover_text("Pan right").clicked() {
                self.sync_manager.pan(PanelId::Reference, pan_delta);
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
        });
    }

    // -- Panel rendering ----------------------------------------------------

    /// Render a single panel with pileup reads drawn on a canvas.
    fn show_panel(ui: &mut egui::Ui, params: &PanelRenderParams<'_>) {
        let header_color = params.panel.color();

        // Panel header
        let header_rect = ui.allocate_space(egui::vec2(ui.available_width(), 24.0)).1;
        ui.painter().rect_filled(header_rect, 0.0, header_color);

        let mut header_text = format!("{}  |  {}", params.panel.label(), params.region_text);
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

        // Panel body: dark background with pileup rendering
        let body = egui::Frame::new()
            .fill(egui::Color32::from_gray(32))
            .inner_margin(egui::Margin::same(4));

        body.show(ui, |ui| {
            ui.set_min_height(80.0);

            let has_reads = !params.data.rows.is_empty();
            let has_seq = params.data.sequence.is_some();

            if !has_reads && !has_seq {
                // Placeholder when no data is loaded
                ui.centered_and_justified(|ui| {
                    ui.label(
                        egui::RichText::new(format!("{} – no data loaded", params.panel.label()))
                            .color(egui::Color32::from_gray(120))
                            .italics(),
                    );
                });
            } else {
                // Render FASTA sequence summary if available
                if let Some(seq) = &params.data.sequence {
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

                // Render pileup reads
                if has_reads {
                    let panel_width = ui.available_width();
                    let rects = pileup::layout_read_rects(
                        &params.data.rows,
                        params.config,
                        params.view_start,
                        params.view_end,
                        panel_width,
                    );
                    Self::paint_pileup(ui, &rects);
                }
            }
        });
    }

    /// Paint pre-computed read rectangles onto the UI using egui painter.
    fn paint_pileup(ui: &mut egui::Ui, rects: &[ReadRect]) {
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

            // Show coordinate mapper info if available.
            if let Some(index) = &self.coord_index {
                let results = genome::coordinate_mapper::query(
                    index,
                    &entry.ref_region.chrom,
                    entry.ref_region.start,
                    entry.ref_region.end,
                    0,
                );
                if !results.is_empty() {
                    let events: Vec<String> = results
                        .iter()
                        .map(|r| {
                            format!(
                                "{}: {}:{}-{}",
                                r.event_type, r.asm_chrom, r.asm_start, r.asm_end
                            )
                        })
                        .collect();
                    msg.push_str(&format!(" | mapped: {}", events.join(", ")));
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
        }
        if zoom_out {
            self.sync_manager.zoom(PanelId::Reference, 0.5);
        }
        if toggle_sync {
            self.sync_manager.sync_enabled = !self.sync_manager.sync_enabled;
        }
    }
}

impl eframe::App for ViewerApp {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        self.handle_keyboard(ctx);

        // Top toolbar
        egui::TopBottomPanel::top("toolbar").show(ctx, |ui| {
            self.show_toolbar(ui);
        });

        // Bottom status bar
        egui::TopBottomPanel::bottom("status_bar").show(ctx, |ui| {
            ui.horizontal(|ui| {
                ui.label(
                    egui::RichText::new(&self.status_message)
                        .small()
                        .color(egui::Color32::from_gray(180)),
                );
            });
        });

        // Clone display_config for immutable borrow inside closure
        let config = self.display_config.clone();

        // Central area: 3 panels + optional dot plot
        egui::CentralPanel::default().show(ctx, |ui| {
            let current = self.navigator.current();
            let ref_text = current
                .map(|e| e.ref_region.to_string())
                .unwrap_or_else(|| "—".to_string());
            let hap1_text = current
                .and_then(|e| e.hap1_region.as_ref())
                .map(|r| r.to_string())
                .unwrap_or_else(|| "—".to_string());
            let hap2_text = current
                .and_then(|e| e.hap2_region.as_ref())
                .map(|r| r.to_string())
                .unwrap_or_else(|| "—".to_string());

            // Use per-panel view ranges from sync manager
            let ref_view = self.sync_manager.view(PanelId::Reference);
            let h1_view = self.sync_manager.view(PanelId::Haplotype1);
            let h2_view = self.sync_manager.view(PanelId::Haplotype2);

            let (ref_start, ref_end) = (ref_view.view_start, ref_view.view_end);
            let (h1_start, h1_end) = (h1_view.view_start, h1_view.view_end);
            let (h2_start, h2_end) = (h2_view.view_start, h2_view.view_end);

            // Compute coordinate mapper info for header display
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

            // Calculate available space
            let dot_plot_height = if self.dot_plot.show { 280.0 } else { 0.0 };
            let available = ui.available_height() - dot_plot_height;
            let panel_height = (available - 16.0) / 3.0; // 16px for spacing

            ui.allocate_ui(egui::vec2(ui.available_width(), panel_height), |ui| {
                Self::show_panel(
                    ui,
                    &PanelRenderParams {
                        panel: Panel::Reference,
                        region_text: &ref_text,
                        data: &self.ref_data,
                        config: &config,
                        view_start: ref_start,
                        view_end: ref_end,
                        coord_info: ref_coord_info.as_deref(),
                    },
                );
            });

            ui.add_space(4.0);

            ui.allocate_ui(egui::vec2(ui.available_width(), panel_height), |ui| {
                Self::show_panel(
                    ui,
                    &PanelRenderParams {
                        panel: Panel::Haplotype1,
                        region_text: &hap1_text,
                        data: &self.hap1_data,
                        config: &config,
                        view_start: h1_start,
                        view_end: h1_end,
                        coord_info: None,
                    },
                );
            });

            ui.add_space(4.0);

            ui.allocate_ui(egui::vec2(ui.available_width(), panel_height), |ui| {
                Self::show_panel(
                    ui,
                    &PanelRenderParams {
                        panel: Panel::Haplotype2,
                        region_text: &hap2_text,
                        data: &self.hap2_data,
                        config: &config,
                        view_start: h2_start,
                        view_end: h2_end,
                        coord_info: None,
                    },
                );
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
    }
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

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
}
