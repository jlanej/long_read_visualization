use crate::config::AppConfig;
use crate::filters::DisplayFilters;
use crate::genome::region::GenomicRegion;
use crate::pileup::PileupData;
use log::info;

/// Top-level viewer state shared across all panels.
pub struct ViewerState {
    pub config: AppConfig,
    pub current_sample_idx: usize,
    pub current_region_idx: Option<usize>,
    pub ref_region: Option<GenomicRegion>,
    pub hap1_region: Option<GenomicRegion>,
    pub hap2_region: Option<GenomicRegion>,
    pub filters: DisplayFilters,
    pub sync_enabled: bool,
    pub status_message: String,
    // Pileup data for each panel
    pub ref_pileup: Option<PileupData>,
    pub hap1_pileup: Option<PileupData>,
    pub hap2_pileup: Option<PileupData>,
}

impl ViewerState {
    pub fn new(
        config: AppConfig,
        initial_locus: Option<String>,
        indel_threshold: u32,
        squished: bool,
    ) -> Self {
        let ref_region = initial_locus.and_then(|s| GenomicRegion::parse(&s));

        Self {
            config,
            current_sample_idx: 0,
            current_region_idx: None,
            ref_region,
            hap1_region: None,
            hap2_region: None,
            filters: DisplayFilters::new(indel_threshold, squished),
            sync_enabled: true,
            status_message: "Ready".to_string(),
            ref_pileup: None,
            hap1_pileup: None,
            hap2_pileup: None,
        }
    }

    /// Navigate to the next region of interest.
    pub fn next_region(&mut self) {
        let n = self.config.region_count();
        if n == 0 {
            return;
        }
        let idx = match self.current_region_idx {
            Some(i) if i + 1 < n => i + 1,
            _ => 0,
        };
        self.navigate_to_region(idx);
    }

    /// Navigate to the previous region of interest.
    pub fn prev_region(&mut self) {
        let n = self.config.region_count();
        if n == 0 {
            return;
        }
        let idx = match self.current_region_idx {
            Some(0) | None => n.saturating_sub(1),
            Some(i) => i - 1,
        };
        self.navigate_to_region(idx);
    }

    /// Navigate to a specific region by index.
    pub fn navigate_to_region(&mut self, idx: usize) {
        if let Some(region) = self.config.get_region(idx) {
            info!("Navigating to region {idx}: {region}");
            self.current_region_idx = Some(idx);
            self.ref_region = Some(region.clone());
            self.status_message =
                format!("Region {}/{}: {region}", idx + 1, self.config.region_count());
            // Clear pileup data to trigger reload
            self.ref_pileup = None;
            self.hap1_pileup = None;
            self.hap2_pileup = None;
        }
    }

    /// Pan the current view by a fraction of the visible span.
    pub fn pan(&mut self, fraction: f64) {
        if let Some(region) = &mut self.ref_region {
            region.pan(fraction);
            // Clear pileup data to trigger reload
            self.ref_pileup = None;
            self.hap1_pileup = None;
            self.hap2_pileup = None;
        }
    }

    /// Zoom in/out by a factor (>1 = zoom out, <1 = zoom in).
    pub fn zoom(&mut self, factor: f64) {
        if let Some(region) = &mut self.ref_region {
            region.zoom(factor);
            // Clear pileup data to trigger reload
            self.ref_pileup = None;
            self.hap1_pileup = None;
            self.hap2_pileup = None;
        }
    }
}

/// Main eframe application.
pub struct ViewerApp {
    state: ViewerState,
}

impl ViewerApp {
    pub fn new(_cc: &eframe::CreationContext<'_>, state: ViewerState) -> Self {
        Self { state }
    }
}

impl eframe::App for ViewerApp {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        // Handle keyboard navigation
        handle_keyboard(ctx, &mut self.state);

        // Top toolbar
        egui::TopBottomPanel::top("toolbar").show(ctx, |ui| {
            crate::ui::toolbar::render_toolbar(ui, &mut self.state);
        });

        // Bottom status bar
        egui::TopBottomPanel::bottom("status").show(ctx, |ui| {
            ui.horizontal(|ui| {
                ui.label(&self.state.status_message);
                ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                    if let Some(region) = &self.state.ref_region {
                        ui.label(format!("{region}"));
                    }
                });
            });
        });

        // Three genome panels side by side
        egui::CentralPanel::default().show(ctx, |ui| {
            let panel_width = ui.available_width() / 3.0 - 8.0;
            ui.horizontal(|ui| {
                // Reference panel
                ui.vertical(|ui| {
                    ui.set_width(panel_width);
                    crate::ui::panel::render_panel(
                        ui,
                        "Reference",
                        &self.state.ref_region,
                        &self.state.ref_pileup,
                        &self.state.filters,
                    );
                });

                ui.separator();

                // Haplotype 1 panel
                ui.vertical(|ui| {
                    ui.set_width(panel_width);
                    crate::ui::panel::render_panel(
                        ui,
                        "Haplotype 1",
                        &self.state.hap1_region,
                        &self.state.hap1_pileup,
                        &self.state.filters,
                    );
                });

                ui.separator();

                // Haplotype 2 panel
                ui.vertical(|ui| {
                    ui.set_width(panel_width);
                    crate::ui::panel::render_panel(
                        ui,
                        "Haplotype 2",
                        &self.state.hap2_region,
                        &self.state.hap2_pileup,
                        &self.state.filters,
                    );
                });
            });
        });
    }
}

/// Process keyboard shortcuts for navigation.
fn handle_keyboard(ctx: &egui::Context, state: &mut ViewerState) {
    ctx.input(|input| {
        // Region stepping
        if input.key_pressed(egui::Key::N)
            || input.key_pressed(egui::Key::ArrowRight)
                && input.modifiers.shift
        {
            state.next_region();
        }
        if input.key_pressed(egui::Key::P)
            || input.key_pressed(egui::Key::ArrowLeft)
                && input.modifiers.shift
        {
            state.prev_region();
        }

        // Panning (arrow keys)
        if input.key_pressed(egui::Key::ArrowRight) && !input.modifiers.shift {
            state.pan(0.25);
        }
        if input.key_pressed(egui::Key::ArrowLeft) && !input.modifiers.shift {
            state.pan(-0.25);
        }

        // Zooming
        if input.key_pressed(egui::Key::Minus)
            || input.key_pressed(egui::Key::OpenBracket)
        {
            state.zoom(2.0); // zoom out
        }
        if input.key_pressed(egui::Key::Plus)
            || input.key_pressed(egui::Key::Equals)
            || input.key_pressed(egui::Key::CloseBracket)
        {
            state.zoom(0.5); // zoom in
        }

        // Toggle squished/expanded
        if input.key_pressed(egui::Key::S) {
            state.filters.squished = !state.filters.squished;
        }

        // Toggle indel filter
        if input.key_pressed(egui::Key::I) {
            state.filters.hide_small_indels = !state.filters.hide_small_indels;
        }
    });
}
