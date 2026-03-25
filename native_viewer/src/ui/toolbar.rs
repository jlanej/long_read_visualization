//! Top toolbar with navigation controls and display toggles.

use egui::Ui;

use crate::app::ViewerState;

/// Render the toolbar at the top of the application.
pub fn render_toolbar(ui: &mut Ui, state: &mut ViewerState) {
    ui.horizontal(|ui| {
        ui.spacing_mut().item_spacing.x = 8.0;

        // Application title
        ui.heading("Long Read Viewer");
        ui.separator();

        // Sample selector (if multiple samples)
        if state.config.samples.len() > 1 {
            let current_name = state
                .config
                .samples
                .get(state.current_sample_idx)
                .map(|s| s.sample_id.as_str())
                .unwrap_or("(none)");

            egui::ComboBox::from_label("Sample")
                .selected_text(current_name)
                .show_ui(ui, |ui| {
                    for (i, sample) in state.config.samples.iter().enumerate() {
                        if ui
                            .selectable_label(
                                i == state.current_sample_idx,
                                &sample.sample_id,
                            )
                            .clicked()
                        {
                            state.current_sample_idx = i;
                            state.ref_pileup = None;
                            state.hap1_pileup = None;
                            state.hap2_pileup = None;
                        }
                    }
                });
            ui.separator();
        }

        // Region navigation
        let region_count = state.config.region_count();
        if region_count > 0 {
            if ui.button("◀ Prev").clicked() {
                state.prev_region();
            }

            let region_label = match state.current_region_idx {
                Some(idx) => format!("{}/{}", idx + 1, region_count),
                None => format!("—/{}", region_count),
            };
            ui.label(region_label);

            if ui.button("Next ▶").clicked() {
                state.next_region();
            }
            ui.separator();
        }

        // Pan/zoom controls
        if ui.button("← Pan").clicked() {
            state.pan(-0.25);
        }
        if ui.button("Pan →").clicked() {
            state.pan(0.25);
        }
        if ui.button("➖ Zoom Out").clicked() {
            state.zoom(2.0);
        }
        if ui.button("➕ Zoom In").clicked() {
            state.zoom(0.5);
        }
        ui.separator();

        // Display toggles
        let squish_label = if state.filters.squished {
            "📐 Squished"
        } else {
            "📏 Expanded"
        };
        if ui.button(squish_label).clicked() {
            state.filters.squished = !state.filters.squished;
        }

        let indel_label = if state.filters.hide_small_indels {
            format!("🔍 Indels ≤{}bp hidden", state.filters.indel_threshold)
        } else {
            "🔍 All indels shown".to_string()
        };
        if ui.button(indel_label).clicked() {
            state.filters.hide_small_indels = !state.filters.hide_small_indels;
        }

        // Sync toggle
        let sync_label = if state.sync_enabled {
            "🔗 Sync ON"
        } else {
            "🔗 Sync OFF"
        };
        if ui.button(sync_label).clicked() {
            state.sync_enabled = !state.sync_enabled;
        }

        // Right-aligned region input
        ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
            ui.label("Keyboard: ←/→ pan, +/- zoom, N/P regions, S squish, I indels");
        });
    });
}
