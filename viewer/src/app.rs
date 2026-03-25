use eframe::egui;

use crate::region::RegionNavigator;

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
}

/// Main application state.
pub struct ViewerApp {
    pub navigator: RegionNavigator,
    pub status_message: String,
}

impl Default for ViewerApp {
    fn default() -> Self {
        Self {
            navigator: RegionNavigator::new(),
            status_message: "No regions loaded. Use File > Load Manifest to open a region manifest JSON."
                .to_string(),
        }
    }
}

impl ViewerApp {
    pub fn new(cc: &eframe::CreationContext<'_>) -> Self {
        configure_fonts(&cc.egui_ctx);
        Self::default()
    }

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
                            }
                            Err(e) => {
                                self.status_message = format!("Error: {e}");
                            }
                        }
                    }
                    ui.close_menu();
                }
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
                self.update_status_for_current_region();
            }

            ui.label(self.navigator.counter_text());
        });
    }

    /// Render a single panel placeholder.
    fn show_panel(ui: &mut egui::Ui, panel: Panel, region_text: &str) {
        let header_color = panel.color();

        // Panel header
        let header_rect = ui.allocate_space(egui::vec2(ui.available_width(), 24.0)).1;
        ui.painter().rect_filled(header_rect, 0.0, header_color);
        ui.painter().text(
            header_rect.left_center() + egui::vec2(8.0, 0.0),
            egui::Align2::LEFT_CENTER,
            format!("{}  |  {}", panel.label(), region_text),
            egui::FontId::proportional(13.0),
            egui::Color32::WHITE,
        );

        // Placeholder body
        let body = egui::Frame::new()
            .fill(egui::Color32::from_gray(32))
            .inner_margin(egui::Margin::same(12))
            .stroke(egui::Stroke::new(1.0, egui::Color32::from_gray(60)));

        body.show(ui, |ui| {
            ui.set_min_height(80.0);
            ui.centered_and_justified(|ui| {
                ui.label(
                    egui::RichText::new(format!("{} – pileup placeholder", panel.label()))
                        .color(egui::Color32::from_gray(120))
                        .italics(),
                );
            });
        });
    }

    /// Update status message to reflect the current region.
    fn update_status_for_current_region(&mut self) {
        if let Some(entry) = self.navigator.current() {
            self.status_message = format!(
                "Region {}: {} → ref {}",
                self.navigator.counter_text(),
                entry.label,
                entry.ref_region
            );
        }
    }

    /// Handle keyboard shortcuts for region navigation.
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

        if prev && !self.navigator.is_empty() {
            self.navigator.prev();
            self.update_status_for_current_region();
        }
        if next && !self.navigator.is_empty() {
            self.navigator.next();
            self.update_status_for_current_region();
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

        // Central area: 3 panels stacked vertically
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

            // Use vertical layout with equal panel sizes
            let available = ui.available_height();
            let panel_height = (available - 16.0) / 3.0; // 16px for spacing

            ui.allocate_ui(egui::vec2(ui.available_width(), panel_height), |ui| {
                Self::show_panel(ui, Panel::Reference, &ref_text);
            });

            ui.add_space(4.0);

            ui.allocate_ui(egui::vec2(ui.available_width(), panel_height), |ui| {
                Self::show_panel(ui, Panel::Haplotype1, &hap1_text);
            });

            ui.add_space(4.0);

            ui.allocate_ui(egui::vec2(ui.available_width(), panel_height), |ui| {
                Self::show_panel(ui, Panel::Haplotype2, &hap2_text);
            });
        });
    }
}

/// Configure default fonts/styles.
fn configure_fonts(ctx: &egui::Context) {
    let mut style = (*ctx.style()).clone();
    style.text_styles.insert(
        egui::TextStyle::Body,
        egui::FontId::proportional(14.0),
    );
    style.text_styles.insert(
        egui::TextStyle::Button,
        egui::FontId::proportional(13.0),
    );
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
    fn test_viewer_app_default() {
        let app = ViewerApp::default();
        assert!(app.navigator.is_empty());
        assert!(app.status_message.contains("No regions loaded"));
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
        app.update_status_for_current_region();
        assert!(app.status_message.contains("1/1"));
        assert!(app.status_message.contains("chr1:100-150"));
    }
}
