// Allow dead code: many modules expose APIs that are tested but not yet wired
// into the GUI entry-point.  This will be removed as features are connected.
#![allow(dead_code, unused_imports)]

mod app;
mod dot_plot;
mod genome;
mod panel_sync;
mod region;

use std::path::PathBuf;

use app::ViewerApp;

fn main() -> eframe::Result<()> {
    // Accept an optional manifest path as the first CLI argument.
    let manifest_path: Option<PathBuf> = std::env::args().nth(1).map(PathBuf::from);

    let options = eframe::NativeOptions {
        viewport: eframe::egui::ViewportBuilder::default()
            .with_inner_size([1200.0, 800.0])
            .with_min_inner_size([600.0, 400.0])
            .with_title("Long Read Viewer"),
        ..Default::default()
    };

    eframe::run_native(
        "Long Read Viewer",
        options,
        Box::new(move |cc| Ok(Box::new(ViewerApp::new(cc, manifest_path.as_deref())))),
    )
}
