mod app;
mod region;

use app::ViewerApp;

fn main() -> eframe::Result<()> {
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
        Box::new(|cc| Ok(Box::new(ViewerApp::new(cc)))),
    )
}
