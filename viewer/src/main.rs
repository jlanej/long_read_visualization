mod app;
mod dot_plot;
mod genome;
mod panel_sync;
mod region;

use std::path::PathBuf;

use app::{DataPaths, ViewerApp};

fn main() -> eframe::Result<()> {
    let args = parse_args();

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
        Box::new(move |cc| {
            Ok(Box::new(ViewerApp::new(
                cc,
                args.manifest.as_deref(),
                args.data_paths,
            )))
        }),
    )
}

/// Parsed command-line arguments.
struct CliArgs {
    manifest: Option<PathBuf>,
    data_paths: DataPaths,
}

fn parse_args() -> CliArgs {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let mut manifest: Option<PathBuf> = None;
    let mut data_paths = DataPaths::default();
    let mut config_tsv: Option<PathBuf> = None;

    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "--config" => {
                i += 1;
                if i < args.len() {
                    config_tsv = Some(PathBuf::from(&args[i]));
                }
            }
            "--bam" | "--reads" => {
                i += 1;
                if i < args.len() {
                    data_paths.reads_bam = Some(PathBuf::from(&args[i]));
                }
            }
            "--reference" | "--ref" => {
                i += 1;
                if i < args.len() {
                    data_paths.reference_fasta = Some(PathBuf::from(&args[i]));
                }
            }
            "--hap1" => {
                i += 1;
                if i < args.len() {
                    data_paths.hap1_fasta = Some(PathBuf::from(&args[i]));
                }
            }
            "--hap2" => {
                i += 1;
                if i < args.len() {
                    data_paths.hap2_fasta = Some(PathBuf::from(&args[i]));
                }
            }
            "--index" => {
                i += 1;
                if i < args.len() {
                    data_paths.coordinate_index = Some(PathBuf::from(&args[i]));
                }
            }
            other => {
                // Positional arg: treat as manifest path if not set.
                if manifest.is_none() {
                    manifest = Some(PathBuf::from(other));
                }
            }
        }
        i += 1;
    }

    // If a config TSV was provided, load data paths from it.
    // Explicit CLI flags override TSV values.
    if let Some(tsv_path) = &config_tsv {
        match DataPaths::from_tsv(tsv_path) {
            Ok(tsv_paths) => {
                if data_paths.reference_fasta.is_none() {
                    data_paths.reference_fasta = tsv_paths.reference_fasta;
                }
                if data_paths.hap1_fasta.is_none() {
                    data_paths.hap1_fasta = tsv_paths.hap1_fasta;
                }
                if data_paths.hap2_fasta.is_none() {
                    data_paths.hap2_fasta = tsv_paths.hap2_fasta;
                }
                if data_paths.reads_bam.is_none() {
                    data_paths.reads_bam = tsv_paths.reads_bam;
                }
                // Use regions from config if no manifest given on CLI
                if manifest.is_none() {
                    manifest = tsv_paths.regions;
                }
            }
            Err(e) => {
                eprintln!("Warning: failed to load config TSV: {e}");
            }
        }
    }

    CliArgs {
        manifest,
        data_paths,
    }
}
