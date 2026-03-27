mod app;
mod data_loader;
mod dot_plot;
mod genome;
mod panel_sync;
mod region;
mod region_cache;
mod ruler;
mod verbose;

use std::path::PathBuf;

use app::{DataPaths, ViewerApp, parse_all_samples};
use verbose::LogBuffer;

fn main() -> eframe::Result<()> {
    let args = parse_args();

    // Apply verbose/quiet setting before anything else.
    verbose::set_verbose(args.verbose);

    let log_buf = LogBuffer::new();

    if args.verbose {
        log_buf.log("Verbose logging enabled (use --quiet to suppress)");
    }

    let options = eframe::NativeOptions {
        viewport: eframe::egui::ViewportBuilder::default()
            .with_inner_size([1200.0, 800.0])
            .with_min_inner_size([600.0, 400.0])
            .with_title("Long Read Viewer"),
        ..Default::default()
    };

    let result = eframe::run_native(
        "Long Read Viewer",
        options,
        Box::new(move |cc| {
            Ok(Box::new(ViewerApp::new(
                cc,
                args.manifest.as_deref(),
                args.data_paths,
                args.samples,
                log_buf.clone(),
            )))
        }),
    );

    if let Err(ref e) = result {
        let msg = e.to_string();
        if msg.contains("NoGlutinConfigs") || msg.contains("NotFound") {
            eprintln!("Error: {e}\n");
            eprintln!("The viewer could not find a suitable OpenGL/EGL configuration.");
            eprintln!("This usually means one or more of the following:");
            eprintln!("  1. Required GPU/OpenGL libraries are not installed.");
            eprintln!("     On Ubuntu/Debian, install: libegl1 libgl1 libgles2 libgl1-mesa-dri");
            eprintln!("  2. No display server is available (DISPLAY or WAYLAND_DISPLAY not set).");
            eprintln!("     For Apptainer: apptainer exec --env DISPLAY=$DISPLAY ...");
            eprintln!("     For SSH:       ssh -X user@host  (enable X11 forwarding)");
            eprintln!("  3. Software rendering is not enabled.");
            eprintln!("     Set:  export LIBGL_ALWAYS_SOFTWARE=1");
            eprintln!("  4. Running inside a container without GPU passthrough or X11 forwarding.");
        }
    }

    result
}

/// Parsed command-line arguments.
struct CliArgs {
    manifest: Option<PathBuf>,
    data_paths: DataPaths,
    samples: Vec<(String, DataPaths)>,
    /// Verbose logging (default: true, disabled by --quiet).
    verbose: bool,
}

fn parse_args() -> CliArgs {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let mut manifest: Option<PathBuf> = None;
    let mut data_paths = DataPaths::default();
    let mut config_tsv: Option<PathBuf> = None;
    let mut verbose = true;

    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "--quiet" | "-q" => {
                verbose = false;
            }
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
                merge_tsv_defaults(&mut data_paths, &mut manifest, tsv_paths);
            }
            Err(e) => {
                eprintln!("Warning: failed to load config TSV: {e}");
            }
        }
    }

    // Parse all samples from TSV if available.
    let samples = config_tsv
        .as_ref()
        .and_then(|path| parse_all_samples(path).ok())
        .unwrap_or_default();

    CliArgs {
        manifest,
        data_paths,
        samples,
        verbose,
    }
}

fn merge_tsv_defaults(
    data_paths: &mut DataPaths,
    manifest: &mut Option<PathBuf>,
    tsv_paths: DataPaths,
) {
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
    // Discovered preprocessing output (always from TSV output_dir)
    if data_paths.hap1_coord_index.is_none() {
        data_paths.hap1_coord_index = tsv_paths.hap1_coord_index;
    }
    if data_paths.hap2_coord_index.is_none() {
        data_paths.hap2_coord_index = tsv_paths.hap2_coord_index;
    }
    if data_paths.reads_to_hap1_bam.is_none() {
        data_paths.reads_to_hap1_bam = tsv_paths.reads_to_hap1_bam;
    }
    if data_paths.reads_to_hap2_bam.is_none() {
        data_paths.reads_to_hap2_bam = tsv_paths.reads_to_hap2_bam;
    }
    if data_paths.cram_ref.is_none() {
        data_paths.cram_ref = tsv_paths.cram_ref;
    }
    if data_paths.hap1_to_ref_bam.is_none() {
        data_paths.hap1_to_ref_bam = tsv_paths.hap1_to_ref_bam;
    }
    if data_paths.hap2_to_ref_bam.is_none() {
        data_paths.hap2_to_ref_bam = tsv_paths.hap2_to_ref_bam;
    }
    if data_paths.ref_to_hap1_bam.is_none() {
        data_paths.ref_to_hap1_bam = tsv_paths.ref_to_hap1_bam;
    }
    if data_paths.ref_to_hap2_bam.is_none() {
        data_paths.ref_to_hap2_bam = tsv_paths.ref_to_hap2_bam;
    }
    if data_paths.hap1_to_hap2_bam.is_none() {
        data_paths.hap1_to_hap2_bam = tsv_paths.hap1_to_hap2_bam;
    }
    if data_paths.hap2_to_hap1_bam.is_none() {
        data_paths.hap2_to_hap1_bam = tsv_paths.hap2_to_hap1_bam;
    }
    if data_paths.output_dir.is_none() {
        data_paths.output_dir = tsv_paths.output_dir;
    }
    if data_paths.sample_id.is_none() {
        data_paths.sample_id = tsv_paths.sample_id;
    }
    // Use regions from config if no manifest given on CLI
    if manifest.is_none() {
        *manifest = tsv_paths.regions;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_merge_tsv_defaults_includes_discovery_context_and_cross_bams() {
        let mut data_paths = DataPaths::default();
        let mut manifest = None;

        let tsv_paths = DataPaths {
            output_dir: Some(PathBuf::from("/work/output/NA21110")),
            sample_id: Some("NA21110".to_string()),
            hap1_to_ref_bam: Some(PathBuf::from(
                "/work/output/NA21110/NA21110_hap1_to_ref.bam",
            )),
            hap2_to_ref_bam: Some(PathBuf::from(
                "/work/output/NA21110/NA21110_hap2_to_ref.bam",
            )),
            ref_to_hap1_bam: Some(PathBuf::from(
                "/work/output/NA21110/NA21110_ref_to_hap1.bam",
            )),
            ref_to_hap2_bam: Some(PathBuf::from(
                "/work/output/NA21110/NA21110_ref_to_hap2.bam",
            )),
            hap1_to_hap2_bam: Some(PathBuf::from(
                "/work/output/NA21110/NA21110_hap1_to_hap2.bam",
            )),
            hap2_to_hap1_bam: Some(PathBuf::from(
                "/work/output/NA21110/NA21110_hap2_to_hap1.bam",
            )),
            ..DataPaths::default()
        };

        merge_tsv_defaults(&mut data_paths, &mut manifest, tsv_paths);

        assert_eq!(
            data_paths.output_dir.as_deref(),
            Some(std::path::Path::new("/work/output/NA21110"))
        );
        assert_eq!(data_paths.sample_id.as_deref(), Some("NA21110"));
        assert!(data_paths.hap1_to_ref_bam.is_some());
        assert!(data_paths.hap2_to_ref_bam.is_some());
        assert!(data_paths.ref_to_hap1_bam.is_some());
        assert!(data_paths.ref_to_hap2_bam.is_some());
        assert!(data_paths.hap1_to_hap2_bam.is_some());
        assert!(data_paths.hap2_to_hap1_bam.is_some());
    }
}
