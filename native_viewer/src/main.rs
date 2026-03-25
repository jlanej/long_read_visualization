use clap::Parser;
use log::info;

/// Native multi-panel long-read genome visualization viewer.
///
/// Displays three synchronized genome browser panels (Reference, Haplotype 1,
/// Haplotype 2) backed by indexed BAM/CRAM and FASTA files. Reads are colored
/// and sorted by haplotype phase (HP tag), with configurable indel filtering
/// and squished/expanded display modes.
#[derive(Parser, Debug)]
#[command(name = "lrv-viewer", version, about)]
struct Cli {
    /// Path to the sample configuration TSV file.
    ///
    /// Each row defines a sample with columns:
    /// sample_id, output_dir, reference, hap1_assembly, hap2_assembly,
    /// reads_bam, regions, cram_ref
    #[arg(short, long)]
    config: Option<String>,

    /// Path to reference FASTA (.fa.gz with .fai and .gzi indices).
    #[arg(long)]
    reference: Option<String>,

    /// Path to reads BAM file (.bam with .bai index).
    #[arg(long)]
    bam: Option<String>,

    /// Path to manifest JSON or VCF file with regions of interest.
    #[arg(long)]
    regions: Option<String>,

    /// Initial genomic region to display (e.g. "chr1:1000-2000").
    #[arg(long)]
    locus: Option<String>,

    /// Default indel filter threshold in bp (indels ≤ this size are hidden).
    #[arg(long, default_value = "3")]
    indel_threshold: u32,

    /// Start in expanded (non-squished) display mode.
    #[arg(long)]
    expanded: bool,
}

fn main() -> anyhow::Result<()> {
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info")).init();

    let cli = Cli::parse();
    info!("Starting lrv-viewer");

    // Load configuration
    let app_config = if let Some(config_path) = &cli.config {
        lrv_native_viewer::config::AppConfig::from_tsv(config_path)?
    } else {
        lrv_native_viewer::config::AppConfig::from_cli(
            cli.reference.as_deref(),
            cli.bam.as_deref(),
            cli.regions.as_deref(),
        )
    };

    let viewer_state = lrv_native_viewer::app::ViewerState::new(
        app_config,
        cli.locus,
        cli.indel_threshold,
        !cli.expanded,
    );

    let native_options = eframe::NativeOptions {
        viewport: egui::ViewportBuilder::default()
            .with_inner_size([1600.0, 900.0])
            .with_min_inner_size([800.0, 400.0])
            .with_title("Long Read Viewer"),
        ..Default::default()
    };

    eframe::run_native(
        "Long Read Viewer",
        native_options,
        Box::new(|cc| Ok(Box::new(lrv_native_viewer::app::ViewerApp::new(cc, viewer_state)))),
    )
    .map_err(|e| anyhow::anyhow!("eframe error: {e}"))?;

    Ok(())
}
