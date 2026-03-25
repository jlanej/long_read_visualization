# Native Viewer (`lrv-viewer`)

A native Rust GUI for multi-panel long-read genome visualization, built with
[egui](https://github.com/emilk/egui) (via eframe) and
[noodles](https://github.com/zaeleus/noodles).

This is a companion to the existing IGV.js-based web viewer — it replicates the
core three-panel visualization (Reference, Haplotype 1, Haplotype 2) while
eliminating the need for a web server and browser.

## Quick start

```bash
cd native_viewer
cargo build --release
./target/release/lrv-viewer \
  --reference ../resources/toy_dataset/toy_reference.fa.gz \
  --bam ../resources/toy_dataset/toy_reads.bam \
  --regions ../resources/toy_dataset/toy_manifest.json
```

## Features

- **Three synchronized genome panels** (Reference, Hap1, Hap2)
- **Indexed BAM/CRAM reading** via noodles with random access
- **Haplotype coloring** — reads colored by HP tag (blue/red/grey)
- **Indel filtering** — small indels (≤3bp) hidden by default
- **Squished/expanded** display mode toggle
- **Keyboard navigation** — arrow keys to pan, +/− to zoom, N/P for regions
- **CIGAR-aware coordinate translation** — base-level ref↔haplotype mapping
- **Same config format** — reads the existing pipeline TSV and manifest JSON

## Documentation

- [Building and Running](docs/BUILD.md)
- [Migration Guide](docs/MIGRATION.md) — how to use existing pipeline outputs

## Tests

```bash
cargo test
```

49 unit tests covering coordinate translation, pileup layout, genomic region
operations, BAM reading, config parsing, display filters, and UI colors.
