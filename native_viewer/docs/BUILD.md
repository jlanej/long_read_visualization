# Building and Running the Native Viewer

## Prerequisites

- **Rust** 1.88 or later (install via [rustup](https://rustup.rs/))
- **System libraries** (Linux): OpenGL, X11, and Wayland development libraries

### Linux (Ubuntu/Debian)

```bash
sudo apt-get install -y \
  libgl1-mesa-dev libx11-dev libxcb1-dev libxrandr-dev \
  libxcursor-dev libxi-dev libxkbcommon-dev libwayland-dev
```

### macOS

No additional system libraries needed — Xcode command-line tools provide
everything.

```bash
xcode-select --install
```

## Building

From the repository root:

```bash
cd native_viewer
cargo build --release
```

The compiled binary will be at `native_viewer/target/release/lrv-viewer`.

## Running

### With a config file (TSV)

Use the same TSV configuration as the Python server:

```bash
./target/release/lrv-viewer --config /path/to/config.tsv
```

### With individual file arguments

```bash
./target/release/lrv-viewer \
  --reference /path/to/reference.fa.gz \
  --bam /path/to/reads.bam \
  --regions /path/to/manifest.json
```

### With the toy dataset

```bash
./target/release/lrv-viewer \
  --reference ../resources/toy_dataset/toy_reference.fa.gz \
  --bam ../resources/toy_dataset/toy_reads.bam \
  --regions ../resources/toy_dataset/toy_manifest.json
```

### Starting at a specific locus

```bash
./target/release/lrv-viewer \
  --config config.tsv \
  --locus "chr1:112164095-112177007"
```

## Command-line options

| Flag | Description | Default |
|------|-------------|---------|
| `--config <TSV>` | Sample configuration TSV (same format as Python server) | — |
| `--reference <FASTA>` | Reference FASTA (.fa.gz with .fai and .gzi) | — |
| `--bam <BAM>` | Reads BAM (.bam with .bai) | — |
| `--regions <JSON/VCF>` | Manifest JSON or VCF file | — |
| `--locus <REGION>` | Initial region, e.g. `chr1:1000-2000` | — |
| `--indel-threshold <N>` | Hide indels ≤ N bp | 3 |
| `--expanded` | Start in expanded (non-squished) mode | false (squished) |

## Keyboard shortcuts

| Key | Action |
|-----|--------|
| `←` / `→` | Pan left / right |
| `Shift+←` / `Shift+→` | Previous / next region |
| `N` / `P` | Next / previous region |
| `+` / `=` / `]` | Zoom in |
| `-` / `[` | Zoom out |
| `S` | Toggle squished / expanded |
| `I` | Toggle indel filtering |

## Running tests

```bash
cd native_viewer
cargo test
```

## Architecture

```
native_viewer/
├── Cargo.toml               # Dependencies and project metadata
├── src/
│   ├── main.rs               # CLI entry point (clap)
│   ├── lib.rs                 # Library root
│   ├── app.rs                 # eframe::App — viewer state & update loop
│   ├── config.rs              # Config loading (TSV, manifest JSON)
│   ├── filters.rs             # Display filters (indel threshold, squish)
│   ├── genome/
│   │   ├── bam.rs             # BAM reading via noodles (indexed queries)
│   │   ├── fasta.rs           # FASTA reading via noodles (indexed)
│   │   └── region.rs          # GenomicRegion type (parse, pan, zoom)
│   ├── pileup/
│   │   ├── mod.rs             # Pileup data types
│   │   └── layout.rs          # Read packing (greedy first-fit)
│   ├── coordinate/
│   │   └── translator.rs      # Ref↔hap coordinate translation (CIGAR-aware)
│   └── ui/
│       ├── colors.rs           # HP tag coloring, base colors
│       ├── panel.rs            # Genome panel rendering (pileup + axis)
│       └── toolbar.rs          # Navigation toolbar
├── tests/                      # Integration tests
└── docs/
    ├── BUILD.md                # This file
    └── MIGRATION.md            # Migration from Python/IGV.js config
```
