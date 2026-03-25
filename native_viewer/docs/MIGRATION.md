# Migration Guide: IGV.js → Native Rust Viewer

This document describes how to use existing pipeline outputs (from the Python
`preprocess.sh` pipeline) with the new native Rust viewer.

## Config TSV format

The native viewer accepts the **same TSV config file** used by the Python
server (`server/app.py`). No changes are required.

### TSV columns

| Column | Description | Example |
|--------|-------------|---------|
| `sample_id` | Unique sample identifier | `NA21110` |
| `output_dir` | Pipeline output directory | `/work/output/NA21110` |
| `reference` | Reference FASTA (.fa.gz) | `/references/chm13v2.0.fa.gz` |
| `hap1_assembly` | Haplotype 1 assembly FASTA | `/assemblies/NA21110_hap1.fa.gz` |
| `hap2_assembly` | Haplotype 2 assembly FASTA | `/assemblies/NA21110_hap2.fa.gz` |
| `reads_bam` | Reads BAM/CRAM | `/data/NA21110_reads.bam` |
| `regions` | Manifest JSON or VCF | `/work/output/NA21110/manifest.json` |
| `cram_ref` | CRAM encoding reference | `/references/encoding_ref.fa.gz` |

### Example

```tsv
#sample_id	output_dir	reference	hap1_assembly	hap2_assembly	reads_bam	regions	cram_ref
NA21110	/work/output/NA21110	/refs/chm13v2.0.fa.gz	/asm/NA21110_hap1.fa.gz	/asm/NA21110_hap2.fa.gz	/data/reads.bam	/work/output/NA21110/manifest.json	
```

## Pipeline output file discovery

The native viewer auto-discovers these pipeline outputs from `output_dir`:

| Pattern | Purpose |
|---------|---------|
| `*_hap1_to_ref.bam` | Assembly-to-reference alignment (hap1) |
| `*_hap2_to_ref.bam` | Assembly-to-reference alignment (hap2) |
| `*_reads_to_hap1.bam` / `*.cram` | Reads aligned to hap1 assembly |
| `*_reads_to_hap2.bam` / `*.cram` | Reads aligned to hap2 assembly |
| `*_hap1_to_ref.mapping.json.gz` | Coordinate mapping index (hap1) |
| `*_hap2_to_ref.mapping.json.gz` | Coordinate mapping index (hap2) |

## Manifest JSON format

The viewer reads the same `manifest.json` format produced by the pipeline:

```json
{
  "description": "Sample structural variants",
  "variants": [
    {
      "chrom": "chr1",
      "pos": 112164095,
      "size": 12912,
      "genotype": "0|1",
      "ref_region": "chr1:112164095-112177007",
      "fasta_region": "chr1:112064095-112277008",
      "hap1_regions": ["NA21110#1#CM089663.1:111967298-112167366"],
      "hap2_regions": ["NA21110#2#CM089678.1:111935960-112149228"]
    }
  ]
}
```

## Mapping index format

Coordinate translation uses `.mapping.json.gz` files produced by the Python
pipeline. These are gzip-compressed JSON files mapping reference chromosomes
to alignment blocks:

```json
{
  "chr1": [
    {
      "rs": 1000,    // ref_start
      "re": 2000,    // ref_end
      "ac": "ctg1",  // asm_chrom
      "as": 5000,    // asm_start
      "ae": 6000,    // asm_end
      "st": "+",     // strand
      "mq": 60,      // mapq
      "cg": "1000M", // CIGAR (optional)
      "tp": "P"      // alignment type (optional)
    }
  ]
}
```

## Feature mapping

| IGV.js Feature | Native Viewer | Notes |
|---------------|---------------|-------|
| Three synchronized panels | ✅ Reference + Hap1 + Hap2 | Same layout |
| Squished/expanded toggle | ✅ `S` key or toolbar | Same behavior |
| Small indel hiding | ✅ `I` key or toolbar | Default threshold: 3bp |
| HP tag coloring | ✅ Blue (HP1) / Red (HP2) / Grey | Matches IGV.js palette |
| Region navigation | ✅ `N`/`P` keys or toolbar | Prev/Next region |
| Pan/zoom | ✅ Arrow keys, `+`/`-` | Keyboard and toolbar |
| BAM byte-range HTTP | ❌ Local files only | No HTTP serving needed |
| Dot plots | ❌ Not yet implemented | Future enhancement |
| VCF variant tracks | ❌ Not yet implemented | Future enhancement |

## Key differences from IGV.js version

1. **No web server required** — The native viewer reads files directly from
   disk; no HTTP server or browser needed.

2. **Performance** — Direct file I/O via noodles is more efficient than HTTP
   byte-range requests. The egui immediate-mode rendering avoids DOM overhead.

3. **Coordinate synchronization** — The native viewer uses the same CIGAR-aware
   projection algorithm as the Python server, but runs it in-process rather
   than via HTTP API calls.

4. **No container dependency** — The native viewer is a single static binary
   that runs on any Linux or macOS system with a display.
