# Toy Dataset Generation

This document describes how the minimal toy dataset used for integration
testing is created.  The toy dataset is a small, self-contained set of files
that exercises the full `preprocess.sh` pipeline without requiring large
genomic inputs.

---

## Overview

The toy dataset is derived from real data for sample **NA21110** using the
[HPRC](https://humanpangenome.org/) assemblies and the T2T-CHM13 v2.0
reference.  It contains only the genomic regions surrounding a handful of
large structural-variant deletions, making it small enough to commit to
the repository and fast enough to run in CI.

### Data sources and citations

The curated training set is sourced from the following publicly available
resources.  **Please cite them if you use this dataset in your work.**

| Source | Description | Reference |
|---|---|---|
| **shapeit5-phased callset** | Phased, sequence-resolved SV VCF used as the primary truth set; included in `resources/` | Ebler *et al.* (2022). Pangenome-based genome inference. *Nat Genet* — and the 1KG ONT Vienna consortium (see below) |
| **1KG ONT Vienna** | Population-scale long-read SV calls from 1,019 samples across the 1000 Genomes Project, sequenced with Oxford Nanopore | Liao *et al.* (2025). Structural variation in 1,019 diverse humans based on long-read sequencing. *Nature* [https://doi.org/10.1038/s41586-025-09290-7](https://doi.org/10.1038/s41586-025-09290-7). Data: [EBI FTP](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1KG_ONT_VIENNA/) |
| **PMC12350158** | Curated structural variant truth sets used to validate the callset | [https://www.ncbi.nlm.nih.gov/pmc/articles/PMC12350158/](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC12350158/) |

### Source files

| File | Description |
|---|---|
| `resources/NA21110.shapeit5-phased-callset_final-vcf.phased.vcf.gz` | Phased SV callset for NA21110 — shapeit5-phased callset (bundled in repo; see citation above) |
| Full NA21110 haplotype assemblies (hap1 & hap2) | HPRC r2 v1.0.1, sourced from the 1KG ONT Vienna project |
| Full NA21110 long-read CRAM | ONT reads from the 1KG ONT Vienna project, aligned to T2T-CHM13 v2.0 |
| T2T-CHM13 v2.0 reference genome | Downloaded via `--genome chm13v2.0` |

### Generated toy files

| File | Description |
|---|---|
| `toy_reference.fa.gz` | Reference FASTA — only the selected regions |
| `toy_hap1.fa.gz` | Haplotype 1 assembly — only the mapped regions |
| `toy_hap2.fa.gz` | Haplotype 2 assembly — only the mapped regions |
| `toy_reads.bam` | Reads extracted from the CRAM for the selected regions |
| `toy_manifest.json` | JSON manifest listing every variant and region |

---

## Step-by-step procedure

### Prerequisites

The generation script requires a **completed pipeline run**.  Specifically,
it needs the coordinate-mapping JSON indices
(`*_hap{1,2}_to_ref.mapping.json.gz`) produced by `preprocess.sh`, because
these indices translate reference coordinates to assembly coordinates.

Run the full pipeline first:

```bash
apptainer run \
    --bind "${PWD}:/work" \
    docker://ghcr.io/jlanej/long_read_visualization:main \
    --hap1       /work/NA21110/NA21110_hap1_hprc_r2_v1.0.1.fa.gz \
    --hap2       /work/NA21110/NA21110_hap2_hprc_r2_v1.0.1.fa.gz \
    --cram       /work/NA21110/NA21110.t2t.cram \
    --genome      chm13v2.0 \
    --output-dir  /work/output/NA21110 \
    --threads     32 \
    --ont
```

### Step 1 — Select large deletions from the VCF

The script parses the NA21110 SV callset VCF and identifies deletions where:

1. `len(REF) - len(ALT) >= --min-size` (default 1000 bp).
2. The NA21110 genotype contains at least one alternate allele (`1|0`,
   `0|1`, or `1|1`).

Deletions are sorted by size (largest first).  The script then selects
`--num-variants` (default 10) deletions, preferring diversity across
chromosomes: it picks the single largest deletion per chromosome first,
then fills remaining slots with the next-largest deletions regardless of
chromosome.

### Step 2 — Compute padded reference regions

For each selected deletion at `chrom:pos` with reference allele length
`ref_len`, a padded region is computed:

```
start = max(0, pos - padding)
end   = pos + ref_len + padding
```

The default `--padding` is 50000 bp, giving roughly 100 kb of context
around each deletion.

### Step 3 — Map reference regions to assembly coordinates

The coordinate-mapping JSON indices built by the pipeline are loaded.
For each padded reference region the script queries both hap1 and hap2
indices to find the corresponding assembly contig(s) and coordinate
ranges.

Overlapping assembly hits on the same contig are merged into
non-overlapping intervals.

### Step 4 — Extract reference regions

Reference subsequences are extracted with `samtools faidx`, then
bgzip-compressed and re-indexed:

```bash
samtools faidx reference.fa.gz chr1:1000-2000 chr3:5000-6000 ... > toy_reference.fa
bgzip toy_reference.fa
samtools faidx toy_reference.fa.gz
```

### Step 5 — Extract assembly regions

The same approach is used for each haplotype assembly.  The assembly
regions identified in Step 3 are extracted with `samtools faidx`:

```bash
samtools faidx hap1.fa.gz contig1:100-5000 contig2:200-8000 ... > toy_hap1.fa
bgzip toy_hap1.fa
samtools faidx toy_hap1.fa.gz
```

### Step 6 — Extract reads

Reads overlapping the padded reference regions are extracted from the
CRAM, sorted, and indexed into a BAM:

```bash
samtools view -b -h --reference ref.fa cram chr1:1000-2000 ... \
    | samtools sort -o toy_reads.bam
samtools index toy_reads.bam
```

### Step 7 — Write manifest

A JSON manifest (`toy_manifest.json`) is written listing each variant's
chromosome, position, size, genotype, padded reference region, and the
corresponding hap1/hap2 assembly regions.  This manifest can be used by
tests to verify that the correct regions were extracted and that the
pipeline produces the expected outputs.

---

## Running the generation script

### With Apptainer (recommended)

```bash
apptainer run \
    --bind "${PWD}:/work" \
    docker://ghcr.io/jlanej/long_read_visualization:main \
    python3 /opt/long_read_visualization/src/generate_toy_dataset.py \
    --vcf /opt/long_read_visualization/resources/NA21110.shapeit5-phased-callset_final-vcf.phased.vcf.gz \
    --hap1-index /work/output/NA21110/NA21110_hap1_to_ref.mapping.json.gz \
    --hap2-index /work/output/NA21110/NA21110_hap2_to_ref.mapping.json.gz \
    --hap1       /work/NA21110/NA21110_hap1_hprc_r2_v1.0.1.fa.gz \
    --hap2       /work/NA21110/NA21110_hap2_hprc_r2_v1.0.1.fa.gz \
    --reference  /work/references/chm13v2.0.fa.gz \
    --cram       /work/NA21110/NA21110.t2t.cram \
    --output-dir /work/toy_dataset \
    --num-variants 10 \
    --padding 50000
```

### Using the shell wrapper (Apptainer)

```bash
cd /path/to/projects

apptainer run \
    --bind "${PWD}:/work" \
    --entrypoint generate_toy_dataset.sh \
    docker://ghcr.io/jlanej/long_read_visualization:main \
    --pipeline-output /work/output/NA21110 \
    --hap1       /work/NA21110/NA21110_hap1_hprc_r2_v1.0.1.fa.gz \
    --hap2       /work/NA21110/NA21110_hap2_hprc_r2_v1.0.1.fa.gz \
    --reference  /work/references/chm13v2.0.fa.gz \
    --cram       /work/NA21110/NA21110.t2t.cram \
    --output-dir /work/toy_dataset
```

---

## Running the pipeline on the toy dataset

Once generated, the toy dataset can be used as input to the main pipeline:

```bash
apptainer run \
    --bind "${PWD}:/work" \
    docker://ghcr.io/jlanej/long_read_visualization:main \
    --hap1       /work/toy_dataset/toy_hap1.fa.gz \
    --hap2       /work/toy_dataset/toy_hap2.fa.gz \
    --reference  /work/toy_dataset/toy_reference.fa.gz \
    --cram       /work/toy_dataset/toy_reads.bam \
    --output-dir /work/toy_dataset/pipeline_output \
    --threads     4 \
    --ont
```

---

## Customisation

| Flag | Default | Description |
|---|---|---|
| `--num-variants` | 10 | Number of deletions to include |
| `--min-size` | 1000 | Minimum deletion size (bp) |
| `--padding` | 50000 | Context bases on each side of a deletion |

Larger `--padding` values produce a bigger dataset but provide more
flanking context for alignment and visualization.
