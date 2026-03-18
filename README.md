# long_read_visualization

Pre-processing pipeline for multi-panel long-read visualization. Given
haplotype assemblies and long-read alignments, produces all files needed for
linked IGV.js browsing across a reference genome and both haplotypes.

## Overview

The pipeline takes three inputs for a sample:

| Input | Description |
|---|---|
| **Hap1 / Hap2 assemblies** | FASTA (`.fa.gz`) with index (`.fa.gz.fai`) |
| **Long-read CRAM** | ONT or PacBio HiFi reads aligned to any reference |
| **Target reference genome** | FASTA the assemblies will be aligned to |

It then:

1. **Aligns each haplotype assembly to the target reference** using
   `minimap2 -x asm5` (BAM + PAF output).
2. **Builds coordinate-mapping indices** from the PAF alignments — a
   tabix-indexed BED and a compact JSON index that support fast
   *reference → assembly* coordinate lookups.
3. **Extracts reads** from the CRAM and **re-aligns them to each haplotype
   assembly** using `minimap2` (`-x map-ont` or `-x map-hifi`), producing
   sorted, indexed BAM files.

### Expected sample directory layout

```
NA21110/
  NA21110.t2t.cram
  NA21110.t2t.cram.crai
  NA21110_hap1_hprc_r2_v1.0.1.fa.gz
  NA21110_hap1_hprc_r2_v1.0.1.fa.gz.fai
  NA21110_hap2_hprc_r2_v1.0.1.fa.gz
  NA21110_hap2_hprc_r2_v1.0.1.fa.gz.fai
```

### Pipeline outputs

| File | Description |
|---|---|
| `*_hap{1,2}_to_ref.bam(.bai)` | Assembly-to-reference alignment |
| `*_hap{1,2}_to_ref.paf` | PAF alignment for coordinate mapping |
| `*_hap{1,2}_to_ref.mapping.bed.gz(.tbi)` | Tabix-indexed coordinate map |
| `*_hap{1,2}_to_ref.mapping.json.gz` | JSON coordinate map for programmatic use |
| `*_reads_to_hap{1,2}.bam(.bai)` | Reads aligned to each haplotype |

---

## Quick start

### With Docker

```bash
docker run --rm -v /data:/data ghcr.io/jlanej/long_read_visualization:main \
    -s /data/NA21110 \
    -r /data/ref/chm13v2.0.fa \
    -o /data/output/NA21110 \
    -t 16
```

### With Apptainer / Singularity (HPC)

```bash
apptainer run \
    docker://ghcr.io/jlanej/long_read_visualization:main \
    -s /data/NA21110 \
    -r /data/ref/chm13v2.0.fa \
    -o /data/output/NA21110 \
    -t 32 \
    --read-type ont
```

### Locally (requires minimap2, samtools, htslib, python3)

```bash
bash scripts/preprocess.sh \
    -s /data/NA21110 \
    -r /data/ref/chm13v2.0.fa \
    -o /data/output/NA21110 \
    -t 16 \
    --read-type ont
```

---

## Command-line options

```
Usage: preprocess.sh [options]

Required:
  -s, --sample-dir DIR    Sample directory with CRAM + assembly FASTAs
  -r, --reference  FILE   Target reference genome FASTA
  -o, --output-dir DIR    Output directory

Optional:
  -t, --threads    INT    Number of threads  [4]
  --read-type      STR    ont | hifi         [ont]
  --cram-ref       FILE   Reference FASTA used to encode the CRAM
                          (needed only when different from --reference)
```

---

## Coordinate mapping

After pre-processing, you can query the mapping to translate reference
coordinates to assembly coordinates:

```bash
# Using the Python tool (loads JSON index into memory)
python3 src/coordinate_mapper.py query \
    -i output/NA21110_hap1_to_ref.mapping.json.gz \
    -r chr1:1000000-2000000

# Using tabix (standard bioinformatics tool)
tabix output/NA21110_hap1_to_ref.mapping.bed.gz chr1:1000000-2000000
```

### Coordinate mapper CLI

```bash
# Build index from PAF
python3 src/coordinate_mapper.py build \
    -p hap1_to_ref.paf \
    -o output_prefix

# Query index
python3 src/coordinate_mapper.py query \
    -i output_prefix.mapping.json.gz \
    -r chr1:1000000-2000000
```

The JSON index uses sorted intervals with binary search for O(log n) overlap
queries. Each query returns assembly coordinates corresponding to the
requested reference region, correctly handling both forward and reverse strand
alignments.

---

## Building the Docker image

```bash
docker build -t long_read_visualization .
```

The image is also built and published automatically via GitHub Actions on
pushes to `main` and on version tags.

---

## Running tests

```bash
python3 -m unittest discover -s tests -v
```

---

## Tools and versions

| Tool | Version | Purpose |
|---|---|---|
| minimap2 | 2.28 | Assembly and read alignment |
| samtools | 1.21 | BAM/CRAM manipulation |
| htslib | 1.21 | bgzip / tabix indexing |
| Python 3 | 3.10+ | Coordinate mapping tool |
