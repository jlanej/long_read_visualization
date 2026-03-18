# long_read_visualization

Pre-processing pipeline for multi-panel long-read visualization. Given
haplotype assemblies and long-read alignments, produces all files needed for
linked IGV.js browsing across a reference genome and both haplotypes.

## Overview

The pipeline takes three inputs for a sample:

| Input | Description |
|---|---|
| **Hap1 / Hap2 assemblies** | FASTA (`.fa` or `.fa.gz`) with index (`.fai`) |
| **Long-read CRAM** | ONT or PacBio HiFi reads aligned to any reference |
| **Target reference genome** | FASTA, or a named genome auto-downloaded with `--genome` |

It then:

1. **Aligns each haplotype assembly to the target reference** using
   `minimap2 -x asm5` (BAM + PAF output).
2. **Builds coordinate-mapping indices** from the PAF alignments — a
   tabix-indexed BED and a compact JSON index that support fast
   *reference → assembly* coordinate lookups.
3. **Extracts reads** from the CRAM and **re-aligns them to each haplotype
   assembly** using `minimap2` (`-x map-ont` or `-x map-hifi`), producing
   sorted, indexed BAM files.

### Example input files

Files can be in any location — pass their paths directly to `--hap1`,
`--hap2`, and `--cram`. A common layout when using HPRC assemblies is:

```
NA21110/
  NA21110.t2t.cram
  NA21110.t2t.cram.crai
  NA21110_hap1_hprc_r2_v1.0.1.fa.gz
  NA21110_hap1_hprc_r2_v1.0.1.fa.gz.fai
  NA21110_hap2_hprc_r2_v1.0.1.fa.gz
  NA21110_hap2_hprc_r2_v1.0.1.fa.gz.fai
```

The **sample name** used in output filenames is derived from the CRAM
filename (everything before the first `.`): `NA21110.t2t.cram` → `NA21110`.

### Pipeline outputs

| File | Description |
|---|---|
| `*_hap{1,2}_to_ref.bam(.bai)` | Assembly-to-reference alignment |
| `*_hap{1,2}_to_ref.paf` | PAF alignment for coordinate mapping |
| `*_hap{1,2}_to_ref.mapping.bed.gz(.tbi)` | Tabix-indexed coordinate map |
| `*_hap{1,2}_to_ref.mapping.json.gz` | JSON coordinate map for programmatic use |
| `*_reads.fastq.gz` | Reads extracted from the input CRAM |
| `*_reads_to_hap{1,2}.bam(.bai)` | Reads aligned to each haplotype |

---

## Quick start

### With Apptainer (recommended for HPC)

Run from the **parent directory** that contains the sample folder. The
`--bind "${PWD}:/work"` flag makes the current directory available inside
the container at `/work`. Downloaded references are cached inside the
bind-mounted working directory at `/work/references/` (i.e. `${PWD}/references/`
on the host) so they persist across runs.

**ONT reads:**

```bash
# cd to the parent directory of your sample folder first
cd /path/to/projects   # e.g. this directory contains NA21110/

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

**PacBio HiFi reads:**

```bash
cd /path/to/projects

apptainer run \
    --bind "${PWD}:/work" \
    docker://ghcr.io/jlanej/long_read_visualization:main \
    --hap1       /work/NA21110/NA21110_hap1_hprc_r2_v1.0.1.fa.gz \
    --hap2       /work/NA21110/NA21110_hap2_hprc_r2_v1.0.1.fa.gz \
    --cram       /work/NA21110/NA21110.hifi.cram \
    --genome      chm13v2.0 \
    --output-dir  /work/output/NA21110 \
    --threads     32 \
    --hifi
```

The `--genome chm13v2.0` flag downloads the reference automatically on first
run and caches it at `${PWD}/references/` (or `/work/references/` inside the
container). Because this directory is inside the bind mount, the cache
persists across container runs — no extra mount is needed.

---

## Command-line options

```
Usage: preprocess.sh [options]

Required inputs:
  --hap1        FILE   Haplotype 1 assembly FASTA (.fa or .fa.gz)
  --hap2        FILE   Haplotype 2 assembly FASTA (.fa or .fa.gz)
  --cram        FILE   Long-read CRAM file
  -o, --output-dir DIR Output directory

Read type (controls minimap2 alignment preset):
  --ont                Oxford Nanopore reads  (minimap2 -x map-ont)  [default]
  --hifi               PacBio HiFi reads      (minimap2 -x map-hifi)

Reference (one of):
  -r, --reference FILE Target reference genome FASTA (local file)
  --genome        STR  Download a known reference if not already cached.
                       Supported: hg38 hg19 chm13v2.0 grch38
  --ref-dir       DIR  Cache directory for downloaded references
                       [./references]

Optional:
  -t, --threads   INT  Number of threads [4]
  --cram-ref      FILE Reference FASTA used to encode the CRAM
                       (needed only when different from --reference / --genome)
```

### Supported `--genome` values

| Name | Source |
|---|---|
| `chm13v2.0` | T2T-CHM13 v2.0 (human-pangenomics S3) |
| `hg38` | UCSC hg38 |
| `hg19` | UCSC hg19 |
| `grch38` | NCBI GRCh38 no-alt analysis set |

---

## Coordinate mapping

### How it works

Each haplotype assembly is aligned to the target reference genome with
`minimap2 -x asm5`, which produces both a BAM file (for browsing) and a PAF
file (for coordinate translation). The PAF records exact base-level
correspondences between reference and assembly positions via extended CIGAR
strings.

Those alignments are parsed into a sorted interval index (compact JSON +
tabix-indexed BED). For any reference coordinate range, a binary search
locates every overlapping alignment record and recomputes the corresponding
assembly position by walking the CIGAR — accounting for insertions, deletions,
and reverse-complement (inversion) orientations. This gives O(log n) lookups
over the full genome.

### Why it works

`asm5` is minimap2's preset for sequences that are **≥99% identical** — the
expected divergence between a high-quality human assembly (e.g., HPRC) and a
population reference. Within well-assembled, collinear regions the mapping is
unambiguous and coordinates translate accurately even across small indels and
SNPs.

### Known limitations

| Situation | Effect |
|---|---|
| Unmapped / highly-divergent regions (centromeres, segmental duplications, novel insertions) | No alignment → coordinate query returns no result |
| Large structural variants (inversions, translocations) | Alignment breaks at SV boundaries; coordinates spanning a breakpoint cannot be translated as a single interval |
| Segmental duplications / copy-number variants | The same reference region may map to multiple assembly loci, producing ambiguous results |
| Assembly gaps or low-quality sequence | Alignments may be clipped short, leaving bases near the gap untranslatable |
| Supplementary / chimeric alignments | Multiple partial alignments for one contig can overlap in reference space; the index retains all of them, so queries in those regions may return multiple hits |

These limitations are inherent to any alignment-based liftover approach. For
the primary use-case — visualising long reads across well-assembled diploid
genomes at specific variant sites — the vast majority of query coordinates fall
in high-confidence, uniquely-mappable regions where the method works reliably.

### Querying coordinate mappings

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

## Container image

The image is built and published automatically via GitHub Actions on pushes to
`main` and on version tags. Pull it with Apptainer using the
`docker://ghcr.io/jlanej/long_read_visualization:main` URI shown in the
[Quick start](#quick-start) examples above.

---

## Toy dataset for integration testing

A minimal toy dataset can be generated from the bundled NA21110 SV callset.
The generation script selects large deletions from the VCF, uses the
coordinate-mapping indices to find the corresponding assembly regions, and
extracts small subsets of the reference, assemblies, and reads.

**Note:** this requires a completed pipeline run (it uses the mapping indices).

### Data sources

The curated training data is sourced from the following resources.  Please
cite them if you use this dataset in your work:

- **PMC12350158** — curated structural variant truth sets:
  <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC12350158/>
- **1KG ONT Vienna** — population-scale long-read SV calls from 1,019 samples
  across the 1000 Genomes Project (Oxford Nanopore, aligned to T2T-CHM13 v2.0):
  Liao *et al.* (2025). *Nature* <https://doi.org/10.1038/s41586-025-09290-7>.
  Data: <https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1KG_ONT_VIENNA/>
- **shapeit5-phased callset** — phased, sequence-resolved SV VCF used as the
  primary training truth set; included in `resources/`.

```bash
cd /path/to/projects   # same parent directory used for the pipeline run

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

See [docs/toy_dataset.md](docs/toy_dataset.md) for the full procedure.

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
