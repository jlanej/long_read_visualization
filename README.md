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
4. **Aligns the reference genome back to each haplotype assembly** using
   `minimap2 -x asm5 --eqx` (query and target swapped relative to item 1),
   producing assembly-coordinate-sorted BAM files that can be loaded as a
   reference track in any assembly-space genome browser.  These reciprocal
   alignments serve as a cross-check against the hap→ref alignments: every
   block visible in the hap→ref BAM should have a complementary block in
   the ref→hap BAM, making spurious or missed mappings immediately apparent.

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
| `*_hap{1,2}_to_ref.bam(.bai)` | Assembly aligned to reference (reference-coordinate BAM) |
| `*_hap{1,2}_to_ref.paf` | PAF alignment for coordinate mapping |
| `*_hap{1,2}_to_ref.mapping.bed.gz(.tbi)` | Tabix-indexed coordinate map |
| `*_hap{1,2}_to_ref.mapping.json.gz` | JSON coordinate map for programmatic use |
| `*_reads.fastq.gz` | Reads extracted from the input CRAM |
| `*_reads_to_hap{1,2}.bam(.bai)` | Reads aligned to each haplotype assembly |
| `*_ref_to_hap{1,2}.bam(.bai)` | Reference genome aligned to each haplotype assembly (assembly-coordinate BAM, for cross-checking and assembly-panel display) |

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

## Reference-to-assembly cross-check BAMs

### What they are

In addition to the assembly-to-reference BAMs produced in Step 2,
the pipeline generates two **reference-to-assembly** BAMs (Step 6):

```
*_ref_to_hap1.bam(.bai)
*_ref_to_hap2.bam(.bai)
```

These files contain the same alignment information as the hap→ref BAMs
(Step 2) but with query and target swapped: the *reference genome* is the
query and each *haplotype assembly* is the target.  As a result the BAM is
sorted in **assembly coordinate space**, making it directly loadable as a
track in any assembly-panel genome browser (e.g., an IGV.js panel whose
reference sequence is `hap1.fa`).

### Why they are useful

| Use case | Detail |
|---|---|
| **Assembly-panel reference track** | Load `*_ref_to_hap1.bam` in an IGV.js panel whose reference is `hap1.fa`. The reference sequence appears as an aligned read track, giving instant visual context for what the reference looks like at any assembly locus. |
| **Reciprocal alignment cross-check** | Every alignment block in the hap→ref BAM should have a complementary block in the ref→hap BAM. Discordant blocks (present in one direction but not the other) flag potentially spurious mappings or missed alignments in low-complexity / segmental-duplication regions. |
| **SV breakpoint orientation** | Viewing both directions simultaneously in linked panels lets you visually confirm whether a structural variant is supported by both the assembly orientation and the reference orientation, reducing false-positive calls. |

### Alignment parameters

| Parameter | Value | Rationale |
|---|---|---|
| `minimap2 -x asm5` | asm5 preset | Optimised for ≥ 99% identity — the expected divergence between a high-quality human assembly (e.g., HPRC) and a population reference. Same preset used for the hap→ref BAMs for consistency. |
| `--eqx` | Extended CIGAR (`=` / `X`) | Explicit match/mismatch encoding; enables per-base mismatch visualisation in IGV.js and simplifies downstream CIGAR parsing. |
| `samtools sort` | By coordinate | BAM is coordinate-sorted in assembly space so that any assembly-coordinate range can be fetched with a standard random-access query. |

### Idempotency

Step 6 checks for the existence of both the BAM file **and** its index
before running.  Re-running the pipeline after the BAMs are already present
skips the alignment and indexing without overwriting or corrupting the
existing files:

```
[...] ref_to_hap1.bam + index exist, skipping
[...] ref_to_hap2.bam + index exist, skipping
```

### Loading in IGV.js

To display the reference-on-assembly track alongside reads:

```json
{
  "reference": { "fastaURL": "hap1.fa", "indexURL": "hap1.fa.fai" },
  "tracks": [
    { "name": "Reads → Hap1",     "url": "sample_reads_to_hap1.bam",  "type": "alignment" },
    { "name": "Reference → Hap1", "url": "sample_ref_to_hap1.bam",    "type": "alignment" }
  ]
}
```

---

## Coordinate mapping

### How it works

Each haplotype assembly is aligned to the target reference genome in two
separate minimap2 invocations:

1. **BAM alignment** (`minimap2 -a --eqx -x asm5`): produces a sorted,
   indexed BAM for genome-browser display.
2. **PAF alignment** (`minimap2 --eqx -c -x asm5`): produces a PAF file
   with `cg:Z:` CIGAR tags used for coordinate translation.  Passing `-c`
   causes minimap2 to emit the full CIGAR string in the PAF `cg:Z:` tag,
   capturing every match, mismatch, insertion, and deletion at base-level
   resolution.

Those alignments are parsed into a sorted interval index (compact JSON +
tabix-indexed BED).  For any reference coordinate range a binary search
locates every overlapping alignment record and **walks the CIGAR string** to
project the query interval onto assembly space.  CIGAR-aware projection gives
correct coordinates even when the block contains insertions or deletions —
critical for structural-variant regions.  Blocks without a CIGAR tag fall
back to linear interpolation for compatibility with older PAF files.

### Structural-variant awareness

Because SV breakpoints typically interrupt alignment blocks, a reference query
that spans a breakpoint will overlap *two or more* blocks separated by a gap.
The mapper detects these gaps and classifies them:

| `event_type` | Condition | Meaning |
|---|---|---|
| `alignment` | within an alignment block | Reference bases have direct assembly equivalents |
| `deletion` | same contig, same strand; assembly gap < reference gap | Assembly deleted sequence relative to reference |
| `insertion` | same contig, same strand; assembly gap > reference gap | Assembly has novel sequence not in reference |
| `inversion` | same contig, opposite strands at junction | Orientation flip — likely an inversion breakpoint |
| `translocation` | different assembly contigs at junction | Inter-contig event |
| `complex` | same contig, same strand; assembly gap is negative (overlap) | Duplication or other complex rearrangement |

Every result dict returned by `query()` carries an `event_type` field.
Callers can filter to `"alignment"` records for strict coordinate lookups or
inspect all records to understand the full SV landscape of a region.

### Why it works

`asm5` is minimap2's preset for sequences that are **≥99% identical** — the
expected divergence between a high-quality human assembly (e.g., HPRC) and a
population reference.  Within well-assembled, collinear regions the mapping is
unambiguous and coordinates translate accurately even across small indels and
SNPs.  At SV boundaries, alignment blocks break at the structural rearrangement
and the gap-classification logic exposes what kind of event is present, making
those regions *visible* rather than silently absent.

### Known limitations and how they are addressed

| Situation | Status | Detail |
|---|---|---|
| **Large SVs (deletions, insertions)** | ✅ Addressed | CIGAR-aware projection handles indels within blocks; gaps between blocks are classified as `deletion` / `insertion` with accurate size estimates |
| **Inversions** | ✅ Addressed | Opposite-strand blocks at a junction are detected and reported as `inversion` events |
| **Translocations** | ✅ Addressed | Blocks from different assembly contigs at a junction are reported as `translocation` events |
| **Complex rearrangements** | ✅ Addressed | Overlapping assembly ranges (negative asm gap) are classified as `complex` events |
| **Segmental duplications / CNVs** | ⚠️ Annotated | The same reference region may map to multiple assembly loci; all hits are returned so the caller sees the full picture.  Use `min_mapq` to prefer primary alignments |
| **Supplementary / chimeric alignments** | ⚠️ Handled | Multiple partial alignments are retained in the index; overlapping blocks all appear in results.  Pass `min_mapq` (e.g. `--min-mapq 5`) to exclude low-confidence supplementary hits |
| **PAF without `-c` flag** | ⚠️ Fallback | No CIGAR → coordinate projection falls back to linear interpolation, which is less accurate across indels.  The pipeline uses `-c` by default |
| **Unmapped / highly-divergent regions** | ❌ Inherent | No alignment → query returns no result.  Centromeres, acrocentric arms, and novel insertions without flanking alignment are invisible |
| **Assembly gaps or low-quality sequence** | ❌ Inherent | Alignments may be clipped short, leaving bases near the gap untranslatable |

### Quality filtering

Low-confidence alignments (supplementary hits with `mapq=0`, multi-mapped
blocks in segmental duplications) can be excluded with the `min_mapq`
parameter:

```bash
# CLI: exclude supplementary alignments with mapq < 5
python3 src/coordinate_mapper.py query \
    -i output_prefix.mapping.json.gz \
    -r chr1:1000000-2000000 \
    --min-mapq 5

# Python API
results = coordinate_mapper.query(index, "chr1", 1000000, 2000000, min_mapq=5)
```

When blocks are filtered by quality, gap events are computed between the
*remaining* blocks only — ensuring that gap classifications reflect the
high-confidence alignment landscape.

### Querying coordinate mappings

```bash
# Using the Python tool (loads JSON index into memory)
python3 src/coordinate_mapper.py query \
    -i output/NA21110_hap1_to_ref.mapping.json.gz \
    -r chr1:1000000-2000000

# Using tabix (standard bioinformatics tool)
tabix output/NA21110_hap1_to_ref.mapping.bed.gz chr1:1000000-2000000
```

The TSV output includes an `event_type` column: `alignment` for regions with
direct assembly equivalents, or `deletion` / `insertion` / `inversion` /
`translocation` / `complex` for gap events at SV breakpoints.

### Coordinate mapper CLI

```bash
# Build index from PAF (use -c to include CIGAR for precise SV mapping)
minimap2 -x asm5 -c ref.fa asm.fa > asm_to_ref.paf

python3 src/coordinate_mapper.py build \
    -p asm_to_ref.paf \
    -o output_prefix

# Query index
python3 src/coordinate_mapper.py query \
    -i output_prefix.mapping.json.gz \
    -r chr1:1000000-2000000
```

The JSON index uses sorted intervals with binary search for O(log n) overlap
queries.  When a `cg:Z:` CIGAR is present, coordinate projection walks the
CIGAR for base-exact translation.  Gaps between alignment blocks are
classified by SV type and returned alongside alignment results.

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

apptainer exec \
    --bind "${PWD}:/work" \
    docker://ghcr.io/jlanej/long_read_visualization:main \
    generate_toy_dataset.sh \
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
