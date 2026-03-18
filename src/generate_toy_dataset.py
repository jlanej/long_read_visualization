#!/usr/bin/env python3
"""
generate_toy_dataset.py - Generate a minimal toy dataset for integration testing.

Selects large deletions from an SV callset VCF, identifies corresponding
assembly regions via the coordinate-mapping index produced by the pipeline,
and extracts small subsets of the reference, assemblies, and reads.

This script is intended to be run AFTER a full pipeline run, because it
requires the coordinate-mapping JSON indices to translate reference
coordinates into assembly coordinates.

Usage:
    python3 generate_toy_dataset.py \
        --vcf resources/NA21110.shapeit5-phased-callset_final-vcf.phased.vcf.gz \
        --hap1-index output/NA21110_hap1_to_ref.mapping.json.gz \
        --hap2-index output/NA21110_hap2_to_ref.mapping.json.gz \
        --hap1 NA21110/NA21110_hap1_hprc_r2_v1.0.1.fa.gz \
        --hap2 NA21110/NA21110_hap2_hprc_r2_v1.0.1.fa.gz \
        --reference ref/chm13v2.0.fa.gz \
        --cram NA21110/NA21110.t2t.cram \
        --output-dir toy_dataset \
        --num-variants 10 \
        --padding 50000
"""

import argparse
import gzip
import json
import os
import shlex
import subprocess
import sys

# Allow importing coordinate_mapper from the same src/ directory.
sys.path.insert(0, os.path.dirname(__file__))
import coordinate_mapper  # noqa: E402


# ---------------------------------------------------------------------------
# VCF helpers
# ---------------------------------------------------------------------------

def parse_vcf_deletions(vcf_path, min_size=1000):
    """Parse a VCF and return deletions carried by the sample.

    A deletion is identified when len(REF) - len(ALT) >= *min_size* and
    the sample genotype contains at least one alternate allele.

    Args:
        vcf_path: Path to bgzipped or plain-text VCF.
        min_size: Minimum deletion size in bp.

    Returns:
        List of dicts with chrom, pos, end, size, genotype for each
        qualifying deletion, sorted by size descending.
    """
    opener = gzip.open if vcf_path.endswith(".gz") else open
    deletions = []
    with opener(vcf_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if len(fields) < 10:
                continue

            ref_len = len(fields[3])
            alt_len = len(fields[4])
            if ref_len <= alt_len:
                continue

            size = ref_len - alt_len
            if size < min_size:
                continue

            gt = fields[9].split(":")[0]
            if "1" not in gt:
                continue

            chrom = fields[0]
            pos = int(fields[1])
            deletions.append(
                {
                    "chrom": chrom,
                    "pos": pos,
                    "end": pos + ref_len,
                    "size": size,
                    "genotype": gt,
                }
            )

    deletions.sort(key=lambda d: -d["size"])
    return deletions


def select_variants(deletions, num_variants=10):
    """Select *num_variants* deletions spread across different chromosomes.

    Prefers the largest deletions while ensuring chromosome diversity.

    Args:
        deletions: List of deletion dicts (sorted by size desc).
        num_variants: Number of deletions to select.

    Returns:
        List of selected deletion dicts.
    """
    if len(deletions) <= num_variants:
        return deletions

    selected = []
    seen_chroms = set()

    # First pass: pick the largest deletion per chromosome.
    for d in deletions:
        if d["chrom"] not in seen_chroms:
            selected.append(d)
            seen_chroms.add(d["chrom"])
            if len(selected) >= num_variants:
                break

    # Second pass: fill remaining slots with the largest unselected deletions.
    if len(selected) < num_variants:
        selected_set = {(d["chrom"], d["pos"]) for d in selected}
        for d in deletions:
            if (d["chrom"], d["pos"]) not in selected_set:
                selected.append(d)
                if len(selected) >= num_variants:
                    break

    selected.sort(key=lambda d: (d["chrom"], d["pos"]))
    return selected


# ---------------------------------------------------------------------------
# Region helpers
# ---------------------------------------------------------------------------

def compute_regions(variants, hap1_index, hap2_index, padding=50000):
    """Compute reference and assembly regions for each selected variant.

    For each variant a padded reference region is created.  The coordinate
    mapping indices are then queried to find the corresponding assembly
    regions on each haplotype.

    Args:
        variants: List of selected variant dicts.
        hap1_index: Loaded coordinate-mapping index for haplotype 1.
        hap2_index: Loaded coordinate-mapping index for haplotype 2.
        padding: Bases to add on each side of the deletion.

    Returns:
        List of dicts, each containing:
          - ref_region: "chrom:start-end"
          - hap1_regions: list of "contig:start-end" strings
          - hap2_regions: list of "contig:start-end" strings
          - variant: the original variant dict
    """
    regions = []
    for v in variants:
        ref_start = max(0, v["pos"] - padding)
        ref_end = v["end"] + padding
        ref_region = f"{v['chrom']}:{ref_start}-{ref_end}"

        hap1_hits = coordinate_mapper.query(
            hap1_index, v["chrom"], ref_start, ref_end
        )
        hap2_hits = coordinate_mapper.query(
            hap2_index, v["chrom"], ref_start, ref_end
        )

        hap1_regions = _collapse_asm_regions(hap1_hits)
        hap2_regions = _collapse_asm_regions(hap2_hits)

        regions.append(
            {
                "ref_region": ref_region,
                "hap1_regions": hap1_regions,
                "hap2_regions": hap2_regions,
                "variant": v,
            }
        )
    return regions


def _collapse_asm_regions(hits):
    """Merge assembly hits per contig into non-overlapping intervals.

    Only ``"alignment"`` events are considered; structural-variant gap events
    (deletions, insertions, inversions, translocations) are skipped so that
    the returned regions correspond strictly to aligned assembly sequence.

    Returns a list of "contig:start-end" strings.
    """
    by_contig = {}
    for h in hits:
        if h.get("event_type", "alignment") != "alignment":
            continue
        contig = h["asm_chrom"]
        s, e = h["asm_start"], h["asm_end"]
        if s > e:
            s, e = e, s
        by_contig.setdefault(contig, []).append((s, e))

    regions = []
    for contig in sorted(by_contig):
        intervals = sorted(by_contig[contig])
        merged = [intervals[0]]
        for s, e in intervals[1:]:
            if s <= merged[-1][1]:
                merged[-1] = (merged[-1][0], max(merged[-1][1], e))
            else:
                merged.append((s, e))
        for s, e in merged:
            regions.append(f"{contig}:{s}-{e}")
    return regions


# ---------------------------------------------------------------------------
# Extraction helpers (call out to samtools)
# ---------------------------------------------------------------------------

def _run(cmd, description=""):
    """Run a shell command, raising on failure.

    Always prints the exact command being executed so the invocation is
    fully visible in the log.  If *description* is also provided it is
    printed on a separate line before the command.
    """
    cmd_str = " ".join(shlex.quote(c) for c in cmd)
    if description:
        print(f"  {description}", file=sys.stderr)
    print(f"  CMD: {cmd_str}", file=sys.stderr)
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"ERROR: {cmd_str}\n{result.stderr}", file=sys.stderr)
        sys.exit(1)
    return result


def extract_fasta_regions(input_fasta, regions, output_fasta):
    """Extract regions from a FASTA using ``samtools faidx``.

    Args:
        input_fasta: Path to the indexed input FASTA.
        regions: List of "chrom:start-end" strings.
        output_fasta: Path for the output FASTA.
    """
    if not regions:
        return
    cmd = ["samtools", "faidx", input_fasta] + regions
    result = _run(cmd, f"Extracting {len(regions)} region(s) → {output_fasta}")
    with open(output_fasta, "w") as fh:
        fh.write(result.stdout)


def extract_reads(cram_path, ref_regions, output_bam, reference=None,
                  cram_ref=None):
    """Extract reads overlapping *ref_regions* from a CRAM/BAM.

    Produces a sorted, indexed BAM file containing only the reads that
    overlap the requested reference regions.

    Args:
        cram_path: Path to the input CRAM/BAM.
        ref_regions: List of "chrom:start-end" strings.
        output_bam: Output BAM path.
        reference: Reference FASTA (for CRAM decoding, optional).
        cram_ref: Explicit CRAM reference (overrides *reference*).
    """
    if not ref_regions:
        return

    ref_opt = []
    if cram_ref:
        ref_opt = ["--reference", cram_ref]
    elif reference:
        ref_opt = ["--reference", reference]

    cmd = ["samtools", "view", "-b", "-h"] + ref_opt + [cram_path] + ref_regions
    unsorted = output_bam + ".unsorted.bam"
    cmd_str = " ".join(shlex.quote(c) for c in cmd) + " > " + shlex.quote(unsorted)
    print(f"  Extracting reads for {len(ref_regions)} region(s)", file=sys.stderr)
    print(f"  CMD: {cmd_str}", file=sys.stderr)
    with open(unsorted, "wb") as fh:
        result = subprocess.run(cmd, stdout=fh, stderr=subprocess.PIPE)
        if result.returncode != 0:
            print(f"ERROR: {' '.join(cmd)}\n{result.stderr.decode()}", file=sys.stderr)
            sys.exit(1)

    _run(
        ["samtools", "sort", "-o", output_bam, unsorted],
        f"Sorting → {output_bam}",
    )
    _run(["samtools", "index", output_bam], f"Indexing {output_bam}")
    os.remove(unsorted)


def index_and_compress_fasta(fasta_path):
    """bgzip-compress and index a FASTA file.

    Args:
        fasta_path: Path to the uncompressed FASTA.

    Returns:
        Path to the compressed FASTA (.fa.gz).
    """
    gz_path = fasta_path + ".gz"
    _run(["bgzip", "-f", fasta_path], f"Compressing {fasta_path}")
    _run(["samtools", "faidx", gz_path], f"Indexing {gz_path}")
    return gz_path


# ---------------------------------------------------------------------------
# Manifest / summary
# ---------------------------------------------------------------------------

def write_manifest(regions, output_dir):
    """Write a JSON manifest describing the toy dataset.

    Args:
        regions: List of region dicts from compute_regions().
        output_dir: Output directory.
    """
    manifest = {
        "description": "Toy dataset for long_read_visualization integration testing",
        "variants": [],
    }
    for r in regions:
        manifest["variants"].append(
            {
                "chrom": r["variant"]["chrom"],
                "pos": r["variant"]["pos"],
                "size": r["variant"]["size"],
                "genotype": r["variant"]["genotype"],
                "ref_region": r["ref_region"],
                "hap1_regions": r["hap1_regions"],
                "hap2_regions": r["hap2_regions"],
            }
        )

    manifest_path = os.path.join(output_dir, "toy_manifest.json")
    with open(manifest_path, "w") as fh:
        json.dump(manifest, fh, indent=2)
    print(f"Manifest written to {manifest_path}", file=sys.stderr)
    return manifest_path


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Generate a toy dataset for integration testing"
    )
    parser.add_argument(
        "--vcf",
        required=True,
        help="SV callset VCF (bgzipped or plain text)",
    )
    parser.add_argument(
        "--hap1-index",
        required=True,
        help="Hap1 coordinate-mapping JSON index (.mapping.json.gz)",
    )
    parser.add_argument(
        "--hap2-index",
        required=True,
        help="Hap2 coordinate-mapping JSON index (.mapping.json.gz)",
    )
    parser.add_argument(
        "--hap1",
        required=True,
        help="Haplotype 1 assembly FASTA (.fa.gz)",
    )
    parser.add_argument(
        "--hap2",
        required=True,
        help="Haplotype 2 assembly FASTA (.fa.gz)",
    )
    parser.add_argument(
        "--reference",
        required=True,
        help="Reference genome FASTA",
    )
    parser.add_argument(
        "--cram",
        required=True,
        help="Long-read CRAM file",
    )
    parser.add_argument(
        "--cram-ref",
        default=None,
        help="Reference FASTA used to encode the CRAM (if different from --reference)",
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        required=True,
        help="Output directory for the toy dataset",
    )
    parser.add_argument(
        "-n",
        "--num-variants",
        type=int,
        default=10,
        help="Number of large deletions to include [10]",
    )
    parser.add_argument(
        "--min-size",
        type=int,
        default=1000,
        help="Minimum deletion size in bp [1000]",
    )
    parser.add_argument(
        "--padding",
        type=int,
        default=50000,
        help="Padding around each deletion in bp [50000]",
    )

    args = parser.parse_args()

    # ── Step 1: Select variants ─────────────────────────────────────────
    print("Step 1: Selecting large deletions from VCF", file=sys.stderr)
    deletions = parse_vcf_deletions(args.vcf, min_size=args.min_size)
    print(
        f"  Found {len(deletions)} deletions >= {args.min_size} bp",
        file=sys.stderr,
    )
    if len(deletions) == 0:
        print("ERROR: no qualifying deletions found in VCF", file=sys.stderr)
        sys.exit(1)

    selected = select_variants(deletions, num_variants=args.num_variants)
    print(
        f"  Selected {len(selected)} variant(s) across "
        f"{len({v['chrom'] for v in selected})} chromosome(s)",
        file=sys.stderr,
    )

    # ── Step 2: Load mapping indices and compute regions ────────────────
    print("Step 2: Computing assembly regions via coordinate mapping",
          file=sys.stderr)
    hap1_index = coordinate_mapper.load_index(args.hap1_index)
    hap2_index = coordinate_mapper.load_index(args.hap2_index)
    regions = compute_regions(selected, hap1_index, hap2_index,
                              padding=args.padding)

    # ── Step 3: Create output directory and extract reference regions ───
    os.makedirs(args.output_dir, exist_ok=True)
    print("Step 3: Extracting reference regions", file=sys.stderr)
    ref_regions = [r["ref_region"] for r in regions]
    toy_ref = os.path.join(args.output_dir, "toy_reference.fa")
    extract_fasta_regions(args.reference, ref_regions, toy_ref)
    index_and_compress_fasta(toy_ref)

    # ── Step 4: Extract assembly regions ────────────────────────────────
    print("Step 4: Extracting haplotype assembly regions", file=sys.stderr)
    hap1_asm_regions = []
    hap2_asm_regions = []
    for r in regions:
        hap1_asm_regions.extend(r["hap1_regions"])
        hap2_asm_regions.extend(r["hap2_regions"])

    # Deduplicate assembly regions — when multiple variants map to overlapping
    # assembly regions, _collapse_asm_regions returns the same region string
    # for each variant.  Passing duplicates to samtools faidx produces a FASTA
    # with duplicate sequence names, which causes downstream tools (minimap2,
    # samtools sort) to reject the file.
    # dict.fromkeys() preserves insertion order (important for reproducible
    # FASTA output) while removing duplicates in O(n) time.
    hap1_asm_regions = list(dict.fromkeys(hap1_asm_regions))
    hap2_asm_regions = list(dict.fromkeys(hap2_asm_regions))

    toy_hap1 = os.path.join(args.output_dir, "toy_hap1.fa")
    toy_hap2 = os.path.join(args.output_dir, "toy_hap2.fa")
    extract_fasta_regions(args.hap1, hap1_asm_regions, toy_hap1)
    extract_fasta_regions(args.hap2, hap2_asm_regions, toy_hap2)
    index_and_compress_fasta(toy_hap1)
    index_and_compress_fasta(toy_hap2)

    # ── Step 5: Extract reads ───────────────────────────────────────────
    print("Step 5: Extracting reads for selected regions", file=sys.stderr)
    toy_reads_bam = os.path.join(args.output_dir, "toy_reads.bam")
    extract_reads(
        args.cram,
        ref_regions,
        toy_reads_bam,
        reference=args.reference,
        cram_ref=args.cram_ref,
    )
    _run(["samtools", "index", toy_reads_bam], f"Indexing {toy_reads_bam}")

    # ── Step 6: Write manifest ──────────────────────────────────────────
    print("Step 6: Writing manifest", file=sys.stderr)
    write_manifest(regions, args.output_dir)

    # ── Summary ─────────────────────────────────────────────────────────
    print("", file=sys.stderr)
    print("Toy dataset generated successfully!", file=sys.stderr)
    print(f"  Output directory: {args.output_dir}", file=sys.stderr)
    print("  Files:", file=sys.stderr)
    for f in sorted(os.listdir(args.output_dir)):
        fpath = os.path.join(args.output_dir, f)
        size_kb = os.path.getsize(fpath) / 1024
        print(f"    {f}  ({size_kb:.1f} KB)", file=sys.stderr)


if __name__ == "__main__":
    main()
