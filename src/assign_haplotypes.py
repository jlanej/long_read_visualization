#!/usr/bin/env python3
"""Assign haplotype phase tags (HP) to reads via competitive alignment scoring.

Given two BAM files — reads aligned to haplotype 1 and reads aligned to
haplotype 2 — this script compares the minimap2 Alignment Score (``AS:i:``)
for every read and assigns an HP tag:

* ``HP:i:1`` — read aligns better to hap1
* ``HP:i:2`` — read aligns better to hap2
* ``HP:i:0`` — ambiguous (scores within *tolerance* of each other)

The script performs two operations:

1. **In-place tagging** of the pipeline-generated ``reads_to_hap{1,2}.bam``
   files (haplotype-space BAMs used in the assembly panels).  The original
   ``--cram`` input supplied by the user is *never* modified.

2. **HP-tagged reads file for reference panel**:
   either write a new file (``--reads-out``) derived from ``--reads-in``,
   or tag ``--reads-in`` in place with ``--reads-in-place``.  This file is
   used as the reads track in the reference panel IGV browser.

IGV natively understands the ``HP`` tag and can sort, group, and colour
reads by haplotype phase.

The script is **idempotent**: HP tags are always recomputed from the AS
scores so re-running produces the same result.

Usage
-----
    python3 assign_haplotypes.py \\
        --hap1-bam  sample_reads_to_hap1.bam \\
        --hap2-bam  sample_reads_to_hap2.bam \\
        --reads-in  sample.cram \\
        [--reads-out sample_reads.hp.cram | --reads-in-place] \\
        [--reference ref.fa] \\
        [--tolerance 0]
"""

import argparse
import os
import sys
import tempfile
import shutil

try:
    import pysam
except ImportError:
    sys.exit(
        "ERROR: pysam is required for haplotype assignment.\n"
        "Install with: pip install pysam"
    )


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def get_best_scores(bam_path):
    """Return a dict mapping read name → best Alignment Score (AS).

    When a read has multiple alignments (supplementary / secondary), the
    highest AS value is kept.
    """
    scores = {}
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam.fetch(until_eof=True):
            if read.is_unmapped:
                continue
            qname = read.query_name
            try:
                als = read.get_tag("AS")
            except KeyError:
                continue
            if qname not in scores or als > scores[qname]:
                scores[qname] = als
    return scores


def assign_haplotypes(hap1_scores, hap2_scores, tolerance=0):
    """Compare alignment scores and return read → HP assignment.

    Returns:
        dict mapping read_name → 1 (hap1), 2 (hap2), or 0 (ambiguous).
    """
    all_reads = set(hap1_scores) | set(hap2_scores)
    assignments = {}
    for qname in all_reads:
        s1 = hap1_scores.get(qname)
        s2 = hap2_scores.get(qname)

        if s1 is not None and s2 is not None:
            diff = s1 - s2
            if diff > tolerance:
                assignments[qname] = 1
            elif diff < -tolerance:
                assignments[qname] = 2
            else:
                assignments[qname] = 0
        elif s1 is not None:
            assignments[qname] = 1
        elif s2 is not None:
            assignments[qname] = 2
        # Both None should not happen (at least one is set).
    return assignments


def tag_bam(input_bam, output_bam, assignments):
    """Write *output_bam* with HP:i: tags injected/updated from *assignments*.

    Reads whose query name is not in *assignments* are written through
    unchanged (they keep any pre-existing HP tag).
    """
    with pysam.AlignmentFile(input_bam, "rb") as infile:
        with pysam.AlignmentFile(output_bam, "wb", header=infile.header) as outfile:
            for read in infile.fetch(until_eof=True):
                hp = assignments.get(read.query_name)
                if hp is not None:
                    read.set_tag("HP", hp, value_type="i")
                outfile.write(read)


def tag_reads_file(input_path, output_path, assignments, reference=None):
    """Write a new *output_path* with HP:i: tags from *assignments*.

    Supports both BAM (``.bam``) and CRAM (``.cram``) output formats.  When
    writing CRAM, a *reference* FASTA path must be supplied.

    The original *input_path* is never modified.  Reads whose query name is
    not in *assignments* pass through unchanged.
    """
    is_cram_out = output_path.endswith(".cram")
    if is_cram_out and not reference:
        raise ValueError(
            f"--reference is required when writing CRAM output: {output_path}"
        )

    open_kw_in = {}
    if input_path.endswith(".cram") and reference:
        open_kw_in["reference_filename"] = reference

    if is_cram_out:
        mode_out = "wc"
        open_kw_out = {"reference_filename": reference}
    else:
        mode_out = "wb"
        open_kw_out = {}

    with pysam.AlignmentFile(input_path, "rb", **open_kw_in) as infile:
        with pysam.AlignmentFile(
            output_path, mode_out, header=infile.header, **open_kw_out
        ) as outfile:
            for read in infile.fetch(until_eof=True):
                hp = assignments.get(read.query_name)
                if hp is not None:
                    read.set_tag("HP", hp, value_type="i")
                outfile.write(read)


def tag_reads_file_in_place(path, assignments, reference=None):
    """Tag a BAM/CRAM with HP values, replacing the file atomically."""
    reads_dir = os.path.dirname(path) or "."
    suffix = ".cram" if path.endswith(".cram") else ".bam"
    fd, tmp_path = tempfile.mkstemp(suffix=suffix, dir=reads_dir)
    os.close(fd)
    try:
        tag_reads_file(path, tmp_path, assignments, reference=reference)
        shutil.move(tmp_path, path)
        index_bam(path)
    except Exception:
        if os.path.exists(tmp_path):
            os.remove(tmp_path)
        raise


def index_bam(bam_path):
    """Create a BAI/CRAI index for a BAM or CRAM file."""
    pysam.index(bam_path)


def tag_bam_in_place(bam_path, assignments):
    """Tag a BAM with HP values, replacing the file atomically.

    Writes to a temporary file in the same directory and renames on
    success so that a crash never leaves a half-written BAM.  Only
    pipeline-generated BAMs should be passed here; the original user
    CRAM input is never touched by this script.
    """
    bam_dir = os.path.dirname(bam_path) or "."
    fd, tmp_path = tempfile.mkstemp(suffix=".bam", dir=bam_dir)
    os.close(fd)
    try:
        tag_bam(bam_path, tmp_path, assignments)
        shutil.move(tmp_path, bam_path)
        index_bam(bam_path)
    except Exception:
        if os.path.exists(tmp_path):
            os.remove(tmp_path)
        raise


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main(args=None):
    parser = argparse.ArgumentParser(
        description="Assign HP (haplotype phase) tags to reads via "
                    "competitive alignment scoring."
    )
    parser.add_argument(
        "--hap1-bam", required=True,
        help="Pipeline-generated BAM of reads aligned to haplotype 1",
    )
    parser.add_argument(
        "--hap2-bam", required=True,
        help="Pipeline-generated BAM of reads aligned to haplotype 2",
    )
    parser.add_argument(
        "--reads-in",
        help="Original CRAM/BAM input (e.g. the --cram passed to preprocess.sh). "
             "Used together with --reads-out to create an HP-tagged copy in "
             "reference coordinate space for IGV visualisation.",
    )
    parser.add_argument(
        "--reads-out",
        help="Path for the new HP-tagged output file (CRAM or BAM). "
             "The original --reads-in file is not modified.",
    )
    parser.add_argument(
        "--reads-in-place",
        action="store_true",
        help="Tag --reads-in in place (atomic replace) instead of writing "
             "--reads-out. Useful after optional remapping workflows.",
    )
    parser.add_argument(
        "--reference",
        help="Reference FASTA (required when --reads-out is a CRAM file).",
    )
    parser.add_argument(
        "--tolerance", type=int, default=0,
        help="Alignment-score difference within which a read is called "
             "ambiguous (HP:i:0). Default: 0 (strict equality).",
    )
    opts = parser.parse_args(args)

    hap1_bam = opts.hap1_bam
    hap2_bam = opts.hap2_bam

    if not os.path.isfile(hap1_bam):
        sys.exit(f"ERROR: hap1 BAM not found: {hap1_bam}")
    if not os.path.isfile(hap2_bam):
        sys.exit(f"ERROR: hap2 BAM not found: {hap2_bam}")

    # Validate reads-in/reads-out options
    if opts.reads_in_place:
        if not opts.reads_in:
            sys.exit("ERROR: --reads-in-place requires --reads-in")
        if opts.reads_out:
            sys.exit("ERROR: --reads-out cannot be used with --reads-in-place")
    elif bool(opts.reads_in) != bool(opts.reads_out):
        sys.exit("ERROR: --reads-in and --reads-out must be used together")
    if opts.reads_in and not os.path.isfile(opts.reads_in):
        sys.exit(f"ERROR: reads-in file not found: {opts.reads_in}")

    print(f"[assign_haplotypes] Reading alignment scores from {hap1_bam}",
          file=sys.stderr)
    hap1_scores = get_best_scores(hap1_bam)
    print(f"[assign_haplotypes] Reading alignment scores from {hap2_bam}",
          file=sys.stderr)
    hap2_scores = get_best_scores(hap2_bam)

    print(f"[assign_haplotypes] Reads with AS in hap1: {len(hap1_scores)}, "
          f"hap2: {len(hap2_scores)}", file=sys.stderr)

    assignments = assign_haplotypes(hap1_scores, hap2_scores,
                                    tolerance=opts.tolerance)

    counts = {0: 0, 1: 0, 2: 0}
    for hp in assignments.values():
        counts[hp] += 1
    print(f"[assign_haplotypes] Assignments — "
          f"hap1: {counts[1]}, hap2: {counts[2]}, "
          f"ambiguous: {counts[0]}", file=sys.stderr)

    # ── 1. In-place tagging of the haplotype-space BAMs ─────────────────────
    print(f"[assign_haplotypes] Tagging {hap1_bam}", file=sys.stderr)
    tag_bam_in_place(hap1_bam, assignments)
    print(f"[assign_haplotypes] Tagging {hap2_bam}", file=sys.stderr)
    tag_bam_in_place(hap2_bam, assignments)

    # ── 2. HP-tagging of reads file for the reference panel ──────────────────
    if opts.reads_in and opts.reads_in_place:
        print(f"[assign_haplotypes] Tagging reads file in place: {opts.reads_in}",
              file=sys.stderr)
        tag_reads_file_in_place(opts.reads_in, assignments,
                                reference=opts.reference)
    elif opts.reads_in and opts.reads_out:
        print(
            f"[assign_haplotypes] Writing HP-tagged reads file: "
            f"{opts.reads_in} → {opts.reads_out}",
            file=sys.stderr,
        )
        tag_reads_file(opts.reads_in, opts.reads_out, assignments,
                       reference=opts.reference)
        index_bam(opts.reads_out)

    print("[assign_haplotypes] Done.", file=sys.stderr)


if __name__ == "__main__":
    main()
