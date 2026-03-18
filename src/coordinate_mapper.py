#!/usr/bin/env python3
"""
coordinate_mapper.py - Build and query assembly-to-reference coordinate mappings.

Creates an efficient interval-based lookup from minimap2 PAF alignments,
enabling fast translation between reference genome coordinates and
haplotype assembly coordinates.

Structural-variation aware: when the `cg:Z:` CIGAR tag is present in the PAF
(produced by `minimap2 -c` or `minimap2 --cs`), coordinate projection is done
by walking the CIGAR base-by-base rather than by linear interpolation.  Gaps
between alignment blocks are detected and classified as deletions, insertions,
inversions, or translocations.

Usage:
    # Build mapping index from PAF
    coordinate_mapper.py build -p hap1_to_ref.paf -o output_prefix

    # Query mapping by reference region
    coordinate_mapper.py query -i output_prefix.mapping.json.gz -r chr1:1000000-2000000
"""

import argparse
import bisect
import gzip
import json
import re
import sys
from collections import defaultdict


def parse_paf(paf_file):
    """Parse a minimap2 PAF file and extract alignment blocks.

    Optional PAF tags extracted when present:
        cg:Z:  CIGAR string (produced by ``minimap2 -c`` / ``--cs``).
               Used for precise base-level coordinate translation.
        tp:A:  Alignment type (P=primary, S=secondary, I=inversion).

    Args:
        paf_file: Path to PAF file (plain or gzip-compressed).

    Returns:
        List of alignment block dicts with reference and assembly coordinates.
    """
    blocks = []
    opener = gzip.open if paf_file.endswith(".gz") else open
    with opener(paf_file, "rt") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 12:
                continue

            # Extract optional type:value tags (fields 12+)
            tag_dict = {}
            for field in fields[12:]:
                parts = field.split(":", 2)
                if len(parts) == 3:
                    tag_dict[parts[0]] = parts[2]

            blocks.append(
                {
                    "ref_chrom": fields[5],
                    "ref_start": int(fields[7]),
                    "ref_end": int(fields[8]),
                    "asm_chrom": fields[0],
                    "asm_start": int(fields[2]),
                    "asm_end": int(fields[3]),
                    "strand": fields[4],
                    "mapq": int(fields[11]),
                    "matches": int(fields[9]),
                    "block_len": int(fields[10]),
                    "cigar": tag_dict.get("cg", ""),
                    "tp": tag_dict.get("tp", ""),
                }
            )
    return blocks


# ---------------------------------------------------------------------------
# CIGAR helpers
# ---------------------------------------------------------------------------


def parse_cigar(cigar_str):
    """Parse a CIGAR string into a list of (length, op) tuples.

    Supports the CIGAR operations used in PAF / SAM:
        M (match/mismatch), = (sequence match), X (sequence mismatch),
        D (deletion – reference bases absent from assembly),
        N (reference skip, treated like D),
        I (insertion in assembly), S (soft clip, treated like I),
        H / P (hard clip / padding – consume nothing).

    Args:
        cigar_str: CIGAR string, e.g. ``"100M5D50M10I200M"``.
                   Returns an empty list for ``None`` or ``""``.

    Returns:
        List of ``(length, op)`` tuples.
    """
    if not cigar_str:
        return []
    return [
        (int(m.group(1)), m.group(2))
        for m in re.finditer(r"(\d+)([MIDNSHPX=])", cigar_str)
    ]


def project_cigar(cigar_ops, ref_offset_start, ref_offset_end, strand,
                  asm_block_start, asm_block_end):
    """Walk a CIGAR string to project a reference sub-interval to assembly space.

    Gives precise base-level coordinate translation that correctly handles
    insertions and deletions within an alignment block — essential for
    accurately mapping across structural-variant breakpoints.

    CIGAR operation semantics (PAF / SAM conventions):
        M / = / X  – match/mismatch: advance both ref and asm cursors.
        D / N      – deletion from assembly: advance ref cursor only
                     (reference has bases; assembly does not).
        I / S      – insertion in assembly / soft-clip: advance asm cursor only
                     (assembly has bases; reference does not).
        H / P      – consume nothing.

    For deletions (D/N), the reference bases at the deletion site are
    mapped to the adjacent assembly position (i.e. both edges of the
    deletion clamp to the same assembly coordinate).

    For insertions (I/S), the reference position at the insertion site
    maps to the start of the inserted assembly sequence, so queries
    that straddle an insertion boundary include the full inserted range.

    Args:
        cigar_ops:        List of ``(length, op)`` tuples from
                          :func:`parse_cigar`.
        ref_offset_start: Start offset from the block's ``ref_start`` (≥ 0).
        ref_offset_end:   End offset from the block's ``ref_start``
                          (≥ ``ref_offset_start``).
        strand:           ``'+'`` or ``'-'``.
        asm_block_start:  Assembly start coordinate of the alignment block
                          (0-based, inclusive).
        asm_block_end:    Assembly end coordinate of the alignment block
                          (0-based, exclusive).

    Returns:
        ``(asm_start, asm_end)`` – absolute assembly coordinates
        (0-based, half-open).
    """
    # Sets of operations that advance each cursor
    REF_CONSUMING = frozenset("MDNX=")
    ASM_CONSUMING = frozenset("MISX=")

    ref_cur = 0  # position relative to block ref_start
    asm_cur = 0  # position relative to the start of asm traversal

    asm_q_start = None
    asm_q_end = None

    for length, op in cigar_ops:
        consumes_ref = op in REF_CONSUMING
        consumes_asm = op in ASM_CONSUMING

        if consumes_ref:
            # Does ref_offset_start fall inside this operation?
            if asm_q_start is None and ref_cur + length > ref_offset_start:
                inner = ref_offset_start - ref_cur
                # Deletions clamp asm to current cursor (no asm advance)
                asm_q_start = asm_cur + (inner if consumes_asm else 0)

            # Does ref_offset_end fall inside (or at the end of) this op?
            if asm_q_start is not None and asm_q_end is None:
                if ref_cur + length >= ref_offset_end:
                    inner = ref_offset_end - ref_cur
                    asm_q_end = asm_cur + (inner if consumes_asm else 0)
                    break

        if consumes_ref:
            ref_cur += length
        if consumes_asm:
            asm_cur += length

    # Fallback: query end was beyond the CIGAR; clamp to asm_cur reached.
    if asm_q_start is None:
        asm_q_start = asm_cur
    if asm_q_end is None:
        asm_q_end = asm_cur

    # Apply strand orientation to obtain absolute assembly coordinates.
    if strand == "+":
        return asm_block_start + asm_q_start, asm_block_start + asm_q_end
    else:
        # Minus strand: CIGAR walks asm from asm_block_end downward.
        return asm_block_end - asm_q_end, asm_block_end - asm_q_start


# ---------------------------------------------------------------------------
# SV gap classification
# ---------------------------------------------------------------------------


def classify_sv_gap(prev_block, next_block, gap_ref_start, gap_ref_end, chrom):
    """Classify the structural-variation event represented by a gap between blocks.

    Examines the reference-space gap and the assembly context of the flanking
    blocks to infer whether the gap corresponds to a deletion, insertion,
    inversion, translocation, or a more complex rearrangement.

    Classification rules (in priority order):

    * **translocation** – flanking blocks lie on different assembly contigs.
    * **inversion** – same contig, opposite strands.
    * **complex** – same contig, same strand, but the assembly gap is negative
      (blocks overlap in assembly space, indicating a duplication or other
      complex event).
    * **deletion** – same contig, same strand; assembly gap < reference gap
      (assembly has less sequence than reference at this locus).
    * **insertion** – same contig, same strand; assembly gap ≥ reference gap
      (assembly has more sequence than reference).

    Args:
        prev_block:    Block dict preceding the gap (sorted by ``ref_start``).
        next_block:    Block dict following the gap.
        gap_ref_start: Start of the gap in reference coordinates (0-based).
        gap_ref_end:   End of the gap in reference coordinates (0-based).
        chrom:         Reference chromosome name.

    Returns:
        Result dict with the following keys:

        * ``ref_chrom``, ``ref_start``, ``ref_end`` – reference gap coordinates.
        * ``asm_chrom``, ``asm_start``, ``asm_end`` – assembly coordinates
          flanking or spanning the event.
        * ``strand`` – strand of the preceding block.
        * ``mapq`` – minimum mapq of the two flanking blocks.
        * ``event_type`` – one of ``"deletion"``, ``"insertion"``,
          ``"inversion"``, ``"translocation"``, ``"complex"``.
        * ``ref_gap_size`` – size of the reference gap in bases.
        * ``asm_gap_size`` – estimated size of the assembly gap in bases
          (``None`` for inter-contig and inversion events).
    """
    ref_gap_size = gap_ref_end - gap_ref_start
    same_contig = prev_block["ac"] == next_block["ac"]
    same_strand = prev_block["st"] == next_block["st"]

    if not same_contig:
        event_type = "translocation"
        asm_gap_size = None
        # Report assembly position at the end of the preceding block
        asm_start = prev_block["ae"] if prev_block["st"] == "+" else prev_block["as"]
        asm_end = asm_start

    elif not same_strand:
        event_type = "inversion"
        asm_gap_size = None
        asm_start = prev_block["ae"] if prev_block["st"] == "+" else prev_block["as"]
        asm_end = asm_start

    else:
        # Same contig, same strand – compute signed assembly gap
        if prev_block["st"] == "+":
            asm_gap_size = next_block["as"] - prev_block["ae"]
            asm_start = prev_block["ae"]
            asm_end = next_block["as"]
        else:
            # On − strand, assembly runs high → low as ref runs low → high
            asm_gap_size = prev_block["as"] - next_block["ae"]
            asm_start = next_block["ae"]
            asm_end = prev_block["as"]

        if asm_gap_size < 0:
            event_type = "complex"
        elif asm_gap_size < ref_gap_size:
            event_type = "deletion"
        else:
            event_type = "insertion"

    return {
        "ref_chrom": chrom,
        "ref_start": gap_ref_start,
        "ref_end": gap_ref_end,
        "asm_chrom": prev_block["ac"],
        "asm_start": asm_start,
        "asm_end": asm_end,
        "strand": prev_block["st"],
        "mapq": min(prev_block["mq"], next_block["mq"]),
        "event_type": event_type,
        "ref_gap_size": ref_gap_size,
        "asm_gap_size": asm_gap_size,
    }


def build_mapping_bed(blocks, bed_path):
    """Write alignment blocks as a sorted BED file.

    Columns: ref_chrom, ref_start, ref_end, asm_chrom, asm_start, asm_end,
             strand, mapq, block_len, matches

    Args:
        blocks: List of alignment block dicts.
        bed_path: Output BED file path.
    """
    blocks.sort(key=lambda b: (b["ref_chrom"], b["ref_start"]))
    with open(bed_path, "w") as fh:
        for b in blocks:
            fh.write(
                f"{b['ref_chrom']}\t{b['ref_start']}\t{b['ref_end']}\t"
                f"{b['asm_chrom']}\t{b['asm_start']}\t{b['asm_end']}\t"
                f"{b['strand']}\t{b['mapq']}\t{b['block_len']}\t{b['matches']}\n"
            )
    return bed_path


def build_json_index(blocks, index_path):
    """Build a compact gzipped JSON index for fast programmatic queries.

    Structure: { "chrom": [ { rs, re, ac, as, ae, st, mq[, cg, tp] }, ... ] }
    Blocks within each chromosome are sorted by ref_start.
    The ``cg`` (CIGAR) and ``tp`` (alignment type) fields are stored only when
    non-empty so that indexes built from PAF files without ``-c`` remain compact.

    Args:
        blocks: List of alignment block dicts.
        index_path: Output .json.gz path.
    """
    index = defaultdict(list)
    for b in sorted(blocks, key=lambda b: (b["ref_chrom"], b["ref_start"])):
        entry = {
            "rs": b["ref_start"],
            "re": b["ref_end"],
            "ac": b["asm_chrom"],
            "as": b["asm_start"],
            "ae": b["asm_end"],
            "st": b["strand"],
            "mq": b["mapq"],
        }
        if b.get("cigar"):
            entry["cg"] = b["cigar"]
        if b.get("tp"):
            entry["tp"] = b["tp"]
        index[b["ref_chrom"]].append(entry)
    with gzip.open(index_path, "wt") as fh:
        json.dump(dict(index), fh, separators=(",", ":"))
    return index_path


def load_index(index_path):
    """Load a JSON index and prepare sorted arrays for binary search.

    Returns:
        Dict keyed by chromosome, each value containing 'blocks', 'starts',
        and 'ends' lists for efficient overlap queries.
    """
    with gzip.open(index_path, "rt") as fh:
        data = json.load(fh)

    index = {}
    for chrom, raw_blocks in data.items():
        starts = [b["rs"] for b in raw_blocks]
        ends = [b["re"] for b in raw_blocks]
        max_block_len = max((e - s for s, e in zip(starts, ends)), default=0)
        index[chrom] = {
            "blocks": raw_blocks,
            "starts": starts,
            "ends": ends,
            "max_block_len": max_block_len,
        }
    return index


def _has_sv_at_junction(prev_block, block):
    """Return True if there is a structural event at the junction between blocks.

    Used to detect assembly-level events (insertions, inversions, translocations)
    that occur at ref positions where two alignment blocks are adjacent
    (i.e. there is no reference-space gap).

    Args:
        prev_block: Block dict preceding the junction (sorted by ref_start).
        block:      Block dict following the junction.

    Returns:
        True when the blocks are on different contigs, different strands, or
        have a non-zero assembly gap (positive = insertion, negative = complex).
    """
    if prev_block["ac"] != block["ac"]:
        return True  # different assembly contigs → translocation
    if prev_block["st"] != block["st"]:
        return True  # different strands → inversion
    # Same contig, same strand: check for non-zero assembly gap.
    if prev_block["st"] == "+":
        asm_gap = block["as"] - prev_block["ae"]
    else:
        asm_gap = prev_block["as"] - block["ae"]
    return asm_gap != 0


def query(index, chrom, start, end, min_mapq=0):
    """Find assembly regions and SV events for a reference region.

    Uses binary search on sorted start positions to efficiently find all
    alignment blocks overlapping the query interval [start, end).

    For each overlapping block the assembly coordinates are computed using
    CIGAR-aware projection when a ``cg`` CIGAR string is stored in the index,
    or by linear interpolation otherwise.  Gaps between successive overlapping
    blocks are detected and classified as structural-variation events.

    Args:
        index: Loaded index from :func:`load_index`.
        chrom: Reference chromosome name.
        start: Query start (0-based, inclusive).
        end: Query end (0-based, exclusive).
        min_mapq: Minimum mapping quality.  Blocks with ``mapq < min_mapq``
            are silently excluded.  Useful for filtering supplementary
            or low-confidence alignments (default 0 = keep all).

    Returns:
        List of result dicts, ordered by ``ref_start``.  Each dict contains:

        * ``ref_chrom``, ``ref_start``, ``ref_end`` – reference coordinates.
        * ``asm_chrom``, ``asm_start``, ``asm_end`` – assembly coordinates.
        * ``strand`` – alignment strand (``'+'`` or ``'-'``).
        * ``mapq`` – mapping quality.
        * ``event_type`` – ``"alignment"`` for aligned regions; one of
          ``"deletion"``, ``"insertion"``, ``"inversion"``,
          ``"translocation"``, or ``"complex"`` for gap events.
        * ``ref_gap_size``, ``asm_gap_size`` – present only in gap events.
    """
    if chrom not in index:
        return []

    chrom_data = index[chrom]
    blocks = chrom_data["blocks"]
    starts = chrom_data["starts"]
    ends = chrom_data["ends"]

    # Find candidate blocks: block_start < end AND block_end > start
    right_idx = bisect.bisect_left(starts, end)

    # Use the precomputed maximum block length to skip blocks that are
    # too far left to overlap the query region, turning the linear scan
    # into O(k + log N) where k is the number of overlapping blocks.
    max_len = chrom_data.get("max_block_len", 0)
    left_idx = bisect.bisect_left(starts, start - max_len) if max_len > 0 else 0

    overlapping = []
    for i in range(left_idx, right_idx):
        if ends[i] > start and blocks[i]["mq"] >= min_mapq:
            overlapping.append(blocks[i])

    if not overlapping:
        return []

    results = []
    prev_block = None

    for block in overlapping:
        # Detect and classify the gap between the previous block and this one.
        # A gap event is emitted when:
        #   (a) there is a reference-space gap that intersects the query, OR
        #   (b) blocks are adjacent in reference but have an assembly-level
        #       event (insertion, inversion, translocation) at the junction
        #       and that junction falls inside the query window.
        if prev_block is not None:
            junction = prev_block["re"]   # ref end of the preceding block
            next_start = block["rs"]      # ref start of this block

            gap_start = max(start, junction)
            gap_end = min(end, next_start)

            has_ref_gap = gap_start < gap_end
            has_junction_event = (
                start <= junction <= end
                and junction == next_start
                and _has_sv_at_junction(prev_block, block)
            )

            if has_ref_gap or has_junction_event:
                results.append(
                    classify_sv_gap(prev_block, block, gap_start, gap_end, chrom)
                )

        # Intersect the query interval with this block
        overlap_start = max(start, block["rs"])
        overlap_end = min(end, block["re"])

        if overlap_start >= overlap_end:
            prev_block = block
            continue

        ref_offset_start = overlap_start - block["rs"]
        ref_offset_end = overlap_end - block["rs"]

        # Prefer CIGAR-based projection; fall back to linear interpolation.
        if block.get("cg"):
            cigar_ops = parse_cigar(block["cg"])
            asm_start, asm_end = project_cigar(
                cigar_ops,
                ref_offset_start,
                ref_offset_end,
                block["st"],
                block["as"],
                block["ae"],
            )
        else:
            if block["st"] == "+":
                asm_start = block["as"] + ref_offset_start
                asm_end = block["as"] + ref_offset_end
            else:
                asm_start = block["ae"] - ref_offset_end
                asm_end = block["ae"] - ref_offset_start

        results.append(
            {
                "ref_chrom": chrom,
                "ref_start": overlap_start,
                "ref_end": overlap_end,
                "asm_chrom": block["ac"],
                "asm_start": asm_start,
                "asm_end": asm_end,
                "strand": block["st"],
                "mapq": block["mq"],
                "event_type": "alignment",
            }
        )

        prev_block = block

    return results


def parse_region(region_str):
    """Parse a region string like 'chr1:1000-2000' into (chrom, start, end)."""
    if ":" not in region_str:
        raise ValueError(f"Invalid region format: {region_str} (expected chrom:start-end)")
    chrom, coords = region_str.split(":", 1)
    if "-" not in coords:
        raise ValueError(f"Invalid region format: {region_str} (expected chrom:start-end)")
    start_str, end_str = coords.split("-", 1)
    return chrom, int(start_str), int(end_str)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def cmd_build(args):
    """Build mapping index from one or more PAF files."""
    all_blocks = []
    for paf in args.paf:
        all_blocks.extend(parse_paf(paf))

    if not all_blocks:
        print("WARNING: no alignment blocks found in PAF input", file=sys.stderr)

    bed_path = f"{args.output}.mapping.bed"
    build_mapping_bed(all_blocks, bed_path)
    print(f"BED mapping written to {bed_path}", file=sys.stderr)

    json_path = f"{args.output}.mapping.json.gz"
    build_json_index(all_blocks, json_path)
    print(f"JSON index written to {json_path}", file=sys.stderr)

    print(
        f"Tip: bgzip {bed_path} && tabix -p bed {bed_path}.gz  "
        "for tabix-indexed queries",
        file=sys.stderr,
    )


def cmd_query(args):
    """Query the mapping index for a reference region."""
    index = load_index(args.index)
    chrom, start, end = parse_region(args.region)
    results = query(index, chrom, start, end, min_mapq=args.min_mapq)

    if not results:
        print(f"No mappings found for {args.region}", file=sys.stderr)
        sys.exit(0)

    # Output as TSV – event_type column distinguishes alignments from SV gaps.
    header = (
        "ref_chrom\tref_start\tref_end\t"
        "asm_chrom\tasm_start\tasm_end\t"
        "strand\tmapq\tevent_type"
    )
    print(header)
    for r in results:
        print(
            f"{r['ref_chrom']}\t{r['ref_start']}\t{r['ref_end']}\t"
            f"{r['asm_chrom']}\t{r['asm_start']}\t{r['asm_end']}\t"
            f"{r['strand']}\t{r['mapq']}\t{r['event_type']}"
        )


def main():
    parser = argparse.ArgumentParser(
        description="Build and query assembly-to-reference coordinate mappings"
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    # build
    build_p = subparsers.add_parser(
        "build", help="Build mapping index from PAF alignment(s)"
    )
    build_p.add_argument(
        "-p",
        "--paf",
        nargs="+",
        required=True,
        help="One or more minimap2 PAF files",
    )
    build_p.add_argument(
        "-o",
        "--output",
        required=True,
        help="Output prefix (produces <prefix>.mapping.bed and <prefix>.mapping.json.gz)",
    )

    # query
    query_p = subparsers.add_parser(
        "query", help="Query mapping index for a reference region"
    )
    query_p.add_argument(
        "-i",
        "--index",
        required=True,
        help="Path to .mapping.json.gz index file",
    )
    query_p.add_argument(
        "-r",
        "--region",
        required=True,
        help="Reference region to query (e.g. chr1:1000000-2000000)",
    )
    query_p.add_argument(
        "--min-mapq",
        type=int,
        default=0,
        help="Minimum mapping quality to include (default: 0 = keep all)",
    )

    args = parser.parse_args()
    if args.command == "build":
        cmd_build(args)
    elif args.command == "query":
        cmd_query(args)


if __name__ == "__main__":
    main()
