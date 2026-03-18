#!/usr/bin/env python3
"""
coordinate_mapper.py - Build and query assembly-to-reference coordinate mappings.

Creates an efficient interval-based lookup from minimap2 PAF alignments,
enabling fast translation between reference genome coordinates and
haplotype assembly coordinates.

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
import sys
from collections import defaultdict


def parse_paf(paf_file):
    """Parse a minimap2 PAF file and extract alignment blocks.

    Args:
        paf_file: Path to PAF file.

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
                }
            )
    return blocks


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

    Structure: { "chrom": [ { ref_start, ref_end, asm_chrom, asm_start,
                               asm_end, strand, mapq }, ... ] }
    Blocks within each chromosome are sorted by ref_start.

    Args:
        blocks: List of alignment block dicts.
        index_path: Output .json.gz path.
    """
    index = defaultdict(list)
    for b in sorted(blocks, key=lambda b: (b["ref_chrom"], b["ref_start"])):
        index[b["ref_chrom"]].append(
            {
                "rs": b["ref_start"],
                "re": b["ref_end"],
                "ac": b["asm_chrom"],
                "as": b["asm_start"],
                "ae": b["asm_end"],
                "st": b["strand"],
                "mq": b["mapq"],
            }
        )
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
        index[chrom] = {"blocks": raw_blocks, "starts": starts, "ends": ends}
    return index


def query(index, chrom, start, end):
    """Find assembly regions corresponding to a reference region.

    Uses binary search on sorted start positions to efficiently find all
    alignment blocks overlapping the query interval [start, end).

    For each overlapping block the precise assembly coordinates are computed
    by projecting the overlap offsets through the alignment, accounting for
    strand orientation.

    Args:
        index: Loaded index from load_index().
        chrom: Reference chromosome name.
        start: Query start (0-based, inclusive).
        end: Query end (0-based, exclusive).

    Returns:
        List of dicts with ref and assembly coordinates for each overlap.
    """
    if chrom not in index:
        return []

    chrom_data = index[chrom]
    blocks = chrom_data["blocks"]
    starts = chrom_data["starts"]
    ends = chrom_data["ends"]

    # Find candidate blocks: block_start < end AND block_end > start
    right_idx = bisect.bisect_left(starts, end)

    results = []
    for i in range(right_idx):
        if ends[i] <= start:
            continue

        block = blocks[i]
        # Compute overlap in reference coordinates
        overlap_start = max(start, block["rs"])
        overlap_end = min(end, block["re"])

        # Project to assembly coordinates
        ref_offset_start = overlap_start - block["rs"]
        ref_offset_end = overlap_end - block["rs"]

        if block["st"] == "+":
            asm_overlap_start = block["as"] + ref_offset_start
            asm_overlap_end = block["as"] + ref_offset_end
        else:
            asm_overlap_start = block["ae"] - ref_offset_end
            asm_overlap_end = block["ae"] - ref_offset_start

        results.append(
            {
                "ref_chrom": chrom,
                "ref_start": overlap_start,
                "ref_end": overlap_end,
                "asm_chrom": block["ac"],
                "asm_start": asm_overlap_start,
                "asm_end": asm_overlap_end,
                "strand": block["st"],
                "mapq": block["mq"],
            }
        )
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
    results = query(index, chrom, start, end)

    if not results:
        print(f"No mappings found for {args.region}", file=sys.stderr)
        sys.exit(0)

    # Output as TSV
    header = "ref_chrom\tref_start\tref_end\tasm_chrom\tasm_start\tasm_end\tstrand\tmapq"
    print(header)
    for r in results:
        print(
            f"{r['ref_chrom']}\t{r['ref_start']}\t{r['ref_end']}\t"
            f"{r['asm_chrom']}\t{r['asm_start']}\t{r['asm_end']}\t"
            f"{r['strand']}\t{r['mapq']}"
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

    args = parser.parse_args()
    if args.command == "build":
        cmd_build(args)
    elif args.command == "query":
        cmd_query(args)


if __name__ == "__main__":
    main()
