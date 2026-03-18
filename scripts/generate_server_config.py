#!/usr/bin/env python3
"""Generate a server configuration TSV from pipeline output directories.

Scans one or more pipeline output directories and optional VCF/manifest
files to produce a tab-separated configuration file for the IGV.js
visualization server.

Usage examples
--------------
    # Single sample from pipeline output
    python3 scripts/generate_server_config.py \\
        --output-dir /path/to/output/NA21110 \\
        --reference /path/to/chm13v2.0.fa.gz \\
        --hap1 /path/to/hap1.fa.gz \\
        --hap2 /path/to/hap2.fa.gz \\
        --regions /path/to/manifest.json \\
        -o config.tsv

    # Multiple samples
    python3 scripts/generate_server_config.py \\
        --output-dir /path/to/output/NA21110 /path/to/output/HG002 \\
        --reference /path/to/ref.fa.gz \\
        --hap1 /path/to/NA21110_hap1.fa.gz /path/to/HG002_hap1.fa.gz \\
        --hap2 /path/to/NA21110_hap2.fa.gz /path/to/HG002_hap2.fa.gz \\
        -o config.tsv

    # Auto-discover from a parent directory containing sample subdirs
    python3 scripts/generate_server_config.py \\
        --scan-dir /path/to/output \\
        --reference /path/to/ref.fa.gz \\
        -o config.tsv

    # From the toy dataset (after running the pipeline)
    python3 scripts/generate_server_config.py \\
        --output-dir /tmp/toy_output \\
        --reference resources/toy_dataset/toy_reference.fa.gz \\
        --hap1 resources/toy_dataset/toy_hap1.fa.gz \\
        --hap2 resources/toy_dataset/toy_hap2.fa.gz \\
        --regions resources/toy_dataset/toy_manifest.json \\
        --reads-bam resources/toy_dataset/toy_reads.bam \\
        -o config.tsv
"""

import argparse
import glob
import os
import sys


REQUIRED_SUFFIXES = [
    "_hap1_to_ref.bam",
    "_hap2_to_ref.bam",
    "_reads_to_hap1.bam",
    "_reads_to_hap2.bam",
]


def detect_sample_prefix(output_dir):
    """Detect the sample prefix from a pipeline output directory.

    Returns the prefix string, or None if no pipeline files are found.
    The prefix is the common filename stem before ``_hap1_to_ref.bam``.
    """
    for fname in os.listdir(output_dir):
        if fname.endswith("_hap1_to_ref.bam"):
            return fname.replace("_hap1_to_ref.bam", "")
    return None


def find_assembly_fasta(output_dir, prefix, hap_label):
    """Try to find a haplotype assembly FASTA near the output directory.

    Searches the output directory and its parent for patterns like:
        {prefix}_{hap_label}*.fa.gz
    """
    patterns = [
        os.path.join(output_dir, f"*{hap_label}*.fa.gz"),
        os.path.join(os.path.dirname(output_dir), f"*{hap_label}*.fa.gz"),
        os.path.join(os.path.dirname(output_dir), prefix,
                     f"*{hap_label}*.fa.gz"),
    ]
    for pat in patterns:
        matches = glob.glob(pat)
        # Filter out coordinate mapping index files
        matches = [m for m in matches if not m.endswith(".mapping.json.gz")]
        if matches:
            return os.path.abspath(matches[0])
    return ""


def find_regions_file(output_dir, prefix):
    """Try to find a regions file (manifest JSON or VCF) near the output dir."""
    candidates = [
        os.path.join(output_dir, "manifest.json"),
        os.path.join(output_dir, f"{prefix}_manifest.json"),
        os.path.join(os.path.dirname(output_dir), f"{prefix}_manifest.json"),
    ]
    for c in candidates:
        if os.path.isfile(c):
            return os.path.abspath(c)
    # Look for VCF files
    for ext in ("*.vcf.gz", "*.vcf"):
        matches = glob.glob(os.path.join(output_dir, ext))
        if not matches:
            matches = glob.glob(
                os.path.join(os.path.dirname(output_dir), ext))
        if matches:
            return os.path.abspath(matches[0])
    return ""


def scan_parent_dir(scan_dir, reference):
    """Scan a parent directory for subdirectories containing pipeline outputs.

    Returns a list of (output_dir, prefix) tuples.
    """
    results = []
    for entry in sorted(os.listdir(scan_dir)):
        subdir = os.path.join(scan_dir, entry)
        if not os.path.isdir(subdir):
            continue
        prefix = detect_sample_prefix(subdir)
        if prefix:
            results.append((subdir, prefix))
    return results


def build_sample_row(output_dir, prefix, reference="", hap1="", hap2="",
                     regions="", reads_bam=""):
    """Build a sample configuration dict.

    Resolves all paths to absolute paths and fills in missing values
    by auto-discovery when possible.
    """
    output_dir = os.path.abspath(output_dir)
    sample_id = prefix

    # Resolve assemblies
    if not hap1:
        hap1 = find_assembly_fasta(output_dir, prefix, "hap1")
    if not hap2:
        hap2 = find_assembly_fasta(output_dir, prefix, "hap2")

    # Resolve regions
    if not regions:
        regions = find_regions_file(output_dir, prefix)

    return {
        "sample_id": sample_id,
        "output_dir": output_dir,
        "reference": os.path.abspath(reference) if reference else "",
        "hap1_assembly": os.path.abspath(hap1) if hap1 else "",
        "hap2_assembly": os.path.abspath(hap2) if hap2 else "",
        "reads_bam": os.path.abspath(reads_bam) if reads_bam else "",
        "regions": os.path.abspath(regions) if regions else "",
    }


def validate_sample(row):
    """Validate that essential files exist. Print warnings for missing files."""
    warnings = []
    for key in ("reference", "hap1_assembly", "hap2_assembly"):
        path = row.get(key, "")
        if not path:
            warnings.append(f"  WARNING: {key} not specified")
        elif not os.path.isfile(path):
            warnings.append(f"  WARNING: {key} file not found: {path}")

    # Check pipeline output files
    output_dir = row["output_dir"]
    prefix = row["sample_id"]
    for suffix in REQUIRED_SUFFIXES:
        fpath = os.path.join(output_dir, f"{prefix}{suffix}")
        if not os.path.isfile(fpath):
            warnings.append(f"  WARNING: Missing pipeline output: "
                            f"{prefix}{suffix}")

    return warnings


def write_config(samples, output_path):
    """Write samples to a TSV configuration file."""
    header = [
        "sample_id", "output_dir", "reference", "hap1_assembly",
        "hap2_assembly", "reads_bam", "regions",
    ]
    with open(output_path, "w") as fh:
        fh.write("#" + "\t".join(header) + "\n")
        for row in samples:
            fields = [row.get(col, "") for col in header]
            fh.write("\t".join(fields) + "\n")


def main():
    parser = argparse.ArgumentParser(
        description="Generate server configuration TSV from pipeline outputs",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "--output-dir", "-d", nargs="+", default=[],
        help="Pipeline output directory (one per sample)")
    parser.add_argument(
        "--scan-dir", "-s",
        help="Parent directory to scan for sample subdirectories")
    parser.add_argument(
        "--reference", "-r", default="",
        help="Reference FASTA (shared across all samples, or per-sample "
             "if multiple values match --output-dir count)")
    parser.add_argument(
        "--hap1", nargs="*", default=[],
        help="Hap1 assembly FASTA(s), one per sample")
    parser.add_argument(
        "--hap2", nargs="*", default=[],
        help="Hap2 assembly FASTA(s), one per sample")
    parser.add_argument(
        "--reads-bam", nargs="*", default=[],
        help="Reads BAM(s) aligned to reference, one per sample (optional)")
    parser.add_argument(
        "--regions", nargs="*", default=[],
        help="Region file(s) — manifest JSON or VCF (optional)")
    parser.add_argument(
        "--output", "-o", default="samples.tsv",
        help="Output TSV file path (default: samples.tsv)")
    args = parser.parse_args()

    # Collect sample directories
    dirs_and_prefixes = []

    if args.scan_dir:
        if not os.path.isdir(args.scan_dir):
            print(f"ERROR: --scan-dir not found: {args.scan_dir}",
                  file=sys.stderr)
            sys.exit(1)
        dirs_and_prefixes = scan_parent_dir(args.scan_dir, args.reference)
        if not dirs_and_prefixes:
            print(f"WARNING: No pipeline outputs found in {args.scan_dir}",
                  file=sys.stderr)

    for d in args.output_dir:
        d = os.path.abspath(d)
        if not os.path.isdir(d):
            print(f"ERROR: Output directory not found: {d}", file=sys.stderr)
            sys.exit(1)
        prefix = detect_sample_prefix(d)
        if not prefix:
            print(f"WARNING: No pipeline outputs found in {d}", file=sys.stderr)
            continue
        dirs_and_prefixes.append((d, prefix))

    if not dirs_and_prefixes:
        print("ERROR: No sample directories specified or found. "
              "Use --output-dir or --scan-dir.", file=sys.stderr)
        sys.exit(1)

    # Build sample rows
    samples = []
    for i, (out_dir, prefix) in enumerate(dirs_and_prefixes):
        hap1 = args.hap1[i] if i < len(args.hap1) else ""
        hap2 = args.hap2[i] if i < len(args.hap2) else ""
        reads = args.reads_bam[i] if i < len(args.reads_bam) else ""
        rgn = args.regions[i] if i < len(args.regions) else ""

        row = build_sample_row(
            out_dir, prefix,
            reference=args.reference,
            hap1=hap1, hap2=hap2,
            regions=rgn, reads_bam=reads,
        )
        samples.append(row)

    # Validate and report
    print(f"Found {len(samples)} sample(s):\n")
    all_ok = True
    for row in samples:
        print(f"  {row['sample_id']}:")
        print(f"    output_dir:    {row['output_dir']}")
        print(f"    reference:     {row['reference']}")
        print(f"    hap1_assembly: {row['hap1_assembly']}")
        print(f"    hap2_assembly: {row['hap2_assembly']}")
        if row["reads_bam"]:
            print(f"    reads_bam:     {row['reads_bam']}")
        if row["regions"]:
            print(f"    regions:       {row['regions']}")
        warnings = validate_sample(row)
        for w in warnings:
            print(w)
            all_ok = False
        print()

    # Write config
    write_config(samples, args.output)
    print(f"Configuration written to {args.output}")
    if not all_ok:
        print("\nSome warnings were generated. Check paths above.",
              file=sys.stderr)


if __name__ == "__main__":
    main()
