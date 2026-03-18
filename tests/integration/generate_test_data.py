#!/usr/bin/env python3
"""Generate minimal synthetic test data for integration testing.

Creates a small reference FASTA (two 5 kb chromosomes), identical hap1/hap2
assembly FASTAs, and a SAM file with one read taken directly from the reference
so minimap2 can align it end-to-end.

Usage: python3 generate_test_data.py <output_dir>
"""
import os
import random
import sys


def main() -> None:
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} <output_dir>", file=sys.stderr)
        sys.exit(1)

    outdir = sys.argv[1]
    os.makedirs(outdir, exist_ok=True)

    random.seed(42)
    bases = "ACGT"
    chrom_len = 5000

    # Build two chromosomes of deterministic random sequence
    ref_seqs: dict[str, str] = {}
    for chrom in ("chr1", "chr2"):
        ref_seqs[chrom] = "".join(random.choice(bases) for _ in range(chrom_len))

    # Reference and both assembly FASTAs are identical (perfect alignment)
    for fname in ("ref.fa", "hap1.fa", "hap2.fa"):
        with open(os.path.join(outdir, fname), "w") as f:
            for chrom, seq in ref_seqs.items():
                f.write(f">{chrom}\n")
                for i in range(0, len(seq), 80):
                    f.write(seq[i : i + 80] + "\n")

    # One synthetic read: exact substring from chr1 positions 100–149 (0-based)
    read_seq = ref_seqs["chr1"][100:150]
    read_qual = "I" * len(read_seq)

    with open(os.path.join(outdir, "reads.sam"), "w") as f:
        f.write("@HD\tVN:1.6\tSO:coordinate\n")
        f.write(f"@SQ\tSN:chr1\tLN:{chrom_len}\n")
        f.write(f"@SQ\tSN:chr2\tLN:{chrom_len}\n")
        # SAM positions are 1-based; the read starts at reference base 101
        f.write(
            f"read1\t0\tchr1\t101\t60\t{len(read_seq)}M"
            f"\t*\t0\t0\t{read_seq}\t{read_qual}\n"
        )


if __name__ == "__main__":
    main()
