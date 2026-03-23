#!/usr/bin/env python3
"""Unit tests for the bundled toy dataset.

Validates that the committed toy dataset in ``resources/toy_dataset/`` is
complete, internally consistent, and ready to be fed into the pipeline.
These checks run without Docker or bioinformatics tools — only the Python
standard library is needed.
"""

import gzip
import json
import os
import struct
import sys
import unittest

# repo root → resources/toy_dataset
_REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
_TOY_DIR = os.path.join(_REPO_ROOT, "resources", "toy_dataset")
_MANIFEST = os.path.join(_TOY_DIR, "toy_manifest.json")


class TestToyDatasetFiles(unittest.TestCase):
    """Verify that every required file exists and is non-empty."""

    REQUIRED_FILES = [
        "toy_reference.fa.gz",
        "toy_reference.fa.gz.fai",
        "toy_reference.fa.gz.gzi",
        "toy_hap1.fa.gz",
        "toy_hap1.fa.gz.fai",
        "toy_hap1.fa.gz.gzi",
        "toy_hap2.fa.gz",
        "toy_hap2.fa.gz.fai",
        "toy_hap2.fa.gz.gzi",
        "toy_reads.bam",
        "toy_reads.bam.bai",
        "toy_manifest.json",
    ]

    def test_all_required_files_exist(self):
        """Every required toy dataset file exists."""
        for fname in self.REQUIRED_FILES:
            path = os.path.join(_TOY_DIR, fname)
            self.assertTrue(
                os.path.isfile(path), f"Missing required file: {fname}"
            )

    def test_files_are_non_empty(self):
        """Every required file has non-zero size."""
        for fname in self.REQUIRED_FILES:
            path = os.path.join(_TOY_DIR, fname)
            if os.path.isfile(path):
                self.assertGreater(
                    os.path.getsize(path), 0, f"Empty file: {fname}"
                )


class TestToyManifest(unittest.TestCase):
    """Validate the toy_manifest.json structure and content."""

    @classmethod
    def setUpClass(cls):
        with open(_MANIFEST) as fh:
            cls.manifest = json.load(fh)

    def test_has_description(self):
        """Manifest has a description field."""
        self.assertIn("description", self.manifest)

    def test_has_variants(self):
        """Manifest contains a non-empty variants list."""
        self.assertIn("variants", self.manifest)
        self.assertIsInstance(self.manifest["variants"], list)
        self.assertGreater(len(self.manifest["variants"]), 0)

    def test_variant_count_is_10(self):
        """Manifest contains exactly 10 variants (deletions)."""
        self.assertEqual(len(self.manifest["variants"]), 10)

    def test_variant_fields(self):
        """Every variant has the expected fields."""
        required_keys = {
            "chrom", "pos", "size", "genotype",
            "ref_region", "fasta_region", "hap1_regions", "hap2_regions",
        }
        for i, v in enumerate(self.manifest["variants"]):
            with self.subTest(variant=i):
                self.assertTrue(
                    required_keys.issubset(v.keys()),
                    f"Variant {i} missing keys: {required_keys - v.keys()}",
                )

    def test_variant_regions_format(self):
        """Ref, fasta, and assembly regions follow the 'chrom:start-end' format."""
        import re
        region_re = re.compile(r"^[^\s:]+:\d+-\d+$")
        for i, v in enumerate(self.manifest["variants"]):
            with self.subTest(variant=i, field="ref_region"):
                self.assertRegex(v["ref_region"], region_re)
            with self.subTest(variant=i, field="fasta_region"):
                self.assertRegex(v["fasta_region"], region_re)
            for field in ("hap1_regions", "hap2_regions"):
                for region in v[field]:
                    with self.subTest(variant=i, field=field, region=region):
                        self.assertRegex(region, region_re)

    def test_deletion_sizes_positive(self):
        """All variant sizes are positive integers (deletions)."""
        for v in self.manifest["variants"]:
            self.assertGreater(v["size"], 0)

    def test_genotypes_valid(self):
        """All genotypes are phased and contain at least one alt allele."""
        for v in self.manifest["variants"]:
            gt = v["genotype"]
            self.assertIn("|", gt, f"Genotype not phased: {gt}")
            self.assertIn("1", gt, f"No alt allele in genotype: {gt}")

    def test_chromosome_diversity(self):
        """Variants span multiple chromosomes."""
        chroms = {v["chrom"] for v in self.manifest["variants"]}
        self.assertGreaterEqual(len(chroms), 10)

    def test_ref_regions_match_chromosome(self):
        """ref_region chromosome matches the variant chrom field."""
        for v in self.manifest["variants"]:
            region_chrom = v["ref_region"].split(":")[0]
            self.assertEqual(region_chrom, v["chrom"])

    def test_ref_region_is_sv_coords(self):
        """ref_region should be actual SV coords (not padded) and match pos/size."""
        for v in self.manifest["variants"]:
            parts = v["ref_region"].split(":")[1].split("-")
            ref_start = int(parts[0])
            ref_end = int(parts[1])
            self.assertEqual(ref_start, v["pos"],
                             f"ref_region start should equal pos for {v['chrom']}")
            self.assertEqual(ref_end, v["pos"] + v["size"],
                             f"ref_region end should equal pos+size for {v['chrom']}")

    def test_fasta_region_is_padded(self):
        """fasta_region should be the padded FASTA sequence name (larger than SV)."""
        for v in self.manifest["variants"]:
            fasta = v["fasta_region"]
            fasta_chrom = fasta.split(":")[0]
            self.assertEqual(fasta_chrom, v["chrom"],
                             "fasta_region chromosome should match variant chrom")
            parts = fasta.split(":")[1].split("-")
            fasta_start = int(parts[0])
            fasta_end = int(parts[1])
            fasta_span = fasta_end - fasta_start
            self.assertGreater(fasta_span, v["size"],
                               "fasta_region should be larger than the SV size")


class TestToyFastaIndices(unittest.TestCase):
    """Validate that FASTA index (.fai) files have reasonable content."""

    FAI_FILES = [
        "toy_reference.fa.gz.fai",
        "toy_hap1.fa.gz.fai",
        "toy_hap2.fa.gz.fai",
    ]

    def test_fai_has_entries(self):
        """Each .fai file has at least one contig entry."""
        for fname in self.FAI_FILES:
            path = os.path.join(_TOY_DIR, fname)
            if not os.path.isfile(path):
                self.skipTest(f"{fname} not found")
            with open(path) as fh:
                lines = [l.strip() for l in fh if l.strip()]
            self.assertGreater(
                len(lines), 0, f"No entries in {fname}"
            )

    def test_fai_tab_separated(self):
        """FAI entries are tab-separated with at least 5 columns."""
        for fname in self.FAI_FILES:
            path = os.path.join(_TOY_DIR, fname)
            if not os.path.isfile(path):
                continue
            with open(path) as fh:
                for lineno, line in enumerate(fh, 1):
                    cols = line.strip().split("\t")
                    self.assertGreaterEqual(
                        len(cols), 5,
                        f"{fname}:{lineno} has only {len(cols)} columns",
                    )


class TestToyReadsIndexing(unittest.TestCase):
    """Validate toy reads BAM/BAI are non-empty and indexable by range."""

    def test_toy_reads_bam_contains_alignments(self):
        """toy_reads.bam should include at least one mapped alignment."""
        bam_path = os.path.join(_TOY_DIR, "toy_reads.bam")
        with gzip.open(bam_path, "rb") as fh:
            data = fh.read()

        self.assertGreaterEqual(len(data), 8, "BAM payload too small")
        self.assertEqual(data[:4], b"BAM\x01")

        p = 4
        l_text, = struct.unpack_from("<i", data, p)
        p += 4 + l_text
        n_ref, = struct.unpack_from("<i", data, p)
        p += 4
        for _ in range(n_ref):
            l_name, = struct.unpack_from("<i", data, p)
            p += 4 + l_name + 4  # name bytes + l_ref

        mapped = 0
        while p + 4 <= len(data):
            block_size, = struct.unpack_from("<i", data, p)
            p += 4
            if block_size <= 0 or p + block_size > len(data):
                break

            rec = data[p:p + block_size]
            p += block_size

            _, _, _, _, _, _, flag, _, _, _, _ = struct.unpack_from(
                "<iiBBHHHiiii", rec, 0
            )
            if (flag & 0x4) == 0:
                mapped += 1

        self.assertGreater(
            mapped, 0, "toy_reads.bam has no mapped alignments"
        )

    def test_toy_reads_bai_has_nonzero_virtual_offsets(self):
        """toy_reads.bam.bai should not be a degenerate all-zero index."""
        bam_path = os.path.join(_TOY_DIR, "toy_reads.bam")
        bai_path = os.path.join(_TOY_DIR, "toy_reads.bam.bai")
        with open(bai_path, "rb") as fh:
            data = fh.read()

        self.assertGreaterEqual(len(data), 8, "BAI payload too small")
        self.assertEqual(data[:4], b"BAI\x01")

        p = 4
        n_ref, = struct.unpack_from("<i", data, p)
        p += 4
        max_vo = 0
        for _ in range(n_ref):
            n_bin, = struct.unpack_from("<i", data, p)
            p += 4
            for _ in range(n_bin):
                _, n_chunk = struct.unpack_from("<II", data, p)
                p += 8
                for _ in range(n_chunk):
                    cs, ce = struct.unpack_from("<QQ", data, p)
                    p += 16
                    max_vo = max(max_vo, cs, ce)

            n_intv, = struct.unpack_from("<i", data, p)
            p += 4
            for _ in range(n_intv):
                io, = struct.unpack_from("<Q", data, p)
                p += 8
                max_vo = max(max_vo, io)

        self.assertGreater(
            max_vo, 0, "toy_reads.bam.bai has only zero virtual offsets"
        )
        self.assertLess(
            max_vo >> 16,
            os.path.getsize(bam_path),
            "toy_reads.bam.bai virtual offsets exceed BAM file size",
        )


if __name__ == "__main__":
    unittest.main()
