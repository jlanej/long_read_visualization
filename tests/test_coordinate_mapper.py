#!/usr/bin/env python3
"""Tests for coordinate_mapper.py"""

import gzip
import json
import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import coordinate_mapper  # noqa: E402

# Minimal PAF data (tab-separated, 12+ columns):
# query_name query_len query_start query_end strand target_name target_len
# target_start target_end matches block_len mapq
SAMPLE_PAF = """\
asm_chr1\t5000000\t1000\t4000\t+\tchr1\t248956422\t10000\t13000\t2900\t3000\t60
asm_chr1\t5000000\t5000\t8000\t-\tchr1\t248956422\t20000\t23000\t2800\t3000\t50
asm_chr2\t3000000\t0\t2000\t+\tchr2\t242193529\t50000\t52000\t1900\t2000\t40
"""


class TestParsePaf(unittest.TestCase):
    def test_parse_basic(self):
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".paf", delete=False
        ) as fh:
            fh.write(SAMPLE_PAF)
            paf_path = fh.name

        try:
            blocks = coordinate_mapper.parse_paf(paf_path)
            self.assertEqual(len(blocks), 3)

            b0 = blocks[0]
            self.assertEqual(b0["ref_chrom"], "chr1")
            self.assertEqual(b0["ref_start"], 10000)
            self.assertEqual(b0["ref_end"], 13000)
            self.assertEqual(b0["asm_chrom"], "asm_chr1")
            self.assertEqual(b0["asm_start"], 1000)
            self.assertEqual(b0["asm_end"], 4000)
            self.assertEqual(b0["strand"], "+")
            self.assertEqual(b0["mapq"], 60)
        finally:
            os.unlink(paf_path)

    def test_parse_gzipped(self):
        with tempfile.NamedTemporaryFile(
            suffix=".paf.gz", delete=False
        ) as fh:
            fh.write(gzip.compress(SAMPLE_PAF.encode()))
            paf_path = fh.name

        try:
            blocks = coordinate_mapper.parse_paf(paf_path)
            self.assertEqual(len(blocks), 3)
        finally:
            os.unlink(paf_path)

    def test_parse_empty(self):
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".paf", delete=False
        ) as fh:
            fh.write("")
            paf_path = fh.name

        try:
            blocks = coordinate_mapper.parse_paf(paf_path)
            self.assertEqual(len(blocks), 0)
        finally:
            os.unlink(paf_path)

    def test_parse_comments_and_short_lines(self):
        content = "# comment line\nshort\tline\n" + SAMPLE_PAF
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".paf", delete=False
        ) as fh:
            fh.write(content)
            paf_path = fh.name

        try:
            blocks = coordinate_mapper.parse_paf(paf_path)
            self.assertEqual(len(blocks), 3)
        finally:
            os.unlink(paf_path)


class TestBuildIndex(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        paf_path = os.path.join(self.tmpdir, "test.paf")
        with open(paf_path, "w") as fh:
            fh.write(SAMPLE_PAF)

        self.blocks = coordinate_mapper.parse_paf(paf_path)
        self.prefix = os.path.join(self.tmpdir, "test")

    def tearDown(self):
        import shutil

        shutil.rmtree(self.tmpdir)

    def test_build_bed(self):
        bed_path = coordinate_mapper.build_mapping_bed(
            self.blocks, f"{self.prefix}.mapping.bed"
        )
        self.assertTrue(os.path.exists(bed_path))
        with open(bed_path) as fh:
            lines = fh.readlines()
        self.assertEqual(len(lines), 3)
        # Verify sorted by ref_chrom, ref_start
        chroms_starts = []
        for line in lines:
            fields = line.strip().split("\t")
            chroms_starts.append((fields[0], int(fields[1])))
        self.assertEqual(
            chroms_starts,
            sorted(chroms_starts, key=lambda x: (x[0], x[1])),
        )

    def test_build_json_index(self):
        json_path = coordinate_mapper.build_json_index(
            self.blocks, f"{self.prefix}.mapping.json.gz"
        )
        self.assertTrue(os.path.exists(json_path))

        with gzip.open(json_path, "rt") as fh:
            data = json.load(fh)
        self.assertIn("chr1", data)
        self.assertIn("chr2", data)
        self.assertEqual(len(data["chr1"]), 2)
        self.assertEqual(len(data["chr2"]), 1)


class TestQuery(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        paf_path = os.path.join(self.tmpdir, "test.paf")
        with open(paf_path, "w") as fh:
            fh.write(SAMPLE_PAF)

        blocks = coordinate_mapper.parse_paf(paf_path)
        json_path = os.path.join(self.tmpdir, "test.mapping.json.gz")
        coordinate_mapper.build_json_index(blocks, json_path)
        self.index = coordinate_mapper.load_index(json_path)

    def tearDown(self):
        import shutil

        shutil.rmtree(self.tmpdir)

    def test_exact_overlap_plus_strand(self):
        """Query exactly matching a + strand alignment block."""
        results = coordinate_mapper.query(self.index, "chr1", 10000, 13000)
        self.assertEqual(len(results), 1)
        r = results[0]
        self.assertEqual(r["asm_chrom"], "asm_chr1")
        self.assertEqual(r["asm_start"], 1000)
        self.assertEqual(r["asm_end"], 4000)
        self.assertEqual(r["strand"], "+")

    def test_partial_overlap_plus_strand(self):
        """Query partially overlapping a + strand block."""
        results = coordinate_mapper.query(self.index, "chr1", 11000, 12000)
        self.assertEqual(len(results), 1)
        r = results[0]
        # ref 11000-12000 is offset 1000-2000 into the block starting at ref 10000
        self.assertEqual(r["asm_start"], 2000)  # 1000 + 1000
        self.assertEqual(r["asm_end"], 3000)  # 1000 + 2000

    def test_minus_strand(self):
        """Query overlapping a - strand alignment block."""
        # Block: ref 20000-23000, asm 5000-8000 (-)
        results = coordinate_mapper.query(self.index, "chr1", 20000, 23000)
        self.assertEqual(len(results), 1)
        r = results[0]
        self.assertEqual(r["strand"], "-")
        # Full block: asm_end - offset_end to asm_end - offset_start
        # offset_start=0, offset_end=3000, so asm: 8000-3000=5000 to 8000-0=8000
        self.assertEqual(r["asm_start"], 5000)
        self.assertEqual(r["asm_end"], 8000)

    def test_partial_minus_strand(self):
        """Query partially overlapping a - strand block."""
        # Block: ref 20000-23000, asm 5000-8000 (-)
        # Query ref 21000-22000 → offset 1000-2000
        # On - strand: asm_end - 2000 = 6000, asm_end - 1000 = 7000
        results = coordinate_mapper.query(self.index, "chr1", 21000, 22000)
        self.assertEqual(len(results), 1)
        r = results[0]
        self.assertEqual(r["asm_start"], 6000)
        self.assertEqual(r["asm_end"], 7000)

    def test_no_overlap(self):
        """Query a region with no alignments."""
        results = coordinate_mapper.query(self.index, "chr1", 0, 5000)
        self.assertEqual(len(results), 0)

    def test_unknown_chrom(self):
        """Query a chromosome not in the index."""
        results = coordinate_mapper.query(self.index, "chrX", 0, 1000000)
        self.assertEqual(len(results), 0)

    def test_multi_block_overlap(self):
        """Query spanning two blocks on the same chromosome."""
        results = coordinate_mapper.query(self.index, "chr1", 10000, 23000)
        self.assertEqual(len(results), 2)

    def test_chr2_query(self):
        """Query chr2 alignment."""
        results = coordinate_mapper.query(self.index, "chr2", 50000, 52000)
        self.assertEqual(len(results), 1)
        r = results[0]
        self.assertEqual(r["asm_chrom"], "asm_chr2")
        self.assertEqual(r["asm_start"], 0)
        self.assertEqual(r["asm_end"], 2000)


class TestParseRegion(unittest.TestCase):
    def test_valid_region(self):
        chrom, start, end = coordinate_mapper.parse_region("chr1:1000-2000")
        self.assertEqual(chrom, "chr1")
        self.assertEqual(start, 1000)
        self.assertEqual(end, 2000)

    def test_invalid_no_colon(self):
        with self.assertRaises(ValueError):
            coordinate_mapper.parse_region("chr1")

    def test_invalid_no_dash(self):
        with self.assertRaises(ValueError):
            coordinate_mapper.parse_region("chr1:1000")


class TestRoundTrip(unittest.TestCase):
    """Test full build → load → query round trip."""

    def test_roundtrip(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            paf_path = os.path.join(tmpdir, "test.paf")
            with open(paf_path, "w") as fh:
                fh.write(SAMPLE_PAF)

            blocks = coordinate_mapper.parse_paf(paf_path)
            prefix = os.path.join(tmpdir, "out")
            coordinate_mapper.build_mapping_bed(blocks, f"{prefix}.mapping.bed")
            coordinate_mapper.build_json_index(
                blocks, f"{prefix}.mapping.json.gz"
            )

            index = coordinate_mapper.load_index(f"{prefix}.mapping.json.gz")
            results = coordinate_mapper.query(index, "chr1", 10000, 13000)
            self.assertEqual(len(results), 1)
            self.assertEqual(results[0]["asm_chrom"], "asm_chr1")


if __name__ == "__main__":
    unittest.main()
