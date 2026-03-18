#!/usr/bin/env python3
"""Tests for coordinate_mapper.py"""

import gzip
import json
import os
import shutil
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
        self.assertEqual(r["event_type"], "alignment")

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
        """Query spanning two blocks on the same chromosome.

        The gap between the blocks (ref 13000-20000) falls within the query
        window (ref 10000-23000). The two blocks are on opposite strands, so
        the gap is classified as an inversion event.  Results: two alignment
        records flanking one inversion gap record = 3 total.
        """
        results = coordinate_mapper.query(self.index, "chr1", 10000, 23000)
        self.assertEqual(len(results), 3)
        event_types = [r["event_type"] for r in results]
        self.assertEqual(event_types, ["alignment", "inversion", "alignment"])

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


# ---------------------------------------------------------------------------
# CIGAR parsing tests
# ---------------------------------------------------------------------------

class TestParseCigar(unittest.TestCase):
    def test_empty(self):
        self.assertEqual(coordinate_mapper.parse_cigar(""), [])
        self.assertEqual(coordinate_mapper.parse_cigar(None), [])

    def test_simple_match(self):
        ops = coordinate_mapper.parse_cigar("100M")
        self.assertEqual(ops, [(100, "M")])

    def test_complex_cigar(self):
        ops = coordinate_mapper.parse_cigar("50M5D30M10I20M")
        self.assertEqual(ops, [(50, "M"), (5, "D"), (30, "M"), (10, "I"), (20, "M")])

    def test_extended_cigar(self):
        ops = coordinate_mapper.parse_cigar("10=5X3=")
        self.assertEqual(ops, [(10, "="), (5, "X"), (3, "=")])


# ---------------------------------------------------------------------------
# CIGAR projection tests
# ---------------------------------------------------------------------------

# PAF with CIGAR strings for projection tests:
# Block A: ref chr1 [10000, 10200), asm_chr1 [1000, 1195), +
#   CIGAR: 100M5D95M
#   - ref[10000:10100] → asm[1000:1100]  (100M)
#   - ref[10100:10105] deleted in asm    (5D)
#   - ref[10105:10200] → asm[1100:1195] (95M)
#
# Block B: ref chr1 [20000, 20200), asm_chr1 [2000, 2210), +
#   CIGAR: 100M10I100M
#   - ref[20000:20100] → asm[2000:2100]  (100M)
#   - asm[2100:2110] inserted             (10I)
#   - ref[20100:20200] → asm[2110:2210] (100M)
#
# Block C: ref chr2 [50000, 50200), asm_chr2 [3000, 3195), -
#   CIGAR: 100M5D95M  (same ops, but minus strand)
#   - ref[50000:50100] → asm[3195:3095] (reversed, 100M on - strand)
#   - ref[50100:50105] deleted in asm (5D, no asm advance)
#   - ref[50105:50200] → asm[3095:3000] (95M)

CIGAR_PAF = """\
asm_chr1\t5000000\t1000\t1195\t+\tchr1\t248956422\t10000\t10200\t195\t200\t60\tcg:Z:100M5D95M
asm_chr1\t5000000\t2000\t2210\t+\tchr1\t248956422\t20000\t20200\t200\t210\t60\tcg:Z:100M10I100M
asm_chr2\t3000000\t3000\t3195\t-\tchr2\t242193529\t50000\t50200\t195\t200\t50\tcg:Z:100M5D95M
"""


class TestProjectCigar(unittest.TestCase):
    """Unit tests for project_cigar()."""

    def test_pure_match_plus_full(self):
        """100M block: full block maps 1-1."""
        ops = coordinate_mapper.parse_cigar("100M")
        s, e = coordinate_mapper.project_cigar(ops, 0, 100, "+", 1000, 1100)
        self.assertEqual((s, e), (1000, 1100))

    def test_pure_match_plus_partial(self):
        """100M block: partial query."""
        ops = coordinate_mapper.parse_cigar("100M")
        s, e = coordinate_mapper.project_cigar(ops, 10, 40, "+", 1000, 1100)
        self.assertEqual((s, e), (1010, 1040))

    def test_deletion_query_before(self):
        """100M5D95M: query entirely before the deletion."""
        ops = coordinate_mapper.parse_cigar("100M5D95M")
        s, e = coordinate_mapper.project_cigar(ops, 0, 50, "+", 1000, 1195)
        self.assertEqual((s, e), (1000, 1050))

    def test_deletion_query_spanning(self):
        """100M5D95M: query spanning the 5-base deletion (ref 90-110)."""
        ops = coordinate_mapper.parse_cigar("100M5D95M")
        # ref_offset 90-100 → asm 90-100 (in 100M)
        # ref_offset 100-105 → deletion, asm clamps to 100
        # ref_offset 105-110 → asm 100-105 (in 95M)
        # Total: asm_q_start=90 (from 90 in 100M), asm_q_end=105 (5 into 95M)
        s, e = coordinate_mapper.project_cigar(ops, 90, 110, "+", 1000, 1195)
        self.assertEqual((s, e), (1090, 1105))

    def test_deletion_query_after(self):
        """100M5D95M: query entirely after the deletion."""
        ops = coordinate_mapper.parse_cigar("100M5D95M")
        # ref_offset 110-150 → asm 105-145
        s, e = coordinate_mapper.project_cigar(ops, 110, 150, "+", 1000, 1195)
        self.assertEqual((s, e), (1105, 1145))

    def test_insertion_query_before(self):
        """100M10I100M: query entirely before the insertion."""
        ops = coordinate_mapper.parse_cigar("100M10I100M")
        s, e = coordinate_mapper.project_cigar(ops, 0, 50, "+", 2000, 2210)
        self.assertEqual((s, e), (2000, 2050))

    def test_insertion_query_straddling(self):
        """100M10I100M: query ref[90:110] straddles the insertion point."""
        # ref 0-100 → asm 0-100 (100M), then 10I advances asm to 110,
        # then ref 100-200 → asm 110-210 (100M).
        # ref_offset_start=90 → asm_q_start=90 (in 100M)
        # ref_offset_end=110 → after 100M (ref_cur=100, asm_cur=110 after I),
        #   10 into the second 100M → asm_q_end=110+10=120
        ops = coordinate_mapper.parse_cigar("100M10I100M")
        s, e = coordinate_mapper.project_cigar(ops, 90, 110, "+", 2000, 2210)
        self.assertEqual((s, e), (2090, 2120))

    def test_insertion_query_after(self):
        """100M10I100M: query entirely after the insertion."""
        ops = coordinate_mapper.parse_cigar("100M10I100M")
        # ref_offset 110-150 → asm_offset 120-160
        s, e = coordinate_mapper.project_cigar(ops, 110, 150, "+", 2000, 2210)
        self.assertEqual((s, e), (2120, 2160))

    def test_minus_strand_full(self):
        """100M5D95M on − strand: full block."""
        ops = coordinate_mapper.parse_cigar("100M5D95M")
        # Block asm [3000, 3195) on − strand; total asm consumed = 195.
        # Full ref range [0, 200): asm_q_start=0, asm_q_end=195
        # → asm: [3195-195, 3195-0] = [3000, 3195]
        s, e = coordinate_mapper.project_cigar(ops, 0, 200, "-", 3000, 3195)
        self.assertEqual((s, e), (3000, 3195))

    def test_minus_strand_before_deletion(self):
        """100M5D95M on − strand: query before deletion."""
        ops = coordinate_mapper.parse_cigar("100M5D95M")
        # ref_offset [0, 50) → asm_q_start=0, asm_q_end=50
        # On − strand: [3195-50, 3195-0] = [3145, 3195]
        s, e = coordinate_mapper.project_cigar(ops, 0, 50, "-", 3000, 3195)
        self.assertEqual((s, e), (3145, 3195))

    def test_minus_strand_spanning_deletion(self):
        """100M5D95M on − strand: query spanning the deletion."""
        # ref_offset [90, 110): same asm offsets as + strand test → [90, 105]
        # On − strand: [3195-105, 3195-90] = [3090, 3105]
        ops = coordinate_mapper.parse_cigar("100M5D95M")
        s, e = coordinate_mapper.project_cigar(ops, 90, 110, "-", 3000, 3195)
        self.assertEqual((s, e), (3090, 3105))


class TestCigarRoundTrip(unittest.TestCase):
    """Build index from PAF with CIGAR and query it."""

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        paf_path = os.path.join(self.tmpdir, "cigar.paf")
        with open(paf_path, "w") as fh:
            fh.write(CIGAR_PAF)

        blocks = coordinate_mapper.parse_paf(paf_path)
        json_path = os.path.join(self.tmpdir, "cigar.mapping.json.gz")
        coordinate_mapper.build_json_index(blocks, json_path)
        self.index = coordinate_mapper.load_index(json_path)

    def tearDown(self):
        shutil.rmtree(self.tmpdir)

    def test_cigar_stored_in_index(self):
        """CIGAR strings must be stored in the JSON index."""
        block = self.index["chr1"]["blocks"][0]
        self.assertIn("cg", block)
        self.assertEqual(block["cg"], "100M5D95M")

    def test_plus_strand_deletion_cigar_query(self):
        """Query spanning deletion: asm range correctly contracts."""
        # ref[10090:10110] spans the 5D; asm range should be [1090, 1105]
        results = coordinate_mapper.query(self.index, "chr1", 10090, 10110)
        self.assertEqual(len(results), 1)
        r = results[0]
        self.assertEqual(r["event_type"], "alignment")
        self.assertEqual(r["asm_start"], 1090)
        self.assertEqual(r["asm_end"], 1105)

    def test_plus_strand_insertion_cigar_query(self):
        """Query straddling insertion: asm range includes inserted sequence."""
        # ref[20090:20110] straddles the 10I; asm range should be [2090, 2120]
        results = coordinate_mapper.query(self.index, "chr1", 20090, 20110)
        self.assertEqual(len(results), 1)
        r = results[0]
        self.assertEqual(r["event_type"], "alignment")
        self.assertEqual(r["asm_start"], 2090)
        self.assertEqual(r["asm_end"], 2120)

    def test_minus_strand_deletion_cigar_query(self):
        """− strand block with deletion: assembly coordinates correctly reversed."""
        # ref[50090:50110] spans deletion; asm range should be [3090, 3105]
        results = coordinate_mapper.query(self.index, "chr2", 50090, 50110)
        self.assertEqual(len(results), 1)
        r = results[0]
        self.assertEqual(r["event_type"], "alignment")
        self.assertEqual(r["asm_start"], 3090)
        self.assertEqual(r["asm_end"], 3105)


# ---------------------------------------------------------------------------
# SV gap classification tests
# ---------------------------------------------------------------------------

# PAF designed to produce recognisable inter-block gaps:
#
#  chr3 blocks for deletion / insertion / inversion / translocation tests:
#
#  del_block_1:  ref chr3 [1000, 2000), asm_chr3 [10000, 11000), +
#  del_block_2:  ref chr3 [5000, 6000), asm_chr3 [11000, 12000), +
#    → ref gap 3000 bp, asm gap 0 bp → DELETION
#
#  ins_block_1:  ref chr3 [10000, 11000), asm_chr3 [20000, 21000), +
#  ins_block_2:  ref chr3 [11000, 12000), asm_chr3 [31000, 32000), +
#    → ref gap 0 bp, asm gap 10000 bp → INSERTION
#
#  inv_block_1:  ref chr3 [20000, 21000), asm_chr3 [40000, 41000), +
#  inv_block_2:  ref chr3 [21000, 22000), asm_chr3 [41000, 42000), -
#    → opposite strands → INVERSION
#
#  trl_block_1:  ref chr3 [30000, 31000), asm_chrA [50000, 51000), +
#  trl_block_2:  ref chr3 [32000, 33000), asm_chrB [60000, 61000), +
#    → different asm contigs → TRANSLOCATION

SV_GAP_PAF = """\
asm_chr3\t5000000\t10000\t11000\t+\tchr3\t200000000\t1000\t2000\t1000\t1000\t60
asm_chr3\t5000000\t11000\t12000\t+\tchr3\t200000000\t5000\t6000\t1000\t1000\t60
asm_chr3\t5000000\t20000\t21000\t+\tchr3\t200000000\t10000\t11000\t1000\t1000\t60
asm_chr3\t5000000\t31000\t32000\t+\tchr3\t200000000\t11000\t12000\t1000\t1000\t60
asm_chr3\t5000000\t40000\t41000\t+\tchr3\t200000000\t20000\t21000\t1000\t1000\t60
asm_chr3\t5000000\t41000\t42000\t-\tchr3\t200000000\t21000\t22000\t1000\t1000\t60
asm_chrA\t5000000\t50000\t51000\t+\tchr3\t200000000\t30000\t31000\t1000\t1000\t60
asm_chrB\t5000000\t60000\t61000\t+\tchr3\t200000000\t32000\t33000\t1000\t1000\t60
"""


class TestSvGapClassification(unittest.TestCase):
    """Tests for gap detection and SV event classification in query()."""

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        paf_path = os.path.join(self.tmpdir, "sv.paf")
        with open(paf_path, "w") as fh:
            fh.write(SV_GAP_PAF)

        blocks = coordinate_mapper.parse_paf(paf_path)
        json_path = os.path.join(self.tmpdir, "sv.mapping.json.gz")
        coordinate_mapper.build_json_index(blocks, json_path)
        self.index = coordinate_mapper.load_index(json_path)

    def tearDown(self):
        shutil.rmtree(self.tmpdir)

    def test_deletion_event(self):
        """Gap where ref > asm → deletion in assembly."""
        results = coordinate_mapper.query(self.index, "chr3", 500, 6500)
        event_types = [r["event_type"] for r in results]
        self.assertIn("deletion", event_types)
        gap = next(r for r in results if r["event_type"] == "deletion")
        # ref gap: [2000, 5000] = 3000 bp; asm gap: 0 bp
        self.assertEqual(gap["ref_gap_size"], 3000)
        self.assertEqual(gap["asm_gap_size"], 0)

    def test_insertion_event(self):
        """Gap where asm > ref → insertion in assembly."""
        results = coordinate_mapper.query(self.index, "chr3", 9500, 12500)
        event_types = [r["event_type"] for r in results]
        self.assertIn("insertion", event_types)
        gap = next(r for r in results if r["event_type"] == "insertion")
        # ref gap: 0 bp; asm gap: 10000 bp
        self.assertEqual(gap["ref_gap_size"], 0)
        self.assertEqual(gap["asm_gap_size"], 10000)

    def test_inversion_event(self):
        """Opposite-strand flanking blocks → inversion."""
        results = coordinate_mapper.query(self.index, "chr3", 19500, 22500)
        event_types = [r["event_type"] for r in results]
        self.assertIn("inversion", event_types)

    def test_translocation_event(self):
        """Different assembly contigs → translocation."""
        results = coordinate_mapper.query(self.index, "chr3", 29500, 33500)
        event_types = [r["event_type"] for r in results]
        self.assertIn("translocation", event_types)

    def test_alignment_events_present(self):
        """Alignment events are still returned alongside gap events."""
        results = coordinate_mapper.query(self.index, "chr3", 500, 6500)
        alignment_results = [r for r in results if r["event_type"] == "alignment"]
        self.assertEqual(len(alignment_results), 2)

    def test_gap_event_fields(self):
        """Gap events carry all expected fields."""
        results = coordinate_mapper.query(self.index, "chr3", 500, 6500)
        gap = next(r for r in results if r["event_type"] != "alignment")
        for field in ("ref_chrom", "ref_start", "ref_end",
                      "asm_chrom", "asm_start", "asm_end",
                      "strand", "mapq", "event_type",
                      "ref_gap_size", "asm_gap_size"):
            self.assertIn(field, gap)

    def test_deletion_asm_coordinates(self):
        """Deletion gap assembly coordinates span zero bases (asm_start == asm_end)."""
        results = coordinate_mapper.query(self.index, "chr3", 500, 6500)
        gap = next(r for r in results if r["event_type"] == "deletion")
        self.assertEqual(gap["asm_start"], 11000)
        self.assertEqual(gap["asm_end"], 11000)

    def test_insertion_asm_coordinates(self):
        """Insertion gap assembly coordinates span the inserted sequence."""
        results = coordinate_mapper.query(self.index, "chr3", 9500, 12500)
        gap = next(r for r in results if r["event_type"] == "insertion")
        self.assertEqual(gap["asm_start"], 21000)
        self.assertEqual(gap["asm_end"], 31000)

    def test_translocation_ref_gap_and_mapq(self):
        """Translocation gap carries correct ref coordinates and min mapq."""
        results = coordinate_mapper.query(self.index, "chr3", 29500, 33500)
        gap = next(r for r in results if r["event_type"] == "translocation")
        self.assertEqual(gap["ref_start"], 31000)
        self.assertEqual(gap["ref_end"], 32000)
        self.assertEqual(gap["mapq"], 60)


# =========================================================================
# Biologically realistic SV scenarios
# =========================================================================
#
# These PAF records model real structural variant events as produced by
# minimap2 -x asm5 --eqx -c against diploid assemblies.  Each scenario
# uses realistic chromosome coordinates, contig names, and CIGAR strings
# that would be observed in HPRC-quality assemblies.
# =========================================================================

# ---------------------------------------------------------------------------
# Scenario 1 – Heterozygous ~5 kb deletion
# ---------------------------------------------------------------------------
# The assembly is missing 5 kb relative to the reference.  minimap2 produces
# two collinear blocks flanking the deletion.
#
#   Block A: ref chr1 [1000000, 1005000), asm ctg1_hap1 [500000, 505000), +
#     CIGAR: 5000=  (perfect match)
#   Block B: ref chr1 [1010000, 1015000), asm ctg1_hap1 [505000, 510000), +
#     CIGAR: 5000=
#   → ref gap = 5000 bp, asm gap = 0 bp → deletion
#
# Scenario 2 – Heterozygous ~3 kb insertion
# ---------------------------------------------------------------------------
# The assembly has 3 kb of novel sequence not in reference.  Adjacent blocks
# in ref space, but a 3 kb gap in assembly space.
#
#   Block C: ref chr1 [2000000, 2005000), asm ctg1_hap1 [600000, 605000), +
#     CIGAR: 5000=
#   Block D: ref chr1 [2005000, 2010000), asm ctg1_hap1 [608000, 613000), +
#     CIGAR: 5000=
#   → ref gap = 0 bp, asm gap = 3000 bp → insertion
#
# Scenario 3 – 10 kb inversion
# ---------------------------------------------------------------------------
#   Block E: ref chr2 [3000000, 3005000), asm ctg2_hap1 [700000, 705000), +
#     CIGAR: 5000=
#   Block F: ref chr2 [3005000, 3015000), asm ctg2_hap1 [705000, 715000), -
#     CIGAR: 10000=  (inverted)
#   Block G: ref chr2 [3015000, 3020000), asm ctg2_hap1 [715000, 720000), +
#     CIGAR: 5000=
#
# Scenario 4 – Tandem duplication (overlapping blocks in ref space)
# ---------------------------------------------------------------------------
#   Block H: ref chr3 [4000000, 4002000), asm ctg3_hap1 [800000, 802000), +
#     CIGAR: 2000=
#   Block I: ref chr3 [4000000, 4002000), asm ctg3_hap1 [802000, 804000), +
#     CIGAR: 2000=
#   Both map the same ref region to different asm positions → tandem dup
#
# Scenario 5 – Minus-strand deletion
# ---------------------------------------------------------------------------
#   Block J: ref chr5 [5000000, 5003000), asm ctg5_hap2 [900000, 903000), -
#     CIGAR: 3000=
#   Block K: ref chr5 [5008000, 5011000), asm ctg5_hap2 [895000, 898000), -
#     CIGAR: 3000=
#   ref gap=5000bp, asm gap = prev["as"]-next["ae"] = 900000-898000 = 2000bp
#   → deletion (asm_gap < ref_gap)
#
# Scenario 6 – Minus-strand insertion (ref-adjacent blocks)
# ---------------------------------------------------------------------------
#   Block L: ref chr5 [6000000, 6003000), asm ctg5_hap2 [1000000, 1003000), -
#     CIGAR: 3000=
#   Block M: ref chr5 [6003000, 6006000), asm ctg5_hap2 [993000, 996000), -
#     CIGAR: 3000=
#   ref gap=0bp, asm gap = 1000000-996000 = 4000bp → insertion
#
# Scenario 7 – Inter-contig translocation with CIGAR
# ---------------------------------------------------------------------------
#   Block N: ref chr6 [7000000, 7005000), asm ctg6a [100000, 105000), +
#     CIGAR: 5000=
#   Block O: ref chr6 [7006000, 7011000), asm ctg6b [200000, 205000), +
#     CIGAR: 5000=
#   Different contigs → translocation

BIOLOGICAL_SV_PAF = """\
ctg1_hap1\t10000000\t500000\t505000\t+\tchr1\t248956422\t1000000\t1005000\t5000\t5000\t60\tcg:Z:5000=
ctg1_hap1\t10000000\t505000\t510000\t+\tchr1\t248956422\t1010000\t1015000\t5000\t5000\t60\tcg:Z:5000=
ctg1_hap1\t10000000\t600000\t605000\t+\tchr1\t248956422\t2000000\t2005000\t5000\t5000\t60\tcg:Z:5000=
ctg1_hap1\t10000000\t608000\t613000\t+\tchr1\t248956422\t2005000\t2010000\t5000\t5000\t60\tcg:Z:5000=
ctg2_hap1\t10000000\t700000\t705000\t+\tchr2\t242193529\t3000000\t3005000\t5000\t5000\t60\tcg:Z:5000=
ctg2_hap1\t10000000\t705000\t715000\t-\tchr2\t242193529\t3005000\t3015000\t10000\t10000\t60\tcg:Z:10000=
ctg2_hap1\t10000000\t715000\t720000\t+\tchr2\t242193529\t3015000\t3020000\t5000\t5000\t60\tcg:Z:5000=
ctg3_hap1\t10000000\t800000\t802000\t+\tchr3\t198295559\t4000000\t4002000\t2000\t2000\t60\tcg:Z:2000=
ctg3_hap1\t10000000\t802000\t804000\t+\tchr3\t198295559\t4000000\t4002000\t2000\t2000\t60\tcg:Z:2000=
ctg5_hap2\t10000000\t900000\t903000\t-\tchr5\t181538259\t5000000\t5003000\t3000\t3000\t60\tcg:Z:3000=
ctg5_hap2\t10000000\t895000\t898000\t-\tchr5\t181538259\t5008000\t5011000\t3000\t3000\t60\tcg:Z:3000=
ctg5_hap2\t10000000\t1000000\t1003000\t-\tchr5\t181538259\t6000000\t6003000\t3000\t3000\t60\tcg:Z:3000=
ctg5_hap2\t10000000\t993000\t996000\t-\tchr5\t181538259\t6003000\t6006000\t3000\t3000\t60\tcg:Z:3000=
ctg6a\t10000000\t100000\t105000\t+\tchr6\t170805979\t7000000\t7005000\t5000\t5000\t60\tcg:Z:5000=
ctg6b\t10000000\t200000\t205000\t+\tchr6\t170805979\t7006000\t7011000\t5000\t5000\t60\tcg:Z:5000=
"""


class TestBiologicalSVScenarios(unittest.TestCase):
    """Tests modelling biologically realistic SV events.

    Each test verifies that the mapper correctly identifies the SV event
    type, reports accurate reference and assembly coordinates, and computes
    correct gap sizes — ensuring the mapping is suitable for real genome
    visualization at SV loci.
    """

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        paf_path = os.path.join(self.tmpdir, "bio.paf")
        with open(paf_path, "w") as fh:
            fh.write(BIOLOGICAL_SV_PAF)

        blocks = coordinate_mapper.parse_paf(paf_path)
        json_path = os.path.join(self.tmpdir, "bio.mapping.json.gz")
        coordinate_mapper.build_json_index(blocks, json_path)
        self.index = coordinate_mapper.load_index(json_path)

    def tearDown(self):
        shutil.rmtree(self.tmpdir)

    # -- Scenario 1: Heterozygous deletion (5 kb) --------------------------

    def test_het_deletion_event_detected(self):
        """5 kb deletion detected as deletion event between flanking blocks."""
        results = coordinate_mapper.query(
            self.index, "chr1", 999000, 1016000
        )
        types = [r["event_type"] for r in results]
        self.assertEqual(types, ["alignment", "deletion", "alignment"])

    def test_het_deletion_gap_sizes(self):
        """Deletion: 5 kb ref gap, 0 asm gap."""
        results = coordinate_mapper.query(
            self.index, "chr1", 999000, 1016000
        )
        gap = next(r for r in results if r["event_type"] == "deletion")
        self.assertEqual(gap["ref_gap_size"], 5000)
        self.assertEqual(gap["asm_gap_size"], 0)

    def test_het_deletion_asm_coords(self):
        """Deletion: assembly coords are a zero-width point at the junction."""
        results = coordinate_mapper.query(
            self.index, "chr1", 999000, 1016000
        )
        gap = next(r for r in results if r["event_type"] == "deletion")
        self.assertEqual(gap["asm_start"], 505000)
        self.assertEqual(gap["asm_end"], 505000)

    def test_het_deletion_flanking_alignments(self):
        """Deletion: flanking alignments map correctly with CIGAR."""
        results = coordinate_mapper.query(
            self.index, "chr1", 999000, 1016000
        )
        alns = [r for r in results if r["event_type"] == "alignment"]
        self.assertEqual(len(alns), 2)
        # Block A: ref [1000000, 1005000) → asm [500000, 505000) +
        self.assertEqual(alns[0]["ref_start"], 1000000)
        self.assertEqual(alns[0]["ref_end"], 1005000)
        self.assertEqual(alns[0]["asm_start"], 500000)
        self.assertEqual(alns[0]["asm_end"], 505000)
        # Block B: ref [1010000, 1015000) → asm [505000, 510000) +
        self.assertEqual(alns[1]["ref_start"], 1010000)
        self.assertEqual(alns[1]["ref_end"], 1015000)
        self.assertEqual(alns[1]["asm_start"], 505000)
        self.assertEqual(alns[1]["asm_end"], 510000)

    # -- Scenario 2: Heterozygous insertion (3 kb) -------------------------

    def test_het_insertion_event_detected(self):
        """3 kb insertion detected as insertion event at ref-adjacent junction."""
        results = coordinate_mapper.query(
            self.index, "chr1", 1999000, 2011000
        )
        types = [r["event_type"] for r in results]
        self.assertEqual(types, ["alignment", "insertion", "alignment"])

    def test_het_insertion_gap_sizes(self):
        """Insertion: 0 ref gap, 3000 asm gap."""
        results = coordinate_mapper.query(
            self.index, "chr1", 1999000, 2011000
        )
        gap = next(r for r in results if r["event_type"] == "insertion")
        self.assertEqual(gap["ref_gap_size"], 0)
        self.assertEqual(gap["asm_gap_size"], 3000)

    def test_het_insertion_asm_coords(self):
        """Insertion: assembly coordinates span the 3 kb novel sequence."""
        results = coordinate_mapper.query(
            self.index, "chr1", 1999000, 2011000
        )
        gap = next(r for r in results if r["event_type"] == "insertion")
        self.assertEqual(gap["asm_start"], 605000)
        self.assertEqual(gap["asm_end"], 608000)

    # -- Scenario 3: Inversion (10 kb inverted segment) --------------------

    def test_inversion_three_blocks(self):
        """Inversion: + block, inverted - block, + block with 2 gap events."""
        results = coordinate_mapper.query(
            self.index, "chr2", 2999000, 3021000
        )
        types = [r["event_type"] for r in results]
        self.assertEqual(
            types,
            ["alignment", "inversion", "alignment", "inversion", "alignment"],
        )

    def test_inversion_flanking_strands(self):
        """Inversion: flanking blocks are +, inverted block is -."""
        results = coordinate_mapper.query(
            self.index, "chr2", 2999000, 3021000
        )
        alns = [r for r in results if r["event_type"] == "alignment"]
        self.assertEqual(len(alns), 3)
        self.assertEqual(alns[0]["strand"], "+")
        self.assertEqual(alns[1]["strand"], "-")
        self.assertEqual(alns[2]["strand"], "+")

    def test_inversion_inverted_block_asm_coords(self):
        """Inversion: - strand block correctly projects assembly coordinates.

        Block F is ref [3005000, 3015000) → asm [705000, 715000) minus.
        Querying the full block should yield the full asm range.
        """
        results = coordinate_mapper.query(
            self.index, "chr2", 3005000, 3015000
        )
        self.assertEqual(len(results), 1)
        r = results[0]
        self.assertEqual(r["strand"], "-")
        self.assertEqual(r["asm_start"], 705000)
        self.assertEqual(r["asm_end"], 715000)

    def test_inversion_partial_query_of_inverted_block(self):
        """Inversion: partial query of inverted block reverses correctly.

        ref [3007000, 3012000) is offset [2000, 7000) into block F.
        On − strand: asm [715000-7000, 715000-2000] = [708000, 713000).
        """
        results = coordinate_mapper.query(
            self.index, "chr2", 3007000, 3012000
        )
        self.assertEqual(len(results), 1)
        r = results[0]
        self.assertEqual(r["asm_start"], 708000)
        self.assertEqual(r["asm_end"], 713000)

    # -- Scenario 4: Tandem duplication ------------------------------------

    def test_tandem_dup_two_hits(self):
        """Tandem dup: same ref region maps to two different assembly positions."""
        results = coordinate_mapper.query(
            self.index, "chr3", 4000000, 4002000
        )
        self.assertEqual(len(results), 2)
        alns = [r for r in results if r["event_type"] == "alignment"]
        self.assertEqual(len(alns), 2)
        # Both point to same ref but different asm loci
        asm_starts = sorted(r["asm_start"] for r in alns)
        self.assertEqual(asm_starts, [800000, 802000])

    def test_tandem_dup_partial_query(self):
        """Tandem dup: partial query still returns two hits."""
        results = coordinate_mapper.query(
            self.index, "chr3", 4000500, 4001500
        )
        self.assertEqual(len(results), 2)

    # -- Scenario 5: Minus-strand deletion ---------------------------------

    def test_minus_strand_deletion_event(self):
        """Deletion between two − strand blocks."""
        results = coordinate_mapper.query(
            self.index, "chr5", 4999000, 5012000
        )
        types = [r["event_type"] for r in results]
        self.assertIn("deletion", types)

    def test_minus_strand_deletion_gap_sizes(self):
        """Minus-strand deletion: ref gap=5000, asm gap=2000."""
        results = coordinate_mapper.query(
            self.index, "chr5", 4999000, 5012000
        )
        gap = next(r for r in results if r["event_type"] == "deletion")
        self.assertEqual(gap["ref_gap_size"], 5000)
        self.assertEqual(gap["asm_gap_size"], 2000)

    def test_minus_strand_deletion_asm_coords(self):
        """Minus-strand deletion: assembly coords span the shortened gap.

        prev block asm_start=900000, next block asm_end=898000.
        asm_start = next["ae"] = 898000, asm_end = prev["as"] = 900000.
        """
        results = coordinate_mapper.query(
            self.index, "chr5", 4999000, 5012000
        )
        gap = next(r for r in results if r["event_type"] == "deletion")
        self.assertEqual(gap["asm_start"], 898000)
        self.assertEqual(gap["asm_end"], 900000)

    # -- Scenario 6: Minus-strand insertion --------------------------------

    def test_minus_strand_insertion_event(self):
        """Insertion at ref-adjacent − strand blocks."""
        results = coordinate_mapper.query(
            self.index, "chr5", 5999000, 6007000
        )
        types = [r["event_type"] for r in results]
        self.assertIn("insertion", types)

    def test_minus_strand_insertion_gap_sizes(self):
        """Minus-strand insertion: ref gap=0, asm gap=4000."""
        results = coordinate_mapper.query(
            self.index, "chr5", 5999000, 6007000
        )
        gap = next(r for r in results if r["event_type"] == "insertion")
        self.assertEqual(gap["ref_gap_size"], 0)
        self.assertEqual(gap["asm_gap_size"], 4000)

    def test_minus_strand_insertion_asm_coords(self):
        """Minus-strand insertion: asm coords span the inserted sequence.

        prev["as"]=1000000, next["ae"]=996000.
        asm_start = 996000, asm_end = 1000000.
        """
        results = coordinate_mapper.query(
            self.index, "chr5", 5999000, 6007000
        )
        gap = next(r for r in results if r["event_type"] == "insertion")
        self.assertEqual(gap["asm_start"], 996000)
        self.assertEqual(gap["asm_end"], 1000000)

    # -- Scenario 7: Inter-contig translocation ----------------------------

    def test_translocation_with_cigar(self):
        """Translocation: blocks on different contigs with CIGAR alignment."""
        results = coordinate_mapper.query(
            self.index, "chr6", 6999000, 7012000
        )
        types = [r["event_type"] for r in results]
        self.assertEqual(types, ["alignment", "translocation", "alignment"])

    def test_translocation_different_contigs(self):
        """Translocation: flanking alignments point to different contigs."""
        results = coordinate_mapper.query(
            self.index, "chr6", 6999000, 7012000
        )
        alns = [r for r in results if r["event_type"] == "alignment"]
        self.assertEqual(alns[0]["asm_chrom"], "ctg6a")
        self.assertEqual(alns[1]["asm_chrom"], "ctg6b")


# =========================================================================
# CIGAR projection edge cases
# =========================================================================

# PAF records with complex CIGARs modelling various within-block events:
#
# Block 1: Large deletion inside a single block (1000=5000D1000=)
#   ref chr7 [100000, 107000), asm ctg7 [200000, 202000), +
#   ref consumed: 1000 + 5000 + 1000 = 7000
#   asm consumed: 1000 + 1000 = 2000
#
# Block 2: Multiple small indels – microsatellite expansion
#   100=3I50=2D100=5I50=100=
#   ref chr7 [200000, 200402), asm ctg7 [300000, 300408), +
#   ref consumed: 100 + 50 + 2 + 100 + 50 + 100 = 402
#   asm consumed: 100 + 3 + 50 + 100 + 5 + 50 + 100 = 408
#
# Block 3: Extended CIGAR (=/X) with SNPs and indels
#   500=10X490=5D500=10I500=
#   ref chr7 [300000, 301505), asm ctg7 [400000, 402010), +
#   ref consumed: 500 + 10 + 490 + 5 + 500 + 500 = 2005
#   asm consumed: 500 + 10 + 490 + 500 + 10 + 500 = 2010
#
# Block 4: Minus strand with insertion
#   2000=500I3000=
#   ref chr8 [400000, 405000), asm ctg8 [500000, 505500), -
#   ref consumed: 2000 + 3000 = 5000
#   asm consumed: 2000 + 500 + 3000 = 5500

CIGAR_EDGE_PAF = """\
ctg7\t10000000\t200000\t202000\t+\tchr7\t159345973\t100000\t107000\t2000\t7000\t60\tcg:Z:1000=5000D1000=
ctg7\t10000000\t300000\t300408\t+\tchr7\t159345973\t200000\t200402\t400\t410\t60\tcg:Z:100=3I50=2D100=5I50=100=
ctg7\t10000000\t400000\t402010\t+\tchr7\t159345973\t300000\t302005\t2000\t2010\t60\tcg:Z:500=10X490=5D500=10I500=
ctg8\t10000000\t500000\t505500\t-\tchr8\t145138636\t400000\t405000\t5000\t5500\t60\tcg:Z:2000=500I3000=
"""


class TestCigarProjectionEdgeCases(unittest.TestCase):
    """Edge cases for CIGAR-aware coordinate projection."""

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        paf_path = os.path.join(self.tmpdir, "edge.paf")
        with open(paf_path, "w") as fh:
            fh.write(CIGAR_EDGE_PAF)

        blocks = coordinate_mapper.parse_paf(paf_path)
        json_path = os.path.join(self.tmpdir, "edge.mapping.json.gz")
        coordinate_mapper.build_json_index(blocks, json_path)
        self.index = coordinate_mapper.load_index(json_path)

    def tearDown(self):
        shutil.rmtree(self.tmpdir)

    # -- Block 1: Large deletion (5 kb D inside one CIGAR) --

    def test_large_del_query_before(self):
        """Query before the 5 kb deletion: 1:1 mapping."""
        # ref [100000, 100500) → offset [0, 500) in 1000=
        # asm [200000, 200500)
        results = coordinate_mapper.query(self.index, "chr7", 100000, 100500)
        self.assertEqual(len(results), 1)
        r = results[0]
        self.assertEqual(r["asm_start"], 200000)
        self.assertEqual(r["asm_end"], 200500)

    def test_large_del_query_spanning(self):
        """Query spanning the 5 kb deletion: asm range contracts.

        ref [100500, 106500) → offsets [500, 6500).
        Within 1000= (offset 500-1000): asm 500-1000
        Within 5000D (offset 1000-6000): asm stays at 1000
        Within 1000= (offset 6000-6500): asm 1000-1500
        Result: asm [200500, 201500)
        """
        results = coordinate_mapper.query(self.index, "chr7", 100500, 106500)
        self.assertEqual(len(results), 1)
        r = results[0]
        self.assertEqual(r["asm_start"], 200500)
        self.assertEqual(r["asm_end"], 201500)

    def test_large_del_query_exactly_in_deletion(self):
        """Query entirely within the 5 kb deletion: zero-width asm range.

        ref [101000, 106000) → offsets [1000, 6000) = exactly the 5000D
        asm_q_start = 1000, asm_q_end = 1000 (D doesn't advance asm)
        Result: asm [201000, 201000)
        """
        results = coordinate_mapper.query(self.index, "chr7", 101000, 106000)
        self.assertEqual(len(results), 1)
        r = results[0]
        self.assertEqual(r["asm_start"], 201000)
        self.assertEqual(r["asm_end"], 201000)

    def test_large_del_query_after(self):
        """Query after the 5 kb deletion: offset accounts for the D."""
        # ref [106500, 107000) → offsets [6500, 7000) in 1000= after 5000D
        # asm offset: 1000 + 500 = 1500 to 1000 + 1000 = 2000
        # Result: asm [201500, 202000)
        results = coordinate_mapper.query(self.index, "chr7", 106500, 107000)
        self.assertEqual(len(results), 1)
        r = results[0]
        self.assertEqual(r["asm_start"], 201500)
        self.assertEqual(r["asm_end"], 202000)

    # -- Block 2: Multiple small indels (microsatellite-like) --

    def test_multi_indel_full_block(self):
        """Full block query with many small indels: correct total asm range."""
        results = coordinate_mapper.query(self.index, "chr7", 200000, 200402)
        self.assertEqual(len(results), 1)
        r = results[0]
        self.assertEqual(r["asm_start"], 300000)
        self.assertEqual(r["asm_end"], 300408)

    def test_multi_indel_query_before_first_insertion(self):
        """Query before first I: 1:1 mapping in initial = region."""
        # ref [200000, 200050) → offset [0, 50) in 100=
        results = coordinate_mapper.query(self.index, "chr7", 200000, 200050)
        self.assertEqual(len(results), 1)
        r = results[0]
        self.assertEqual(r["asm_start"], 300000)
        self.assertEqual(r["asm_end"], 300050)

    def test_multi_indel_query_spanning_first_insertion(self):
        """Query straddling the first 3I: asm expands by 3."""
        # ref [200090, 200110) → offsets [90, 110)
        # 100= consumed ref 0-100, asm 0-100; 3I: asm 100-103
        # 50= starts at ref 100, asm 103
        # offset 90 in 100= → asm 90
        # offset 110 = 100 (past 100=) + 10 into 50= → asm 103+10 = 113
        results = coordinate_mapper.query(self.index, "chr7", 200090, 200110)
        self.assertEqual(len(results), 1)
        r = results[0]
        self.assertEqual(r["asm_start"], 300090)
        self.assertEqual(r["asm_end"], 300113)

    def test_multi_indel_query_spanning_deletion(self):
        """Query spanning the 2D: asm contracts by 2.

        CIGAR walk to 2D: 100= 3I 50= → ref_cur=150, asm_cur=153
        2D: ref [150, 152), no asm advance
        Next 100=: ref [152, 252), asm [153, 253)

        Query ref [200145, 200160) → offset [145, 160)
        offset 145 in 50= (starts at ref 100, asm 103): 45 into 50= → asm 103+45=148
        offset 152 starts 100= at asm 153.
        offset 160: 8 into 100= after 2D → asm 153+8=161

        Result: asm [300148, 300161)  — 15 ref bases → 13 asm bases (2D removed)
        """
        results = coordinate_mapper.query(self.index, "chr7", 200145, 200160)
        self.assertEqual(len(results), 1)
        r = results[0]
        self.assertEqual(r["asm_start"], 300148)
        self.assertEqual(r["asm_end"], 300161)

    # -- Block 3: Extended CIGAR with =/X --

    def test_eqx_cigar_mismatch_region(self):
        """Query spanning X (mismatch) region: maps 1:1 like M.

        500= 10X 490= ...
        ref [300495, 300515) → offset [495, 515)
        offset 495 in 500= → asm 495
        offset 500 enters 10X → asm 500
        offset 510 enters 490= → asm 510
        offset 515 → asm 515
        Result: asm [400495, 400515)
        """
        results = coordinate_mapper.query(self.index, "chr7", 300495, 300515)
        self.assertEqual(len(results), 1)
        r = results[0]
        self.assertEqual(r["asm_start"], 400495)
        self.assertEqual(r["asm_end"], 400515)

    def test_eqx_cigar_spanning_del_and_ins(self):
        """Query spanning both the 5D and 10I in extended CIGAR.

        500= 10X 490= 5D 500= 10I 500=
        ref_cur / asm_cur at each boundary:
          500=:  ref 500,   asm 500
          10X:   ref 510,   asm 510
          490=:  ref 1000,  asm 1000
          5D:    ref 1005,  asm 1000
          500=:  ref 1505,  asm 1500
          10I:   ref 1505,  asm 1510
          500=:  ref 2005,  asm 2010

        Query ref [300990, 301510) → offset [990, 1510)
        offset 990 in 490= (starts ref 510, asm 510): 480 in → asm 990
        offset 1000 → end of 490=: asm 1000
        5D: ref [1000,1005), asm stays at 1000
        500=: ref [1005,1505), asm [1000,1500)
        offset 1505: end of 500=, asm 1500; then 10I → asm 1510
        next 500=: ref 1505, asm 1510
        offset 1510: 5 into last 500= → asm 1510+5 = 1515

        Result: asm [400990, 401515) — 520 ref → 525 asm (gained 10I, lost 5D)
        """
        results = coordinate_mapper.query(self.index, "chr7", 300990, 301510)
        self.assertEqual(len(results), 1)
        r = results[0]
        self.assertEqual(r["asm_start"], 400990)
        self.assertEqual(r["asm_end"], 401515)

    # -- Block 4: Minus strand with insertion --

    def test_minus_strand_insertion_in_cigar_full(self):
        """Full block query on − strand with 500I: correct total range.

        Block 4: ref [400000, 405000) → asm [500000, 505500) −
        CIGAR 2000=500I3000=: ref consumed = 5000, asm consumed = 5500
        Full query → asm [500000, 505500).
        """
        results = coordinate_mapper.query(self.index, "chr8", 400000, 405000)
        self.assertEqual(len(results), 1)
        r = results[0]
        self.assertEqual(r["asm_start"], 500000)
        self.assertEqual(r["asm_end"], 505500)

    def test_minus_strand_insertion_in_cigar_partial_before(self):
        """Minus strand: query ref [400000, 401000) before the I.

        CIGAR 2000= 500I 3000=; ref offset [0, 1000)
        asm_q_start=0, asm_q_end=1000
        On − strand: [505500-1000, 505500-0] = [504500, 505500)
        """
        results = coordinate_mapper.query(self.index, "chr8", 400000, 401000)
        self.assertEqual(len(results), 1)
        r = results[0]
        self.assertEqual(r["asm_start"], 504500)
        self.assertEqual(r["asm_end"], 505500)

    def test_minus_strand_insertion_in_cigar_straddling(self):
        """Minus strand: query straddling the 500I.

        ref offset [1500, 2500):
        - in 2000= offset 1500 → asm 1500
        - 2000= ends at ref_cur=2000, asm_cur=2000; then 500I → asm_cur=2500
        - in 3000= offset 2000 → asm 2500; offset 2500 → asm 2500+(2500-2000)=3000
        asm_q_start=1500, asm_q_end=3000
        On − strand: [505500-3000, 505500-1500] = [502500, 504000)
        """
        results = coordinate_mapper.query(self.index, "chr8", 401500, 402500)
        self.assertEqual(len(results), 1)
        r = results[0]
        self.assertEqual(r["asm_start"], 502500)
        self.assertEqual(r["asm_end"], 504000)

    def test_minus_strand_insertion_in_cigar_after(self):
        """Minus strand: query after the I.

        ref offset [3000, 4000):
        After 2000= (ref=2000, asm=2000), 500I (asm=2500), 3000= starts.
        offset 3000 = 1000 into 3000= → asm 2500+1000 = 3500
        offset 4000 = 2000 into 3000= → asm 2500+2000 = 4500
        On − strand: [505500-4500, 505500-3500] = [501000, 502000)
        """
        results = coordinate_mapper.query(self.index, "chr8", 403000, 404000)
        self.assertEqual(len(results), 1)
        r = results[0]
        self.assertEqual(r["asm_start"], 501000)
        self.assertEqual(r["asm_end"], 502000)


# =========================================================================
# Query edge cases
# =========================================================================

class TestQueryEdgeCases(unittest.TestCase):
    """Edge cases for query(): boundaries, single-base, zero-width, etc."""

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        paf_path = os.path.join(self.tmpdir, "edge.paf")
        with open(paf_path, "w") as fh:
            fh.write(BIOLOGICAL_SV_PAF)

        blocks = coordinate_mapper.parse_paf(paf_path)
        json_path = os.path.join(self.tmpdir, "edge.mapping.json.gz")
        coordinate_mapper.build_json_index(blocks, json_path)
        self.index = coordinate_mapper.load_index(json_path)

    def tearDown(self):
        shutil.rmtree(self.tmpdir)

    def test_single_base_query(self):
        """Single-base query [pos, pos+1) returns correct 1-base mapping."""
        results = coordinate_mapper.query(
            self.index, "chr1", 1002000, 1002001
        )
        self.assertEqual(len(results), 1)
        r = results[0]
        self.assertEqual(r["ref_start"], 1002000)
        self.assertEqual(r["ref_end"], 1002001)
        self.assertEqual(r["asm_start"], 502000)
        self.assertEqual(r["asm_end"], 502001)

    def test_zero_width_query(self):
        """Zero-width query [pos, pos) returns nothing."""
        results = coordinate_mapper.query(
            self.index, "chr1", 1002000, 1002000
        )
        self.assertEqual(len(results), 0)

    def test_query_at_block_start(self):
        """Query starting exactly at block start."""
        results = coordinate_mapper.query(
            self.index, "chr1", 1000000, 1000100
        )
        self.assertEqual(len(results), 1)
        r = results[0]
        self.assertEqual(r["ref_start"], 1000000)
        self.assertEqual(r["asm_start"], 500000)

    def test_query_at_block_end(self):
        """Query ending exactly at block end boundary.

        Block ends at ref 1005000 (exclusive).  Query [1004900, 1005000)
        should produce a result with ref_end = 1005000.
        """
        results = coordinate_mapper.query(
            self.index, "chr1", 1004900, 1005000
        )
        self.assertEqual(len(results), 1)
        r = results[0]
        self.assertEqual(r["ref_end"], 1005000)
        self.assertEqual(r["asm_end"], 505000)

    def test_query_just_past_block_end(self):
        """Query starting at block end (exclusive boundary) returns no hit."""
        results = coordinate_mapper.query(
            self.index, "chr1", 1005000, 1005001
        )
        # This ref position falls in the deletion gap (ref 1005000-1010000),
        # but only one flanking block (before the gap) is available.
        # Without a second block adjacent, no gap event is generated.
        alns = [r for r in results if r["event_type"] == "alignment"]
        self.assertEqual(len(alns), 0)

    def test_query_spanning_unmapped_gap(self):
        """Query in region with no alignment blocks returns nothing."""
        # Region well away from any aligned blocks on chr1
        results = coordinate_mapper.query(
            self.index, "chr1", 50000000, 50001000
        )
        self.assertEqual(len(results), 0)

    def test_query_exactly_matching_block(self):
        """Query coordinates exactly matching block boundaries."""
        results = coordinate_mapper.query(
            self.index, "chr1", 1000000, 1005000
        )
        self.assertEqual(len(results), 1)
        r = results[0]
        self.assertEqual(r["ref_start"], 1000000)
        self.assertEqual(r["ref_end"], 1005000)
        self.assertEqual(r["asm_start"], 500000)
        self.assertEqual(r["asm_end"], 505000)

    def test_wide_query_spans_multiple_chroms_blocks(self):
        """Wide query on chr1 spans deletion and insertion scenarios."""
        results = coordinate_mapper.query(
            self.index, "chr1", 999000, 2011000
        )
        # Should find: 2 blocks + 1 deletion gap (scenario 1)
        # + 2 blocks + 1 insertion gap (scenario 2)
        # Total: 4 alignments + 2 gap events = 6
        types = [r["event_type"] for r in results]
        self.assertEqual(types.count("alignment"), 4)
        self.assertIn("deletion", types)
        self.assertIn("insertion", types)


# =========================================================================
# Mapping quality filtering
# =========================================================================

# PAF with varying mapq values for testing min_mapq filtering:
#   Block 1: mapq 60 (high confidence primary)
#   Block 2: mapq 0  (supplementary / multi-mapped)
#   Block 3: mapq 20 (moderate confidence)
#   Block 4: mapq 60 (high confidence primary)

MAPQ_PAF = """\
ctg1\t5000000\t10000\t15000\t+\tchr1\t248956422\t100000\t105000\t5000\t5000\t60
ctg1\t5000000\t10000\t15000\t+\tchr1\t248956422\t100000\t105000\t5000\t5000\t0
ctg1\t5000000\t20000\t25000\t+\tchr1\t248956422\t110000\t115000\t5000\t5000\t20
ctg1\t5000000\t30000\t35000\t+\tchr1\t248956422\t120000\t125000\t5000\t5000\t60
"""


class TestMapqFiltering(unittest.TestCase):
    """Tests for min_mapq quality-based filtering."""

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        paf_path = os.path.join(self.tmpdir, "mapq.paf")
        with open(paf_path, "w") as fh:
            fh.write(MAPQ_PAF)

        blocks = coordinate_mapper.parse_paf(paf_path)
        json_path = os.path.join(self.tmpdir, "mapq.mapping.json.gz")
        coordinate_mapper.build_json_index(blocks, json_path)
        self.index = coordinate_mapper.load_index(json_path)

    def tearDown(self):
        shutil.rmtree(self.tmpdir)

    def test_default_no_filtering(self):
        """Default min_mapq=0 returns all blocks."""
        results = coordinate_mapper.query(
            self.index, "chr1", 99000, 126000
        )
        alns = [r for r in results if r["event_type"] == "alignment"]
        self.assertEqual(len(alns), 4)

    def test_min_mapq_1_excludes_zero(self):
        """min_mapq=1 excludes the mapq=0 supplementary block."""
        results = coordinate_mapper.query(
            self.index, "chr1", 99000, 126000, min_mapq=1
        )
        alns = [r for r in results if r["event_type"] == "alignment"]
        self.assertEqual(len(alns), 3)
        mapqs = [r["mapq"] for r in alns]
        self.assertNotIn(0, mapqs)

    def test_min_mapq_30_excludes_low_quality(self):
        """min_mapq=30 excludes mapq=0 and mapq=20 blocks."""
        results = coordinate_mapper.query(
            self.index, "chr1", 99000, 126000, min_mapq=30
        )
        alns = [r for r in results if r["event_type"] == "alignment"]
        self.assertEqual(len(alns), 2)
        for r in alns:
            self.assertGreaterEqual(r["mapq"], 30)

    def test_min_mapq_60_only_primary(self):
        """min_mapq=60 returns only high-confidence primary blocks."""
        results = coordinate_mapper.query(
            self.index, "chr1", 99000, 126000, min_mapq=60
        )
        alns = [r for r in results if r["event_type"] == "alignment"]
        self.assertEqual(len(alns), 2)
        for r in alns:
            self.assertEqual(r["mapq"], 60)

    def test_min_mapq_too_high_returns_empty(self):
        """min_mapq higher than any block returns nothing."""
        results = coordinate_mapper.query(
            self.index, "chr1", 99000, 126000, min_mapq=61
        )
        self.assertEqual(len(results), 0)

    def test_mapq_filtering_with_gaps(self):
        """Gap events between filtered blocks reflect the filtered set.

        With min_mapq=30, blocks are: mapq60 at [100000-105000),
        mapq20 at [110000-115000) excluded, mapq60 at [120000-125000).
        The gap is between the two mapq60 blocks: ref [105000-120000).
        """
        results = coordinate_mapper.query(
            self.index, "chr1", 99000, 126000, min_mapq=30
        )
        types = [r["event_type"] for r in results]
        # Two alignment blocks + one gap between them
        gap_events = [r for r in results if r["event_type"] != "alignment"]
        self.assertEqual(len(gap_events), 1)
        gap = gap_events[0]
        self.assertEqual(gap["ref_start"], 105000)
        self.assertEqual(gap["ref_end"], 120000)


# =========================================================================
# Overlapping / supplementary alignment blocks
# =========================================================================

# PAF modelling overlapping blocks from supplementary alignments:
#   Block 1: ref chr9 [1000, 3000), asm ctg9 [10000, 12000), + mapq=60
#   Block 2: ref chr9 [2500, 5000), asm ctg9 [12000, 14500), + mapq=5
#     Blocks overlap in ref space [2500, 3000)
#
# PAF modelling segdup: same ref region maps to two different contigs
#   Block 3: ref chr10 [100000, 102000), asm ctg10a [200000, 202000), +
#   Block 4: ref chr10 [100000, 102000), asm ctg10b [300000, 302000), +

OVERLAP_PAF = """\
ctg9\t5000000\t10000\t12000\t+\tchr9\t138394717\t1000\t3000\t2000\t2000\t60
ctg9\t5000000\t12000\t14500\t+\tchr9\t138394717\t2500\t5000\t2500\t2500\t5
ctg10a\t5000000\t200000\t202000\t+\tchr10\t133797422\t100000\t102000\t2000\t2000\t60
ctg10b\t5000000\t300000\t302000\t+\tchr10\t133797422\t100000\t102000\t2000\t2000\t60
"""


class TestOverlappingBlocks(unittest.TestCase):
    """Supplementary alignments and segdups producing overlapping blocks."""

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        paf_path = os.path.join(self.tmpdir, "overlap.paf")
        with open(paf_path, "w") as fh:
            fh.write(OVERLAP_PAF)

        blocks = coordinate_mapper.parse_paf(paf_path)
        json_path = os.path.join(self.tmpdir, "overlap.mapping.json.gz")
        coordinate_mapper.build_json_index(blocks, json_path)
        self.index = coordinate_mapper.load_index(json_path)

    def tearDown(self):
        shutil.rmtree(self.tmpdir)

    def test_overlapping_blocks_both_returned(self):
        """Overlapping blocks in ref space both appear in results."""
        results = coordinate_mapper.query(self.index, "chr9", 1000, 5000)
        alns = [r for r in results if r["event_type"] == "alignment"]
        # Both blocks overlap the query
        self.assertEqual(len(alns), 2)

    def test_overlapping_blocks_ref_overlap_region(self):
        """Query in the overlap zone returns both blocks."""
        results = coordinate_mapper.query(self.index, "chr9", 2600, 2900)
        alns = [r for r in results if r["event_type"] == "alignment"]
        self.assertEqual(len(alns), 2)
        # Both should report coordinates in [2600, 2900)
        for r in alns:
            self.assertEqual(r["ref_start"], 2600)
            self.assertEqual(r["ref_end"], 2900)

    def test_overlapping_blocks_different_asm_coords(self):
        """Overlapping blocks map to different assembly positions."""
        results = coordinate_mapper.query(self.index, "chr9", 2600, 2900)
        alns = [r for r in results if r["event_type"] == "alignment"]
        asm_starts = [r["asm_start"] for r in alns]
        # Different assembly positions
        self.assertEqual(len(set(asm_starts)), 2)

    def test_overlapping_filtered_by_mapq(self):
        """min_mapq can remove the low-quality supplementary alignment."""
        results = coordinate_mapper.query(
            self.index, "chr9", 2600, 2900, min_mapq=10
        )
        alns = [r for r in results if r["event_type"] == "alignment"]
        self.assertEqual(len(alns), 1)
        self.assertEqual(alns[0]["mapq"], 60)

    def test_segdup_two_contigs_same_ref(self):
        """Segdup: same ref region maps to two different contigs."""
        results = coordinate_mapper.query(
            self.index, "chr10", 100000, 102000
        )
        alns = [r for r in results if r["event_type"] == "alignment"]
        self.assertEqual(len(alns), 2)
        contigs = sorted(r["asm_chrom"] for r in alns)
        self.assertEqual(contigs, ["ctg10a", "ctg10b"])

    def test_segdup_different_asm_ranges(self):
        """Segdup: the two hits map to different assembly loci."""
        results = coordinate_mapper.query(
            self.index, "chr10", 100500, 101500
        )
        alns = [r for r in results if r["event_type"] == "alignment"]
        self.assertEqual(len(alns), 2)
        asm_ranges = sorted((r["asm_start"], r["asm_end"]) for r in alns)
        self.assertEqual(asm_ranges, [(200500, 201500), (300500, 301500)])


# =========================================================================
# PAF parsing edge cases
# =========================================================================


class TestParsePafWithTags(unittest.TestCase):
    """Test PAF parsing extracts optional tags correctly."""

    def test_cigar_tag_extracted(self):
        """cg:Z: tag is extracted from PAF."""
        paf = "q\t1000\t0\t100\t+\tr\t2000\t0\t100\t100\t100\t60\tcg:Z:100M\n"
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".paf", delete=False
        ) as fh:
            fh.write(paf)
            path = fh.name
        try:
            blocks = coordinate_mapper.parse_paf(path)
            self.assertEqual(blocks[0]["cigar"], "100M")
        finally:
            os.unlink(path)

    def test_tp_tag_extracted(self):
        """tp:A: tag is extracted from PAF."""
        paf = "q\t1000\t0\t100\t+\tr\t2000\t0\t100\t100\t100\t60\ttp:A:P\n"
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".paf", delete=False
        ) as fh:
            fh.write(paf)
            path = fh.name
        try:
            blocks = coordinate_mapper.parse_paf(path)
            self.assertEqual(blocks[0]["tp"], "P")
        finally:
            os.unlink(path)

    def test_multiple_tags(self):
        """Multiple optional tags are extracted."""
        paf = (
            "q\t1000\t0\t100\t+\tr\t2000\t0\t100\t100\t100\t60\t"
            "NM:i:5\tcg:Z:95M5I\ttp:A:S\tms:i:190\n"
        )
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".paf", delete=False
        ) as fh:
            fh.write(paf)
            path = fh.name
        try:
            blocks = coordinate_mapper.parse_paf(path)
            self.assertEqual(blocks[0]["cigar"], "95M5I")
            self.assertEqual(blocks[0]["tp"], "S")
        finally:
            os.unlink(path)

    def test_no_optional_tags(self):
        """Blocks without optional tags get empty strings."""
        paf = "q\t1000\t0\t100\t+\tr\t2000\t0\t100\t100\t100\t60\n"
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".paf", delete=False
        ) as fh:
            fh.write(paf)
            path = fh.name
        try:
            blocks = coordinate_mapper.parse_paf(path)
            self.assertEqual(blocks[0]["cigar"], "")
            self.assertEqual(blocks[0]["tp"], "")
        finally:
            os.unlink(path)


# =========================================================================
# project_cigar unit-level edge cases
# =========================================================================


class TestProjectCigarUnitEdges(unittest.TestCase):
    """Unit tests for project_cigar edge cases at the function level."""

    def test_query_exactly_at_cigar_end(self):
        """Query end coincides exactly with CIGAR ref length."""
        ops = coordinate_mapper.parse_cigar("100M50D100M")
        # Total ref = 250.  Query [0, 250).
        s, e = coordinate_mapper.project_cigar(ops, 0, 250, "+", 0, 200)
        self.assertEqual((s, e), (0, 200))

    def test_query_beyond_cigar_end(self):
        """Query extends past CIGAR: result clamped to CIGAR end."""
        ops = coordinate_mapper.parse_cigar("100M")
        s, e = coordinate_mapper.project_cigar(ops, 0, 200, "+", 0, 100)
        self.assertEqual((s, e), (0, 100))

    def test_single_base_at_deletion_start(self):
        """Single-base query at the exact start of a deletion.

        50M10D50M: ref_offset [50, 51) = first D base.
        D doesn't advance asm → asm_q_start = asm_q_end = 50.
        """
        ops = coordinate_mapper.parse_cigar("50M10D50M")
        s, e = coordinate_mapper.project_cigar(ops, 50, 51, "+", 0, 100)
        self.assertEqual(s, 50)
        self.assertEqual(e, 50)

    def test_single_base_at_deletion_end(self):
        """Single-base query at the last D base.

        50M10D50M: ref_offset [59, 60) = last D base.
        Still in D → asm stays at 50.
        """
        ops = coordinate_mapper.parse_cigar("50M10D50M")
        s, e = coordinate_mapper.project_cigar(ops, 59, 60, "+", 0, 100)
        self.assertEqual(s, 50)
        self.assertEqual(e, 50)

    def test_single_base_right_after_deletion(self):
        """Single-base query at first base after deletion.

        50M10D50M: ref_offset [60, 61) = first base of trailing 50M.
        asm offset = 50 + (60-60) = 50 → asm [50, 51).
        """
        ops = coordinate_mapper.parse_cigar("50M10D50M")
        s, e = coordinate_mapper.project_cigar(ops, 60, 61, "+", 0, 100)
        self.assertEqual(s, 50)
        self.assertEqual(e, 51)

    def test_insertion_at_cigar_start(self):
        """CIGAR starting with an insertion.

        10I90M: the I doesn't consume ref, so ref starts at 0.
        Query [0, 10): in the 90M → asm [10, 20).
        """
        ops = coordinate_mapper.parse_cigar("10I90M")
        s, e = coordinate_mapper.project_cigar(ops, 0, 10, "+", 0, 100)
        self.assertEqual(s, 10)
        self.assertEqual(e, 20)

    def test_consecutive_deletions(self):
        """Consecutive D ops: 50M5D10D50M.

        ref consumed: 50 + 5 + 10 + 50 = 115
        asm consumed: 50 + 50 = 100
        Query [45, 70): spans both D regions.
        offset 45 in 50M → asm 45
        offsets 50-65 in Ds → asm stays at 50
        offset 65 = 0 into last 50M → asm 50
        offset 70 = 5 into last 50M → asm 55
        """
        ops = coordinate_mapper.parse_cigar("50M5D10D50M")
        s, e = coordinate_mapper.project_cigar(ops, 45, 70, "+", 0, 100)
        self.assertEqual(s, 45)
        self.assertEqual(e, 55)

    def test_consecutive_insertions(self):
        """Consecutive I ops: 50M10I5I50M.

        ref consumed: 50 + 50 = 100
        asm consumed: 50 + 10 + 5 + 50 = 115
        Query [45, 55): straddles the I region.
        offset 45 in 50M → asm 45
        50M ends at ref_cur=50, asm_cur=50; 10I → asm=60; 5I → asm=65
        55 = 5 into last 50M → asm 65+5=70
        """
        ops = coordinate_mapper.parse_cigar("50M10I5I50M")
        s, e = coordinate_mapper.project_cigar(ops, 45, 55, "+", 0, 115)
        self.assertEqual(s, 45)
        self.assertEqual(e, 70)

    def test_minus_strand_with_insertion_at_cigar_start(self):
        """Minus strand: I at start of CIGAR.

        10I100M: ref consumed = 100, asm consumed = 110.
        Block: asm [200, 310), minus.
        Query full [0, 100): asm_q=[10, 110).
        On − strand: [310-110, 310-10] = [200, 300).
        """
        ops = coordinate_mapper.parse_cigar("10I100M")
        s, e = coordinate_mapper.project_cigar(ops, 0, 100, "-", 200, 310)
        self.assertEqual((s, e), (200, 300))


# =========================================================================
# Backward compatibility: no-CIGAR fallback
# =========================================================================


class TestNoCigarFallback(unittest.TestCase):
    """Verify that blocks without CIGAR fall back to linear interpolation."""

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        paf_path = os.path.join(self.tmpdir, "no_cigar.paf")
        # PAF without cg:Z: tag
        with open(paf_path, "w") as fh:
            fh.write(
                "ctg1\t5000000\t1000\t4000\t+\tchr1\t248956422\t10000\t13000"
                "\t2900\t3000\t60\n"
                "ctg1\t5000000\t5000\t8000\t-\tchr1\t248956422\t20000\t23000"
                "\t2800\t3000\t50\n"
            )

        blocks = coordinate_mapper.parse_paf(paf_path)
        json_path = os.path.join(self.tmpdir, "no_cigar.mapping.json.gz")
        coordinate_mapper.build_json_index(blocks, json_path)
        self.index = coordinate_mapper.load_index(json_path)

    def tearDown(self):
        shutil.rmtree(self.tmpdir)

    def test_no_cigar_in_index(self):
        """Blocks without CIGAR should not have cg field in index."""
        block = self.index["chr1"]["blocks"][0]
        self.assertNotIn("cg", block)

    def test_plus_strand_linear(self):
        """Without CIGAR: plus-strand uses linear interpolation."""
        results = coordinate_mapper.query(self.index, "chr1", 11000, 12000)
        self.assertEqual(len(results), 1)
        r = results[0]
        # Linear: asm_start + offset = 1000 + 1000 = 2000
        self.assertEqual(r["asm_start"], 2000)
        self.assertEqual(r["asm_end"], 3000)

    def test_minus_strand_linear(self):
        """Without CIGAR: minus-strand uses linear interpolation."""
        results = coordinate_mapper.query(self.index, "chr1", 21000, 22000)
        self.assertEqual(len(results), 1)
        r = results[0]
        # Minus: ae - offset_end, ae - offset_start = 8000-2000, 8000-1000
        self.assertEqual(r["asm_start"], 6000)
        self.assertEqual(r["asm_end"], 7000)

    def test_event_type_present(self):
        """Even without CIGAR, results carry event_type='alignment'."""
        results = coordinate_mapper.query(self.index, "chr1", 11000, 12000)
        self.assertEqual(results[0]["event_type"], "alignment")


# =========================================================================
# classify_sv_gap unit tests
# =========================================================================


class TestClassifySvGapUnit(unittest.TestCase):
    """Direct unit tests for classify_sv_gap()."""

    def _block(self, ac="ctg1", asm_start=0, asm_end=1000,
               strand="+", mapq=60, ref_start=0, ref_end=1000):
        return {
            "ac": ac, "as": asm_start, "ae": asm_end,
            "st": strand, "mq": mapq, "rs": ref_start, "re": ref_end,
        }

    def test_deletion_plus_strand(self):
        prev = self._block(asm_end=1000)
        nxt = self._block(asm_start=1000, asm_end=2000)
        result = coordinate_mapper.classify_sv_gap(
            prev, nxt, 1000, 6000, "chr1"
        )
        self.assertEqual(result["event_type"], "deletion")
        self.assertEqual(result["ref_gap_size"], 5000)
        self.assertEqual(result["asm_gap_size"], 0)

    def test_insertion_plus_strand(self):
        prev = self._block(asm_end=1000)
        nxt = self._block(asm_start=6000, asm_end=7000)
        result = coordinate_mapper.classify_sv_gap(
            prev, nxt, 1000, 1000, "chr1"
        )
        self.assertEqual(result["event_type"], "insertion")
        self.assertEqual(result["ref_gap_size"], 0)
        self.assertEqual(result["asm_gap_size"], 5000)

    def test_complex_plus_strand(self):
        """Overlapping assembly ranges (negative gap) → complex."""
        prev = self._block(asm_end=5000)
        nxt = self._block(asm_start=3000, asm_end=6000)
        result = coordinate_mapper.classify_sv_gap(
            prev, nxt, 1000, 2000, "chr1"
        )
        self.assertEqual(result["event_type"], "complex")
        self.assertEqual(result["asm_gap_size"], -2000)

    def test_inversion(self):
        prev = self._block(strand="+")
        nxt = self._block(strand="-")
        result = coordinate_mapper.classify_sv_gap(
            prev, nxt, 1000, 2000, "chr1"
        )
        self.assertEqual(result["event_type"], "inversion")
        self.assertIsNone(result["asm_gap_size"])

    def test_translocation(self):
        prev = self._block(ac="ctgA")
        nxt = self._block(ac="ctgB")
        result = coordinate_mapper.classify_sv_gap(
            prev, nxt, 1000, 2000, "chr1"
        )
        self.assertEqual(result["event_type"], "translocation")
        self.assertIsNone(result["asm_gap_size"])

    def test_deletion_minus_strand(self):
        """Deletion on minus strand: asm runs in reverse."""
        prev = self._block(asm_start=5000, asm_end=8000, strand="-")
        nxt = self._block(asm_start=3000, asm_end=4000, strand="-")
        # asm_gap = prev["as"] - nxt["ae"] = 5000 - 4000 = 1000
        result = coordinate_mapper.classify_sv_gap(
            prev, nxt, 8000, 13000, "chr1"
        )
        self.assertEqual(result["event_type"], "deletion")
        self.assertEqual(result["ref_gap_size"], 5000)
        self.assertEqual(result["asm_gap_size"], 1000)

    def test_insertion_minus_strand(self):
        """Insertion on minus strand."""
        prev = self._block(asm_start=10000, asm_end=13000, strand="-")
        nxt = self._block(asm_start=3000, asm_end=6000, strand="-")
        # asm_gap = prev["as"] - nxt["ae"] = 10000 - 6000 = 4000
        # ref_gap = 0
        result = coordinate_mapper.classify_sv_gap(
            prev, nxt, 13000, 13000, "chr1"
        )
        self.assertEqual(result["event_type"], "insertion")
        self.assertEqual(result["asm_gap_size"], 4000)

    def test_mapq_is_minimum(self):
        prev = self._block(mapq=60)
        nxt = self._block(mapq=20)
        result = coordinate_mapper.classify_sv_gap(
            prev, nxt, 1000, 2000, "chr1"
        )
        self.assertEqual(result["mapq"], 20)


if __name__ == '__main__':
    unittest.main()
