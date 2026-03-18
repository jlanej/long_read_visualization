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

if __name__ == '__main__':
    unittest.main()
