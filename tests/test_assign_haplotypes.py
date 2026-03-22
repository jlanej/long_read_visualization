#!/usr/bin/env python3
"""Tests for assign_haplotypes.py"""

import os
import shutil
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

try:
    import pysam
    HAS_PYSAM = True
except ImportError:
    HAS_PYSAM = False

if HAS_PYSAM:
    import assign_haplotypes  # noqa: E402


def _make_bam(path, reads, contigs=None):
    """Create a minimal BAM with the given reads.

    *reads* is a list of dicts with keys: qname, contig, pos, seq, qual, tags.
    *contigs* defaults to [("chr1", 100000)].
    """
    if contigs is None:
        contigs = [("chr1", 100000)]
    header = pysam.AlignmentHeader.from_dict({
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": c, "LN": l} for c, l in contigs],
    })
    with pysam.AlignmentFile(path, "wb", header=header) as outf:
        for r in reads:
            a = pysam.AlignedSegment(header)
            a.query_name = r["qname"]
            a.reference_id = header.get_tid(r.get("contig", "chr1"))
            a.reference_start = r.get("pos", 0)
            seq = r.get("seq", "ACGT")
            a.query_sequence = seq
            a.query_qualities = pysam.qualitystring_to_array(
                r.get("qual", "I" * len(seq))
            )
            a.cigar = [(0, len(seq))]  # M
            a.mapping_quality = 60
            for tag, val in r.get("tags", []):
                a.set_tag(tag, val)
            outf.write(a)
    pysam.sort("-o", path, path)
    pysam.index(path)


def _read_hp_tags(bam_path):
    """Return {qname: HP_value} for all reads carrying an HP tag."""
    result = {}
    with pysam.AlignmentFile(bam_path, "rb") as f:
        for read in f.fetch(until_eof=True):
            try:
                result[read.query_name] = read.get_tag("HP")
            except KeyError:
                pass
    return result


@unittest.skipUnless(HAS_PYSAM, "pysam not installed")
class TestGetBestScores(unittest.TestCase):
    """Tests for get_best_scores()."""

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmpdir)

    def test_single_alignment_per_read(self):
        """Each read appears once with a known AS tag."""
        bam = os.path.join(self.tmpdir, "test.bam")
        _make_bam(bam, [
            {"qname": "r1", "tags": [("AS", 100)]},
            {"qname": "r2", "tags": [("AS", 200)]},
        ])
        scores = assign_haplotypes.get_best_scores(bam)
        self.assertEqual(scores, {"r1": 100, "r2": 200})

    def test_best_of_multiple_alignments(self):
        """When a read has multiple alignments, keep the highest AS."""
        bam = os.path.join(self.tmpdir, "test.bam")
        _make_bam(bam, [
            {"qname": "r1", "pos": 0, "tags": [("AS", 50)]},
            {"qname": "r1", "pos": 10, "tags": [("AS", 200)]},
            {"qname": "r1", "pos": 20, "tags": [("AS", 150)]},
        ])
        scores = assign_haplotypes.get_best_scores(bam)
        self.assertEqual(scores["r1"], 200)

    def test_missing_as_tag_skipped(self):
        """Reads without an AS tag are ignored."""
        bam = os.path.join(self.tmpdir, "test.bam")
        _make_bam(bam, [
            {"qname": "r1", "tags": [("AS", 100)]},
            {"qname": "r2", "tags": []},
        ])
        scores = assign_haplotypes.get_best_scores(bam)
        self.assertIn("r1", scores)
        self.assertNotIn("r2", scores)

    def test_unmapped_reads_skipped(self):
        """Unmapped reads are not included in scores."""
        bam = os.path.join(self.tmpdir, "test.bam")
        header = pysam.AlignmentHeader.from_dict({
            "HD": {"VN": "1.6", "SO": "coordinate"},
            "SQ": [{"SN": "chr1", "LN": 100000}],
        })
        with pysam.AlignmentFile(bam, "wb", header=header) as outf:
            # Mapped read
            a = pysam.AlignedSegment(header)
            a.query_name = "r1"
            a.reference_id = 0
            a.reference_start = 0
            a.query_sequence = "ACGT"
            a.query_qualities = pysam.qualitystring_to_array("IIII")
            a.cigar = [(0, 4)]
            a.mapping_quality = 60
            a.set_tag("AS", 100)
            outf.write(a)
            # Unmapped read
            b = pysam.AlignedSegment(header)
            b.query_name = "r2"
            b.flag = 4  # unmapped
            b.query_sequence = "ACGT"
            b.query_qualities = pysam.qualitystring_to_array("IIII")
            b.set_tag("AS", 999)
            outf.write(b)
        pysam.index(bam)
        scores = assign_haplotypes.get_best_scores(bam)
        self.assertIn("r1", scores)
        self.assertNotIn("r2", scores)


@unittest.skipUnless(HAS_PYSAM, "pysam not installed")
class TestAssignHaplotypes(unittest.TestCase):
    """Tests for assign_haplotypes()."""

    def test_clear_hap1_winner(self):
        """Read with higher AS in hap1 → HP=1."""
        result = assign_haplotypes.assign_haplotypes(
            {"r1": 500}, {"r1": 300}
        )
        self.assertEqual(result["r1"], 1)

    def test_clear_hap2_winner(self):
        """Read with higher AS in hap2 → HP=2."""
        result = assign_haplotypes.assign_haplotypes(
            {"r1": 300}, {"r1": 500}
        )
        self.assertEqual(result["r1"], 2)

    def test_equal_scores_ambiguous(self):
        """Equal AS → HP=0."""
        result = assign_haplotypes.assign_haplotypes(
            {"r1": 500}, {"r1": 500}
        )
        self.assertEqual(result["r1"], 0)

    def test_tolerance(self):
        """Scores within tolerance → HP=0."""
        result = assign_haplotypes.assign_haplotypes(
            {"r1": 505}, {"r1": 500}, tolerance=10
        )
        self.assertEqual(result["r1"], 0)

    def test_tolerance_exceeded(self):
        """Score diff exceeding tolerance → assigned to winner."""
        result = assign_haplotypes.assign_haplotypes(
            {"r1": 520}, {"r1": 500}, tolerance=10
        )
        self.assertEqual(result["r1"], 1)

    def test_read_only_in_hap1(self):
        """Read present only in hap1 → HP=1."""
        result = assign_haplotypes.assign_haplotypes(
            {"r1": 500}, {}
        )
        self.assertEqual(result["r1"], 1)

    def test_read_only_in_hap2(self):
        """Read present only in hap2 → HP=2."""
        result = assign_haplotypes.assign_haplotypes(
            {}, {"r1": 500}
        )
        self.assertEqual(result["r1"], 2)

    def test_multiple_reads(self):
        """Multiple reads assigned correctly."""
        result = assign_haplotypes.assign_haplotypes(
            {"r1": 500, "r2": 100, "r3": 300},
            {"r1": 300, "r2": 500, "r3": 300},
        )
        self.assertEqual(result["r1"], 1)
        self.assertEqual(result["r2"], 2)
        self.assertEqual(result["r3"], 0)

    def test_empty_inputs(self):
        """Empty inputs → empty assignments."""
        result = assign_haplotypes.assign_haplotypes({}, {})
        self.assertEqual(result, {})


@unittest.skipUnless(HAS_PYSAM, "pysam not installed")
class TestTagBam(unittest.TestCase):
    """Tests for tag_bam() and tag_bam_in_place()."""

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmpdir)

    def test_tag_bam_writes_hp(self):
        """tag_bam() injects HP tags into reads."""
        src = os.path.join(self.tmpdir, "src.bam")
        dst = os.path.join(self.tmpdir, "dst.bam")
        _make_bam(src, [
            {"qname": "r1", "tags": [("AS", 100)]},
            {"qname": "r2", "tags": [("AS", 200)]},
        ])
        assignments = {"r1": 1, "r2": 2}
        assign_haplotypes.tag_bam(src, dst, assignments)
        hp = _read_hp_tags(dst)
        self.assertEqual(hp, {"r1": 1, "r2": 2})

    def test_tag_bam_preserves_unknown_reads(self):
        """Reads not in assignments pass through unchanged."""
        src = os.path.join(self.tmpdir, "src.bam")
        dst = os.path.join(self.tmpdir, "dst.bam")
        _make_bam(src, [
            {"qname": "r1", "tags": [("AS", 100)]},
            {"qname": "rX", "tags": [("AS", 50)]},
        ])
        assign_haplotypes.tag_bam(src, dst, {"r1": 1})
        hp = _read_hp_tags(dst)
        self.assertEqual(hp, {"r1": 1})
        # rX should still be in the file
        with pysam.AlignmentFile(dst, "rb") as f:
            names = [r.query_name for r in f.fetch(until_eof=True)]
        self.assertIn("rX", names)

    def test_tag_bam_in_place(self):
        """tag_bam_in_place() replaces the BAM and re-indexes."""
        bam = os.path.join(self.tmpdir, "test.bam")
        _make_bam(bam, [
            {"qname": "r1", "tags": [("AS", 100)]},
            {"qname": "r2", "tags": [("AS", 200)]},
        ])
        assign_haplotypes.tag_bam_in_place(bam, {"r1": 2, "r2": 1})
        hp = _read_hp_tags(bam)
        self.assertEqual(hp, {"r1": 2, "r2": 1})
        # Index should exist
        self.assertTrue(os.path.isfile(bam + ".bai"))

    def test_overwrites_existing_hp(self):
        """Re-running on already-tagged BAMs produces the same result."""
        bam = os.path.join(self.tmpdir, "test.bam")
        _make_bam(bam, [
            {"qname": "r1", "tags": [("AS", 100), ("HP", 2)]},
            {"qname": "r2", "tags": [("AS", 200), ("HP", 0)]},
        ])
        assign_haplotypes.tag_bam_in_place(bam, {"r1": 1, "r2": 2})
        hp = _read_hp_tags(bam)
        self.assertEqual(hp, {"r1": 1, "r2": 2})


@unittest.skipUnless(HAS_PYSAM, "pysam not installed")
class TestIdempotency(unittest.TestCase):
    """End-to-end idempotency: running assign twice gives same result."""

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmpdir)

    def test_double_run_same_result(self):
        """Running the full pipeline twice yields identical HP assignments."""
        hap1 = os.path.join(self.tmpdir, "h1.bam")
        hap2 = os.path.join(self.tmpdir, "h2.bam")
        reads = [
            {"qname": "r1", "tags": [("AS", 500)]},
            {"qname": "r2", "tags": [("AS", 100)]},
            {"qname": "r3", "tags": [("AS", 300)]},
        ]
        _make_bam(hap1, reads)
        _make_bam(hap2, [
            {"qname": "r1", "tags": [("AS", 300)]},
            {"qname": "r2", "tags": [("AS", 500)]},
            {"qname": "r3", "tags": [("AS", 300)]},
        ])
        # First run
        assign_haplotypes.main([
            "--hap1-bam", hap1, "--hap2-bam", hap2,
        ])
        hp_first = _read_hp_tags(hap1)

        # Second run (idempotency)
        assign_haplotypes.main([
            "--hap1-bam", hap1, "--hap2-bam", hap2,
        ])
        hp_second = _read_hp_tags(hap1)

        self.assertEqual(hp_first, hp_second)


@unittest.skipUnless(HAS_PYSAM, "pysam not installed")
class TestMainCLI(unittest.TestCase):
    """Tests for the main() CLI entry point."""

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmpdir)

    def test_main_tags_both_bams(self):
        """main() tags both hap1 and hap2 BAMs."""
        hap1 = os.path.join(self.tmpdir, "h1.bam")
        hap2 = os.path.join(self.tmpdir, "h2.bam")
        _make_bam(hap1, [
            {"qname": "r1", "tags": [("AS", 500)]},
            {"qname": "r2", "tags": [("AS", 100)]},
        ])
        _make_bam(hap2, [
            {"qname": "r1", "tags": [("AS", 300)]},
            {"qname": "r2", "tags": [("AS", 400)]},
        ])
        assign_haplotypes.main([
            "--hap1-bam", hap1, "--hap2-bam", hap2,
        ])
        hp1 = _read_hp_tags(hap1)
        hp2 = _read_hp_tags(hap2)
        self.assertEqual(hp1, {"r1": 1, "r2": 2})
        self.assertEqual(hp2, {"r1": 1, "r2": 2})

    def test_main_with_tolerance(self):
        """main() respects --tolerance flag."""
        hap1 = os.path.join(self.tmpdir, "h1.bam")
        hap2 = os.path.join(self.tmpdir, "h2.bam")
        _make_bam(hap1, [
            {"qname": "r1", "tags": [("AS", 505)]},
        ])
        _make_bam(hap2, [
            {"qname": "r1", "tags": [("AS", 500)]},
        ])
        assign_haplotypes.main([
            "--hap1-bam", hap1, "--hap2-bam", hap2,
            "--tolerance", "10",
        ])
        hp = _read_hp_tags(hap1)
        self.assertEqual(hp["r1"], 0)  # within tolerance → ambiguous

    def test_main_missing_file(self):
        """main() exits with error for missing BAM."""
        with self.assertRaises(SystemExit):
            assign_haplotypes.main([
                "--hap1-bam", "/nonexistent/h1.bam",
                "--hap2-bam", "/nonexistent/h2.bam",
            ])


if __name__ == "__main__":
    unittest.main()
