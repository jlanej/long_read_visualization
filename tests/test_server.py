"""Tests for the IGV.js visualization server."""

import json
import os
import sys
import tempfile
import unittest

# Allow importing from server/ and src/ directories
_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(_REPO_ROOT, "server"))
sys.path.insert(0, os.path.join(_REPO_ROOT, "src"))

import app as server_app


class TestLoadConfig(unittest.TestCase):
    """Tests for loading the TSV configuration file."""

    def test_basic_config(self):
        """Load a minimal config with one sample."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv",
                                         delete=False) as f:
            f.write("#sample_id\toutput_dir\treference\thap1_assembly\t"
                    "hap2_assembly\treads_bam\tregions\n")
            f.write("S1\t/out\t/ref.fa\t/h1.fa\t/h2.fa\t/reads.bam\t/r.json\n")
            f.flush()
            samples = server_app.load_config(f.name)
        os.unlink(f.name)
        self.assertEqual(len(samples), 1)
        self.assertEqual(samples[0]["sample_id"], "S1")
        self.assertEqual(samples[0]["reference"], "/ref.fa")
        self.assertEqual(samples[0]["hap1_assembly"], "/h1.fa")

    def test_multiple_samples(self):
        """Load a config with two samples."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv",
                                         delete=False) as f:
            f.write("#sample_id\toutput_dir\treference\n")
            f.write("A\t/outA\t/ref.fa\n")
            f.write("B\t/outB\t/ref.fa\n")
            f.flush()
            samples = server_app.load_config(f.name)
        os.unlink(f.name)
        self.assertEqual(len(samples), 2)
        self.assertEqual(samples[0]["sample_id"], "A")
        self.assertEqual(samples[1]["sample_id"], "B")

    def test_skips_comments_and_blank_lines(self):
        """Comments and blank lines are ignored."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv",
                                         delete=False) as f:
            f.write("#sample_id\toutput_dir\n")
            f.write("# This is a comment\n")
            f.write("\n")
            f.write("S1\t/out\n")
            f.flush()
            samples = server_app.load_config(f.name)
        os.unlink(f.name)
        self.assertEqual(len(samples), 1)

    def test_header_without_hash(self):
        """A header line without '#' prefix is also valid."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv",
                                         delete=False) as f:
            f.write("sample_id\toutput_dir\n")
            f.write("S1\t/out\n")
            f.flush()
            samples = server_app.load_config(f.name)
        os.unlink(f.name)
        self.assertEqual(len(samples), 1)
        self.assertEqual(samples[0]["sample_id"], "S1")


class TestDiscoverPipelineFiles(unittest.TestCase):
    """Tests for auto-discovering pipeline output files."""

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        # Create mock pipeline output files
        prefix = "test_sample"
        suffixes = [
            "_hap1_to_ref.bam",
            "_hap1_to_ref.bam.bai",
            "_hap2_to_ref.bam",
            "_hap2_to_ref.bam.bai",
            "_reads_to_hap1.bam",
            "_reads_to_hap1.bam.bai",
            "_reads_to_hap2.bam",
            "_reads_to_hap2.bam.bai",
            "_ref_to_hap1.bam",
            "_ref_to_hap1.bam.bai",
            "_ref_to_hap2.bam",
            "_ref_to_hap2.bam.bai",
            "_hap1_to_ref.mapping.json.gz",
            "_hap2_to_ref.mapping.json.gz",
        ]
        for s in suffixes:
            open(os.path.join(self.tmpdir, prefix + s), "w").close()

    def tearDown(self):
        import shutil
        shutil.rmtree(self.tmpdir)

    def test_discovers_all_files(self):
        """All pipeline files are discovered by convention."""
        sample = {"sample_id": "test", "output_dir": self.tmpdir}
        sample = server_app.discover_pipeline_files(sample)
        self.assertEqual(sample["_prefix"], "test_sample")
        self.assertIn("hap1_to_ref_bam", sample)
        self.assertIn("hap2_to_ref_bam", sample)
        self.assertIn("reads_to_hap1_bam", sample)
        self.assertIn("reads_to_hap2_bam", sample)
        self.assertIn("ref_to_hap1_bam", sample)
        self.assertIn("ref_to_hap2_bam", sample)
        self.assertIn("hap1_mapping_index", sample)
        self.assertIn("hap2_mapping_index", sample)

    def test_missing_output_dir(self):
        """Missing output_dir is handled gracefully."""
        sample = {"sample_id": "test", "output_dir": "/nonexistent"}
        result = server_app.discover_pipeline_files(sample)
        self.assertNotIn("hap1_to_ref_bam", result)


class TestLoadRegions(unittest.TestCase):
    """Tests for loading regions of interest."""

    def test_load_manifest_json(self):
        """Load regions from a manifest JSON file."""
        manifest = {
            "description": "test",
            "variants": [
                {
                    "chrom": "chr1",
                    "pos": 1000,
                    "size": 500,
                    "genotype": "1|0",
                    "ref_region": "chr1:500-2000",
                    "hap1_regions": ["asm1:500-2000"],
                    "hap2_regions": ["asm2:500-2000"],
                }
            ]
        }
        with tempfile.NamedTemporaryFile(mode="w", suffix=".json",
                                         delete=False) as f:
            json.dump(manifest, f)
            f.flush()
            regions = server_app.load_regions(f.name)
        os.unlink(f.name)
        self.assertEqual(len(regions), 1)
        self.assertEqual(regions[0]["chrom"], "chr1")
        self.assertEqual(regions[0]["start"], 500)
        self.assertEqual(regions[0]["end"], 2000)
        self.assertEqual(regions[0]["size"], 500)
        self.assertIn("1|0", regions[0]["label"])

    def test_load_nonexistent_file(self):
        """Non-existent file returns empty list."""
        self.assertEqual(server_app.load_regions("/nonexistent"), [])

    def test_load_empty_path(self):
        """Empty path returns empty list."""
        self.assertEqual(server_app.load_regions(""), [])

    def test_load_toy_manifest(self):
        """Load the actual toy manifest and get 40 regions."""
        manifest_path = os.path.join(
            _REPO_ROOT, "resources", "toy_dataset", "toy_manifest.json")
        if not os.path.isfile(manifest_path):
            self.skipTest("Toy manifest not found")
        regions = server_app.load_regions(manifest_path)
        self.assertEqual(len(regions), 40)
        # Verify regions have expected fields
        for r in regions:
            self.assertIn("chrom", r)
            self.assertIn("start", r)
            self.assertIn("end", r)
            self.assertIn("label", r)
            self.assertIn("ref_region", r)
            self.assertTrue(r["end"] > r["start"])


class TestMergeRegions(unittest.TestCase):
    """Tests for assembly region merging."""

    def test_non_overlapping(self):
        """Non-overlapping regions are kept separate."""
        regions = [
            {"chrom": "chr1", "start": 0, "end": 100, "strand": "+"},
            {"chrom": "chr1", "start": 200, "end": 300, "strand": "+"},
        ]
        merged = server_app._merge_regions(regions)
        self.assertEqual(len(merged), 2)

    def test_overlapping(self):
        """Overlapping regions are merged."""
        regions = [
            {"chrom": "chr1", "start": 0, "end": 150, "strand": "+"},
            {"chrom": "chr1", "start": 100, "end": 300, "strand": "+"},
        ]
        merged = server_app._merge_regions(regions)
        self.assertEqual(len(merged), 1)
        self.assertEqual(merged[0]["start"], 0)
        self.assertEqual(merged[0]["end"], 300)

    def test_different_contigs(self):
        """Regions on different contigs are not merged."""
        regions = [
            {"chrom": "chr1", "start": 0, "end": 100, "strand": "+"},
            {"chrom": "chr2", "start": 0, "end": 100, "strand": "+"},
        ]
        merged = server_app._merge_regions(regions)
        self.assertEqual(len(merged), 2)

    def test_adjacent(self):
        """Adjacent (touching) regions are merged."""
        regions = [
            {"chrom": "chr1", "start": 0, "end": 100, "strand": "+"},
            {"chrom": "chr1", "start": 100, "end": 200, "strand": "+"},
        ]
        merged = server_app._merge_regions(regions)
        self.assertEqual(len(merged), 1)


class TestFormatSize(unittest.TestCase):
    """Tests for human-readable size formatting."""

    def test_bp(self):
        self.assertEqual(server_app._format_size(500), "500 bp")

    def test_kb(self):
        self.assertEqual(server_app._format_size(5000), "5.0 kb")

    def test_mb(self):
        self.assertEqual(server_app._format_size(1500000), "1.5 Mb")


class TestCoordinateTranslator(unittest.TestCase):
    """Tests for the coordinate translator with toy data."""

    @classmethod
    def setUpClass(cls):
        """Load the toy mapping indices if available."""
        toy_output = "/tmp/toy_output"
        cls.indices_available = False

        if not os.path.isdir(toy_output):
            return

        hap1_idx = os.path.join(toy_output,
                                "toy_reads_hap1_to_ref.mapping.json.gz")
        hap2_idx = os.path.join(toy_output,
                                "toy_reads_hap2_to_ref.mapping.json.gz")

        if os.path.isfile(hap1_idx) and os.path.isfile(hap2_idx):
            cls.translator = server_app.CoordinateTranslator()
            cls.translator.load_sample({
                "sample_id": "test",
                "hap1_mapping_index": hap1_idx,
                "hap2_mapping_index": hap2_idx,
            })
            cls.indices_available = True

    def test_translate_known_region(self):
        """Translate the first variant's reference region."""
        if not self.indices_available:
            self.skipTest("Toy output not available")

        # First variant: chr1:9279383-9389426 region
        result = self.translator.translate(
            "test", "chr1:9279383-9389426", 0, 110043)
        self.assertIn("hap1", result)
        self.assertIn("hap2", result)
        self.assertGreater(len(result["hap1"]), 0)
        self.assertGreater(len(result["hap2"]), 0)

    def test_translate_unknown_chrom(self):
        """Unknown chromosome returns empty results."""
        if not self.indices_available:
            self.skipTest("Toy output not available")

        result = self.translator.translate(
            "test", "chrUnknown", 0, 1000)
        self.assertEqual(result["hap1"], [])
        self.assertEqual(result["hap2"], [])

    def test_translate_unknown_sample(self):
        """Unknown sample returns empty results."""
        if not self.indices_available:
            self.skipTest("Toy output not available")

        result = self.translator.translate(
            "nonexistent", "chr1", 0, 1000)
        self.assertEqual(result["hap1"], [])
        self.assertEqual(result["hap2"], [])


if __name__ == "__main__":
    unittest.main()
