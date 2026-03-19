"""Tests for the server configuration generation script."""

import os
import sys
import tempfile
import shutil
import unittest

# Allow importing from scripts/ directory
_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(_REPO_ROOT, "scripts"))

import generate_server_config as gen_config


class TestDetectSamplePrefix(unittest.TestCase):
    """Tests for sample prefix detection."""

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmpdir)

    def test_detects_prefix(self):
        """Prefix is extracted from the *_hap1_to_ref.bam file."""
        open(os.path.join(self.tmpdir, "NA21110_hap1_to_ref.bam"), "w").close()
        prefix = gen_config.detect_sample_prefix(self.tmpdir)
        self.assertEqual(prefix, "NA21110")

    def test_no_pipeline_files(self):
        """Returns None when no pipeline files are found."""
        open(os.path.join(self.tmpdir, "random.txt"), "w").close()
        prefix = gen_config.detect_sample_prefix(self.tmpdir)
        self.assertIsNone(prefix)


class TestBuildSampleRow(unittest.TestCase):
    """Tests for building sample configuration rows."""

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        # Create mock pipeline files
        prefix = "test"
        for suffix in ["_hap1_to_ref.bam", "_hap2_to_ref.bam",
                        "_reads_to_hap1.bam", "_reads_to_hap2.bam"]:
            open(os.path.join(self.tmpdir, prefix + suffix), "w").close()

    def tearDown(self):
        shutil.rmtree(self.tmpdir)

    def test_basic_row(self):
        """Build a sample row with explicitly provided paths."""
        row = gen_config.build_sample_row(
            self.tmpdir, "test",
            reference="/ref.fa",
            hap1="/h1.fa",
            hap2="/h2.fa",
        )
        self.assertEqual(row["sample_id"], "test")
        self.assertEqual(row["reference"], "/ref.fa")
        self.assertEqual(row["hap1_assembly"], "/h1.fa")
        self.assertEqual(row["hap2_assembly"], "/h2.fa")
        self.assertTrue(os.path.isabs(row["output_dir"]))


class TestWriteConfig(unittest.TestCase):
    """Tests for writing the TSV configuration file."""

    def test_writes_valid_tsv(self):
        """Written TSV can be re-loaded."""
        samples = [{
            "sample_id": "S1",
            "output_dir": "/out",
            "reference": "/ref.fa",
            "hap1_assembly": "/h1.fa",
            "hap2_assembly": "/h2.fa",
            "reads_bam": "",
            "regions": "/r.json",
        }]
        with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv",
                                         delete=False) as f:
            gen_config.write_config(samples, f.name)

        # Read back
        with open(f.name) as fh:
            lines = fh.readlines()
        os.unlink(f.name)

        self.assertTrue(lines[0].startswith("#"))
        self.assertEqual(len(lines), 2)  # header + 1 sample
        fields = lines[1].strip().split("\t")
        self.assertEqual(fields[0], "S1")


class TestValidateSample(unittest.TestCase):
    """Tests for sample validation warnings."""

    def test_missing_reference_warns(self):
        """Missing reference generates a warning."""
        row = {
            "sample_id": "test",
            "output_dir": "/tmp",
            "reference": "",
            "hap1_assembly": "",
            "hap2_assembly": "",
        }
        warnings = gen_config.validate_sample(row)
        self.assertTrue(any("reference" in w for w in warnings))


class TestScanParentDir(unittest.TestCase):
    """Tests for scanning parent directories."""

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmpdir)

    def test_finds_subdirectories(self):
        """Finds sample subdirectories with pipeline outputs."""
        subdir = os.path.join(self.tmpdir, "sample1")
        os.makedirs(subdir)
        open(os.path.join(subdir, "S1_hap1_to_ref.bam"), "w").close()

        results = gen_config.scan_parent_dir(self.tmpdir, "/ref.fa")
        self.assertEqual(len(results), 1)
        self.assertEqual(results[0][1], "S1")

    def test_ignores_non_pipeline_dirs(self):
        """Directories without pipeline outputs are skipped."""
        subdir = os.path.join(self.tmpdir, "other")
        os.makedirs(subdir)
        open(os.path.join(subdir, "random.txt"), "w").close()

        results = gen_config.scan_parent_dir(self.tmpdir, "/ref.fa")
        self.assertEqual(len(results), 0)


if __name__ == "__main__":
    unittest.main()
