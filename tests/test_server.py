"""Tests for the IGV.js visualization server."""

import json
import os
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

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

    def test_cram_in_reads_bam_reclassified(self):
        """A .cram path in reads_bam is reclassified as reads_cram."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv",
                                         delete=False) as f:
            f.write("#sample_id\treads_bam\n")
            f.write("S1\t/data/sample.cram\n")
            f.flush()
            samples = server_app.load_config(f.name)
        os.unlink(f.name)
        self.assertEqual(len(samples), 1)
        self.assertEqual(samples[0]["reads_bam"], "")
        self.assertEqual(samples[0]["reads_cram"], "/data/sample.cram")

    def test_bam_in_reads_bam_unchanged(self):
        """A .bam path in reads_bam is left unchanged."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv",
                                         delete=False) as f:
            f.write("#sample_id\treads_bam\n")
            f.write("S1\t/data/sample.bam\n")
            f.flush()
            samples = server_app.load_config(f.name)
        os.unlink(f.name)
        self.assertEqual(len(samples), 1)
        self.assertEqual(samples[0]["reads_bam"], "/data/sample.bam")
        self.assertNotIn("reads_cram", samples[0])

    def test_cram_reclassify_preserves_existing_reads_cram(self):
        """If reads_cram already exists, reads_bam .cram does not overwrite."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv",
                                         delete=False) as f:
            f.write("#sample_id\treads_bam\treads_cram\n")
            f.write("S1\t/data/extra.cram\t/data/primary.cram\n")
            f.flush()
            samples = server_app.load_config(f.name)
        os.unlink(f.name)
        self.assertEqual(len(samples), 1)
        self.assertEqual(samples[0]["reads_bam"], "")
        self.assertEqual(samples[0]["reads_cram"], "/data/primary.cram")

    def test_cram_ref_column_loaded(self):
        """cram_ref column is loaded from TSV."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv",
                                         delete=False) as f:
            f.write("#sample_id\treads_bam\tcram_ref\n")
            f.write("S1\t/data/sample.cram\t/ref/genome.fa.gz\n")
            f.flush()
            samples = server_app.load_config(f.name)
        os.unlink(f.name)
        self.assertEqual(len(samples), 1)
        self.assertEqual(samples[0]["cram_ref"], "/ref/genome.fa.gz")

    def test_cram_ref_column_empty_when_absent(self):
        """Missing cram_ref column results in empty string via .get()."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv",
                                         delete=False) as f:
            f.write("#sample_id\treads_bam\n")
            f.write("S1\t/data/sample.bam\n")
            f.flush()
            samples = server_app.load_config(f.name)
        os.unlink(f.name)
        self.assertEqual(samples[0].get("cram_ref", ""), "")


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
            "_hap2_to_hap1.bam",
            "_hap2_to_hap1.bam.bai",
            "_hap1_to_hap2.bam",
            "_hap1_to_hap2.bam.bai",
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
        self.assertIn("hap2_to_hap1_bam", sample)
        self.assertIn("hap1_to_hap2_bam", sample)
        self.assertIn("hap1_mapping_index", sample)
        self.assertIn("hap2_mapping_index", sample)

    def test_missing_output_dir(self):
        """Missing output_dir is handled gracefully."""
        sample = {"sample_id": "test", "output_dir": "/nonexistent"}
        result = server_app.discover_pipeline_files(sample)
        self.assertNotIn("hap1_to_ref_bam", result)

    def test_discovers_cram_files(self):
        """CRAM files are discovered alongside BAM files."""
        prefix = "test_sample"
        # Create CRAM files
        for suffix in ("_reads_to_hap1.cram", "_reads_to_hap1.cram.crai",
                        "_reads_to_hap2.cram", "_reads_to_hap2.cram.crai",
                        "_reads.cram", "_reads.cram.crai"):
            open(os.path.join(self.tmpdir, prefix + suffix), "w").close()
        sample = {"sample_id": "test", "output_dir": self.tmpdir}
        sample = server_app.discover_pipeline_files(sample)
        self.assertIn("reads_to_hap1_cram", sample)
        self.assertIn("reads_to_hap2_cram", sample)
        self.assertIn("reads_cram", sample)

    def test_discovers_hp_tagged_reads(self):
        """HP-tagged reads files (reads.hp.cram / reads.hp.bam) are discovered."""
        prefix = "test_sample"
        for suffix in ("_reads.hp.cram", "_reads.hp.bam"):
            open(os.path.join(self.tmpdir, prefix + suffix), "w").close()
        sample = {"sample_id": "test", "output_dir": self.tmpdir}
        sample = server_app.discover_pipeline_files(sample)
        self.assertIn("reads_hp_cram", sample)
        self.assertIn("reads_hp_bam", sample)


class TestLoadRegions(unittest.TestCase):
    """Tests for loading regions of interest."""

    def test_load_manifest_json(self):
        """Load regions from a manifest JSON file.

        start/end should now come from pos/size (actual SV region),
        not from the (potentially padded) ref_region field.
        The original ref_region is preserved as fasta_region.
        """
        manifest = {
            "description": "test",
            "variants": [
                {
                    "chrom": "chr1",
                    "pos": 1000,
                    "size": 500,
                    "genotype": "1|0",
                    "ref_region": "chr1:500-2000",
                    "fasta_region": "chr1:900-2100",
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
        # start/end come from pos/size, NOT from the padded ref_region
        self.assertEqual(regions[0]["start"], 1000)
        self.assertEqual(regions[0]["end"], 1500)
        self.assertEqual(regions[0]["size"], 500)
        self.assertIn("1|0", regions[0]["label"])
        # ref_region reflects actual SV; fasta_region uses explicit field
        self.assertEqual(regions[0]["ref_region"], "chr1:1000-1500")
        self.assertEqual(regions[0]["fasta_region"], "chr1:900-2100")

    def test_load_manifest_json_fasta_region_falls_back_to_ref_region(self):
        """When fasta_region is absent, fallback keeps backward compatibility."""
        manifest = {
            "description": "test",
            "variants": [
                {
                    "chrom": "chr2",
                    "pos": 2000,
                    "size": 400,
                    "genotype": "1|1",
                    "ref_region": "chr2:1500-2600",
                    "hap1_regions": [],
                    "hap2_regions": [],
                }
            ]
        }
        with tempfile.NamedTemporaryFile(mode="w", suffix=".json",
                                         delete=False) as f:
            json.dump(manifest, f)
            f.flush()
            regions = server_app.load_regions(f.name)
        os.unlink(f.name)
        self.assertEqual(regions[0]["fasta_region"], "chr2:1500-2600")

    def test_load_manifest_json_label_uses_sv_coords(self):
        """The region label shows actual SV coordinates, not padded ref_region."""
        manifest = {
            "description": "test",
            "variants": [
                {
                    "chrom": "chr1",
                    "pos": 460749,
                    "size": 834,
                    "genotype": "0|1",
                    "ref_region": "chr1:410749-511583",
                    "hap1_regions": [],
                    "hap2_regions": [],
                }
            ]
        }
        with tempfile.NamedTemporaryFile(mode="w", suffix=".json",
                                         delete=False) as f:
            json.dump(manifest, f)
            f.flush()
            regions = server_app.load_regions(f.name)
        os.unlink(f.name)
        r = regions[0]
        self.assertEqual(r["start"], 460749)
        self.assertEqual(r["end"], 461583)
        # The label must show the actual SV region
        self.assertIn("460749", r["label"])
        self.assertIn("461583", r["label"])
        # fasta_region is the original padded ref_region
        self.assertEqual(r["fasta_region"], "chr1:410749-511583")

    def test_load_vcf_regions_no_padding(self):
        """VCF regions store actual SV coords without 50 kb padding."""
        ref_allele = "A" * 835
        vcf_lines = [
            "##fileformat=VCFv4.2",
            "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE",
            # A deletion: REF is 835 bp, ALT is 1 bp -> sv_len = 834
            f"chr1	460749	.	{ref_allele}	A	.	.	SVLEN=834	GT	0|1",
        ]
        vcf_content = chr(10).join(vcf_lines) + chr(10)
        with tempfile.NamedTemporaryFile(mode="w", suffix=".vcf",
                                         delete=False) as f:
            f.write(vcf_content)
            f.flush()
            regions = server_app.load_regions(f.name)
        os.unlink(f.name)
        self.assertEqual(len(regions), 1)
        r = regions[0]
        self.assertEqual(r["chrom"], "chr1")
        # start/end must be actual SV coords — NO 50 kb padding
        self.assertEqual(r["start"], 460749)
        self.assertEqual(r["end"], 460749 + 834)
        self.assertEqual(r["size"], 834)
        self.assertEqual(r["ref_region"], "chr1:460749-461583")
        # VCF regions have no fasta_region (that is for toy-dataset manifests)
        self.assertNotIn("fasta_region", r)

    def test_load_nonexistent_file(self):
        """Non-existent file returns empty list."""
        self.assertEqual(server_app.load_regions("/nonexistent"), [])

    def test_load_empty_path(self):
        """Empty path returns empty list."""
        self.assertEqual(server_app.load_regions(""), [])

    def test_load_toy_manifest(self):
        """Load the actual toy manifest and get 10 regions."""
        manifest_path = os.path.join(
            _REPO_ROOT, "resources", "toy_dataset", "toy_manifest.json")
        if not os.path.isfile(manifest_path):
            self.skipTest("Toy manifest not found")
        regions = server_app.load_regions(manifest_path)
        self.assertEqual(len(regions), 10)
        # Verify regions have expected fields
        for r in regions:
            self.assertIn("chrom", r)
            self.assertIn("start", r)
            self.assertIn("end", r)
            self.assertIn("label", r)
            self.assertIn("ref_region", r)
            self.assertIn("fasta_region", r)
            self.assertTrue(r["end"] > r["start"])
            # ref_region must be the actual SV region "chrom:start-end"
            # (not a padded region); fasta_region holds the padded FASTA name
            chrom = r["chrom"]
            expected_ref = f"{chrom}:{r['start']}-{r['end']}"
            self.assertEqual(r["ref_region"], expected_ref,
                             "ref_region should be 'chrom:start-end' matching pos/size")


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

    def test_different_strands_not_merged(self):
        """Regions on opposite strands of the same contig are NOT merged.

        This is critical for inversions: a + and - strand alignment on the
        same contig represent opposite orientations and merging them would
        produce an incorrect, inflated assembly region.
        """
        regions = [
            {"chrom": "chr1", "start": 0, "end": 150, "strand": "+"},
            {"chrom": "chr1", "start": 100, "end": 300, "strand": "-"},
        ]
        merged = server_app._merge_regions(regions)
        self.assertEqual(len(merged), 2)
        strands = sorted(r["strand"] for r in merged)
        self.assertEqual(strands, ["+", "-"])

    def test_same_strand_merged_across_intervals(self):
        """Multiple overlapping + strand regions on the same contig merge."""
        regions = [
            {"chrom": "chr1", "start": 0, "end": 100, "strand": "+"},
            {"chrom": "chr1", "start": 50, "end": 200, "strand": "+"},
            {"chrom": "chr1", "start": 180, "end": 350, "strand": "+"},
        ]
        merged = server_app._merge_regions(regions)
        self.assertEqual(len(merged), 1)
        self.assertEqual(merged[0]["start"], 0)
        self.assertEqual(merged[0]["end"], 350)

    def test_mixed_contigs_and_strands(self):
        """Merge groups correctly by (contig, strand) tuple."""
        regions = [
            {"chrom": "chr1", "start": 0, "end": 100, "strand": "+"},
            {"chrom": "chr1", "start": 50, "end": 200, "strand": "+"},
            {"chrom": "chr1", "start": 0, "end": 100, "strand": "-"},
            {"chrom": "chr2", "start": 0, "end": 100, "strand": "+"},
        ]
        merged = server_app._merge_regions(regions)
        # chr1/+ → merged to 1, chr1/- → 1, chr2/+ → 1 = total 3
        self.assertEqual(len(merged), 3)


class TestFormatSize(unittest.TestCase):
    """Tests for human-readable size formatting."""

    def test_bp(self):
        self.assertEqual(server_app._format_size(500), "500 bp")

    def test_kb(self):
        self.assertEqual(server_app._format_size(5000), "5.0 kb")

    def test_mb(self):
        self.assertEqual(server_app._format_size(1500000), "1.5 Mb")


class TestGuessType(unittest.TestCase):
    """Tests for MIME type detection including CRAM."""

    def test_cram_type(self):
        self.assertEqual(
            server_app.IGVHandler._guess_type("sample.cram"),
            "application/octet-stream")

    def test_crai_type(self):
        self.assertEqual(
            server_app.IGVHandler._guess_type("sample.cram.crai"),
            "application/octet-stream")

    def test_bam_type(self):
        self.assertEqual(
            server_app.IGVHandler._guess_type("sample.bam"),
            "application/octet-stream")

    def test_fai_type(self):
        self.assertEqual(
            server_app.IGVHandler._guess_type("ref.fa.gz.fai"),
            "text/plain")


class TestValidateCramRef(unittest.TestCase):
    """Tests for CRAM reference validation."""

    def test_empty_path_no_warnings(self):
        """Empty cram_ref produces no warnings."""
        self.assertEqual(server_app._validate_cram_ref(""), [])

    def test_missing_file_warns(self):
        """Non-existent file produces a warning."""
        warnings = server_app._validate_cram_ref("/no/such/file.fa")
        self.assertTrue(any("not found" in w for w in warnings))

    def test_missing_fai_warns(self):
        """A FASTA without .fai index produces a warning."""
        with tempfile.NamedTemporaryFile(suffix=".fa", delete=False) as f:
            f.write(b">chr1\nACGT\n")
        try:
            warnings = server_app._validate_cram_ref(f.name)
            self.assertTrue(any(".fai" in w for w in warnings))
        finally:
            os.unlink(f.name)

    def test_bgzip_missing_gzi_warns(self):
        """A .fa.gz file without .gzi produces a warning."""
        import gzip as gz
        with tempfile.NamedTemporaryFile(suffix=".fa.gz",
                                          delete=False) as tmp:
            tmp_name = tmp.name
        with gz.open(tmp_name, "wb") as f:
            f.write(b">chr1\nACGT\n")
        try:
            warnings = server_app._validate_cram_ref(tmp_name)
            self.assertTrue(any(".gzi" in w for w in warnings))
        finally:
            os.unlink(tmp_name)

    def test_valid_ref_no_warnings(self):
        """A properly indexed reference produces no warnings."""
        with tempfile.NamedTemporaryFile(suffix=".fa", delete=False) as f:
            f.write(b">chr1\nACGT\n")
            fai_path = f.name + ".fai"
        open(fai_path, "w").close()
        try:
            warnings = server_app._validate_cram_ref(f.name)
            self.assertEqual(warnings, [])
        finally:
            os.unlink(f.name)
            os.unlink(fai_path)


class TestReadItf8(unittest.TestCase):
    """Tests for the ITF-8 integer decoding used by the CRAM parser."""

    def _decode(self, byte_values):
        """Helper: decode ITF-8 from a list of byte values."""
        import io
        return server_app._read_itf8(io.BytesIO(bytes(byte_values)))

    def test_single_byte(self):
        """Values 0-127 are encoded as a single byte."""
        self.assertEqual(self._decode([0]), 0)
        self.assertEqual(self._decode([42]), 42)
        self.assertEqual(self._decode([127]), 127)

    def test_two_byte(self):
        """Values 128-16383 use two bytes (high bit set)."""
        # 0x80 | 0x01 = 0x81, second byte 0x00  → (1 << 8) | 0 = 256
        self.assertEqual(self._decode([0x81, 0x00]), 256)

    def test_eof_raises(self):
        """Empty input raises EOFError."""
        import io
        with self.assertRaises(EOFError):
            server_app._read_itf8(io.BytesIO(b""))


class TestCheckCramEmbeddedRef(unittest.TestCase):
    """Tests for CRAM header UR path extraction."""

    def test_non_cram_file_returns_empty(self):
        """A non-CRAM file returns empty results."""
        with tempfile.NamedTemporaryFile(suffix=".cram", delete=False) as f:
            f.write(b"NOT_CRAM_DATA")
        try:
            ur_paths, warnings = server_app._check_cram_embedded_ref(f.name)
            self.assertEqual(ur_paths, set())
            self.assertEqual(warnings, [])
        finally:
            os.unlink(f.name)

    def test_nonexistent_file_returns_empty(self):
        """A non-existent file returns empty results."""
        ur_paths, warnings = server_app._check_cram_embedded_ref(
            "/no/such/file.cram")
        self.assertEqual(ur_paths, set())
        self.assertEqual(warnings, [])


class TestCoordinateTranslator(unittest.TestCase):
    """Tests for the coordinate translator with toy data."""

    @classmethod
    def setUpClass(cls):
        """Load the toy mapping indices if available."""
        toy_output = os.environ.get("TOY_OUTPUT_DIR", "/tmp/toy_output")
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


class TestTranslatorCache(unittest.TestCase):
    """Verify the CoordinateTranslator LRU cache."""

    def test_cache_returns_same_result(self):
        """Repeated identical queries return cached results."""
        translator = server_app.CoordinateTranslator()
        # No indices loaded — every translate returns empty lists.
        r1 = translator.translate("s", "chr1", 100, 200)
        r2 = translator.translate("s", "chr1", 100, 200)
        self.assertEqual(r1, r2)
        # Both should be the *same object* (cache hit).
        self.assertIs(r1, r2)

    def test_different_params_miss_cache(self):
        """Different coordinates are separate cache entries."""
        translator = server_app.CoordinateTranslator()
        r1 = translator.translate("s", "chr1", 100, 200)
        r2 = translator.translate("s", "chr1", 200, 300)
        self.assertIsNot(r1, r2)

    def test_cache_eviction(self):
        """Cache evicts oldest entries when full."""
        translator = server_app.CoordinateTranslator(cache_size=3)

        # Fill the cache
        translator.translate("s", "chr1", 0, 100)
        translator.translate("s", "chr1", 100, 200)
        translator.translate("s", "chr1", 200, 300)
        self.assertEqual(len(translator._cache), 3)

        # One more should evict the oldest
        translator.translate("s", "chr1", 300, 400)
        self.assertEqual(len(translator._cache), 3)

    def test_has_index_no_indices(self):
        """has_index returns False when no mapping indices are loaded."""
        translator = server_app.CoordinateTranslator()
        self.assertFalse(translator.has_index("s"))

    def test_has_index_with_loaded_sample(self):
        """has_index returns True after loading a sample with indices."""
        translator = server_app.CoordinateTranslator()
        # Simulate: only register in _indices dict directly
        translator._indices[("s", "hap1")] = {"chr1": {"blocks": []}}
        self.assertTrue(translator.has_index("s"))
        # Different sample should still be False
        self.assertFalse(translator.has_index("other"))


class TestApiDotplot(unittest.TestCase):
    """Tests for the _api_dotplot handler logic."""

    def _make_handler_class(self, samples, translator=None):
        """Create a handler class with the given samples list."""
        # We don't instantiate IGVHandler directly (it needs socket args),
        # so we test the method by constructing a minimal wrapper.
        if translator is None:
            translator = type("_NoopTranslator", (), {
                "has_index": lambda self, sid: False,
                "translate": lambda self, sid, c, s, e, **kw: {"hap1": [], "hap2": []},
            })()

        handler = type("FakeHandler", (), {
            "samples": samples,
            "MAX_DOTPLOT_REGION_BP": server_app.IGVHandler.MAX_DOTPLOT_REGION_BP,
            "_find_sample": lambda self, sid: next(
                (s for s in self.samples if s["sample_id"] == sid), None
            ),
        })()
        handler.translator = translator
        # Bind the unbound _api_dotplot method from IGVHandler
        import types
        handler._api_dotplot = types.MethodType(
            server_app.IGVHandler._api_dotplot, handler
        )
        return handler

    def test_missing_sample(self):
        """Unknown sample_id returns an error."""
        handler = self._make_handler_class([])
        result = handler._api_dotplot({"sample": ["no_such"]})
        self.assertIn("error", result)

    def test_empty_regions(self):
        """No regions produces empty dot plot results."""
        handler = self._make_handler_class([{
            "sample_id": "S1",
            "reference": "",
            "hap1_assembly": "",
            "hap2_assembly": "",
        }])
        result = handler._api_dotplot({"sample": ["S1"]})
        # Should succeed (no error) even without regions.
        self.assertNotIn("error", result)
        self.assertIn("hap1_vs_ref", result)
        self.assertIn("hap2_vs_ref", result)
        self.assertIn("hap1_vs_hap2", result)
        # All dot plots empty because no sequences
        self.assertEqual(len(result["hap1_vs_ref"]["forward"]), 0)

    def test_k_clamped(self):
        """K parameter is clamped to valid range."""
        handler = self._make_handler_class([{
            "sample_id": "S1",
            "reference": "",
            "hap1_assembly": "",
            "hap2_assembly": "",
        }])
        # k=200 should be clamped to 101
        result = handler._api_dotplot({"sample": ["S1"], "k": ["200"]})
        self.assertEqual(result["k"], 101)
        # k=-5 should be clamped to 1
        result = handler._api_dotplot({"sample": ["S1"], "k": ["-5"]})
        self.assertEqual(result["k"], 1)

    def test_buffer_parameter_expands_ref_region(self):
        """The buffer parameter expands the reference region symmetrically."""
        handler = self._make_handler_class([{
            "sample_id": "S1",
            "reference": "",
            "hap1_assembly": "",
            "hap2_assembly": "",
        }])
        # Viewport is chr1:10000-20000 (10 kb), buffer=1.0 → chr1:0-30000
        result = handler._api_dotplot({
            "sample": ["S1"],
            "ref": ["chr1:10000-20000"],
            "buffer": ["1.0"],
        })
        self.assertEqual(result["buffered_ref_region"], "chr1:0-30000")

    def test_buffer_zero_no_expansion(self):
        """buffer=0 uses the exact viewport region."""
        handler = self._make_handler_class([{
            "sample_id": "S1",
            "reference": "",
            "hap1_assembly": "",
            "hap2_assembly": "",
        }])
        result = handler._api_dotplot({
            "sample": ["S1"],
            "ref": ["chr1:10000-20000"],
            "buffer": ["0"],
        })
        self.assertEqual(result["buffered_ref_region"], "chr1:10000-20000")

    def test_buffer_clamped_to_max(self):
        """Buffer multiplier is clamped to [0, 5]."""
        handler = self._make_handler_class([{
            "sample_id": "S1",
            "reference": "",
            "hap1_assembly": "",
            "hap2_assembly": "",
        }])
        # buffer=10 → clamped to 5; viewport 10 kb, 5× on each side → 110 kb total
        result = handler._api_dotplot({
            "sample": ["S1"],
            "ref": ["chr1:100000-110000"],
            "buffer": ["10"],
        })
        self.assertEqual(result["buffered_ref_region"], "chr1:50000-160000")

    def test_event_highlight_offsets(self):
        """Event start/end produce correct offsets relative to buffered region."""
        handler = self._make_handler_class([{
            "sample_id": "S1",
            "reference": "",
            "hap1_assembly": "",
            "hap2_assembly": "",
        }])
        # Viewport chr1:10000-20000, buffer=1.0 → buffered chr1:0-30000
        # Event at chr1:12000-15000 → offsets 12000, 15000
        result = handler._api_dotplot({
            "sample": ["S1"],
            "ref": ["chr1:10000-20000"],
            "buffer": ["1.0"],
            "event_start": ["12000"],
            "event_end": ["15000"],
        })
        hl = result["event_highlight"]
        self.assertIsNotNone(hl)
        self.assertEqual(hl["start_offset"], 12000)
        self.assertEqual(hl["end_offset"], 15000)

    def test_event_highlight_absent_when_no_event(self):
        """event_highlight is None when no event_start/event_end provided."""
        handler = self._make_handler_class([{
            "sample_id": "S1",
            "reference": "",
            "hap1_assembly": "",
            "hap2_assembly": "",
        }])
        result = handler._api_dotplot({
            "sample": ["S1"],
            "ref": ["chr1:10000-20000"],
        })
        self.assertIsNone(result["event_highlight"])

    def test_viewport_highlights_present_for_buffered_ref_and_haps(self):
        """Viewport highlights map raw loci into extracted sequence offsets."""
        handler = self._make_handler_class([{
            "sample_id": "S1",
            "reference": "",
            "hap1_assembly": "",
            "hap2_assembly": "",
        }])
        result = handler._api_dotplot({
            "sample": ["S1"],
            "ref": ["chr1:10000-20000"],
            "hap1": ["asm1:500-1500"],
            "hap2": ["asm2:600-1600"],
            "buffer": ["1.0"],
        })
        hl = result["viewport_highlights"]
        self.assertIsNotNone(hl["ref"])
        self.assertEqual(hl["ref"]["start_offset"], 9999)
        self.assertEqual(hl["ref"]["end_offset"], 19999)
        self.assertIsNotNone(hl["hap1"])
        self.assertEqual(hl["hap1"]["start_offset"], 0)
        self.assertEqual(hl["hap1"]["end_offset"], 1000)
        self.assertIsNotNone(hl["hap2"])
        self.assertEqual(hl["hap2"]["start_offset"], 0)
        self.assertEqual(hl["hap2"]["end_offset"], 1000)

    def test_viewport_highlights_use_translated_hap_regions(self):
        """When translator is used, hap highlights use translated extraction coords."""
        class _FakeTranslator:
            def has_index(self, sid):
                return sid == "S1"

            def translate(self, sid, chrom, start, end, **kw):
                return {
                    "hap1": [{"chrom": "asm1", "start": 1000, "end": 3000, "strand": "+"}],
                    "hap2": [{"chrom": "asm2", "start": 2000, "end": 5000, "strand": "+"}],
                }

        handler = self._make_handler_class([{
            "sample_id": "S1",
            "reference": "",
            "hap1_assembly": "",
            "hap2_assembly": "",
        }], translator=_FakeTranslator())
        result = handler._api_dotplot({
            "sample": ["S1"],
            "ref": ["chr1:10000-20000"],
            "hap1": ["asm1:1200-2200"],
            "hap2": ["asm2:2500-3500"],
            "buffer": ["0"],
        })
        hl = result["viewport_highlights"]
        self.assertEqual(hl["hap1"]["start_offset"], 200)
        self.assertEqual(hl["hap1"]["end_offset"], 1200)
        self.assertEqual(hl["hap2"]["start_offset"], 500)
        self.assertEqual(hl["hap2"]["end_offset"], 1500)

    def test_translator_used_when_available(self):
        """When a mapping index exists, translator.translate() is called
        to derive hap regions from the buffered reference region."""
        class _FakeTranslator:
            def __init__(self):
                self.calls = []

            def has_index(self, sid):
                return sid == "S1"

            def translate(self, sid, chrom, start, end, **kw):
                self.calls.append((sid, chrom, start, end))
                return {
                    "hap1": [{"chrom": "asm1", "start": 500, "end": 1500, "strand": "+"}],
                    "hap2": [{"chrom": "asm2", "start": 600, "end": 1600, "strand": "+"}],
                }

        tr = _FakeTranslator()
        handler = self._make_handler_class([{
            "sample_id": "S1",
            "reference": "",
            "hap1_assembly": "",
            "hap2_assembly": "",
        }], translator=tr)
        result = handler._api_dotplot({
            "sample": ["S1"],
            "ref": ["chr1:10000-20000"],
            "buffer": ["0"],
        })
        # Translator should have been called with the (un-buffered) region
        self.assertEqual(len(tr.calls), 1)
        self.assertEqual(tr.calls[0], ("S1", "chr1", 10000, 20000))
        # Labels should reflect the translated hap regions
        self.assertEqual(result["labels"]["hap1"], "asm1:500-1500")
        self.assertEqual(result["labels"]["hap2"], "asm2:600-1600")

    def test_max_dotplot_region_cap(self):
        """Regions exceeding MAX_DOTPLOT_REGION_BP are clamped."""
        handler = self._make_handler_class([{
            "sample_id": "S1",
            "reference": "",
            "hap1_assembly": "",
            "hap2_assembly": "",
        }])
        # Viewport 200 kb, buffer=2.0 → would be 1 Mb total → exceeds 500 kb cap
        result = handler._api_dotplot({
            "sample": ["S1"],
            "ref": ["chr1:100000-300000"],
            "buffer": ["2.0"],
        })
        # The buffered region should be clamped to ~500 kb
        m = __import__("re").match(r"chr1:(\d+)-(\d+)", result["buffered_ref_region"])
        self.assertIsNotNone(m)
        actual_span = int(m.group(2)) - int(m.group(1))
        self.assertLessEqual(actual_span, 500_000)

    def test_translator_disjoint_regions_use_envelope(self):
        """Disjoint translated hits on same contig are merged for extraction."""
        class _FakeTranslator:
            def has_index(self, sid):
                return sid == "S1"

            def translate(self, sid, chrom, start, end, **kw):
                return {
                    "hap1": [
                        {"chrom": "asm1", "start": 100, "end": 200, "strand": "+"},
                        {"chrom": "asm1", "start": 500, "end": 800, "strand": "+"},
                    ],
                    "hap2": [
                        {"chrom": "asm2", "start": 1000, "end": 1100, "strand": "+"},
                    ],
                }

        tr = _FakeTranslator()
        handler = self._make_handler_class([{
            "sample_id": "S1",
            "reference": "",
            "hap1_assembly": "",
            "hap2_assembly": "",
        }], translator=tr)
        result = handler._api_dotplot({
            "sample": ["S1"],
            "ref": ["chr1:10000-20000"],
            "buffer": ["0"],
        })
        self.assertEqual(result["labels"]["hap1"], "asm1:100-800")
        self.assertEqual(result["labels"]["hap2"], "asm2:1000-1100")

    def test_dotplot_extraction_regions_are_normalized_for_faidx(self):
        """Dotplot extraction uses 1-based starts even if loci include 0."""
        handler = self._make_handler_class([{
            "sample_id": "S1",
            "reference": "/ref.fa.gz",
            "hap1_assembly": "/hap1.fa.gz",
            "hap2_assembly": "/hap2.fa.gz",
        }])

        calls = []

        def _fake_extract_sequence(path, region):
            calls.append((path, region))
            return "ACGT"

        with patch.object(server_app.dot_plot, "extract_sequence", side_effect=_fake_extract_sequence):
            result = handler._api_dotplot({
                "sample": ["S1"],
                "ref": ["chr1:0-100"],
                "hap1": ["ctg1:0-200"],
                "hap2": ["ctg2:0-300"],
                "buffer": ["0"],
            })

        self.assertEqual(calls, [
            ("/ref.fa.gz", "chr1:1-100"),
            ("/hap1.fa.gz", "ctg1:1-200"),
            ("/hap2.fa.gz", "ctg2:1-300"),
        ])
        self.assertEqual(result["labels"]["ref"], "chr1:0-100")
        self.assertEqual(result["labels"]["hap1"], "ctg1:0-200")
        self.assertEqual(result["labels"]["hap2"], "ctg2:0-300")
        self.assertIsNotNone(result["viewport_highlights"]["ref"])


class TestSelectDotplotRegion(unittest.TestCase):
    """Unit tests for dot-plot region selection helper."""

    def test_empty(self):
        self.assertIsNone(server_app._select_dotplot_region([]))

    def test_prefers_most_covered_bp_group(self):
        regions = [
            {"chrom": "ctgA", "start": 100, "end": 150, "strand": "+"},   # 50 bp
            {"chrom": "ctgA", "start": 300, "end": 350, "strand": "+"},   # 50 bp
            {"chrom": "ctgB", "start": 1000, "end": 1200, "strand": "+"},  # 200 bp
        ]
        best = server_app._select_dotplot_region(regions)
        self.assertIsNotNone(best)
        self.assertEqual(best["chrom"], "ctgB")
        self.assertEqual(best["start"], 1000)
        self.assertEqual(best["end"], 1200)

    def test_returns_group_envelope_for_disjoint_intervals(self):
        regions = [
            {"chrom": "ctgA", "start": 100, "end": 200, "strand": "+"},
            {"chrom": "ctgA", "start": 500, "end": 800, "strand": "+"},
            {"chrom": "ctgB", "start": 50, "end": 120, "strand": "+"},
        ]
        best = server_app._select_dotplot_region(regions)
        self.assertIsNotNone(best)
        self.assertEqual(best["chrom"], "ctgA")
        self.assertEqual(best["start"], 100)
        self.assertEqual(best["end"], 800)

    def test_ignores_invalid_interval_records(self):
        regions = [
            {"chrom": "ctgA", "start": 100, "end": 200, "strand": "+"},
            {"chrom": "ctgA", "start": "x", "end": 300, "strand": "+"},
            {"chrom": "ctgA", "start": False, "end": 300, "strand": "+"},
            {"chrom": "", "start": 0, "end": 100, "strand": "+"},
            {"chrom": "ctgA", "start": 500, "end": 500, "strand": "+"},
        ]
        best = server_app._select_dotplot_region(regions)
        self.assertIsNotNone(best)
        self.assertEqual(best["chrom"], "ctgA")
        self.assertEqual(best["start"], 100)
        self.assertEqual(best["end"], 200)

    def test_tie_prefers_tighter_envelope(self):
        regions = [
            {"chrom": "ctgA", "start": 0, "end": 100, "strand": "+"},
            {"chrom": "ctgA", "start": 200, "end": 300, "strand": "+"},
            {"chrom": "ctgB", "start": 0, "end": 100, "strand": "+"},
            {"chrom": "ctgB", "start": 100, "end": 200, "strand": "+"},
        ]
        best = server_app._select_dotplot_region(regions)
        self.assertIsNotNone(best)
        self.assertEqual(best["chrom"], "ctgB")
        self.assertEqual(best["start"], 0)
        self.assertEqual(best["end"], 200)

    def test_inversion_keeps_both_strands_together(self):
        """An inversion has flanking + hits and an inverted - hit on the
        same contig.  All three should be grouped together so the full
        envelope is returned."""
        regions = [
            {"chrom": "ctgA", "start": 1000, "end": 4000, "strand": "+"},
            {"chrom": "ctgA", "start": 4000, "end": 7000, "strand": "-"},
            {"chrom": "ctgA", "start": 7000, "end": 10000, "strand": "+"},
        ]
        best = server_app._select_dotplot_region(regions)
        self.assertIsNotNone(best)
        self.assertEqual(best["chrom"], "ctgA")
        self.assertEqual(best["start"], 1000)
        self.assertEqual(best["end"], 10000)

    def test_no_strand_in_result(self):
        """Result dict should not carry a strand key (grouping is by
        chrom only)."""
        regions = [
            {"chrom": "ctgA", "start": 0, "end": 100, "strand": "+"},
        ]
        best = server_app._select_dotplot_region(regions)
        self.assertNotIn("strand", best)

    def test_runaway_envelope_isolated_by_clustering(self):
        """A spurious distant hit on the same contig must not blow up the
        envelope.  The nearby cluster with more covered bp should win."""
        regions = [
            # True locus: 1 kb–11 kb (10 kb covered)
            {"chrom": "ctgA", "start": 1000, "end": 6000, "strand": "+"},
            {"chrom": "ctgA", "start": 6000, "end": 11000, "strand": "+"},
            # Spurious ALU hit far away (1 kb covered, 50 Mb away)
            {"chrom": "ctgA", "start": 50000000, "end": 50001000, "strand": "+"},
        ]
        best = server_app._select_dotplot_region(regions)
        self.assertIsNotNone(best)
        self.assertEqual(best["chrom"], "ctgA")
        self.assertEqual(best["start"], 1000)
        self.assertEqual(best["end"], 11000)
        # Envelope must NOT span 50 Mb
        self.assertLess(best["end"] - best["start"], 100000)

    def test_max_gap_parameter(self):
        """Intervals within max_gap should stay in one cluster; intervals
        beyond max_gap should split."""
        regions = [
            {"chrom": "ctgA", "start": 0, "end": 100, "strand": "+"},
            {"chrom": "ctgA", "start": 200, "end": 300, "strand": "+"},
        ]
        # With default max_gap (50000), both in one cluster
        best = server_app._select_dotplot_region(regions)
        self.assertEqual(best["start"], 0)
        self.assertEqual(best["end"], 300)

        # With tiny max_gap, they split into separate clusters
        best = server_app._select_dotplot_region(regions, max_gap=50)
        self.assertEqual(best["start"], 0)
        self.assertEqual(best["end"], 100)
        # Both clusters have 100 bp but first has tighter envelope


class TestTranslatorSelectionRobustness(unittest.TestCase):
    """Tests for robust translation hit filtering in CoordinateTranslator."""

    def test_translate_ignores_invalid_alignment_hits(self):
        class _StubCM:
            @staticmethod
            def query(idx, chrom, start, end, min_mapq=0):
                return [
                    {"event_type": "alignment", "asm_chrom": "ctgA", "asm_start": 100, "asm_end": 150, "strand": "+"},
                    {"event_type": "alignment", "asm_chrom": "", "asm_start": 10, "asm_end": 20, "strand": "+"},
                    {"event_type": "alignment", "asm_chrom": "ctgA", "asm_start": "100", "asm_end": 200, "strand": "+"},
                    {"event_type": "alignment", "asm_chrom": "ctgA", "asm_start": 220, "asm_end": "260", "strand": "+"},
                    {"event_type": "alignment", "asm_chrom": "ctgA", "asm_start": False, "asm_end": 260, "strand": "+"},
                    {"event_type": "alignment", "asm_chrom": "ctgA", "asm_start": 260, "asm_end": True, "strand": "+"},
                    {"event_type": "alignment", "asm_chrom": "ctgA", "asm_start": 300, "asm_end": 300, "strand": "+"},
                    {"event_type": "deletion", "asm_chrom": "ctgA", "asm_start": 500, "asm_end": 600, "strand": "+"},
                ]

        translator = server_app.CoordinateTranslator()
        translator._indices[("S1", "hap1")] = {"dummy": True}
        translator._indices[("S1", "hap2")] = {"dummy": True}

        from unittest.mock import patch
        with patch.object(server_app, "coordinate_mapper", _StubCM):
            result = translator.translate("S1", "chr1", 0, 1000)

        self.assertEqual(result["hap1"], [{"chrom": "ctgA", "start": 100, "end": 150, "strand": "+"}])
        self.assertEqual(result["hap2"], [{"chrom": "ctgA", "start": 100, "end": 150, "strand": "+"}])


class TestApiTranslate(unittest.TestCase):
    """Tests for the _api_translate handler logic."""

    def _make_handler(self):
        class _Translator:
            MAX_TRANSLATE_SPAN_BP = 2_000_000

            def __init__(self):
                self.calls = []

            def translate(self, sample_id, chrom, start, end, **kwargs):
                self.calls.append((
                    sample_id, chrom, start, end, kwargs.get("min_mapq", 0)
                ))
                return {"hap1": [], "hap2": []}

        handler = type("FakeHandler", (), {})()
        handler.translator = _Translator()
        import types
        handler._api_translate = types.MethodType(
            server_app.IGVHandler._api_translate, handler
        )
        return handler

    def test_invalid_min_mapq_returns_error(self):
        handler = self._make_handler()
        result = handler._api_translate({
            "sample": ["S1"],
            "chrom": ["chr1"],
            "start": ["10"],
            "end": ["20"],
            "min_mapq": ["abc"],
        })
        self.assertIn("error", result)
        self.assertIn("min_mapq", result["error"])

    def test_invalid_coordinate_range_returns_error(self):
        handler = self._make_handler()
        result = handler._api_translate({
            "sample": ["S1"],
            "chrom": ["chr1"],
            "start": ["50"],
            "end": ["10"],
        })
        self.assertEqual(result, {"error": "Invalid coordinate range"})

    def test_oversized_translate_is_guarded(self):
        handler = self._make_handler()
        result = handler._api_translate({
            "sample": ["S1"],
            "chrom": ["chr1"],
            "start": ["0"],
            "end": ["3000000"],
        })
        self.assertEqual(result["hap1"], [])
        self.assertEqual(result["hap2"], [])
        self.assertIn("warning", result)
        self.assertEqual(handler.translator.calls, [])

    def test_valid_translate_calls_translator(self):
        handler = self._make_handler()
        result = handler._api_translate({
            "sample": ["S1"],
            "chrom": ["chr1"],
            "start": ["100"],
            "end": ["200"],
            "min_mapq": ["5"],
        })
        self.assertEqual(result, {"hap1": [], "hap2": []})
        self.assertEqual(
            handler.translator.calls,
            [("S1", "chr1", 100, 200, 5)]
        )

    def test_zero_width_query_returns_error(self):
        """start == end is an empty half-open interval and should return an error."""
        handler = self._make_handler()
        result = handler._api_translate({
            "sample": ["S1"],
            "chrom": ["chr1"],
            "start": ["1000"],
            "end": ["1000"],
        })
        self.assertIn("error", result)
        self.assertIn("Empty", result["error"])
        # Translator must NOT be called for empty intervals
        self.assertEqual(handler.translator.calls, [])

    def test_span_guard_uses_half_open_convention(self):
        """The span guard uses half-open math (end - start), not inclusive (end - start + 1).

        With MAX_TRANSLATE_SPAN_BP = 2_000_000:
        - span = 2_000_000 - 0 = 2_000_000 (exactly at limit) → should pass through
        - span = 2_000_001 - 0 = 2_000_001 → should be guarded
        """
        handler = self._make_handler()
        # Exactly at limit: should call translator
        result = handler._api_translate({
            "sample": ["S1"],
            "chrom": ["chr1"],
            "start": ["0"],
            "end": ["2000000"],
        })
        self.assertNotIn("warning", result)
        self.assertEqual(len(handler.translator.calls), 1)

        # One past the limit: should be guarded
        handler.translator.calls.clear()
        result = handler._api_translate({
            "sample": ["S1"],
            "chrom": ["chr1"],
            "start": ["0"],
            "end": ["2000001"],
        })
        self.assertIn("warning", result)
        self.assertEqual(handler.translator.calls, [])


class TestApiSampleCramRef(unittest.TestCase):
    """Tests for explicit cram_ref URLs returned by _api_sample."""

    def _make_handler(self, samples):
        class _Translator:
            def has_index(self, _sample_id):
                return False

        handler = type("FakeHandler", (), {})()
        handler.samples = samples
        handler.translator = _Translator()
        handler.file_registry = {}
        handler._find_sample = lambda sid: next(
            (s for s in samples if s["sample_id"] == sid), None
        )
        import types
        handler._api_sample = types.MethodType(
            server_app.IGVHandler._api_sample, handler
        )
        return handler

    def test_cram_ref_gzi_url_included_when_present(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            ref = os.path.join(tmpdir, "ref.fa.gz")
            h1 = os.path.join(tmpdir, "h1.fa.gz")
            h2 = os.path.join(tmpdir, "h2.fa.gz")
            cram_ref = os.path.join(tmpdir, "cram_ref.fa.gz")
            for p in (ref, h1, h2, cram_ref):
                Path(p).touch()
                Path(p + ".fai").touch()
                Path(p + ".gzi").touch()

            samples = [{
                "sample_id": "S1",
                "reference": ref,
                "hap1_assembly": h1,
                "hap2_assembly": h2,
                "cram_ref": cram_ref,
            }]
            handler = self._make_handler(samples)
            result = handler._api_sample("S1")
            self.assertEqual(result["cram_ref"], "/data/S1/cram_ref.fa.gz")
            self.assertEqual(result["cram_ref_gzi"], "/data/S1/cram_ref.fa.gz.gzi")

    def test_cram_ref_gzi_url_none_when_missing(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            ref = os.path.join(tmpdir, "ref.fa")
            h1 = os.path.join(tmpdir, "h1.fa")
            h2 = os.path.join(tmpdir, "h2.fa")
            cram_ref = os.path.join(tmpdir, "cram_ref.fa.gz")
            for p in (ref, h1, h2, cram_ref):
                Path(p).touch()
                Path(p + ".fai").touch()

            samples = [{
                "sample_id": "S1",
                "reference": ref,
                "hap1_assembly": h1,
                "hap2_assembly": h2,
                "cram_ref": cram_ref,
            }]
            handler = self._make_handler(samples)
            result = handler._api_sample("S1")
            self.assertEqual(result["cram_ref"], "/data/S1/cram_ref.fa.gz")
            self.assertIsNone(result["cram_ref_gzi"])


class TestApiSampleCrossHaplotypeTracks(unittest.TestCase):
    """Tests for cross-haplotype BAM URLs returned by _api_sample."""

    def _make_handler(self, samples):
        class _Translator:
            def has_index(self, _sample_id):
                return False

        handler = type("FakeHandler", (), {})()
        handler.samples = samples
        handler.translator = _Translator()
        handler.file_registry = {}
        handler._find_sample = lambda sid: next(
            (s for s in samples if s["sample_id"] == sid), None
        )
        import types
        handler._api_sample = types.MethodType(
            server_app.IGVHandler._api_sample, handler
        )
        return handler

    def test_cross_haplotype_bams_in_tracks(self):
        """Cross-haplotype BAM URLs are included when files exist."""
        with tempfile.TemporaryDirectory() as tmpdir:
            ref = os.path.join(tmpdir, "ref.fa")
            h1 = os.path.join(tmpdir, "h1.fa")
            h2 = os.path.join(tmpdir, "h2.fa")
            hap2_to_hap1 = os.path.join(tmpdir, "hap2_to_hap1.bam")
            hap1_to_hap2 = os.path.join(tmpdir, "hap1_to_hap2.bam")
            for p in (ref, h1, h2, hap2_to_hap1, hap1_to_hap2):
                Path(p).touch()

            samples = [{
                "sample_id": "S1",
                "reference": ref,
                "hap1_assembly": h1,
                "hap2_assembly": h2,
                "hap2_to_hap1_bam": hap2_to_hap1,
                "hap1_to_hap2_bam": hap1_to_hap2,
            }]
            handler = self._make_handler(samples)
            result = handler._api_sample("S1")
            self.assertIn("hap2_to_hap1_bam", result["tracks"])
            self.assertIn("hap1_to_hap2_bam", result["tracks"])
            self.assertIsNotNone(result["tracks"]["hap2_to_hap1_bam"])
            self.assertIsNotNone(result["tracks"]["hap1_to_hap2_bam"])

    def test_cross_haplotype_bams_absent_gracefully(self):
        """When cross-haplotype BAMs don't exist, tracks are None."""
        with tempfile.TemporaryDirectory() as tmpdir:
            ref = os.path.join(tmpdir, "ref.fa")
            h1 = os.path.join(tmpdir, "h1.fa")
            h2 = os.path.join(tmpdir, "h2.fa")
            for p in (ref, h1, h2):
                Path(p).touch()

            samples = [{
                "sample_id": "S1",
                "reference": ref,
                "hap1_assembly": h1,
                "hap2_assembly": h2,
            }]
            handler = self._make_handler(samples)
            result = handler._api_sample("S1")
            self.assertIn("hap2_to_hap1_bam", result["tracks"])
            self.assertIn("hap1_to_hap2_bam", result["tracks"])
            self.assertIsNone(result["tracks"]["hap2_to_hap1_bam"])
            self.assertIsNone(result["tracks"]["hap1_to_hap2_bam"])


class TestFrontendMemoryGuards(unittest.TestCase):
    """Tests for memory-safety guards in the region navigation frontend."""

    def test_buffer_and_visibility_window_caps_present(self):
        index_path = os.path.join(_REPO_ROOT, "server", "static", "index.html")
        with open(index_path, encoding="utf-8") as fh:
            html = fh.read()

        self.assertIn("const MAX_REGION_BUFFER = 250000;", html)
        self.assertIn("const MAX_VISIBILITY_WINDOW = 500000;", html)
        self.assertIn("const buffer = Math.min(rawBuffer, MAX_REGION_BUFFER);", html)
        self.assertIn("width = Math.min(width, MAX_VISIBILITY_WINDOW);", html)

    def test_frontend_uses_stable_region_selection(self):
        index_path = os.path.join(_REPO_ROOT, "server", "static", "index.html")
        with open(index_path, encoding="utf-8") as fh:
            html = fh.read()

        self.assertIn("function selectBestRegion(regions, coordType = \"half-open\", maxGap = 50000)", html)
        self.assertIn("function selectBestLocus(loci)", html)
        self.assertIn("selectBestRegion(regions, \"inclusive\")", html)
        self.assertIn("selectBestRegion(result.hap1)", html)
        self.assertIn("selectBestRegion(result.hap2)", html)
        self.assertIn("selectBestLocus(region.hap1_regions)", html)
        self.assertIn("selectBestLocus(region.hap2_regions)", html)

    def test_frontend_has_cross_haplotype_tracks(self):
        """Frontend includes cross-haplotype track definitions."""
        index_path = os.path.join(_REPO_ROOT, "server", "static", "index.html")
        with open(index_path, encoding="utf-8") as fh:
            html = fh.read()
        # Hap2 → Hap1 track in hap1 panel
        self.assertIn("Hap2 → Hap1", html)
        self.assertIn("config.tracks.hap2_to_hap1_bam", html)
        # Hap1 → Hap2 track in hap2 panel
        self.assertIn("Hap1 → Hap2", html)
        self.assertIn("config.tracks.hap1_to_hap2_bam", html)


class TestFrontendDisplayToggles(unittest.TestCase):
    """Tests for alignment display toggles: soft clips, mismatches, display mode, and 3rd gen view.

    These verify the implementation pattern rather than just code presence:
    - State variables are declared with correct defaults.
    - Track configs use the state variables (not hardcoded values).
    - Update functions apply the double-set pattern: both the live track
      property (t.showSoftClips / t.showMismatches / etc.) AND the config
      object property (t.config.*) are updated so the setting persists across
      pan/zoom re-renders.
    - Buttons have the correct initial active/inactive CSS class.
    - The indel threshold is set to LONG_READ_INDEL_THRESHOLD (50 bp) —
      the 3-bp default is meaningless for PacBio/ONT, which have many
      sequencing-error indels in the 1–50 bp range.
    """

    @classmethod
    def setUpClass(cls):
        index_path = os.path.join(_REPO_ROOT, "server", "static", "index.html")
        with open(index_path, encoding="utf-8") as fh:
            cls.html = fh.read()

    # ── State variables ────────────────────────────────────────────────────

    def test_soft_clips_state_variable_initialised_false(self):
        """showSoftClips must default to false (hide clipped bases by default)."""
        self.assertIn("let showSoftClips = false;", self.html)

    def test_mismatches_state_variable_initialised_true(self):
        """showMismatches must default to true (show SNVs by default)."""
        self.assertIn("let showMismatches = true;", self.html)

    # ── Track configs use state variables ─────────────────────────────────

    def test_read_track_uses_soft_clip_state_variable(self):
        """makeReadTrack must use showSoftClips variable, not hardcoded false."""
        self.assertIn("showSoftClips: showSoftClips,", self.html)
        self.assertNotIn("showSoftClips: false,", self.html)

    def test_read_track_uses_mismatch_state_variable(self):
        """makeReadTrack must use showMismatches variable."""
        self.assertIn("showMismatches: showMismatches,", self.html)

    # ── Double-set pattern in update functions ─────────────────────────────

    def test_update_soft_clips_sets_track_property(self):
        """updateAllTrackSoftClipDisplay must set t.showSoftClips on the track."""
        self.assertIn("t.showSoftClips = showSoftClips;", self.html)

    def test_update_soft_clips_sets_config_property(self):
        """updateAllTrackSoftClipDisplay must also update t.config.showSoftClips.

        Setting only the live property causes toggles to revert when IGV.js
        re-reads the config object after a pan/zoom fetch.  The config must
        be kept in sync to make the change durable.
        """
        self.assertIn("t.config.showSoftClips = showSoftClips;", self.html)

    def test_update_mismatches_sets_track_property(self):
        """updateAllTrackMismatchDisplay must set t.showMismatches on the track."""
        self.assertIn("t.showMismatches = showMismatches;", self.html)

    def test_update_mismatches_sets_config_property(self):
        """updateAllTrackMismatchDisplay must also update t.config.showMismatches."""
        self.assertIn("t.config.showMismatches = showMismatches;", self.html)

    # ── Buttons existence and initial state ────────────────────────────────

    def test_soft_clips_button_exists(self):
        """toggleSoftClips button must be present in the toolbar."""
        self.assertIn('id="toggleSoftClips"', self.html)

    def test_soft_clips_button_initially_inactive(self):
        """toggleSoftClips button must NOT carry 'active' class initially.

        The default is showSoftClips=false, so the button should appear
        inactive (not pressed) to match the state.
        """
        # Locate the button element and confirm it has no 'active' class.
        import re
        btn_match = re.search(
            r'<button[^>]+id="toggleSoftClips"[^>]*>', self.html
        )
        self.assertIsNotNone(btn_match, "toggleSoftClips button not found")
        self.assertNotIn("active", btn_match.group())

    def test_mismatches_button_exists(self):
        """toggleMismatches button must be present in the toolbar."""
        self.assertIn('id="toggleMismatches"', self.html)

    def test_mismatches_button_initially_active(self):
        """toggleMismatches button must carry 'active' class initially.

        The default is showMismatches=true, so the button should appear
        active (pressed) to reflect the enabled state.
        """
        import re
        btn_match = re.search(
            r'<button[^>]+id="toggleMismatches"[^>]*>', self.html
        )
        self.assertIsNotNone(btn_match, "toggleMismatches button not found")
        self.assertIn("active", btn_match.group())

    # ── Click handlers call update functions ──────────────────────────────

    def test_soft_clips_handler_calls_update_function(self):
        """The soft-clips click handler must call updateAllTrackSoftClipDisplay()."""
        self.assertIn("updateAllTrackSoftClipDisplay();", self.html)

    def test_mismatches_handler_calls_update_function(self):
        """The mismatches click handler must call updateAllTrackMismatchDisplay()."""
        self.assertIn("updateAllTrackMismatchDisplay();", self.html)

    # ── Display mode (squished / expanded) ────────────────────────────────

    def test_display_squished_state_variable_initialised_true(self):
        """displaySquished must default to true (squished is the long-read default)."""
        self.assertIn("let displaySquished = true;", self.html)

    def test_read_track_uses_display_mode_state_variable(self):
        """makeReadTrack must derive displayMode from displaySquished, not hardcode it.

        The local 'readDisplayMode' variable is computed from the state and
        passed into the track config so the initial panel creation respects
        any toggle state set before tracks are loaded.
        """
        self.assertIn(
            'const readDisplayMode = displaySquished ? "SQUISHED" : "EXPANDED";',
            self.html,
        )
        self.assertIn("displayMode: readDisplayMode,", self.html)

    def test_update_display_mode_sets_track_property(self):
        """updateAllTrackDisplayModes must set t.displayMode on the live track."""
        self.assertIn("t.displayMode = mode;", self.html)

    def test_update_display_mode_sets_config_property(self):
        """updateAllTrackDisplayModes must also update t.config.displayMode.

        Without the config update, the display mode reverts to the original
        value whenever IGV.js re-reads the config during a pan/zoom re-render.
        """
        self.assertIn("if (t.config) t.config.displayMode = mode;", self.html)

    def test_display_mode_button_exists(self):
        """toggleDisplayMode button must be present in the toolbar."""
        self.assertIn('id="toggleDisplayMode"', self.html)

    def test_display_mode_button_initially_inactive(self):
        """toggleDisplayMode button must NOT carry 'active' class initially.

        The default is displaySquished=true.  The handler applies 'active'
        only when expanded (i.e. classList.toggle('active', !displaySquished)),
        so on first load the button must be inactive.
        """
        import re
        btn_match = re.search(
            r'<button[^>]+id="toggleDisplayMode"[^>]*>', self.html
        )
        self.assertIsNotNone(btn_match, "toggleDisplayMode button not found")
        self.assertNotIn("active", btn_match.group())

    def test_display_mode_handler_calls_update_function(self):
        """The display-mode click handler must call updateAllTrackDisplayModes()."""
        self.assertIn("updateAllTrackDisplayModes();", self.html)

    def test_display_mode_active_class_tied_to_expanded_not_squished(self):
        """The 'active' class must be applied when EXPANDED, not when SQUISHED.

        The handler uses classList.toggle('active', !displaySquished) so the
        button is highlighted when the user switches to expanded view.
        This is the semantically correct UX (button lit = non-default state).
        """
        self.assertIn(
            "toggleDisplayModeBtn.classList.toggle(\"active\", !displaySquished);",
            self.html,
        )

    # ── 3rd gen view (hide small indels) ──────────────────────────────────

    def test_indel_threshold_constant_is_fifty(self):
        """LONG_READ_INDEL_THRESHOLD must be 50 bp.

        PacBio/ONT reads carry many sequencing-error indels in the 1–50 bp
        range.  A threshold of 3 bp (the short-read default) hides almost
        nothing on long reads, making the toggle visually ineffective.
        50 bp suppresses most noise while preserving true structural variants.
        """
        self.assertIn("const LONG_READ_INDEL_THRESHOLD = 50;", self.html)

    def test_indel_threshold_not_hardcoded_to_three(self):
        """smallIndelThreshold must use the named constant, not the literal '3'.

        Hardcoding '3' obscures intent and makes the value trivially easy to
        overlook when reviewing long-read best practices.
        """
        self.assertNotIn("smallIndelThreshold: 3,", self.html)

    def test_read_track_uses_indel_threshold_constant(self):
        """makeReadTrack must use LONG_READ_INDEL_THRESHOLD for smallIndelThreshold."""
        self.assertIn("smallIndelThreshold: LONG_READ_INDEL_THRESHOLD,", self.html)

    def test_update_indels_sets_track_property_with_constant(self):
        """updateAllTrackIndelDisplay must set t.smallIndelThreshold using the constant."""
        self.assertIn("t.smallIndelThreshold = LONG_READ_INDEL_THRESHOLD;", self.html)

    def test_update_indels_sets_config_property_with_constant(self):
        """updateAllTrackIndelDisplay must also update t.config.smallIndelThreshold.

        Without the config update, the threshold reverts on pan/zoom re-renders
        (same async-overwrite problem as soft clips/mismatches).
        """
        self.assertIn(
            "t.config.smallIndelThreshold = LONG_READ_INDEL_THRESHOLD;", self.html
        )

    def test_indels_button_labeled_3rd_gen_view_when_active(self):
        """When hideSmallIndels=true (active), the button must say '3rd gen view'.

        '3rd gen view' clearly signals that this setting is optimised for
        third-generation (long-read) sequencing noise filtering, which is
        more informative than the generic 'Hide indels' label.
        """
        self.assertIn("3rd gen view", self.html)

    def test_indels_button_labeled_all_indels_when_inactive(self):
        """When hideSmallIndels=false (inactive), the button must say 'All indels'."""
        self.assertIn("All indels", self.html)

    def test_indels_button_initially_active(self):
        """toggleSmallIndels button must carry 'active' class initially.

        The default is hideSmallIndels=true (filter on by default) so the
        button should appear active (pressed) to reflect the enabled state.
        """
        import re
        btn_match = re.search(
            r'<button[^>]+id="toggleSmallIndels"[^>]*>', self.html
        )
        self.assertIsNotNone(btn_match, "toggleSmallIndels button not found")
        self.assertIn("active", btn_match.group())

    def test_indels_handler_calls_update_function(self):
        """The 3rd-gen-view click handler must call updateAllTrackIndelDisplay()."""
        self.assertIn("updateAllTrackIndelDisplay();", self.html)

    def test_hide_small_indels_state_variable_initialised_true(self):
        """hideSmallIndels must default to true (filter on for long-read mode)."""
        self.assertIn("let hideSmallIndels = true;", self.html)

    def test_update_indels_sets_hide_small_indels_track_property(self):
        """updateAllTrackIndelDisplay must set t.hideSmallIndels on the live track."""
        self.assertIn("t.hideSmallIndels = hideSmallIndels;", self.html)

    def test_update_indels_sets_hide_small_indels_config_property(self):
        """updateAllTrackIndelDisplay must also update t.config.hideSmallIndels."""
        self.assertIn("t.config.hideSmallIndels = hideSmallIndels;", self.html)


if __name__ == "__main__":
    unittest.main()
