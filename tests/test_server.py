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

    def _make_handler_class(self, samples):
        """Create a handler class with the given samples list."""
        # We don't instantiate IGVHandler directly (it needs socket args),
        # so we test the method by constructing a minimal wrapper.
        handler = type("FakeHandler", (), {
            "samples": samples,
            "_find_sample": lambda self, sid: next(
                (s for s in self.samples if s["sample_id"] == sid), None
            ),
        })()
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
                open(p, "w").close()
                open(p + ".fai", "w").close()
                open(p + ".gzi", "w").close()

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
                open(p, "w").close()
                open(p + ".fai", "w").close()

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


if __name__ == "__main__":
    unittest.main()
