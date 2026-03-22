#!/usr/bin/env python3
"""Tests for generate_toy_dataset.py"""

import gzip
import json
import os
import shutil
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import generate_toy_dataset  # noqa: E402
import coordinate_mapper  # noqa: E402

# Minimal VCF data for testing.  The REF/ALT lengths encode the deletion
# size; positions and allele sequences are synthetic.
VCF_HEADER = """\
##fileformat=VCFv4.2
##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE
"""

# Deletion sizes: 5000, 3000, 2000, 1500, 500 (too small), 4000, 8000
# Genotypes: some carried (1|0, 0|1, 1|1) some not (0|0)
VCF_RECORDS = [
    # chr1 - 5000bp del, carried (1|0)
    "chr1\t100000\t.\t{ref}\t{alt}\t.\t.\tAC=1\tGT\t1|0".format(
        ref="A" * 5001, alt="A"
    ),
    # chr1 - 3000bp del, carried (0|1)
    "chr1\t500000\t.\t{ref}\t{alt}\t.\t.\tAC=1\tGT\t0|1".format(
        ref="G" * 3001, alt="G"
    ),
    # chr2 - 2000bp del, carried (1|1)
    "chr2\t200000\t.\t{ref}\t{alt}\t.\t.\tAC=2\tGT\t1|1".format(
        ref="C" * 2001, alt="C"
    ),
    # chr3 - 1500bp del, NOT carried (0|0) — should be skipped
    "chr3\t300000\t.\t{ref}\t{alt}\t.\t.\tAC=0\tGT\t0|0".format(
        ref="T" * 1501, alt="T"
    ),
    # chr4 - 500bp del, carried but too small for default min_size=1000
    "chr4\t400000\t.\t{ref}\t{alt}\t.\t.\tAC=1\tGT\t1|0".format(
        ref="A" * 501, alt="A"
    ),
    # chr5 - 4000bp del, carried (1|1)
    "chr5\t600000\t.\t{ref}\t{alt}\t.\t.\tAC=2\tGT\t1|1".format(
        ref="G" * 4001, alt="G"
    ),
    # chr6 - 8000bp del, carried (0|1)
    "chr6\t700000\t.\t{ref}\t{alt}\t.\t.\tAC=1\tGT\t0|1".format(
        ref="C" * 8001, alt="C"
    ),
]


def _write_vcf(path, records=None, compress=False):
    """Write a test VCF file."""
    if records is None:
        records = VCF_RECORDS
    content = VCF_HEADER + "\n".join(records) + "\n"
    if compress:
        with gzip.open(path, "wt") as fh:
            fh.write(content)
    else:
        with open(path, "w") as fh:
            fh.write(content)


class TestParseVcfDeletions(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmpdir)

    def test_basic_parsing(self):
        """Finds carried deletions >= 1000 bp."""
        vcf_path = os.path.join(self.tmpdir, "test.vcf")
        _write_vcf(vcf_path)
        dels = generate_toy_dataset.parse_vcf_deletions(vcf_path, min_size=1000)
        # Expected: chr1:100000 (5000), chr1:500000 (3000), chr2:200000 (2000),
        #           chr5:600000 (4000), chr6:700000 (8000)
        # NOT: chr3 (0|0), chr4 (500bp)
        self.assertEqual(len(dels), 5)
        # Sorted by size descending
        sizes = [d["size"] for d in dels]
        self.assertEqual(sizes, sorted(sizes, reverse=True))

    def test_gzipped_vcf(self):
        """Handles bgzipped VCF files."""
        vcf_path = os.path.join(self.tmpdir, "test.vcf.gz")
        _write_vcf(vcf_path, compress=True)
        dels = generate_toy_dataset.parse_vcf_deletions(vcf_path, min_size=1000)
        self.assertEqual(len(dels), 5)

    def test_min_size_filter(self):
        """Respects min_size parameter."""
        vcf_path = os.path.join(self.tmpdir, "test.vcf")
        _write_vcf(vcf_path)
        # With min_size=100, should also pick up the 500bp deletion
        dels = generate_toy_dataset.parse_vcf_deletions(vcf_path, min_size=100)
        self.assertEqual(len(dels), 6)

    def test_empty_vcf(self):
        """Returns empty list for header-only VCF."""
        vcf_path = os.path.join(self.tmpdir, "test.vcf")
        _write_vcf(vcf_path, records=[])
        dels = generate_toy_dataset.parse_vcf_deletions(vcf_path, min_size=1000)
        self.assertEqual(len(dels), 0)

    def test_deletion_fields(self):
        """Parsed deletion dicts have correct fields."""
        vcf_path = os.path.join(self.tmpdir, "test.vcf")
        _write_vcf(vcf_path)
        dels = generate_toy_dataset.parse_vcf_deletions(vcf_path, min_size=1000)
        for d in dels:
            self.assertIn("chrom", d)
            self.assertIn("pos", d)
            self.assertIn("end", d)
            self.assertIn("size", d)
            self.assertIn("genotype", d)
            self.assertEqual(d["end"], d["pos"] + d["size"] + 1)


class TestSelectVariants(unittest.TestCase):
    def _make_deletions(self, specs):
        """Create deletion dicts from (chrom, pos, size) tuples."""
        return [
            {"chrom": c, "pos": p, "end": p + s, "size": s, "genotype": "1|0"}
            for c, p, s in specs
        ]

    def test_fewer_than_requested(self):
        """Returns all when fewer than num_variants available."""
        dels = self._make_deletions([("chr1", 100, 5000), ("chr2", 200, 3000)])
        selected = generate_toy_dataset.select_variants(dels, num_variants=10)
        self.assertEqual(len(selected), 2)

    def test_exact_count(self):
        """Returns exactly num_variants."""
        dels = self._make_deletions(
            [
                ("chr1", 100, 5000),
                ("chr2", 200, 3000),
                ("chr3", 300, 2000),
                ("chr4", 400, 1500),
                ("chr5", 500, 1000),
            ]
        )
        # Sort by size desc as parse_vcf_deletions does
        dels.sort(key=lambda d: -d["size"])
        selected = generate_toy_dataset.select_variants(dels, num_variants=3)
        self.assertEqual(len(selected), 3)

    def test_chromosome_diversity(self):
        """Prefers different chromosomes."""
        dels = self._make_deletions(
            [
                ("chr1", 100, 9000),
                ("chr1", 200, 8000),
                ("chr2", 300, 7000),
                ("chr3", 400, 6000),
                ("chr4", 500, 5000),
            ]
        )
        dels.sort(key=lambda d: -d["size"])
        selected = generate_toy_dataset.select_variants(dels, num_variants=4)
        self.assertEqual(len(selected), 4)
        chroms = {d["chrom"] for d in selected}
        self.assertEqual(len(chroms), 4)

    def test_sorted_output(self):
        """Selected variants are sorted by chrom, pos."""
        dels = self._make_deletions(
            [
                ("chr5", 500, 5000),
                ("chr1", 100, 9000),
                ("chr3", 300, 7000),
            ]
        )
        dels.sort(key=lambda d: -d["size"])
        selected = generate_toy_dataset.select_variants(dels, num_variants=3)
        positions = [(d["chrom"], d["pos"]) for d in selected]
        self.assertEqual(positions, sorted(positions))


class TestCollapseAsmRegions(unittest.TestCase):
    def test_non_overlapping(self):
        hits = [
            {"asm_chrom": "ctg1", "asm_start": 100, "asm_end": 200},
            {"asm_chrom": "ctg1", "asm_start": 300, "asm_end": 400},
        ]
        result = generate_toy_dataset._collapse_asm_regions(hits)
        self.assertEqual(result, ["ctg1:100-200", "ctg1:300-400"])

    def test_overlapping(self):
        hits = [
            {"asm_chrom": "ctg1", "asm_start": 100, "asm_end": 300},
            {"asm_chrom": "ctg1", "asm_start": 200, "asm_end": 400},
        ]
        result = generate_toy_dataset._collapse_asm_regions(hits)
        self.assertEqual(result, ["ctg1:100-400"])

    def test_multiple_contigs(self):
        hits = [
            {"asm_chrom": "ctg2", "asm_start": 500, "asm_end": 600},
            {"asm_chrom": "ctg1", "asm_start": 100, "asm_end": 200},
        ]
        result = generate_toy_dataset._collapse_asm_regions(hits)
        self.assertEqual(result, ["ctg1:100-200", "ctg2:500-600"])

    def test_reverse_strand_coords(self):
        """Handles reversed start/end from negative strand alignments."""
        hits = [
            {"asm_chrom": "ctg1", "asm_start": 400, "asm_end": 100},
        ]
        result = generate_toy_dataset._collapse_asm_regions(hits)
        self.assertEqual(result, ["ctg1:100-400"])

    def test_empty(self):
        result = generate_toy_dataset._collapse_asm_regions([])
        self.assertEqual(result, [])


class TestComputeRegions(unittest.TestCase):
    """Test compute_regions with a mock coordinate mapping index."""

    def setUp(self):
        # Build a simple index with one block per chromosome per haplotype.
        self.tmpdir = tempfile.mkdtemp()

        # Import coordinate_mapper for index building
        import coordinate_mapper

        blocks_hap1 = [
            {
                "ref_chrom": "chr1",
                "ref_start": 0,
                "ref_end": 1000000,
                "asm_chrom": "h1_ctg1",
                "asm_start": 0,
                "asm_end": 1000000,
                "strand": "+",
                "mapq": 60,
                "matches": 999000,
                "block_len": 1000000,
            },
        ]
        blocks_hap2 = [
            {
                "ref_chrom": "chr1",
                "ref_start": 0,
                "ref_end": 1000000,
                "asm_chrom": "h2_ctg1",
                "asm_start": 5000,
                "asm_end": 1005000,
                "strand": "+",
                "mapq": 60,
                "matches": 999000,
                "block_len": 1000000,
            },
        ]

        hap1_path = os.path.join(self.tmpdir, "hap1.mapping.json.gz")
        hap2_path = os.path.join(self.tmpdir, "hap2.mapping.json.gz")
        coordinate_mapper.build_json_index(blocks_hap1, hap1_path)
        coordinate_mapper.build_json_index(blocks_hap2, hap2_path)

        self.hap1_index = coordinate_mapper.load_index(hap1_path)
        self.hap2_index = coordinate_mapper.load_index(hap2_path)

    def tearDown(self):
        shutil.rmtree(self.tmpdir)

    def test_basic_region_computation(self):
        variants = [
            {"chrom": "chr1", "pos": 500000, "end": 510000, "size": 10000,
             "genotype": "1|0"},
        ]
        regions = generate_toy_dataset.compute_regions(
            variants, self.hap1_index, self.hap2_index, padding=1000
        )
        self.assertEqual(len(regions), 1)
        r = regions[0]
        self.assertEqual(r["ref_region"], "chr1:499000-511000")
        self.assertTrue(len(r["hap1_regions"]) > 0)
        self.assertTrue(len(r["hap2_regions"]) > 0)

    def test_no_mapping_for_unknown_chrom(self):
        variants = [
            {"chrom": "chrX", "pos": 100, "end": 200, "size": 100,
             "genotype": "1|0"},
        ]
        regions = generate_toy_dataset.compute_regions(
            variants, self.hap1_index, self.hap2_index, padding=100
        )
        self.assertEqual(len(regions), 1)
        self.assertEqual(regions[0]["hap1_regions"], [])
        self.assertEqual(regions[0]["hap2_regions"], [])


class TestWriteManifest(unittest.TestCase):
    def test_manifest_json(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            regions = [
                {
                    "ref_region": "chr1:1000-2000",
                    "hap1_regions": ["ctg1:100-200"],
                    "hap2_regions": ["ctg2:300-400"],
                    "variant": {
                        "chrom": "chr1",
                        "pos": 1500,
                        "size": 500,
                        "genotype": "1|0",
                    },
                }
            ]
            path = generate_toy_dataset.write_manifest(regions, tmpdir)
            self.assertTrue(os.path.exists(path))

            with open(path) as fh:
                data = json.load(fh)
            self.assertIn("variants", data)
            self.assertEqual(len(data["variants"]), 1)
            self.assertEqual(data["variants"][0]["chrom"], "chr1")


class TestRealVcf(unittest.TestCase):
    """Test against the real bundled VCF file."""

    VCF_PATH = os.path.join(
        os.path.dirname(__file__),
        "..",
        "resources",
        "NA21110.shapeit5-phased-callset_final-vcf.phased.vcf.gz",
    )

    @unittest.skipUnless(
        os.path.exists(
            os.path.join(
                os.path.dirname(__file__),
                "..",
                "resources",
                "NA21110.shapeit5-phased-callset_final-vcf.phased.vcf.gz",
            )
        ),
        "Bundled VCF not found",
    )
    def test_real_vcf_has_deletions(self):
        """Bundled NA21110 VCF contains large deletions."""
        dels = generate_toy_dataset.parse_vcf_deletions(self.VCF_PATH)
        self.assertGreater(len(dels), 100)

    @unittest.skipUnless(
        os.path.exists(
            os.path.join(
                os.path.dirname(__file__),
                "..",
                "resources",
                "NA21110.shapeit5-phased-callset_final-vcf.phased.vcf.gz",
            )
        ),
        "Bundled VCF not found",
    )
    def test_select_10_from_real_vcf(self):
        """Can select 10 diverse variants from the real VCF."""
        dels = generate_toy_dataset.parse_vcf_deletions(self.VCF_PATH)
        selected = generate_toy_dataset.select_variants(dels, num_variants=10)
        self.assertEqual(len(selected), 10)
        # Should span multiple chromosomes
        chroms = {d["chrom"] for d in selected}
        self.assertGreater(len(chroms), 1)


class TestRemapSaTag(unittest.TestCase):
    """Tests for _remap_sa_tag SA supplementary-alignment tag remapping."""

    def test_basic_remap(self):
        """Remaps SA entry to the correct toy contig."""
        region_map = {"chr1": [("chr1:1000-2000", 1000, 2000)]}
        result = generate_toy_dataset._remap_sa_tag(
            "SA:Z:chr1,1500,+,50M,60,0;", region_map
        )
        self.assertEqual(result, "SA:Z:chr1:1000-2000,501,+,50M,60,0;")

    def test_multiple_entries(self):
        """Remaps multiple SA entries in a single tag."""
        region_map = {
            "chr1": [("chr1:1000-2000", 1000, 2000)],
            "chr2": [("chr2:5000-6000", 5000, 6000)],
        }
        result = generate_toy_dataset._remap_sa_tag(
            "SA:Z:chr1,1200,+,30M,50,1;chr2,5500,-,40M,55,2;",
            region_map,
        )
        self.assertIn("chr1:1000-2000,201,+,30M,50,1", result)
        self.assertIn("chr2:5000-6000,501,-,40M,55,2", result)

    def test_drops_unknown_chromosome(self):
        """Drops SA entries on chromosomes not in the toy reference."""
        region_map = {"chr1": [("chr1:1000-2000", 1000, 2000)]}
        result = generate_toy_dataset._remap_sa_tag(
            "SA:Z:chr1,1500,+,50M,60,0;chrX,100,+,20M,30,0;",
            region_map,
        )
        # chrX entry should be dropped
        self.assertNotIn("chrX", result)
        self.assertIn("chr1:1000-2000,501,+,50M,60,0", result)

    def test_drops_out_of_range_position(self):
        """Drops SA entries whose position is outside any toy contig."""
        region_map = {"chr1": [("chr1:1000-2000", 1000, 2000)]}
        result = generate_toy_dataset._remap_sa_tag(
            "SA:Z:chr1,9999,+,50M,60,0;", region_map
        )
        # Position 9999 is outside chr1:1000-2000 — should be dropped
        self.assertEqual(result, "SA:Z:*")

    def test_clamp_position_to_1(self):
        """Positions at the very start of the region clamp to 1."""
        region_map = {"chr1": [("chr1:1000-2000", 1000, 2000)]}
        result = generate_toy_dataset._remap_sa_tag(
            "SA:Z:chr1,1000,+,50M,60,0;", region_map
        )
        self.assertIn("chr1:1000-2000,1,+,50M,60,0", result)

    def test_not_sa_tag(self):
        """Returns non-SA fields unchanged."""
        result = generate_toy_dataset._remap_sa_tag(
            "XY:Z:something", {"chr1": [("chr1:1000-2000", 1000, 2000)]}
        )
        self.assertEqual(result, "XY:Z:something")

    def test_multiple_regions_same_chrom(self):
        """Correctly picks the right toy contig when multiple exist."""
        region_map = {
            "chr1": [
                ("chr1:1000-2000", 1000, 2000),
                ("chr1:5000-6000", 5000, 6000),
            ]
        }
        result = generate_toy_dataset._remap_sa_tag(
            "SA:Z:chr1,5500,+,50M,60,0;", region_map
        )
        self.assertIn("chr1:5000-6000,501,+,50M,60,0", result)


class TestMinimumPadding(unittest.TestCase):
    """Verify compute_regions enforces a minimum padding."""

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()

        # Simple alignment block spanning chr1 [490000, 520000)
        paf_line = (
            "ctg1\t5000000\t100000\t130000\t+\tchr1\t248956422"
            "\t490000\t520000\t30000\t30000\t60\n"
        )
        paf_path = os.path.join(self.tmpdir, "min_pad.paf")
        with open(paf_path, "w") as fh:
            fh.write(paf_line)

        blocks = coordinate_mapper.parse_paf(paf_path)
        json_path = os.path.join(self.tmpdir, "min_pad.mapping.json.gz")
        coordinate_mapper.build_json_index(blocks, json_path)
        self.index = coordinate_mapper.load_index(json_path)

    def tearDown(self):
        shutil.rmtree(self.tmpdir)

    def test_zero_padding_raised_to_minimum(self):
        """padding=0 should be silently raised to _MIN_PADDING (100)."""
        variants = [
            {"chrom": "chr1", "pos": 500000, "end": 510000, "size": 10000,
             "genotype": "1|0"},
        ]
        regions = generate_toy_dataset.compute_regions(
            variants, self.index, self.index, padding=0
        )
        r = regions[0]
        # With _MIN_PADDING=100, region should be chr1:499900-510100
        self.assertEqual(r["ref_region"], "chr1:499900-510100")

    def test_small_padding_raised_to_minimum(self):
        """padding=50 should be silently raised to _MIN_PADDING (100)."""
        variants = [
            {"chrom": "chr1", "pos": 500000, "end": 510000, "size": 10000,
             "genotype": "1|0"},
        ]
        regions = generate_toy_dataset.compute_regions(
            variants, self.index, self.index, padding=50
        )
        r = regions[0]
        self.assertEqual(r["ref_region"], "chr1:499900-510100")

    def test_large_padding_not_reduced(self):
        """padding=5000 should be used as-is (above minimum)."""
        variants = [
            {"chrom": "chr1", "pos": 500000, "end": 510000, "size": 10000,
             "genotype": "1|0"},
        ]
        regions = generate_toy_dataset.compute_regions(
            variants, self.index, self.index, padding=5000
        )
        r = regions[0]
        self.assertEqual(r["ref_region"], "chr1:495000-515000")


if __name__ == "__main__":
    unittest.main()
