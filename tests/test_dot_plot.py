"""Tests for the k-mer dot plot module."""

import unittest
import sys
from pathlib import Path

# Ensure src/ is importable
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
import dot_plot


class TestReverseComplement(unittest.TestCase):
    """Tests for reverse_complement()."""

    def test_simple_palindrome(self):
        """ACGT reverse-complements to ACGT (it is a palindromic k-mer)."""
        self.assertEqual(dot_plot.reverse_complement("ACGT"), "ACGT")

    def test_poly_a(self):
        self.assertEqual(dot_plot.reverse_complement("AAAA"), "TTTT")

    def test_single_base(self):
        self.assertEqual(dot_plot.reverse_complement("A"), "T")
        self.assertEqual(dot_plot.reverse_complement("C"), "G")

    def test_mixed_case(self):
        self.assertEqual(dot_plot.reverse_complement("AcGt"), "aCgT")

    def test_empty(self):
        self.assertEqual(dot_plot.reverse_complement(""), "")

    def test_palindrome_sequence(self):
        """ATAT is a palindrome (equals its own reverse complement)."""
        seq = "ATAT"
        self.assertEqual(dot_plot.reverse_complement(seq), seq)


class TestComputeDotplot(unittest.TestCase):
    """Tests for compute_dotplot()."""

    def test_identical_sequences(self):
        """Identical sequences produce a diagonal of forward matches."""
        seq = "ACGTACGTACGTACGTACGTACGTACGTACGT"  # 32 bases, k=31 → 2 k-mers
        result = dot_plot.compute_dotplot(seq, seq, k=31)
        self.assertEqual(result["seq1_len"], 32)
        self.assertEqual(result["seq2_len"], 32)
        self.assertEqual(result["k"], 31)
        self.assertFalse(result["truncated"])
        # Should have forward matches along the diagonal
        self.assertGreater(len(result["forward"]) + len(result["palindrome"]), 0)
        # No reverse complement matches for a non-palindromic sequence
        for match in result["forward"]:
            self.assertEqual(match[0], match[1], "Diagonal matches expected")

    def test_small_k(self):
        """k=1 should find individual base matches."""
        result = dot_plot.compute_dotplot("ACGT", "ACGT", k=1)
        # Each of the 4 bases matches exactly once on the diagonal
        fwd = result["forward"]
        pal = result["palindrome"]
        total = len(fwd) + len(pal)
        # At minimum, diagonal matches: (0,0), (1,1), (2,2), (3,3)
        diagonal = [m for m in fwd + pal if m[0] == m[1]]
        self.assertEqual(len(diagonal), 4)

    def test_reverse_complement_detection(self):
        """Reverse-complement sequence produces RC matches."""
        seq1 = "AAACCCGGGTTT"  # 12 bases
        seq2 = dot_plot.reverse_complement(seq1)  # "AAACCCGGGTTT" reversed
        result = dot_plot.compute_dotplot(seq1, seq2, k=3)
        # Should have reverse complement matches
        total_rc = len(result["reverse"])
        self.assertGreater(total_rc, 0, "Should detect reverse-complement matches")

    def test_no_matches(self):
        """Completely different sequences produce no matches."""
        seq1 = "AAAAAAAAAAAA"  # 12 A's
        seq2 = "CCCCCCCCCCCC"  # 12 C's
        result = dot_plot.compute_dotplot(seq1, seq2, k=3)
        total = (len(result["forward"]) + len(result["reverse"]) +
                 len(result["palindrome"]))
        self.assertEqual(total, 0)

    def test_empty_sequence(self):
        result = dot_plot.compute_dotplot("", "ACGT", k=3)
        self.assertEqual(len(result["forward"]), 0)
        self.assertEqual(result["seq1_len"], 0)

    def test_sequence_shorter_than_k(self):
        result = dot_plot.compute_dotplot("ACG", "ACG", k=5)
        self.assertEqual(len(result["forward"]), 0)
        self.assertFalse(result["truncated"])

    def test_k_equals_1(self):
        result = dot_plot.compute_dotplot("A", "A", k=1)
        total = (len(result["forward"]) + len(result["palindrome"]))
        self.assertGreater(total, 0)

    def test_invalid_k(self):
        with self.assertRaises(ValueError):
            dot_plot.compute_dotplot("ACGT", "ACGT", k=0)

    def test_n_bases_skipped(self):
        """K-mers containing N should be skipped."""
        seq1 = "ACGNACGT"
        seq2 = "ACGTACGT"
        result = dot_plot.compute_dotplot(seq1, seq2, k=4)
        # The k-mer at position 1 in seq1 is "CGNA" → contains N → skipped
        # Only k-mers without N should match
        for match in result["forward"] + result["palindrome"]:
            i = match[0]
            kmer = seq1[i:i + 4]
            self.assertNotIn("N", kmer)

    def test_palindrome_detection(self):
        """A palindromic k-mer (equals its own reverse complement) should be classified as palindrome."""
        # "ACGT" is palindromic: reverse_complement("ACGT") == "ACGT"
        seq = "ACGT"
        result = dot_plot.compute_dotplot(seq, seq, k=4)
        # The single 4-mer "ACGT" is palindromic
        self.assertEqual(len(result["palindrome"]), 1)
        self.assertEqual(len(result["forward"]), 0)

    def test_result_structure(self):
        """Result dict has all expected keys."""
        result = dot_plot.compute_dotplot("ACGTACGT", "ACGTACGT", k=3)
        expected_keys = {"forward", "reverse", "palindrome",
                         "seq1_len", "seq2_len", "k", "truncated"}
        self.assertEqual(set(result.keys()), expected_keys)

    def test_match_coordinates_in_bounds(self):
        """All match coordinates are within valid range."""
        seq1 = "ACGTACGTACGTACGT"
        seq2 = "ACGTACGTACGT"
        k = 5
        result = dot_plot.compute_dotplot(seq1, seq2, k=k)
        max_x = len(seq1) - k
        max_y = len(seq2) - k
        for match_list in [result["forward"], result["reverse"],
                           result["palindrome"]]:
            for x, y in match_list:
                self.assertGreaterEqual(x, 0)
                self.assertGreaterEqual(y, 0)
                self.assertLessEqual(x, max_x)
                self.assertLessEqual(y, max_y)

    def test_different_length_sequences(self):
        """Dot plot works with sequences of different lengths."""
        seq1 = "ACGTACGTACGT"
        seq2 = "ACGTACGT"
        result = dot_plot.compute_dotplot(seq1, seq2, k=4)
        self.assertEqual(result["seq1_len"], 12)
        self.assertEqual(result["seq2_len"], 8)
        self.assertGreater(
            len(result["forward"]) + len(result["palindrome"]), 0)


class TestExtractSequence(unittest.TestCase):
    """Tests for extract_sequence() — requires samtools on PATH."""

    def test_missing_file(self):
        """Non-existent file returns empty string."""
        result = dot_plot.extract_sequence("/nonexistent.fa", "chr1:1-100")
        self.assertEqual(result, "")

    def test_empty_region(self):
        """Empty region string handled gracefully."""
        result = dot_plot.extract_sequence("/nonexistent.fa", "")
        self.assertEqual(result, "")


if __name__ == "__main__":
    unittest.main()
