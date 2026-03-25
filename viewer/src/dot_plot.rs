//! K-mer dot plot computation for sequence comparison.
//!
//! Generates dot plot data comparing two DNA sequences by finding shared
//! k-mer matches.  The algorithm follows the approach described in wotplot
//! (Fedarko 2023, <https://github.com/fedarko/wotplot>) — every exact k-mer
//! match between the two input sequences is identified and classified as a
//! forward match, reverse-complement match, or palindromic match.
//!
//! This is a Rust port of `src/dot_plot.py`.

use std::collections::HashMap;

// ---------------------------------------------------------------------------
// Types
// ---------------------------------------------------------------------------

/// Match type codes following wotplot conventions.
#[allow(dead_code)] // Part of the dot plot API; used for classification.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MatchType {
    /// K-mers are identical.
    Forward,
    /// K-mer in seq1 equals the reverse-complement in seq2.
    ReverseComplement,
    /// K-mer equals its own reverse complement.
    Palindrome,
}

/// A single match position: (x, y) in the dot plot coordinate space.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DotMatch {
    /// Position in seq1 (x-axis).
    pub x: usize,
    /// Position in seq2 (y-axis).
    pub y: usize,
}

/// Result of a dot plot computation.
#[derive(Debug, Clone)]
pub struct DotPlotResult {
    /// Forward matches (k-mer identical).
    pub forward: Vec<DotMatch>,
    /// Reverse-complement matches.
    pub reverse: Vec<DotMatch>,
    /// Palindromic matches.
    pub palindrome: Vec<DotMatch>,
    /// Length of seq1.
    pub seq1_len: usize,
    /// Length of seq2.
    pub seq2_len: usize,
    /// K-mer size used.
    #[allow(dead_code)] // Metadata field; will be displayed in dot plot UI.
    pub k: usize,
    /// Whether the result was truncated due to hitting `MAX_MATCHES`.
    pub truncated: bool,
}

impl DotPlotResult {
    fn empty(seq1_len: usize, seq2_len: usize, k: usize) -> Self {
        Self {
            forward: Vec::new(),
            reverse: Vec::new(),
            palindrome: Vec::new(),
            seq1_len,
            seq2_len,
            k,
            truncated: false,
        }
    }

    /// Total number of matches.
    pub fn total_matches(&self) -> usize {
        self.forward.len() + self.reverse.len() + self.palindrome.len()
    }
}

/// Identifies which pair of panels a dot plot compares.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DotPlotComparison {
    /// Reference vs Haplotype 1.
    RefVsHap1,
    /// Reference vs Haplotype 2.
    RefVsHap2,
    /// Haplotype 1 vs Haplotype 2.
    Hap1VsHap2,
}

impl DotPlotComparison {
    /// All possible comparisons.
    pub const ALL: [DotPlotComparison; 3] = [
        DotPlotComparison::RefVsHap1,
        DotPlotComparison::RefVsHap2,
        DotPlotComparison::Hap1VsHap2,
    ];

    /// Human-readable label.
    pub fn label(self) -> &'static str {
        match self {
            DotPlotComparison::RefVsHap1 => "Ref vs Hap1",
            DotPlotComparison::RefVsHap2 => "Ref vs Hap2",
            DotPlotComparison::Hap1VsHap2 => "Hap1 vs Hap2",
        }
    }
}

// ---------------------------------------------------------------------------
// DNA helpers
// ---------------------------------------------------------------------------

/// Safety limit: maximum number of match positions to return.
const MAX_MATCHES: usize = 2_000_000;

/// Default k-mer size.
pub const DEFAULT_K: usize = 31;

/// Return the complement of a single DNA base (uppercase).
#[inline]
fn complement(base: u8) -> u8 {
    match base {
        b'A' | b'a' => b'T',
        b'T' | b't' => b'A',
        b'C' | b'c' => b'G',
        b'G' | b'g' => b'C',
        _ => base,
    }
}

/// Return the reverse complement of a DNA sequence.
pub fn reverse_complement(seq: &str) -> String {
    seq.as_bytes()
        .iter()
        .rev()
        .map(|&b| complement(b) as char)
        .collect()
}

// ---------------------------------------------------------------------------
// Dot plot computation
// ---------------------------------------------------------------------------

/// Compute a k-mer dot plot between two sequences.
///
/// For each position *i* in `seq1` and *j* in `seq2` where the k-mers
/// match, a match record is produced.
///
/// # Panics
///
/// Panics if `k` is 0.
pub fn compute_dotplot(seq1: &str, seq2: &str, k: usize) -> DotPlotResult {
    assert!(k >= 1, "k must be >= 1");

    let n1 = seq1.len();
    let n2 = seq2.len();

    if seq1.is_empty() || seq2.is_empty() || n1 < k || n2 < k {
        return DotPlotResult::empty(n1, n2, k);
    }

    let s1 = seq1.as_bytes();
    let s2 = seq2.as_bytes();

    // Build k-mer index for seq2.
    let mut kmer_positions: HashMap<&[u8], Vec<usize>> = HashMap::new();
    for j in 0..=(n2 - k) {
        let kmer = &s2[j..j + k];
        if kmer.contains(&b'N') {
            continue;
        }
        kmer_positions.entry(kmer).or_default().push(j);
    }

    let mut forward = Vec::new();
    let mut reverse = Vec::new();
    let mut palindrome = Vec::new();
    let mut total: usize = 0;
    let mut truncated = false;

    for i in 0..=(n1 - k) {
        let kmer = &s1[i..i + k];
        if kmer.contains(&b'N') {
            continue;
        }

        let rc = reverse_complement(std::str::from_utf8(kmer).unwrap_or(""));
        let rc_bytes = rc.as_bytes();
        let is_palindrome = kmer == rc_bytes;

        // Forward matches
        if let Some(positions) = kmer_positions.get(kmer) {
            for &j in positions {
                if is_palindrome {
                    palindrome.push(DotMatch { x: i, y: j });
                } else {
                    forward.push(DotMatch { x: i, y: j });
                }
                total += 1;
                if total >= MAX_MATCHES {
                    truncated = true;
                    break;
                }
            }
            if truncated {
                break;
            }
        }

        // Reverse-complement matches (skip if palindrome to avoid double-counting)
        if !is_palindrome
            && !truncated
            && let Some(positions) = kmer_positions.get(rc_bytes)
        {
            for &j in positions {
                reverse.push(DotMatch { x: i, y: j });
                total += 1;
                if total >= MAX_MATCHES {
                    truncated = true;
                    break;
                }
            }
            if truncated {
                break;
            }
        }
    }

    DotPlotResult {
        forward,
        reverse,
        palindrome,
        seq1_len: n1,
        seq2_len: n2,
        k,
        truncated,
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // -- reverse_complement -------------------------------------------------

    #[test]
    fn test_reverse_complement_basic() {
        assert_eq!(reverse_complement("ACGT"), "ACGT"); // palindrome
        assert_eq!(reverse_complement("AAAA"), "TTTT");
        assert_eq!(reverse_complement("CCCC"), "GGGG");
        assert_eq!(reverse_complement("ATCG"), "CGAT");
    }

    #[test]
    fn test_reverse_complement_empty() {
        assert_eq!(reverse_complement(""), "");
    }

    #[test]
    fn test_reverse_complement_single_base() {
        assert_eq!(reverse_complement("A"), "T");
        assert_eq!(reverse_complement("T"), "A");
        assert_eq!(reverse_complement("C"), "G");
        assert_eq!(reverse_complement("G"), "C");
    }

    #[test]
    fn test_reverse_complement_lowercase_input() {
        assert_eq!(reverse_complement("acgt"), "ACGT");
    }

    // -- compute_dotplot: basic cases ---------------------------------------

    #[test]
    fn test_dotplot_empty_seq1() {
        let result = compute_dotplot("", "ACGTACGT", 3);
        assert_eq!(result.total_matches(), 0);
        assert!(!result.truncated);
    }

    #[test]
    fn test_dotplot_empty_seq2() {
        let result = compute_dotplot("ACGTACGT", "", 3);
        assert_eq!(result.total_matches(), 0);
    }

    #[test]
    fn test_dotplot_seq_shorter_than_k() {
        let result = compute_dotplot("AC", "ACGT", 3);
        assert_eq!(result.total_matches(), 0);
        assert_eq!(result.seq1_len, 2);
    }

    #[test]
    #[should_panic(expected = "k must be >= 1")]
    fn test_dotplot_k_zero_panics() {
        compute_dotplot("ACGT", "ACGT", 0);
    }

    // -- compute_dotplot: identity ------------------------------------------

    #[test]
    fn test_dotplot_identical_short() {
        // k=3 on identical 6bp sequences → 4 forward matches on the diagonal
        let seq = "ATCGAT";
        let result = compute_dotplot(seq, seq, 3);
        assert_eq!(result.forward.len(), 4);
        // All matches should lie on the diagonal
        for m in &result.forward {
            assert_eq!(m.x, m.y, "match should be on diagonal");
        }
    }

    #[test]
    fn test_dotplot_k1_identical() {
        let result = compute_dotplot("ACGT", "ACGT", 1);
        // k=1: each position matches at least itself
        assert!(result.total_matches() >= 4);
    }

    // -- compute_dotplot: forward matches -----------------------------------

    #[test]
    fn test_dotplot_forward_match() {
        let seq1 = "AAAATTTTT"; // 9bp
        let seq2 = "XXXXXAAAAXXX"; // 12bp
        // k=4: "AAAT" appears at positions 2..6 in seq1
        //       and we look for "AAAA" in seq2 at position 5..9
        let result = compute_dotplot(seq1, seq2, 4);
        // At least one forward match should be found for the AAAA kmer
        assert!(
            !result.forward.is_empty() || !result.palindrome.is_empty(),
            "expected at least one match"
        );
    }

    // -- compute_dotplot: reverse complement --------------------------------

    #[test]
    fn test_dotplot_reverse_complement_match() {
        // "AACCGG" → RC = "CCGGTT"
        let seq1 = "AACCGG";
        let seq2 = "CCGGTT";
        let result = compute_dotplot(seq1, seq2, 3);
        // There should be at least one reverse-complement match
        assert!(
            !result.reverse.is_empty(),
            "expected reverse-complement matches, got {:?}",
            result
        );
    }

    // -- compute_dotplot: N-containing k-mers skipped -----------------------

    #[test]
    fn test_dotplot_n_skipped() {
        let seq1 = "ANGT";
        let seq2 = "ANGT";
        let result = compute_dotplot(seq1, seq2, 3);
        // k-mers containing N should be skipped
        assert_eq!(result.total_matches(), 0);
    }

    // -- compute_dotplot: palindromic k-mers --------------------------------

    #[test]
    fn test_dotplot_palindrome() {
        // "ACGT" is a DNA palindrome (its RC is also "ACGT")
        let seq1 = "ACGT";
        let seq2 = "ACGT";
        let result = compute_dotplot(seq1, seq2, 4);
        // The single k-mer should be classified as palindrome
        assert_eq!(result.palindrome.len(), 1);
        assert_eq!(result.forward.len(), 0);
        assert_eq!(result.reverse.len(), 0);
    }

    // -- compute_dotplot: result metadata -----------------------------------

    #[test]
    fn test_dotplot_metadata() {
        let result = compute_dotplot("AACCGGTT", "TTGGCCAA", 4);
        assert_eq!(result.seq1_len, 8);
        assert_eq!(result.seq2_len, 8);
        assert_eq!(result.k, 4);
    }

    // -- DotPlotComparison --------------------------------------------------

    #[test]
    fn test_comparison_labels() {
        assert_eq!(DotPlotComparison::RefVsHap1.label(), "Ref vs Hap1");
        assert_eq!(DotPlotComparison::RefVsHap2.label(), "Ref vs Hap2");
        assert_eq!(DotPlotComparison::Hap1VsHap2.label(), "Hap1 vs Hap2");
    }

    #[test]
    fn test_comparison_all() {
        assert_eq!(DotPlotComparison::ALL.len(), 3);
    }

    // -- Dot plot region calculation: trigger / compute routine --------------

    /// Simulate the dot plot trigger/compute routine that would be called
    /// when the user requests a dot plot for the current region.
    #[test]
    fn test_dotplot_trigger_compute_routine() {
        // Simulate sequences loaded from the three panels.
        // Use non-repetitive sequences so the identity case produces a clean diagonal.
        let ref_seq = "ATCGATCGTTGACCAATCGAATCGTTGACCA";
        let hap1_seq = "ATCGATCGTTGACCAATCGAATCGTTGACCA"; // identical to ref
        let hap2_seq = "TGCATGCAAACTGGTTGCATGCAAACTGGTT"; // different

        let k = 5;

        // Compute all three comparisons
        let ref_vs_hap1 = compute_dotplot(ref_seq, hap1_seq, k);
        let ref_vs_hap2 = compute_dotplot(ref_seq, hap2_seq, k);
        let hap1_vs_hap2 = compute_dotplot(hap1_seq, hap2_seq, k);

        // ref vs hap1: identical sequences → strong diagonal.
        // At minimum, some forward/palindrome matches must exist.
        assert!(
            ref_vs_hap1.total_matches() > 0,
            "identical sequences must produce matches"
        );

        // Verify diagonal matches exist (some x == y)
        let on_diagonal = ref_vs_hap1
            .forward
            .iter()
            .chain(ref_vs_hap1.palindrome.iter())
            .filter(|m| m.x == m.y)
            .count();
        assert!(
            on_diagonal > 0,
            "identical sequences should have diagonal matches"
        );

        // ref vs hap2: different sequences → different match pattern
        // (the exact counts will differ from ref_vs_hap1)
        let _ = ref_vs_hap2;

        // hap1 vs hap2: should also compute without panic
        let _ = hap1_vs_hap2;

        // Verify metadata is correct
        assert_eq!(ref_vs_hap1.seq1_len, ref_seq.len());
        assert_eq!(ref_vs_hap1.seq2_len, hap1_seq.len());
        assert_eq!(ref_vs_hap1.k, k);
    }

    // -- Stress / large sequence test ---------------------------------------

    #[test]
    fn test_dotplot_large_sequence() {
        // 1000bp random-ish sequences
        let seq1: String = "ACGT".repeat(250);
        let seq2: String = "ACGT".repeat(250);
        let result = compute_dotplot(&seq1, &seq2, 11);
        assert!(result.total_matches() > 0);
        assert!(!result.truncated);
    }

    #[test]
    fn test_dotplot_no_match() {
        // Sequences with no shared k-mers
        let seq1 = "AAAAAAAAAA";
        let seq2 = "CCCCCCCCCC";
        let result = compute_dotplot(seq1, seq2, 5);
        assert_eq!(result.forward.len(), 0);
        // But reverse complement of AAAAA = TTTTT, not CCCCC, so no RC match either
        assert_eq!(result.reverse.len(), 0);
    }

    #[test]
    fn test_dotplot_repeat_rich() {
        // Repeat-rich: many matches expected
        let seq1 = "AAAA".repeat(10); // 40 A's
        let seq2 = "AAAA".repeat(10);
        let result = compute_dotplot(&seq1, &seq2, 4);
        // "AAAA" matches everywhere against itself, and RC("AAAA")="TTTT"
        // so no RC matches with seq2
        assert!(!result.forward.is_empty() || !result.palindrome.is_empty());
    }
}
