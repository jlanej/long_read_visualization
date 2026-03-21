"""K-mer dot plot computation for sequence comparison.

Generates dot plot data comparing two DNA sequences by finding shared
k-mer matches.  The algorithm follows the approach described in wotplot
(Fedarko 2023, https://github.com/fedarko/wotplot) — every exact k-mer
match between the two input sequences is identified and classified as a
forward match, reverse-complement match, or palindromic match.

The implementation is optimised for the typical window sizes encountered
in structural-variant visualisation (up to ~200 kb) and uses only the
Python standard library so that no additional dependencies are required.

References
----------
Fedarko, M. (2023). wotplot: Visualizing exact dot plots of two
sequences in a web browser. *Journal of Open Source Software*, 8(92),
6012. https://doi.org/10.21105/joss.06012
"""

import logging
import subprocess

logger = logging.getLogger("lrv")

# ── Reverse complement ──────────────────────────────────────────────────────

_COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")


def reverse_complement(seq):
    """Return the reverse complement of a DNA sequence."""
    return seq.translate(_COMPLEMENT)[::-1]


# ── Sequence extraction ─────────────────────────────────────────────────────

def extract_sequence(fasta_path, region):
    """Extract a DNA sequence from an indexed FASTA using *samtools faidx*.

    Parameters
    ----------
    fasta_path : str
        Path to the (optionally bgzip-compressed) FASTA file.  Must have
        a ``.fai`` index (and ``.gzi`` for bgzipped files).
    region : str
        Samtools-style region string, e.g. ``chr1:1000-2000``.

    Returns
    -------
    str
        Upper-cased sequence with all non-ACGT characters removed.
    """
    try:
        proc = subprocess.run(
            ["samtools", "faidx", fasta_path, region],
            capture_output=True, text=True, timeout=60,
        )
        if proc.returncode != 0:
            logger.warning("samtools faidx failed for %s %s: %s",
                           fasta_path, region, proc.stderr.strip())
            return ""
        lines = proc.stdout.strip().split("\n")
        # Skip the FASTA header line
        seq = "".join(line.strip() for line in lines if not line.startswith(">"))
        return seq.upper()
    except FileNotFoundError:
        logger.error("samtools not found on PATH")
        return ""
    except subprocess.TimeoutExpired:
        logger.warning("samtools faidx timed out for %s %s", fasta_path, region)
        return ""


# ── Dot plot computation ────────────────────────────────────────────────────

# Match-type codes (following wotplot conventions)
MATCH_FORWARD = 1
MATCH_REVERSE_COMPLEMENT = -1
MATCH_PALINDROME = 2

# Safety limit: maximum number of match positions to return.
MAX_MATCHES = 2_000_000


def compute_dotplot(seq1, seq2, k=31):
    """Compute a k-mer dot plot between two sequences.

    For each position *i* in *seq1* and *j* in *seq2* where the k-mers
    match, a match record ``(i, j, match_type)`` is produced.

    Match types follow the wotplot convention (Fedarko 2023):

    * ``1``  — forward match (k-mers identical)
    * ``-1`` — reverse-complement match
    * ``2``  — palindromic match (k-mer equals its own reverse complement)

    Parameters
    ----------
    seq1 : str
        Sequence displayed on the **x-axis** (horizontal).
    seq2 : str
        Sequence displayed on the **y-axis** (vertical).
    k : int
        K-mer size (default 31).

    Returns
    -------
    dict
        ``{"forward": [[x, y], ...], "reverse": [[x, y], ...],
        "palindrome": [[x, y], ...], "seq1_len": int, "seq2_len": int,
        "k": int, "truncated": bool}``
    """
    if k < 1:
        raise ValueError("k must be >= 1")
    if not seq1 or not seq2:
        return _empty_result(len(seq1 or ""), len(seq2 or ""), k)

    n1 = len(seq1)
    n2 = len(seq2)

    if n1 < k or n2 < k:
        return _empty_result(n1, n2, k)

    # Build k-mer index for seq2 (hash table: kmer_string → list of positions)
    kmer_positions = {}
    for j in range(n2 - k + 1):
        kmer = seq2[j:j + k]
        if "N" in kmer:
            continue
        if kmer not in kmer_positions:
            kmer_positions[kmer] = []
        kmer_positions[kmer].append(j)

    forward = []
    reverse = []
    palindrome = []
    total = 0
    truncated = False

    for i in range(n1 - k + 1):
        kmer = seq1[i:i + k]
        if "N" in kmer:
            continue

        rc_kmer = reverse_complement(kmer)
        is_palindrome = (kmer == rc_kmer)

        # Forward matches
        if kmer in kmer_positions:
            for j in kmer_positions[kmer]:
                if is_palindrome:
                    palindrome.append([i, j])
                else:
                    forward.append([i, j])
                total += 1
                if total >= MAX_MATCHES:
                    truncated = True
                    break
            if truncated:
                break

        # Reverse-complement matches (skip if palindrome to avoid double-counting)
        if not is_palindrome and not truncated and rc_kmer in kmer_positions:
            for j in kmer_positions[rc_kmer]:
                reverse.append([i, j])
                total += 1
                if total >= MAX_MATCHES:
                    truncated = True
                    break
            if truncated:
                break

    return {
        "forward": forward,
        "reverse": reverse,
        "palindrome": palindrome,
        "seq1_len": n1,
        "seq2_len": n2,
        "k": k,
        "truncated": truncated,
    }


def _empty_result(n1, n2, k):
    return {
        "forward": [],
        "reverse": [],
        "palindrome": [],
        "seq1_len": n1,
        "seq2_len": n2,
        "k": k,
        "truncated": False,
    }
