#!/usr/bin/env python3
"""Multi-panel IGV.js visualization server.

A lightweight Python HTTP server that serves a multi-panel IGV.js interface
for reviewing long-read structural variant data across a reference genome
and both haplotype assemblies.

Features
--------
- Three coordinated IGV.js browser panels (Reference, Hap1, Hap2)
- Server-side coordinate translation between reference and assembly spaces
- Sample selection and region-of-interest browsing (from manifest JSON or VCF)
- Byte-range HTTP support for BAM/FASTA random access
- Zero external Python dependencies (uses only the standard library)

Usage
-----
    python3 server/app.py --config samples.tsv [--port 8080] [--host 0.0.0.0]
"""

import argparse
import gzip
import json
import logging
import mimetypes
import os
import re
import subprocess
import sys
import urllib.parse
import urllib.request
from functools import lru_cache
from http import HTTPStatus
from http.server import ThreadingHTTPServer, SimpleHTTPRequestHandler
from pathlib import Path

# Allow importing coordinate_mapper from the sibling src/ directory.
_REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_REPO_ROOT / "src"))
import coordinate_mapper  # noqa: E402

# ── Logging setup ───────────────────────────────────────────────────────────

logger = logging.getLogger("lrv")


def _setup_logging(level=logging.INFO):
    """Configure logging for the visualization server."""
    handler = logging.StreamHandler(sys.stderr)
    handler.setFormatter(logging.Formatter(
        "%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    ))
    logger.addHandler(handler)
    logger.setLevel(level)


# ── Configuration loading ───────────────────────────────────────────────────

def load_config(tsv_path):
    """Load sample configuration from a TSV file.

    Expected columns (tab-separated, first row is header):
        sample_id       Unique sample identifier
        output_dir      Path to pipeline output directory
        reference       Path to reference FASTA (.fa.gz)
        hap1_assembly   Path to haplotype 1 assembly FASTA (.fa.gz)
        hap2_assembly   Path to haplotype 2 assembly FASTA (.fa.gz)
        reads_bam       (optional) Path to reads-vs-reference BAM or CRAM
        regions         (optional) Path to regions file (manifest JSON or VCF)
        cram_ref        (optional) Path to the reference FASTA used for CRAM
                        encoding.  When set, this reference is served to
                        igv.js for CRAM decoding instead of the panel
                        reference.

    Lines beginning with ``#`` are ignored.  Additional columns are stored
    as extra metadata.  If ``reads_bam`` points to a ``.cram`` file it is
    automatically reclassified as ``reads_cram`` (an explicit ``reads_cram``
    column takes priority).
    """
    samples = []
    with open(tsv_path) as fh:
        header = None
        for line in fh:
            line = line.rstrip("\n\r")
            if not line or line.startswith("#"):
                if line.startswith("#") and header is None:
                    # Allow header line starting with '#'
                    header = line.lstrip("#").strip().split("\t")
                continue
            if header is None:
                header = line.split("\t")
                continue
            fields = line.split("\t")
            row = {}
            for i, col in enumerate(header):
                row[col.strip()] = fields[i].strip() if i < len(fields) else ""
            # Reclassify reads_bam → reads_cram when the path is a CRAM file
            bam_val = row.get("reads_bam", "")
            if bam_val.endswith(".cram"):
                row.setdefault("reads_cram", bam_val)
                row["reads_bam"] = ""
            samples.append(row)

    # Log loaded samples
    for s in samples:
        cram_keys = [k for k in s if "cram" in k and s[k]]
        cram_ref = s.get("cram_ref", "")
        logger.info("Config sample=%s cram_files=%s cram_ref=%s",
                     s.get("sample_id"), cram_keys, cram_ref or "(not set)")
    return samples


def discover_pipeline_files(sample):
    """Discover pipeline output files for a sample.

    The pipeline writes files with a ``{prefix}_`` naming convention
    inside the output directory.  This function scans the output_dir
    and populates the sample dict with resolved file paths.
    """
    output_dir = sample.get("output_dir", "")
    if not output_dir or not os.path.isdir(output_dir):
        logger.debug("discover_pipeline_files: skipping (no output_dir)")
        return sample

    # Auto-detect the sample prefix by looking for *_hap1_to_ref.bam
    prefix = None
    for fname in os.listdir(output_dir):
        if fname.endswith("_hap1_to_ref.bam"):
            prefix = fname.replace("_hap1_to_ref.bam", "")
            break

    if prefix is None:
        logger.debug("discover_pipeline_files: no prefix found in %s",
                      output_dir)
        return sample

    sample["_prefix"] = prefix
    sample["_output_dir"] = output_dir

    # Map of logical names → filename suffixes
    file_map = {
        "hap1_to_ref_bam": f"{prefix}_hap1_to_ref.bam",
        "hap2_to_ref_bam": f"{prefix}_hap2_to_ref.bam",
        "reads_to_hap1_bam": f"{prefix}_reads_to_hap1.bam",
        "reads_to_hap2_bam": f"{prefix}_reads_to_hap2.bam",
        "ref_to_hap1_bam": f"{prefix}_ref_to_hap1.bam",
        "ref_to_hap2_bam": f"{prefix}_ref_to_hap2.bam",
        "hap1_mapping_index": f"{prefix}_hap1_to_ref.mapping.json.gz",
        "hap2_mapping_index": f"{prefix}_hap2_to_ref.mapping.json.gz",
    }
    # Also check for CRAM versions of the read alignment files
    cram_map = {
        "reads_to_hap1_cram": f"{prefix}_reads_to_hap1.cram",
        "reads_to_hap2_cram": f"{prefix}_reads_to_hap2.cram",
        "reads_cram": f"{prefix}_reads.cram",
    }
    for key, fname in file_map.items():
        path = os.path.join(output_dir, fname)
        if os.path.isfile(path):
            sample[key] = path
    for key, fname in cram_map.items():
        path = os.path.join(output_dir, fname)
        if os.path.isfile(path):
            sample[key] = path
            logger.info("Discovered CRAM file %s=%s", key, path)

    return sample


# ── Region loading (manifest JSON or VCF) ──────────────────────────────────

def load_regions(regions_path):
    """Load regions of interest from a manifest JSON or VCF file.

    Returns a list of dicts with at minimum:
        label, chrom, start, end, ref_region
    And optionally:
        hap1_regions, hap2_regions, genotype, size
    """
    if not regions_path or not os.path.isfile(regions_path):
        return []

    if regions_path.endswith(".json"):
        return _load_manifest_regions(regions_path)
    elif regions_path.endswith(".vcf") or regions_path.endswith(".vcf.gz"):
        return _load_vcf_regions(regions_path)
    return []


def _load_manifest_regions(path):
    """Load regions from a toy_manifest.json-style file."""
    with open(path) as fh:
        data = json.load(fh)
    regions = []
    for v in data.get("variants", []):
        chrom = v["chrom"]
        ref_region = v.get("ref_region", "")
        # Parse start/end from ref_region  "chr1:9279383-9389426"
        m = re.match(r"(.+):(\d+)-(\d+)", ref_region)
        start = int(m.group(2)) if m else v.get("pos", 0)
        end = int(m.group(3)) if m else start + v.get("size", 1000)
        size = v.get("size", end - start)
        gt = v.get("genotype", "?")
        label = f"{chrom}:{start}-{end} ({_format_size(size)}, {gt})"
        regions.append({
            "label": label,
            "chrom": chrom,
            "start": start,
            "end": end,
            "ref_region": ref_region,
            "hap1_regions": v.get("hap1_regions", []),
            "hap2_regions": v.get("hap2_regions", []),
            "genotype": gt,
            "size": size,
        })
    return regions


def _load_vcf_regions(path, padding=50000, min_size=500):
    """Load SV deletion regions from a VCF file."""
    open_fn = gzip.open if path.endswith(".gz") else open
    regions = []
    with open_fn(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            cols = line.strip().split("\t")
            if len(cols) < 8:
                continue
            chrom = cols[0]
            pos = int(cols[1])
            ref_allele = cols[3]
            alt_allele = cols[4]
            info = cols[7]

            # Determine SV size
            sv_len = 0
            # Check INFO field for SVLEN
            for field in info.split(";"):
                if field.startswith("SVLEN="):
                    sv_len = abs(int(field.split("=")[1]))
                    break
            if sv_len == 0:
                sv_len = abs(len(ref_allele) - len(alt_allele))
            if sv_len < min_size:
                continue

            end = pos + sv_len
            padded_start = max(0, pos - padding)
            padded_end = end + padding

            # Parse genotype if available
            gt = ""
            if len(cols) >= 10:
                gt_field = cols[9].split(":")[0]
                gt = gt_field.replace("/", "|")

            label = f"{chrom}:{pos}-{end} ({_format_size(sv_len)}"
            if gt:
                label += f", {gt}"
            label += ")"

            regions.append({
                "label": label,
                "chrom": chrom,
                "start": padded_start,
                "end": padded_end,
                "ref_region": f"{chrom}:{padded_start}-{padded_end}",
                "hap1_regions": [],
                "hap2_regions": [],
                "genotype": gt,
                "size": sv_len,
            })
    return regions


def _format_size(bp):
    """Format a base-pair count as a human-readable string."""
    if bp >= 1_000_000:
        return f"{bp / 1_000_000:.1f} Mb"
    if bp >= 1000:
        return f"{bp / 1000:.1f} kb"
    return f"{bp} bp"


# ── Coordinate translation ──────────────────────────────────────────────────

class CoordinateTranslator:
    """Manages coordinate mapping indices for all samples."""

    def __init__(self):
        self._indices = {}  # (sample_id, haplotype) → loaded index

    def load_sample(self, sample):
        """Load mapping indices for a sample."""
        sid = sample["sample_id"]
        for hap in ("hap1", "hap2"):
            key = f"{hap}_mapping_index"
            path = sample.get(key)
            if path and os.path.isfile(path):
                idx = coordinate_mapper.load_index(path)
                self._indices[(sid, hap)] = idx

    def translate(self, sample_id, chrom, start, end, min_mapq=0):
        """Translate a reference region to assembly coordinates.

        Returns a dict with hap1 and hap2 assembly region lists.
        """
        result = {"hap1": [], "hap2": []}
        for hap in ("hap1", "hap2"):
            idx = self._indices.get((sample_id, hap))
            if idx is None:
                continue
            hits = coordinate_mapper.query(idx, chrom, start, end,
                                           min_mapq=min_mapq)
            # Collect assembly regions from alignment hits
            asm_regions = []
            for hit in hits:
                if hit.get("event_type") != "alignment":
                    continue
                asm_chrom = hit["asm_chrom"]
                asm_start = hit["asm_start"]
                asm_end = hit["asm_end"]
                asm_regions.append({
                    "chrom": asm_chrom,
                    "start": asm_start,
                    "end": asm_end,
                    "strand": hit.get("strand", "+"),
                })
            # Merge overlapping regions per contig
            result[hap] = _merge_regions(asm_regions)
        return result


def _merge_regions(regions):
    """Merge overlapping assembly regions, grouped by contig."""
    by_contig = {}
    for r in regions:
        by_contig.setdefault(r["chrom"], []).append(r)
    merged = []
    for contig, intervals in sorted(by_contig.items()):
        intervals.sort(key=lambda x: x["start"])
        cur = intervals[0].copy()
        for nxt in intervals[1:]:
            if nxt["start"] <= cur["end"]:
                cur["end"] = max(cur["end"], nxt["end"])
            else:
                merged.append(cur)
                cur = nxt.copy()
        merged.append(cur)
    return merged


# ── HTTP Request Handler ────────────────────────────────────────────────────

class IGVHandler(SimpleHTTPRequestHandler):
    """Custom HTTP handler with byte-range support and API endpoints."""

    # Class-level state (set before server starts)
    samples = []
    translator = None
    regions_cache = {}
    file_registry = {}    # URL path → filesystem path
    static_dir = None

    def do_GET(self):
        parsed = urllib.parse.urlparse(self.path)
        path = parsed.path
        query = urllib.parse.parse_qs(parsed.query)

        # API endpoints
        if path == "/api/samples":
            return self._json_response(self._api_samples())
        if path.startswith("/api/sample/"):
            sid = path[len("/api/sample/"):]
            return self._json_response(self._api_sample(sid))
        if path == "/api/translate":
            return self._json_response(self._api_translate(query))
        if path.startswith("/api/regions/"):
            sid = path[len("/api/regions/"):]
            return self._json_response(self._api_regions(sid))

        # Data files (BAM, FASTA, etc.) — served with byte-range support
        if path.startswith("/data/"):
            return self._serve_data_file(path)

        # Static files and index
        if path == "/" or path == "/index.html":
            return self._serve_static("index.html")
        if path.startswith("/static/"):
            fname = path[len("/static/"):]
            return self._serve_static(fname)

        self.send_error(HTTPStatus.NOT_FOUND)

    def do_HEAD(self):
        """Handle HEAD requests (needed by igv.js for content-length checks)."""
        parsed = urllib.parse.urlparse(self.path)
        path = parsed.path

        if path.startswith("/data/"):
            return self._serve_data_file(path, head_only=True)

        self.send_error(HTTPStatus.NOT_FOUND)

    def do_OPTIONS(self):
        """Handle CORS preflight."""
        self.send_response(HTTPStatus.NO_CONTENT)
        self._cors_headers()
        self.end_headers()

    # ── API handlers ────────────────────────────────────────────────────────

    def _api_samples(self):
        """Return list of sample IDs and basic info."""
        return [
            {
                "sample_id": s["sample_id"],
                "has_regions": bool(s.get("regions")),
            }
            for s in self.samples
        ]

    def _api_sample(self, sample_id):
        """Return full configuration for a sample, including file URLs."""
        sample = self._find_sample(sample_id)
        if sample is None:
            return {"error": f"Sample '{sample_id}' not found"}

        def data_url(key):
            path = sample.get(key)
            if not path or not os.path.isfile(path):
                return None
            url_path = f"/data/{sample_id}/{os.path.basename(path)}"
            self.file_registry[url_path] = path
            # Also register index files
            for ext in (".bai", ".fai", ".gzi", ".tbi", ".crai"):
                idx = path + ext
                if os.path.isfile(idx):
                    idx_url = url_path + ext
                    self.file_registry[idx_url] = idx
            return url_path

        # Build cram_ref URL if the sample provides one
        cram_ref_url = data_url("cram_ref")
        if cram_ref_url:
            logger.info("Sample %s: using explicit cram_ref=%s",
                         sample_id, sample.get("cram_ref"))
        else:
            cram_keys = [k for k in ("reads_cram", "reads_to_hap1_cram",
                                      "reads_to_hap2_cram")
                         if sample.get(k)]
            if cram_keys:
                logger.info(
                    "Sample %s: CRAM files present (%s) but no cram_ref; "
                    "using panel reference for CRAM decoding",
                    sample_id, ", ".join(cram_keys))

        config = {
            "sample_id": sample_id,
            "reference": {
                "fastaURL": data_url("reference"),
            },
            "hap1_assembly": {
                "fastaURL": data_url("hap1_assembly"),
            },
            "hap2_assembly": {
                "fastaURL": data_url("hap2_assembly"),
            },
            "tracks": {
                "hap1_to_ref_bam": data_url("hap1_to_ref_bam"),
                "hap2_to_ref_bam": data_url("hap2_to_ref_bam"),
                "reads_to_hap1_bam": data_url("reads_to_hap1_bam"),
                "reads_to_hap2_bam": data_url("reads_to_hap2_bam"),
                "ref_to_hap1_bam": data_url("ref_to_hap1_bam"),
                "ref_to_hap2_bam": data_url("ref_to_hap2_bam"),
                "reads_bam": data_url("reads_bam"),
                "reads_cram": data_url("reads_cram"),
                "reads_to_hap1_cram": data_url("reads_to_hap1_cram"),
                "reads_to_hap2_cram": data_url("reads_to_hap2_cram"),
            },
            "cram_ref": cram_ref_url,
        }
        return config

    def _api_translate(self, query):
        """Translate reference coordinates to assembly coordinates."""
        sample_id = query.get("sample", [""])[0]
        chrom = query.get("chrom", [""])[0]
        try:
            start = int(query.get("start", [0])[0])
            end = int(query.get("end", [0])[0])
        except (ValueError, IndexError):
            return {"error": "Invalid start/end coordinates"}
        min_mapq = int(query.get("min_mapq", [0])[0])

        if not sample_id or not chrom:
            return {"error": "Missing sample or chrom parameter"}

        result = self.translator.translate(sample_id, chrom, start, end,
                                           min_mapq=min_mapq)
        return result

    def _api_regions(self, sample_id):
        """Return regions of interest for a sample."""
        if sample_id in self.regions_cache:
            return self.regions_cache[sample_id]

        sample = self._find_sample(sample_id)
        if sample is None:
            return []

        regions = load_regions(sample.get("regions", ""))
        self.regions_cache[sample_id] = regions
        return regions

    # ── File serving ────────────────────────────────────────────────────────

    def _serve_data_file(self, url_path, head_only=False):
        """Serve a data file with HTTP byte-range support."""
        # First check the registry
        fs_path = self.file_registry.get(url_path)

        if fs_path is None:
            # Try to find it by pattern: /data/{sample_id}/{filename}
            parts = url_path.split("/")
            if len(parts) >= 4:
                sample_id = parts[2]
                filename = "/".join(parts[3:])
                sample = self._find_sample(sample_id)
                if sample:
                    # Search in output_dir and assembly/reference directories
                    candidates = [sample.get("_output_dir", "")]
                    for key in ("reference", "hap1_assembly", "hap2_assembly",
                                "reads_bam", "cram_ref"):
                        v = sample.get(key, "")
                        if v:
                            candidates.append(os.path.dirname(v))
                    for d in candidates:
                        check = os.path.join(d, filename)
                        if os.path.isfile(check):
                            fs_path = check
                            self.file_registry[url_path] = fs_path
                            logger.debug("Resolved %s → %s", url_path,
                                          fs_path)
                            break

        if fs_path is None or not os.path.isfile(fs_path):
            logger.warning("File not found: %s (resolved=%s)", url_path,
                            fs_path)
            self.send_error(HTTPStatus.NOT_FOUND)
            return

        file_size = os.path.getsize(fs_path)
        content_type = self._guess_type(fs_path)

        # Log CRAM-related requests at INFO level for diagnostics
        if url_path.endswith((".cram", ".crai")):
            range_hdr = self.headers.get("Range", "none")
            logger.info("CRAM request: %s %s size=%d range=%s",
                         "HEAD" if head_only else "GET", url_path,
                         file_size, range_hdr)

        # Parse Range header
        range_header = self.headers.get("Range")
        if range_header:
            return self._serve_range(fs_path, file_size, content_type,
                                     range_header, head_only)

        # Full file response
        self.send_response(HTTPStatus.OK)
        self._cors_headers()
        self.send_header("Content-Type", content_type)
        self.send_header("Content-Length", str(file_size))
        self.send_header("Accept-Ranges", "bytes")
        self.end_headers()

        if not head_only:
            with open(fs_path, "rb") as f:
                self._copy_file(f, self.wfile, file_size)

    def _serve_range(self, fs_path, file_size, content_type,
                     range_header, head_only=False):
        """Handle an HTTP Range request."""
        m = re.match(r"bytes=(\d*)-(\d*)", range_header)
        if not m:
            self.send_error(HTTPStatus.REQUESTED_RANGE_NOT_SATISFIABLE)
            return

        range_start = m.group(1)
        range_end = m.group(2)

        if range_start:
            start = int(range_start)
            end = int(range_end) if range_end else file_size - 1
        elif range_end:
            # Suffix range: last N bytes
            start = file_size - int(range_end)
            end = file_size - 1
        else:
            self.send_error(HTTPStatus.REQUESTED_RANGE_NOT_SATISFIABLE)
            return

        if start > end or start >= file_size:
            self.send_error(HTTPStatus.REQUESTED_RANGE_NOT_SATISFIABLE)
            return

        end = min(end, file_size - 1)
        length = end - start + 1

        self.send_response(HTTPStatus.PARTIAL_CONTENT)
        self._cors_headers()
        self.send_header("Content-Type", content_type)
        self.send_header("Content-Length", str(length))
        self.send_header("Content-Range",
                         f"bytes {start}-{end}/{file_size}")
        self.send_header("Accept-Ranges", "bytes")
        self.end_headers()

        if not head_only:
            with open(fs_path, "rb") as f:
                f.seek(start)
                self._copy_file(f, self.wfile, length)

    def _serve_static(self, filename):
        """Serve a static file from the static directory."""
        safe = os.path.normpath(filename)
        if safe.startswith(".."):
            self.send_error(HTTPStatus.FORBIDDEN)
            return

        filepath = os.path.join(self.static_dir, safe)
        if not os.path.isfile(filepath):
            self.send_error(HTTPStatus.NOT_FOUND)
            return

        content_type = self._guess_type(filepath)
        with open(filepath, "rb") as f:
            data = f.read()

        self.send_response(HTTPStatus.OK)
        self._cors_headers()
        self.send_header("Content-Type", content_type)
        self.send_header("Content-Length", str(len(data)))
        self.send_header("Cache-Control", "no-cache")
        self.end_headers()
        self.wfile.write(data)

    # ── Helpers ─────────────────────────────────────────────────────────────

    def _find_sample(self, sample_id):
        for s in self.samples:
            if s["sample_id"] == sample_id:
                return s
        return None

    def _json_response(self, data):
        body = json.dumps(data, indent=None, separators=(",", ":"))
        body_bytes = body.encode("utf-8")
        self.send_response(HTTPStatus.OK)
        self._cors_headers()
        self.send_header("Content-Type", "application/json")
        self.send_header("Content-Length", str(len(body_bytes)))
        self.end_headers()
        self.wfile.write(body_bytes)

    def _cors_headers(self):
        self.send_header("Access-Control-Allow-Origin", "*")
        self.send_header("Access-Control-Allow-Headers",
                         "Range, Content-Type")
        self.send_header("Access-Control-Allow-Methods",
                         "GET, HEAD, OPTIONS")
        self.send_header("Access-Control-Expose-Headers",
                         "Content-Range, Content-Length, Accept-Ranges")

    @staticmethod
    def _guess_type(path):
        ext = os.path.splitext(path)[1].lower()
        type_map = {
            ".html": "text/html; charset=utf-8",
            ".js": "application/javascript",
            ".css": "text/css",
            ".json": "application/json",
            ".bam": "application/octet-stream",
            ".bai": "application/octet-stream",
            ".cram": "application/octet-stream",
            ".crai": "application/octet-stream",
            ".fa": "text/plain",
            ".fasta": "text/plain",
            ".gz": "application/octet-stream",
            ".gzi": "application/octet-stream",
            ".tbi": "application/octet-stream",
            ".fai": "text/plain",
            ".bed": "text/plain",
            ".paf": "text/plain",
            ".vcf": "text/plain",
        }
        return type_map.get(ext, "application/octet-stream")

    @staticmethod
    def _copy_file(src, dst, length, chunk_size=65536):
        remaining = length
        while remaining > 0:
            chunk = src.read(min(chunk_size, remaining))
            if not chunk:
                break
            dst.write(chunk)
            remaining -= len(chunk)

    def log_message(self, format, *args):
        """Override default logging.

        Successful /data/ responses are logged at DEBUG level to reduce
        noise while still being available when --verbose is used.
        Everything else is logged at INFO via the standard handler.
        """
        msg = format % args if args else format
        if "/data/" in msg and (" 200 " in msg or " 206 " in msg):
            logger.debug("HTTP %s", msg)
            return
        logger.info("HTTP %s", msg)


# ── IGV.js dependency ───────────────────────────────────────────────────────

# Note: keep in sync with the ARG IGV_JS_VERSION in Dockerfile
IGV_JS_VERSION = "3.1.3"
IGV_JS_URL = (
    f"https://cdn.jsdelivr.net/npm/igv@{IGV_JS_VERSION}/dist/igv.min.js"
)


def _ensure_igv_js(static_dir):
    """Download igv.min.js into static/ if it is not already present."""
    dest = os.path.join(static_dir, "igv.min.js")
    if os.path.isfile(dest):
        return
    print(f"Downloading igv.js v{IGV_JS_VERSION} ...")
    try:
        urllib.request.urlretrieve(IGV_JS_URL, dest)
        print(f"  Saved to {dest}")
    except Exception as exc:
        print(f"  WARNING: Could not download igv.js ({exc})")
        print("  The server will still start, but the UI requires igv.min.js")
        print(f"  in {static_dir}/.")
        print(f"  Download it manually from: {IGV_JS_URL}")


def _validate_cram_ref(path):
    """Check that a CRAM reference FASTA is readable and indexed.

    Returns a list of warning strings (empty when everything is OK).
    """
    warnings = []
    if not path:
        return warnings
    if not os.path.isfile(path):
        warnings.append(f"cram_ref file not found: {path}")
        return warnings

    fai = path + ".fai"
    if not os.path.isfile(fai):
        warnings.append(f"cram_ref missing .fai index: {fai}")

    if path.endswith(".gz"):
        gzi = path + ".gzi"
        if not os.path.isfile(gzi):
            warnings.append(f"cram_ref is bgzipped but missing .gzi index: "
                            f"{gzi}")

    # Quick sanity check: try reading the first few bytes
    try:
        with open(path, "rb") as fh:
            magic = fh.read(2)
        if path.endswith(".gz") and magic != b"\x1f\x8b":
            warnings.append(
                f"cram_ref has .gz extension but does not look bgzipped "
                f"(magic bytes: {magic!r})")
    except OSError as exc:
        warnings.append(f"cram_ref cannot be read: {exc}")

    return warnings


def _check_cram_embedded_ref(cram_path):
    """Read a CRAM file's header and return the set of UR paths found.

    Uses ``samtools view -H`` when *samtools* is on PATH; otherwise
    falls back to a lightweight binary scan of the CRAM container that
    extracts the SAM header text without any external dependency.

    Returns ``(ur_paths, warnings)`` where *ur_paths* is the set of
    unique ``UR:`` values and *warnings* is a list of human-readable
    warning strings for any UR path that is not reachable locally.
    """
    header_text = _read_cram_header(cram_path)
    if header_text is None:
        return set(), []

    ur_paths = set()
    for line in header_text.splitlines():
        if not line.startswith("@SQ"):
            continue
        for field in line.split("\t"):
            if field.startswith("UR:"):
                ur_paths.add(field[3:])

    warnings = []
    for ur in sorted(ur_paths):
        if not os.path.isfile(ur):
            warnings.append(
                f"CRAM header embeds UR:{ur} which does not exist locally. "
                f"Decoding will use the panel reference (or explicit "
                f"cram_ref) instead.")
    return ur_paths, warnings


def _read_cram_header(cram_path):
    """Extract the SAM header text from a CRAM file.

    Tries ``samtools view -H`` first.  When samtools is unavailable,
    does a minimal binary parse of the CRAM v3 container to pull out the
    header without any external dependency.

    Returns the header string, or *None* on failure.
    """
    # ── Try samtools first ──────────────────────────────────────────────
    try:
        proc = subprocess.run(
            ["samtools", "view", "-H", cram_path],
            capture_output=True, text=True, timeout=30,
        )
        if proc.returncode == 0 and proc.stdout:
            return proc.stdout
    except (FileNotFoundError, subprocess.TimeoutExpired):
        pass

    # ── Fallback: lightweight binary parse ──────────────────────────────
    # CRAM v2/v3 layout (first container):
    #   magic  (4 bytes)  "CRAM"
    #   major  (1 byte)
    #   minor  (1 byte)
    #   file_id (20 bytes)
    #   --- first container (the SAM header) ---
    #   container_length   (ITF-8)
    #   ref_seq_id         (ITF-8)
    #   start_pos          (ITF-8)
    #   alignment_span     (ITF-8)
    #   num_records        (ITF-8)
    #   record_counter     (ITF-8 in v3, absent in v2)
    #   num_bases          (ITF-8)
    #   num_blocks         (ITF-8)
    #   landmarks_len      (ITF-8)  + landmark array
    #   crc32              (4 bytes, v3 only)
    #   --- first block (header block) ---
    #   block_method       (1 byte)    0 = raw
    #   block_content_type (1 byte)    0 = FILE_HEADER
    #   block_content_id   (ITF-8)
    #   compressed_size    (ITF-8)
    #   raw_size           (ITF-8)
    #   data               (raw_size bytes — ITF-8 length + SAM header text)
    try:
        with open(cram_path, "rb") as fh:
            magic = fh.read(4)
            if magic != b"CRAM":
                logger.debug("Not a CRAM file (magic=%r): %s",
                              magic, cram_path)
                return None
            major = fh.read(1)[0]
            _minor = fh.read(1)[0]
            _file_id = fh.read(20)

            # Skip the container header — we need to read ITF-8 values
            _cont_len = _read_itf8(fh)
            _ref_id = _read_itf8(fh)
            _start = _read_itf8(fh)
            _span = _read_itf8(fh)
            _n_rec = _read_itf8(fh)
            if major >= 3:
                _rec_ctr = _read_itf8(fh)
            _n_bases = _read_itf8(fh)
            n_blocks = _read_itf8(fh)
            landmarks_len = _read_itf8(fh)
            for _ in range(landmarks_len):
                _read_itf8(fh)
            if major >= 3:
                fh.read(4)  # CRC32

            # First block: the header block
            _blk_method = fh.read(1)[0]
            _blk_ctype = fh.read(1)[0]
            _blk_cid = _read_itf8(fh)
            _comp_sz = _read_itf8(fh)
            raw_sz = _read_itf8(fh)

            # The block data starts with an ITF-8 header length, then
            # the SAM header text.
            hdr_len = _read_itf8(fh)
            hdr_bytes = fh.read(hdr_len)
            return hdr_bytes.decode("utf-8", errors="replace")
    except Exception as exc:
        logger.debug("Could not parse CRAM header from %s: %s",
                      cram_path, exc)
        return None


def _read_itf8(fh):
    """Read an ITF-8 encoded integer from a binary file handle."""
    b0 = fh.read(1)
    if not b0:
        raise EOFError("Unexpected end of CRAM file")
    b0 = b0[0]
    if b0 < 0x80:
        return b0
    if b0 < 0xC0:
        b1 = fh.read(1)[0]
        return ((b0 & 0x3F) << 8) | b1
    if b0 < 0xE0:
        rest = fh.read(2)
        return ((b0 & 0x1F) << 16) | (rest[0] << 8) | rest[1]
    if b0 < 0xF0:
        rest = fh.read(3)
        return (((b0 & 0x0F) << 24) | (rest[0] << 16)
                | (rest[1] << 8) | rest[2])
    rest = fh.read(4)
    return (((b0 & 0x0F) << 28) | (rest[0] << 20)
            | (rest[1] << 12) | (rest[2] << 4) | (rest[3] & 0x0F))


# ── Main ────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Multi-panel IGV.js visualization server")
    parser.add_argument("--config", "-c", required=True,
                        help="Path to samples TSV configuration file")
    parser.add_argument("--port", "-p", type=int, default=8080,
                        help="Server port (default: 8080)")
    parser.add_argument("--host", default="0.0.0.0",
                        help="Server host (default: 0.0.0.0)")
    parser.add_argument("--verbose", "-v", action="store_true",
                        help="Enable verbose (DEBUG) logging")
    args = parser.parse_args()

    # Initialise logging
    _setup_logging(logging.DEBUG if args.verbose else logging.INFO)

    # Load configuration
    config_path = os.path.abspath(args.config)
    logger.info("Loading configuration from %s", config_path)
    samples = load_config(config_path)
    if not samples:
        logger.error("No samples found in configuration file")
        sys.exit(1)

    # Discover pipeline files and load mapping indices
    translator = CoordinateTranslator()
    for sample in samples:
        sample = discover_pipeline_files(sample)
        translator.load_sample(sample)
        logger.info("Loaded sample: %s", sample["sample_id"])

    # Validate CRAM references
    for sample in samples:
        sid = sample["sample_id"]
        cram_ref = sample.get("cram_ref", "")
        cram_keys = [k for k in ("reads_cram", "reads_to_hap1_cram",
                                  "reads_to_hap2_cram")
                     if sample.get(k)]
        if cram_keys and not cram_ref:
            logger.warning(
                "Sample %s has CRAM file(s) (%s) but no cram_ref column. "
                "CRAM decoding will use each panel's reference genome. "
                "If the browser crashes, set cram_ref to the FASTA that "
                "was used when encoding the CRAMs.",
                sid, ", ".join(cram_keys))
        if cram_ref:
            for w in _validate_cram_ref(cram_ref):
                logger.warning("Sample %s: %s", sid, w)

        # Check each CRAM file's embedded UR reference path
        for cram_key in cram_keys:
            cram_path = sample.get(cram_key, "")
            if not cram_path or not os.path.isfile(cram_path):
                continue
            ur_paths, ur_warnings = _check_cram_embedded_ref(cram_path)
            for w in ur_warnings:
                logger.warning("Sample %s (%s): %s", sid, cram_key, w)
            if ur_paths:
                logger.info(
                    "Sample %s (%s): embedded UR path(s): %s",
                    sid, cram_key, ", ".join(sorted(ur_paths)))

    # Pre-register data files for all samples
    file_registry = {}
    for sample in samples:
        sid = sample["sample_id"]
        for key in ("reference", "hap1_assembly", "hap2_assembly",
                     "reads_bam", "hap1_to_ref_bam", "hap2_to_ref_bam",
                     "reads_to_hap1_bam", "reads_to_hap2_bam",
                     "ref_to_hap1_bam", "ref_to_hap2_bam",
                     "reads_cram", "reads_to_hap1_cram",
                     "reads_to_hap2_cram", "cram_ref"):
            path = sample.get(key)
            if path and os.path.isfile(path):
                url = f"/data/{sid}/{os.path.basename(path)}"
                file_registry[url] = path
                # Register associated index files
                for ext in (".bai", ".fai", ".gzi", ".tbi", ".crai"):
                    idx = path + ext
                    if os.path.isfile(idx):
                        file_registry[url + ext] = idx

    logger.info("File registry: %d entries", len(file_registry))
    for url, fspath in sorted(file_registry.items()):
        logger.debug("  %s → %s", url, fspath)

    # Configure the handler
    static_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                              "static")
    _ensure_igv_js(static_dir)
    IGVHandler.samples = samples
    IGVHandler.translator = translator
    IGVHandler.file_registry = file_registry
    IGVHandler.static_dir = static_dir

    # Start server
    server = ThreadingHTTPServer((args.host, args.port), IGVHandler)
    url = f"http://{'localhost' if args.host == '0.0.0.0' else args.host}:{args.port}"
    logger.info("Server running at %s", url)
    logger.info("Serving %d sample(s)", len(samples))
    logger.info("Press Ctrl+C to stop")

    try:
        server.serve_forever()
    except KeyboardInterrupt:
        logger.info("Shutting down...")
        server.shutdown()


if __name__ == "__main__":
    main()
