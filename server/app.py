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
import mimetypes
import os
import re
import sys
import urllib.parse
import urllib.request
from functools import lru_cache
from http import HTTPStatus
from http.server import HTTPServer, SimpleHTTPRequestHandler
from pathlib import Path

# Allow importing coordinate_mapper from the sibling src/ directory.
_REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_REPO_ROOT / "src"))
import coordinate_mapper  # noqa: E402


# ── Configuration loading ───────────────────────────────────────────────────

def load_config(tsv_path):
    """Load sample configuration from a TSV file.

    Expected columns (tab-separated, first row is header):
        sample_id       Unique sample identifier
        output_dir      Path to pipeline output directory
        reference       Path to reference FASTA (.fa.gz)
        hap1_assembly   Path to haplotype 1 assembly FASTA (.fa.gz)
        hap2_assembly   Path to haplotype 2 assembly FASTA (.fa.gz)
        reads_bam       (optional) Path to reads-vs-reference BAM
        regions         (optional) Path to regions file (manifest JSON or VCF)

    Lines beginning with ``#`` are ignored.  Additional columns are stored
    as extra metadata.
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
            samples.append(row)
    return samples


def discover_pipeline_files(sample):
    """Discover pipeline output files for a sample.

    The pipeline writes files with a ``{prefix}_`` naming convention
    inside the output directory.  This function scans the output_dir
    and populates the sample dict with resolved file paths.
    """
    output_dir = sample.get("output_dir", "")
    if not output_dir or not os.path.isdir(output_dir):
        return sample

    # Auto-detect the sample prefix by looking for *_hap1_to_ref.bam
    prefix = None
    for fname in os.listdir(output_dir):
        if fname.endswith("_hap1_to_ref.bam"):
            prefix = fname.replace("_hap1_to_ref.bam", "")
            break

    if prefix is None:
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
                                "reads_bam"):
                        v = sample.get(key, "")
                        if v:
                            candidates.append(os.path.dirname(v))
                    for d in candidates:
                        check = os.path.join(d, filename)
                        if os.path.isfile(check):
                            fs_path = check
                            self.file_registry[url_path] = fs_path
                            break

        if fs_path is None or not os.path.isfile(fs_path):
            self.send_error(HTTPStatus.NOT_FOUND)
            return

        file_size = os.path.getsize(fs_path)
        content_type = self._guess_type(fs_path)

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
        """Quieter logging — skip successful data requests."""
        path = args[0] if args else ""
        if isinstance(path, str) and "/data/" in path and "200" in str(args):
            return
        super().log_message(format, *args)


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
    args = parser.parse_args()

    # Load configuration
    config_path = os.path.abspath(args.config)
    print(f"Loading configuration from {config_path}")
    samples = load_config(config_path)
    if not samples:
        print("ERROR: No samples found in configuration file", file=sys.stderr)
        sys.exit(1)

    # Discover pipeline files and load mapping indices
    translator = CoordinateTranslator()
    for sample in samples:
        sample = discover_pipeline_files(sample)
        translator.load_sample(sample)
        print(f"  Loaded sample: {sample['sample_id']}")

    # Pre-register data files for all samples
    file_registry = {}
    for sample in samples:
        sid = sample["sample_id"]
        for key in ("reference", "hap1_assembly", "hap2_assembly",
                     "reads_bam", "hap1_to_ref_bam", "hap2_to_ref_bam",
                     "reads_to_hap1_bam", "reads_to_hap2_bam",
                     "ref_to_hap1_bam", "ref_to_hap2_bam",
                     "reads_cram", "reads_to_hap1_cram",
                     "reads_to_hap2_cram"):
            path = sample.get(key)
            if path and os.path.isfile(path):
                url = f"/data/{sid}/{os.path.basename(path)}"
                file_registry[url] = path
                # Register associated index files
                for ext in (".bai", ".fai", ".gzi", ".tbi", ".crai"):
                    idx = path + ext
                    if os.path.isfile(idx):
                        file_registry[url + ext] = idx

    # Configure the handler
    static_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                              "static")
    _ensure_igv_js(static_dir)
    IGVHandler.samples = samples
    IGVHandler.translator = translator
    IGVHandler.file_registry = file_registry
    IGVHandler.static_dir = static_dir

    # Start server
    server = HTTPServer((args.host, args.port), IGVHandler)
    url = f"http://{'localhost' if args.host == '0.0.0.0' else args.host}:{args.port}"
    print(f"\nServer running at {url}")
    print(f"Serving {len(samples)} sample(s)")
    print("Press Ctrl+C to stop\n")

    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("\nShutting down...")
        server.shutdown()


if __name__ == "__main__":
    main()
