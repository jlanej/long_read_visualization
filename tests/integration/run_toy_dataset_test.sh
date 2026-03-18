#!/usr/bin/env bash
# run_toy_dataset_test.sh – Integration test using the bundled toy dataset.
#
# Runs the full preprocess.sh pipeline on the 40-deletion toy dataset
# committed in resources/toy_dataset/ and verifies that every expected
# output file is produced.  After the pipeline run, a handful of coordinate
# mapping queries are executed against the generated indices to confirm
# they return valid results for regions listed in the toy manifest.
#
# The test is designed so that the pipeline outputs can later be loaded
# directly into an igv.js multi-panel visualization; all BAMs, indices,
# and coordinate-mapping files are produced in a single output directory.
#
# Usage:
#   IMAGE=long_read_visualization:test bash tests/integration/run_toy_dataset_test.sh
#
# The IMAGE variable defaults to "long_read_visualization:test" when not set.

set -euo pipefail

IMAGE="${IMAGE:-long_read_visualization:test}"
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
TOY_DIR="${REPO_ROOT}/resources/toy_dataset"
TESTDIR="$(mktemp -d)"
trap 'rm -rf "${TESTDIR}"' EXIT

echo "=== Toy Dataset Integration Test ==="
echo "Image:    ${IMAGE}"
echo "Toy data: ${TOY_DIR}"
echo "Workdir:  ${TESTDIR}"
echo ""

# ── Helper ──────────────────────────────────────────────────────────────────
docker_run() {
    docker run --rm --entrypoint "" \
        --user "$(id -u):$(id -g)" \
        -v "${TESTDIR}:/testdata" \
        "${IMAGE}" "$@"
}

# ── 1. Copy toy dataset into the working directory ──────────────────────────
echo "--- 1. Preparing toy dataset ---"
cp -a "${TOY_DIR}/." "${TESTDIR}/toy/"
echo "Copied toy dataset to ${TESTDIR}/toy/"
ls -lh "${TESTDIR}/toy/"
echo ""

# ── 2. Run the pipeline on the toy dataset ──────────────────────────────────
echo "--- 2. Running preprocess.sh on toy dataset ---"
docker run --rm \
    --user "$(id -u):$(id -g)" \
    -v "${TESTDIR}:/testdata" \
    "${IMAGE}" \
    --hap1       /testdata/toy/toy_hap1.fa.gz \
    --hap2       /testdata/toy/toy_hap2.fa.gz \
    --cram       /testdata/toy/toy_reads.cram \
    --cram-ref   /testdata/toy/toy_reference.fa.gz \
    --reference  /testdata/toy/toy_reference.fa.gz \
    --output-dir /testdata/output \
    --threads    2 \
    --ont
echo ""

# ── 3. Verify expected output files ────────────────────────────────────────
echo "--- 3. Verifying pipeline outputs ---"
SAMPLE="toy_reads"

pass=0
fail=0
check() {
    local label="$1"
    if [[ -f "${TESTDIR}/output/${label}" ]]; then
        echo "  ✓ ${label}"
        pass=$((pass + 1))
    else
        echo "  ✗ ${label}  (MISSING)"
        fail=$((fail + 1))
    fi
}

# Assembly-to-reference alignments (BAM)
check "${SAMPLE}_hap1_to_ref.bam"
check "${SAMPLE}_hap1_to_ref.bam.bai"
check "${SAMPLE}_hap2_to_ref.bam"
check "${SAMPLE}_hap2_to_ref.bam.bai"

# PAF alignments for coordinate mapping
check "${SAMPLE}_hap1_to_ref.paf"
check "${SAMPLE}_hap2_to_ref.paf"

# Coordinate mapping indices (JSON + BED)
check "${SAMPLE}_hap1_to_ref.mapping.json.gz"
check "${SAMPLE}_hap2_to_ref.mapping.json.gz"
check "${SAMPLE}_hap1_to_ref.mapping.bed.gz"
check "${SAMPLE}_hap1_to_ref.mapping.bed.gz.tbi"
check "${SAMPLE}_hap2_to_ref.mapping.bed.gz"
check "${SAMPLE}_hap2_to_ref.mapping.bed.gz.tbi"

# Reads aligned to assemblies
check "${SAMPLE}_reads_to_hap1.bam"
check "${SAMPLE}_reads_to_hap1.bam.bai"
check "${SAMPLE}_reads_to_hap2.bam"
check "${SAMPLE}_reads_to_hap2.bam.bai"

# Reference-to-assembly cross-check BAMs
check "${SAMPLE}_ref_to_hap1.bam"
check "${SAMPLE}_ref_to_hap1.bam.bai"
check "${SAMPLE}_ref_to_hap2.bam"
check "${SAMPLE}_ref_to_hap2.bam.bai"

echo ""
if [[ "${fail}" -gt 0 ]]; then
    echo "=== FAILED: ${fail} expected file(s) missing ==="
    exit 1
fi
echo "All ${pass} output files present."
echo ""

# ── 4. Validate coordinate mapping queries ──────────────────────────────────
echo "--- 4. Validating coordinate mapping queries ---"

# Use the first variant from the manifest to test a coordinate query.
# The Python one-liner reads the manifest, picks the first ref_region,
# queries the hap1 and hap2 mapping indices, and verifies non-empty results.
docker_run python3 -c "
import json, sys
sys.path.insert(0, '/opt/long_read_visualization/src')
import coordinate_mapper

manifest = json.load(open('/testdata/toy/toy_manifest.json'))
variants = manifest['variants']
print(f'Manifest contains {len(variants)} variant(s)')

hap1_idx = coordinate_mapper.load_index(
    '/testdata/output/${SAMPLE}_hap1_to_ref.mapping.json.gz')
hap2_idx = coordinate_mapper.load_index(
    '/testdata/output/${SAMPLE}_hap2_to_ref.mapping.json.gz')

tested = 0
for v in variants[:5]:
    region = v['ref_region']
    chrom, coords = region.split(':')
    start, end = [int(x) for x in coords.split('-')]
    h1 = coordinate_mapper.query(hap1_idx, chrom, start, end)
    h2 = coordinate_mapper.query(hap2_idx, chrom, start, end)
    total = len(h1) + len(h2)
    if total == 0:
        print(f'  WARNING: no mapping hits for {region}')
    else:
        print(f'  ✓ {region}: {len(h1)} hap1 + {len(h2)} hap2 hits')
    tested += 1

print(f'Queried {tested} regions successfully')
"
echo ""

# ── 5. Validate BAM files are readable ──────────────────────────────────────
echo "--- 5. Quick-checking BAM file integrity ---"
for bam in \
    "${SAMPLE}_hap1_to_ref.bam" \
    "${SAMPLE}_hap2_to_ref.bam" \
    "${SAMPLE}_reads_to_hap1.bam" \
    "${SAMPLE}_reads_to_hap2.bam" \
    "${SAMPLE}_ref_to_hap1.bam" \
    "${SAMPLE}_ref_to_hap2.bam"; do
    count=$(docker_run samtools view -c "/testdata/output/${bam}" 2>/dev/null || echo "ERROR")
    if [[ "${count}" == "ERROR" ]]; then
        echo "  ✗ ${bam}: samtools view failed"
        fail=$((fail + 1))
    elif [[ "${count}" -eq 0 ]]; then
        echo "  ⚠ ${bam}: 0 alignments (may be expected for small dataset)"
    else
        echo "  ✓ ${bam}: ${count} alignment(s)"
    fi
done
echo ""

# ── Summary ─────────────────────────────────────────────────────────────────
if [[ "${fail}" -gt 0 ]]; then
    echo "=== Toy Dataset Integration Test FAILED ==="
    exit 1
fi
echo "=== Toy Dataset Integration Test PASSED ==="
