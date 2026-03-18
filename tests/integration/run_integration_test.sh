#!/usr/bin/env bash
# run_integration_test.sh – Integration tests for the long_read_visualization
# Docker image.
#
# Tests:
#   1. All packaged tools are available and correctly linked (samtools,
#      minimap2, bgzip, tabix, python3).
#   2. A minimal end-to-end pipeline run with synthetic data succeeds and
#      produces the expected output files.
#
# Usage:
#   IMAGE=long_read_visualization:test bash tests/integration/run_integration_test.sh
#
# The IMAGE variable defaults to "long_read_visualization:test" when not set.

set -euo pipefail

IMAGE="${IMAGE:-long_read_visualization:test}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TESTDIR="$(mktemp -d)"
trap 'rm -rf "${TESTDIR}"' EXIT

echo "=== Integration Tests ==="
echo "Image:   ${IMAGE}"
echo "Testdir: ${TESTDIR}"
echo ""

# Helper: run a command inside the image with /testdata mounted.
# Run as the current user so that output files are owned by the host user
# and can be cleaned up by the EXIT trap without permission errors.
docker_run() {
    docker run --rm --entrypoint "" \
        --user "$(id -u):$(id -g)" \
        -v "${TESTDIR}:/testdata" \
        "${IMAGE}" "$@"
}

# ── 1. Tool availability and shared-library checks ───────────────────────────
echo "--- 1. Tool version checks ---"
docker_run samtools --version | head -1
docker_run minimap2 --version
docker_run bgzip    --version 2>&1 | head -1
docker_run tabix    --version 2>&1 | head -1
docker_run python3  --version
docker_run bash     --version | head -1
echo "All packaged tools are available and linked correctly."
echo ""

# ── 2. Generate minimal synthetic test data ──────────────────────────────────
echo "--- 2. Generating test data ---"
python3 "${SCRIPT_DIR}/generate_test_data.py" "${TESTDIR}"

# Compress + index FASTAs; convert the SAM to a sorted, indexed CRAM.
docker_run bash -c "
    set -euo pipefail
    for fa in /testdata/ref.fa /testdata/hap1.fa /testdata/hap2.fa; do
        bgzip -f \"\${fa}\"
        samtools faidx \"\${fa}.gz\"
    done
    samtools sort -o /testdata/reads.bam /testdata/reads.sam
    samtools index /testdata/reads.bam
    samtools view -C -T /testdata/ref.fa.gz \
        -o /testdata/reads.cram /testdata/reads.bam
    samtools index /testdata/reads.cram
"
echo "Test data ready."
echo ""

# ── 3. End-to-end pipeline run ───────────────────────────────────────────────
echo "--- 3. Running preprocess.sh ---"
docker run --rm \
    --user "$(id -u):$(id -g)" \
    -v "${TESTDIR}:/testdata" \
    "${IMAGE}" \
    --hap1       /testdata/hap1.fa.gz \
    --hap2       /testdata/hap2.fa.gz \
    --cram       /testdata/reads.cram \
    --cram-ref   /testdata/ref.fa.gz \
    --reference  /testdata/ref.fa.gz \
    --output-dir /testdata/output \
    --threads    2 \
    --ont
echo ""

# ── 4. Verify expected output files ─────────────────────────────────────────
echo "--- 4. Verifying outputs ---"
SAMPLE=$(ls "${TESTDIR}/output/" 2>/dev/null \
    | grep "_hap1_to_ref\.bam$" \
    | sed 's/_hap1_to_ref\.bam$//' \
    | head -1)

if [[ -z "${SAMPLE}" ]]; then
    echo "ERROR: No pipeline output files found in ${TESTDIR}/output/" >&2
    exit 1
fi
echo "Detected sample name: ${SAMPLE}"

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

check "${SAMPLE}_hap1_to_ref.bam"
check "${SAMPLE}_hap1_to_ref.bam.bai"
check "${SAMPLE}_hap2_to_ref.bam"
check "${SAMPLE}_hap2_to_ref.bam.bai"
check "${SAMPLE}_hap1_to_ref.paf"
check "${SAMPLE}_hap2_to_ref.paf"
check "${SAMPLE}_hap1_to_ref.mapping.json.gz"
check "${SAMPLE}_hap2_to_ref.mapping.json.gz"
check "${SAMPLE}_reads.fastq.gz"
check "${SAMPLE}_reads_to_hap1.bam"
check "${SAMPLE}_reads_to_hap2.bam"

echo ""
if [[ "${fail}" -gt 0 ]]; then
    echo "=== Integration Test FAILED: ${fail} file(s) missing ==="
    exit 1
fi
echo "=== Integration Test PASSED: ${pass}/${pass} checks OK ==="
