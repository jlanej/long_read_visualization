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

# ── 1b. Verify packaged scripts and viewer binary exist ──────────────────────
echo "--- 1b. Verifying packaged scripts and viewer binary ---"
MISSING=0
check_path() {
    local path="$1"
    if docker_run test -f "${path}"; then
        echo "  ✓ ${path}"
    else
        echo "  ✗ ${path}  (MISSING)"
        MISSING=$((MISSING + 1))
    fi
}
check_path /opt/long_read_visualization/scripts/launch_viewer.sh
check_path /opt/long_read_visualization/scripts/launch_server.sh
check_path /opt/long_read_visualization/scripts/preprocess.sh
check_path /opt/long_read_visualization/scripts/generate_server_config.py
check_path /usr/local/bin/long_read_viewer
if [[ "${MISSING}" -gt 0 ]]; then
    echo "ERROR: ${MISSING} required file(s) are missing from the container." >&2
    exit 1
fi
echo "All required scripts and binaries are present."
echo ""

# ── 1b2. Verify GUI runtime libraries needed by winit/eframe ─────────────────
echo "--- 1b2. Verifying GUI runtime libraries ---"
if docker_run test -e /usr/lib/x86_64-linux-gnu/libX11-xcb.so.1; then
    echo "  ✓ /usr/lib/x86_64-linux-gnu/libX11-xcb.so.1"
else
    echo "ERROR: missing GUI runtime library /usr/lib/x86_64-linux-gnu/libX11-xcb.so.1" >&2
    exit 1
fi
echo ""

# ── 1c. Shared-library and GLIBC compatibility checks ────────────────────────
# Verify that every compiled binary in the image can be loaded by the runtime
# dynamic linker.  This catches glibc version mismatches like:
#   "version `GLIBC_2.39' not found"
# which occur when a binary is compiled on a host with a newer glibc than the
# runtime image provides.  We check all binaries that originate from outside
# the runtime image (Rust build stage, pre-compiled upstream releases).
echo "--- 1c. Shared-library and GLIBC compatibility checks ---"
docker_run bash -c '
    set -euo pipefail
    fail=0

    # Determine the glibc version provided by the runtime image.
    # The first line of "ldd --version" ends with the version number, e.g.:
    #   ldd (Ubuntu GLIBC 2.35-0ubuntu3.8) 2.35
    # Avoid `ldd | head -1` with pipefail enabled: head exits early and can
    # trigger SIGPIPE (141) in ldd, which would abort this script under `set -e`.
    IFS= read -r ldd_first_line < <(ldd --version 2>&1)
    runtime_glibc=$(printf "%s\n" "${ldd_first_line}" \
        | sed -nE "s/.* ([0-9]+(\.[0-9]+)*)$/\1/p")
    if [[ -z "${runtime_glibc}" ]]; then
        echo "ERROR: could not determine runtime glibc version from ldd --version" >&2
        exit 1
    fi
    echo "Runtime glibc: ${runtime_glibc}"
    echo ""

    check_binary() {
        local bin="$1"
        echo "  Checking ${bin}..."

        # ldd: confirm every shared library dependency is resolved.
        # A "not found" line means the linker cannot satisfy a dependency —
        # this is exactly the symptom of the GLIBC_2.39 mismatch.
        if ldd "${bin}" 2>&1 | grep -q "not found"; then
            echo "    ✗ Unresolved shared libraries:"
            ldd "${bin}" 2>&1 | grep "not found" | sed "s/^/      /"
            fail=$((fail + 1))
            return
        fi
        echo "    ✓ All shared libraries resolved"

        # readelf -V lists all versioned symbol requirements from the
        # .gnu.version_r ELF section.  We extract every GLIBC_x.y entry
        # and take the highest — that is the minimum glibc the runtime
        # must provide for the binary to load successfully.
        required=$(readelf -V "${bin}" 2>/dev/null \
            | grep -oE "GLIBC_[0-9]+(\.[0-9]+)*" \
            | sort -V | tail -1 | sed "s/GLIBC_//" || true)
        if [[ -z "${required}" ]]; then
            echo "    ✓ No versioned GLIBC symbols required"
            return
        fi
        echo "    Required GLIBC: ${required}"

        # The binary is compatible only when runtime_glibc >= required.
        # sort -V places the higher version last; if runtime_glibc is last
        # (or equal) then runtime_glibc >= required.
        if [[ "$(printf "%s\n" "${required}" "${runtime_glibc}" \
                  | sort -V | tail -1)" == "${runtime_glibc}" ]]; then
            echo "    ✓ GLIBC version compatible"
        else
            echo "    ✗ Requires GLIBC_${required} but runtime only provides ${runtime_glibc}"
            fail=$((fail + 1))
        fi
    }

    # Rust binary compiled in the builder stage — the primary regression target.
    check_binary /usr/local/bin/long_read_viewer
    # Pre-compiled upstream binaries downloaded during the Docker build.
    check_binary /usr/local/bin/minimap2
    # Binaries compiled inside the runtime image (should always pass).
    check_binary /usr/local/bin/samtools
    check_binary /usr/local/bin/bgzip
    check_binary /usr/local/bin/tabix

    if [[ "${fail}" -gt 0 ]]; then
        echo ""
        echo "ERROR: ${fail} binary check(s) failed." >&2
        exit 1
    fi
    echo ""
    echo "All binaries passed shared-library compatibility checks."
'
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
check "${SAMPLE}_hap1_to_ref.mapping.bed.gz"
check "${SAMPLE}_hap1_to_ref.mapping.bed.gz.tbi"
check "${SAMPLE}_hap2_to_ref.mapping.bed.gz"
check "${SAMPLE}_hap2_to_ref.mapping.bed.gz.tbi"
check "${SAMPLE}_reads_to_hap1.bam"
check "${SAMPLE}_reads_to_hap1.bam.bai"
check "${SAMPLE}_reads_to_hap2.bam"
check "${SAMPLE}_reads_to_hap2.bam.bai"
check "${SAMPLE}_ref_to_hap1.bam"
check "${SAMPLE}_ref_to_hap1.bam.bai"
check "${SAMPLE}_ref_to_hap2.bam"
check "${SAMPLE}_ref_to_hap2.bam.bai"

echo ""
if [[ "${fail}" -gt 0 ]]; then
    echo "=== Integration Test FAILED: ${fail} file(s) missing ==="
    exit 1
fi
echo "=== Integration Test PASSED: ${pass}/${pass} checks OK ==="
