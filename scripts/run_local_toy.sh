#!/usr/bin/env bash
# run_local_toy.sh — Build and launch the native Rust GUI with the toy dataset.
#
# Usage:
#   bash scripts/run_local_toy.sh
#
# Prerequisites:
#   - A working Rust toolchain (rustup / cargo).
#   - On Linux: libxcb-render0-dev, libxcb-shape0-dev, libxcb-xfixes0-dev,
#     libxkbcommon-dev, libgtk-3-dev  (for egui/eframe).
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

TOY="${REPO_ROOT}/resources/toy_dataset"
MANIFEST="${TOY}/toy_manifest.json"

if [[ ! -f "${MANIFEST}" ]]; then
    echo "ERROR: Toy manifest not found at ${MANIFEST}" >&2
    exit 1
fi

echo "=== Building long_read_viewer (release) ==="
cargo build --manifest-path "${REPO_ROOT}/viewer/Cargo.toml" --release

BINARY="${REPO_ROOT}/viewer/target/release/long_read_viewer"
if [[ ! -x "${BINARY}" ]]; then
    echo "ERROR: Binary not found at ${BINARY}" >&2
    exit 1
fi

echo ""
echo "=== Launching viewer with toy dataset ==="
echo "  Manifest:  ${MANIFEST}"
echo "  Reference: ${TOY}/toy_reference.fa.gz"
echo "  Hap1:      ${TOY}/toy_hap1.fa.gz"
echo "  Hap2:      ${TOY}/toy_hap2.fa.gz"
echo "  Reads:     ${TOY}/toy_reads.bam"
echo ""
exec "${BINARY}" "${MANIFEST}" \
    --reference "${TOY}/toy_reference.fa.gz" \
    --hap1      "${TOY}/toy_hap1.fa.gz" \
    --hap2      "${TOY}/toy_hap2.fa.gz" \
    --bam       "${TOY}/toy_reads.bam"
