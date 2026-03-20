#!/usr/bin/env bash
# on_create.sh — One-time setup for the GitHub Codespaces environment.
#
# This script runs automatically after the Codespace container is first created
# (not on subsequent restarts).  It:
#   1. Runs the preprocessing pipeline on the bundled toy dataset.
#   2. Generates the server configuration TSV.
#
# The preprocessing output is stored in toy_example_output/ inside the
# workspace and persists between Codespace restarts.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
TOY_DIR="${REPO_ROOT}/resources/toy_dataset"
TOY_OUTPUT="${REPO_ROOT}/toy_example_output"
CONFIG="${REPO_ROOT}/toy_example_config.tsv"

echo "=== Long Read Visualization — one-time setup ==="
echo ""
echo "Running preprocessing pipeline on the bundled toy dataset."
echo "This takes approximately 1–3 minutes."
echo ""

bash "${REPO_ROOT}/scripts/preprocess.sh" \
    --hap1       "${TOY_DIR}/toy_hap1.fa.gz" \
    --hap2       "${TOY_DIR}/toy_hap2.fa.gz" \
    --cram       "${TOY_DIR}/toy_reads.bam" \
    --reference  "${TOY_DIR}/toy_reference.fa.gz" \
    --output-dir "${TOY_OUTPUT}" \
    --threads    4 \
    --ont

echo ""
echo "=== Generating server configuration ==="

python3 "${REPO_ROOT}/scripts/generate_server_config.py" \
    --output-dir "${TOY_OUTPUT}" \
    --reference  "${TOY_DIR}/toy_reference.fa.gz" \
    --hap1       "${TOY_DIR}/toy_hap1.fa.gz" \
    --hap2       "${TOY_DIR}/toy_hap2.fa.gz" \
    --reads-bam  "${TOY_DIR}/toy_reads.bam" \
    --regions    "${TOY_DIR}/toy_manifest.json" \
    --output     "${CONFIG}"

echo ""
echo "Setup complete.  The visualization server will start automatically."
echo "Open the Ports panel in VS Code and click on port 8080 to view."
