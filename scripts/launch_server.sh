#!/usr/bin/env bash
# launch_server.sh — Launch the multi-panel IGV.js visualization server.
#
# Designed to run inside an Apptainer container, though native Python
# execution is also supported.
#
# Usage (inside Apptainer — recommended):
#   apptainer exec --bind "${PWD}:/work" IMAGE \
#       bash /opt/long_read_visualization/scripts/launch_server.sh \
#       --config /work/samples.tsv --port 8080
#
#   # Quick-start with the bundled toy dataset:
#   apptainer exec --bind "${PWD}:/work" IMAGE \
#       bash /opt/long_read_visualization/scripts/launch_server.sh --toy
#
# Usage (native, for development):
#   bash scripts/launch_server.sh --config samples.tsv --port 8080
#   bash scripts/launch_server.sh --toy
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
PORT=8080
HOST="0.0.0.0"
CONFIG=""
TOY_MODE=false

usage() {
    cat <<EOF
Usage: $(basename "$0") [options]

Start the multi-panel IGV.js visualization server.

Options:
  --config FILE       Path to samples TSV configuration file (required unless --toy)
  --port INT          Server port [${PORT}]
  --host STR          Server host [${HOST}]
  --toy               Quick-start with the bundled toy dataset example
  -h, --help          Show this help message

Apptainer example:
  apptainer exec --bind "\${PWD}:/work" docker://ghcr.io/jlanej/long_read_visualization:main \\
      bash /opt/long_read_visualization/scripts/launch_server.sh \\
      --config /work/samples.tsv --port 8080

Toy dataset example:
  apptainer exec --bind "\${PWD}:/work" docker://ghcr.io/jlanej/long_read_visualization:main \\
      bash /opt/long_read_visualization/scripts/launch_server.sh --toy
EOF
    exit 1
}

[[ $# -eq 0 ]] && usage

while [[ $# -gt 0 ]]; do
    case "$1" in
        --config)     CONFIG="$2";        shift 2 ;;
        --port)       PORT="$2";          shift 2 ;;
        --host)       HOST="$2";          shift 2 ;;
        --toy)        TOY_MODE=true;      shift   ;;
        -h|--help)    usage ;;
        *)            echo "Unknown option: $1" >&2; usage ;;
    esac
done

# ── Toy mode: generate config from toy dataset if needed ────────────────────
if [[ "${TOY_MODE}" = true ]]; then
    TOY_DIR="${REPO_ROOT}/resources/toy_dataset"
    TOY_OUTPUT="${REPO_ROOT}/toy_example_output"

    if [[ ! -d "${TOY_OUTPUT}" ]] || \
       [[ ! -f "${TOY_OUTPUT}/toy_reads_hap1_to_ref.bam" ]]; then
        echo "=== Running pipeline on toy dataset ==="
        echo "(This takes 1-3 minutes depending on your system)"
        echo ""
        bash "${SCRIPT_DIR}/preprocess.sh" \
            --hap1       "${TOY_DIR}/toy_hap1.fa.gz" \
            --hap2       "${TOY_DIR}/toy_hap2.fa.gz" \
            --cram       "${TOY_DIR}/toy_reads.bam" \
            --reference  "${TOY_DIR}/toy_reference.fa.gz" \
            --output-dir "${TOY_OUTPUT}" \
            --threads    4 \
            --ont
        echo ""
    fi

    CONFIG="${REPO_ROOT}/toy_example_config.tsv"
    python3 "${SCRIPT_DIR}/generate_server_config.py" \
        --output-dir "${TOY_OUTPUT}" \
        --reference  "${TOY_DIR}/toy_reference.fa.gz" \
        --hap1       "${TOY_DIR}/toy_hap1.fa.gz" \
        --hap2       "${TOY_DIR}/toy_hap2.fa.gz" \
        --reads-bam  "${TOY_DIR}/toy_reads.bam" \
        --regions    "${TOY_DIR}/toy_manifest.json" \
        --output     "${CONFIG}"
    echo ""
fi

[[ -z "${CONFIG}" ]] && { echo "ERROR: --config is required" >&2; usage; }
[[ -f "${CONFIG}" ]] || { echo "ERROR: Config file not found: ${CONFIG}" >&2; exit 1; }

# ── Launch server ───────────────────────────────────────────────────────────
echo "Starting server..."
echo "  Config: ${CONFIG}"
echo "  URL:    http://localhost:${PORT}"
echo ""
exec python3 "${REPO_ROOT}/server/app.py" \
    --config "${CONFIG}" \
    --port "${PORT}" \
    --host "${HOST}"
