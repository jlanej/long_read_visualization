#!/usr/bin/env bash
# launch_server.sh — Launch the multi-panel IGV.js visualization server.
#
# This script provides a convenient way to start the server either
# natively (with Python 3) or inside an Apptainer container.
#
# Usage:
#   # Native Python
#   bash scripts/launch_server.sh --config samples.tsv [--port 8080]
#
#   # Inside Apptainer container
#   bash scripts/launch_server.sh --config samples.tsv --apptainer IMAGE_URI
#
#   # Quick start with toy dataset (generates config + starts server)
#   bash scripts/launch_server.sh --toy
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
PORT=8080
HOST="0.0.0.0"
CONFIG=""
APPTAINER_URI=""
TOY_MODE=false

usage() {
    cat <<EOF
Usage: $(basename "$0") [options]

Start the multi-panel IGV.js visualization server.

Options:
  --config FILE       Path to samples TSV configuration file (required unless --toy)
  --port INT          Server port [${PORT}]
  --host STR          Server host [${HOST}]
  --apptainer URI     Run inside Apptainer container (e.g. docker://ghcr.io/...)
  --toy               Quick-start with the toy dataset example
  -h, --help          Show this help message
EOF
    exit 1
}

[[ $# -eq 0 ]] && usage

while [[ $# -gt 0 ]]; do
    case "$1" in
        --config)     CONFIG="$2";        shift 2 ;;
        --port)       PORT="$2";          shift 2 ;;
        --host)       HOST="$2";          shift 2 ;;
        --apptainer)  APPTAINER_URI="$2"; shift 2 ;;
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
if [[ -n "${APPTAINER_URI}" ]]; then
    echo "Starting server via Apptainer..."
    echo "  Image: ${APPTAINER_URI}"
    echo "  Config: ${CONFIG}"
    echo "  Port: ${PORT}"
    echo ""

    # Bind the working directory and the config file's directory
    CONFIG_DIR="$(cd "$(dirname "${CONFIG}")" && pwd)"
    BIND_DIRS="${PWD}:/work"
    if [[ "${CONFIG_DIR}" != "${PWD}" ]]; then
        BIND_DIRS="${BIND_DIRS},${CONFIG_DIR}:${CONFIG_DIR}"
    fi

    # Collect unique directories from the config for binding
    while IFS=$'\t' read -r _ output_dir ref hap1 hap2 reads regions; do
        for p in "${output_dir}" "${ref}" "${hap1}" "${hap2}" "${reads}" "${regions}"; do
            if [[ -n "${p}" && -e "${p}" ]]; then
                d="$(cd "$(dirname "${p}")" && pwd)"
                if [[ "${BIND_DIRS}" != *"${d}"* ]]; then
                    BIND_DIRS="${BIND_DIRS},${d}:${d}"
                fi
            fi
        done
    done < <(grep -v '^#' "${CONFIG}")

    apptainer exec \
        --bind "${BIND_DIRS}" \
        "${APPTAINER_URI}" \
        python3 /opt/long_read_visualization/server/app.py \
        --config "${CONFIG}" \
        --port "${PORT}" \
        --host "${HOST}"
else
    echo "Starting server natively..."
    echo "  Config: ${CONFIG}"
    echo "  URL:    http://localhost:${PORT}"
    echo ""
    exec python3 "${REPO_ROOT}/server/app.py" \
        --config "${CONFIG}" \
        --port "${PORT}" \
        --host "${HOST}"
fi
