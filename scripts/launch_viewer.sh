#!/usr/bin/env bash
# launch_viewer.sh — Launch the native Rust long-read viewer.
#
# Designed to run inside an Apptainer container on HPC, though native
# execution is also supported when the binary is on PATH.
#
# Usage (inside Apptainer — recommended for HPC):
#   apptainer exec \
#       --bind "${PWD}:/work" \
#       --env DISPLAY="${DISPLAY}" \
#       docker://ghcr.io/jlanej/long_read_visualization:main \
#       bash /opt/long_read_visualization/scripts/launch_viewer.sh \
#       --config /work/samples.tsv
#
#   # Quick-start with the bundled toy dataset:
#   apptainer exec \
#       --bind "${PWD}:/work" \
#       --env DISPLAY="${DISPLAY}" \
#       docker://ghcr.io/jlanej/long_read_visualization:main \
#       bash /opt/long_read_visualization/scripts/launch_viewer.sh --toy
#
# Usage (native, for development):
#   bash scripts/launch_viewer.sh --config samples.tsv
#   bash scripts/launch_viewer.sh --toy
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
CONFIG=""
TOY_MODE=false

usage() {
    cat <<EOF
Usage: $(basename "$0") [options]

Start the native Rust long-read viewer.

Options:
  --config FILE       Path to samples TSV configuration file (required unless --toy)
  --toy               Quick-start with the bundled toy dataset example
  -h, --help          Show this help message

Apptainer example (HPC):
  apptainer exec \\
      --bind "\${PWD}:/work" \\
      --env DISPLAY="\${DISPLAY}" \\
      docker://ghcr.io/jlanej/long_read_visualization:main \\
      bash /opt/long_read_visualization/scripts/launch_viewer.sh \\
      --config /work/samples.tsv

Toy dataset example:
  apptainer exec \\
      --bind "\${PWD}:/work" \\
      --env DISPLAY="\${DISPLAY}" \\
      docker://ghcr.io/jlanej/long_read_visualization:main \\
      bash /opt/long_read_visualization/scripts/launch_viewer.sh --toy
EOF
    exit 1
}

[[ $# -eq 0 ]] && usage

while [[ $# -gt 0 ]]; do
    case "$1" in
        --config)     CONFIG="$2";   shift 2 ;;
        --toy)        TOY_MODE=true; shift   ;;
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

# ── Locate viewer binary ─────────────────────────────────────────────────────
# Prefer the installed system binary; fall back to a local release build
# when running outside the container.
if command -v long_read_viewer &>/dev/null; then
    VIEWER_BIN="long_read_viewer"
elif [[ -x "/usr/local/bin/long_read_viewer" ]]; then
    VIEWER_BIN="/usr/local/bin/long_read_viewer"
elif [[ -x "${REPO_ROOT}/viewer/target/release/long_read_viewer" ]]; then
    VIEWER_BIN="${REPO_ROOT}/viewer/target/release/long_read_viewer"
else
    echo "ERROR: long_read_viewer binary not found." >&2
    echo "  Build it with: cargo build --manifest-path viewer/Cargo.toml --release" >&2
    exit 1
fi

# ── Launch viewer ───────────────────────────────────────────────────────────
echo "Starting viewer..."
echo "  Config:  ${CONFIG}"
echo "  Binary:  ${VIEWER_BIN}"
echo "  Status:  Initial GUI load can take up to 1 minute while indexes/data warm up"

# Enable Mesa software rendering when no GPU is available (typical in
# containers / Apptainer on HPC login nodes).  Users with a real GPU can
# override by setting LIBGL_ALWAYS_SOFTWARE=0.
if [[ -z "${LIBGL_ALWAYS_SOFTWARE:-}" ]]; then
    export LIBGL_ALWAYS_SOFTWARE=1
    echo "  Note:    LIBGL_ALWAYS_SOFTWARE=1 (set automatically)"
    echo "           Override with: export LIBGL_ALWAYS_SOFTWARE=0"
fi

# Warn early if no display is available — the viewer needs X11 or Wayland.
if [[ -z "${DISPLAY:-}" ]] && [[ -z "${WAYLAND_DISPLAY:-}" ]]; then
    echo ""
    echo "WARNING: Neither DISPLAY nor WAYLAND_DISPLAY is set."
    echo "  The viewer requires an X11 or Wayland display server."
    echo "  For Apptainer on HPC: apptainer exec --env DISPLAY=\$DISPLAY ..."
    echo "  For SSH:              ssh -X user@host  (X11 forwarding)"
fi

echo ""
exec "${VIEWER_BIN}" --config "${CONFIG}"
