#!/usr/bin/env bash
# build_viewer.sh – Build and test the Rust long_read_viewer locally.
#
# Builds the viewer binary and runs all unit tests using either a local Rust
# toolchain (cargo) or Docker, whichever is available.  The toy dataset
# directory is accepted as an optional positional argument (relative path
# from the repository root).
#
# Usage:
#   bash scripts/build_viewer.sh [toy_dataset_dir]
#
# Examples:
#   bash scripts/build_viewer.sh                          # default dataset path
#   bash scripts/build_viewer.sh resources/toy_dataset   # explicit relative path
#   bash scripts/build_viewer.sh /abs/path/to/dataset    # absolute path

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

# ── Resolve toy dataset directory ────────────────────────────────────────────
DATASET_REL="${1:-resources/toy_dataset}"
if [[ "${DATASET_REL}" = /* ]]; then
    DATASET_DIR="${DATASET_REL}"
else
    DATASET_DIR="${REPO_ROOT}/${DATASET_REL}"
fi

if [[ ! -d "${DATASET_DIR}" ]]; then
    echo "ERROR: toy dataset directory not found: ${DATASET_DIR}" >&2
    exit 1
fi

MANIFEST="${DATASET_DIR}/toy_manifest.json"

echo "=== long_read_viewer build ==="
echo "Repository root: ${REPO_ROOT}"
echo "Toy dataset:     ${DATASET_DIR}"
echo ""

# ── Build and test ───────────────────────────────────────────────────────────
if command -v cargo &>/dev/null; then
    # ── Local cargo path ─────────────────────────────────────────────────────
    echo "Using local cargo ($(cargo --version))"
    echo ""

    cd "${REPO_ROOT}/viewer"

    echo "--- Building release binary ---"
    cargo build --release
    echo ""

    echo "--- Running unit tests ---"
    cargo test
    echo ""

    BINARY="${REPO_ROOT}/viewer/target/release/long_read_viewer"
    echo "=== Build complete ==="
    echo ""
    echo "Binary: ${BINARY}"
    if [[ -f "${MANIFEST}" ]]; then
        echo ""
        echo "To launch the viewer with the toy dataset:"
        echo "  ${BINARY}"
        echo ""
        echo "  Then choose  File > Load Manifest  and open:"
        echo "  ${MANIFEST}"
        echo ""
        echo "  (A display/windowing system is required to open the GUI.)"
    fi

elif command -v docker &>/dev/null; then
    # ── Docker path (fallback when cargo is not installed) ───────────────────
    echo "cargo not found – using Docker ($(docker --version | head -1))"
    echo ""

    echo "--- Building Docker image (compile + cargo test) ---"
    docker build \
        -f "${REPO_ROOT}/viewer/Dockerfile" \
        --target tester \
        -t long_read_viewer:local \
        "${REPO_ROOT}"
    echo ""

    echo "--- Building runtime image ---"
    docker build \
        -f "${REPO_ROOT}/viewer/Dockerfile" \
        --target runtime \
        -t long_read_viewer:runtime \
        "${REPO_ROOT}"
    echo ""

    echo "=== Build complete ==="
    echo ""
    echo "The viewer binary is packaged in Docker image 'long_read_viewer:runtime'."
    if [[ -f "${MANIFEST}" ]]; then
        echo ""
        echo "To launch the viewer with the toy dataset (requires a display):"
        echo "  docker run --rm -it \\"
        echo "    -e DISPLAY=\${DISPLAY} \\"
        echo "    -v /tmp/.X11-unix:/tmp/.X11-unix \\"
        echo "    -v \"${DATASET_DIR}:/data\" \\"
        echo "    long_read_viewer:runtime"
        echo ""
        echo "  Then choose  File > Load Manifest  and open:"
        echo "  /data/toy_manifest.json"
    fi

else
    echo "ERROR: neither 'cargo' nor 'docker' found in PATH." >&2
    echo "Install Rust (https://rustup.rs) or Docker to build the viewer." >&2
    exit 1
fi
