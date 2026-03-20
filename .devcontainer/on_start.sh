#!/usr/bin/env bash
# on_start.sh — Start the visualization server on every Codespace start/restart.
#
# This script runs automatically each time the Codespace container starts.
# It launches the Python server in the background and logs output to
# /tmp/lrv-server.log so the terminal remains usable.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
CONFIG="${REPO_ROOT}/toy_example_config.tsv"

# Safety check: if the config is missing (e.g. workspace was reset), re-run
# the one-time setup before starting the server.
if [[ ! -f "${CONFIG}" ]]; then
    echo "Server config not found — running one-time setup first..."
    bash "${SCRIPT_DIR}/on_create.sh"
fi

LOG="/tmp/lrv-server.log"

echo "Starting visualization server on port 8080..."
nohup python3 "${REPO_ROOT}/server/app.py" \
    --config "${CONFIG}" \
    --port 8080 \
    --host 0.0.0.0 \
    > "${LOG}" 2>&1 &
SERVER_PID=$!

echo "Server started (PID ${SERVER_PID})."
echo "Logs: ${LOG}"
echo "Open the Ports panel in VS Code and click on port 8080 to view."
