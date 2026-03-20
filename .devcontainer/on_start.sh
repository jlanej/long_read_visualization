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
TOY_OUTPUT="${REPO_ROOT}/toy_example_output"

# Safety check: if the config is missing (e.g. workspace was reset), re-run
# the one-time setup before starting the server.
if [[ ! -f "${CONFIG}" ]] || [[ ! -d "${TOY_OUTPUT}" ]]; then
    echo "Server config not found — running one-time setup first..."
    bash "${SCRIPT_DIR}/on_create.sh"
fi

LOG="/tmp/lrv-server.log"

# Stop any stale server process from previous starts so port 8080 is clean.
if command -v pgrep >/dev/null 2>&1; then
    for pid in $(pgrep -f "server/app.py.*--port 8080" || true); do
        if [[ -n "${pid}" ]]; then
            kill "${pid}" 2>/dev/null || true
        fi
    done
fi

echo "Starting visualization server on port 8080..."
nohup python3 "${REPO_ROOT}/server/app.py" \
    --config "${CONFIG}" \
    --port 8080 \
    --host 0.0.0.0 \
    > "${LOG}" 2>&1 &
SERVER_PID=$!

# Wait briefly for startup and verify the server is actually serving content.
for _ in {1..20}; do
    if curl -fsS "http://127.0.0.1:8080/" >/dev/null 2>&1; then
        break
    fi
    sleep 1
done

if ! curl -fsS "http://127.0.0.1:8080/" >/dev/null 2>&1; then
    echo "ERROR: Server failed to serve content on port 8080."
    echo "Last server log lines:"
    tail -n 40 "${LOG}" || true
    exit 1
fi

echo "Server started (PID ${SERVER_PID})."
echo "Logs: ${LOG}"
echo "Open the Ports panel in VS Code and click on port 8080 to view."
