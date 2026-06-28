#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PYTHON="${ROOT}/.venv/bin/python"

if [[ ! -x "$PYTHON" ]]; then
    PYTHON="$(command -v python3)"
fi

exec "$PYTHON" "$ROOT/auto_figures/auto_suma_figures.py" "$@"
