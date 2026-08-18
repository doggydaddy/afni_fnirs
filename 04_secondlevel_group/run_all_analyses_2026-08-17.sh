#!/usr/bin/env bash
set -euo pipefail

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
exec "$script_dir/run_group_analyses.sh" amplitude \
    --output-prefix 2026-08-17_ \
    --exclude-symptom-covariates \
    "$@"
