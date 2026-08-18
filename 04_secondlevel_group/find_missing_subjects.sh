#!/usr/bin/env bash
set -euo pipefail

if (( $# != 2 )); then
    echo "usage: $0 SUBJECT_IDS PROCESSED_FILES" >&2
    exit 2
fi

subject_file=$1
processed_file=$2
[[ -f $subject_file ]] || { echo "[ERROR] missing file: $subject_file" >&2; exit 1; }
[[ -f $processed_file ]] || { echo "[ERROR] missing file: $processed_file" >&2; exit 1; }

tmp_dir=$(mktemp -d "${TMPDIR:-/tmp}/fnirs-find-missing.XXXXXX")
trap 'rm -rf -- "$tmp_dir"' EXIT

sed -E 's/\x1B\[[0-9;]*[mK]//g' "$subject_file" \
    | tr -d '\r' \
    | tr ' \t' '\n' \
    | sed '/^$/d' \
    | sort -u > "$tmp_dir/subjects"

sed -E 's/\x1B\[[0-9;]*[mK]//g' "$processed_file" \
    | tr -d '\r' \
    | sed -E -n 's#.*stats\.([^.]+)\..*#\1#p' \
    | sort -u > "$tmp_dir/processed"

comm -23 "$tmp_dir/subjects" "$tmp_dir/processed"
