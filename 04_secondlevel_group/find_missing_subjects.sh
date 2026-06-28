#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 2 ]]; then
  echo "Usage: $0 <subjects.txt> <processed_files.txt>" >&2
  exit 1
fi

SUBJECT_FILE="$1"
PROCESSED_FILE="$2"

# Strip ANSI color codes (like \x1b[01;34m) and CRs
strip_ansi() {
  sed -r 's/\x1B\[[0-9;]*[mK]//g' | tr -d '\r'
}

# ---- Read & normalize subjects: allow space or newline separated ----
subjects_sorted=$(
  cat "$SUBJECT_FILE" \
  | strip_ansi \
  | tr ' \t' '\n' \
  | sed '/^$/d' \
  | sort -u
)

# ---- Extract IDs from processed list ----
# Prefer ID from filename: stats.<ID>.something
# Fallback: ID from directory segment .../<ID>/stats....
processed_sorted=$(
  cat "$PROCESSED_FILE" \
  | strip_ansi \
  | sed -n -E '
      s#.*stats\.([A-Za-z0-9]+)\..*#\1#p
    ' \
  ; cat "$PROCESSED_FILE" \
  | strip_ansi \
  | sed -n -E '
      s#^.*/([A-Za-z0-9]+)/stats\..*#\1#p
    ' \
  | sort -u
)

# Use comm to find subjects not in processed
comm -23 \
  <(printf "%s\n" "$subjects_sorted") \
  <(printf "%s\n" "$processed_sorted")
