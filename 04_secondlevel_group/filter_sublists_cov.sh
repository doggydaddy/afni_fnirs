#!/usr/bin/env bash
set -euo pipefail

if (( $# < 2 || $# > 3 )); then
    echo "usage: $0 SUBLIST COVARIATES.1D [OUTPUT]" >&2
    exit 2
fi

sublist=$1
covariates=$2
output=${3:-$sublist}
[[ -f $sublist ]] || { echo "[ERROR] missing sublist: $sublist" >&2; exit 1; }
[[ -f $covariates ]] || { echo "[ERROR] missing covariates: $covariates" >&2; exit 1; }

tmp=$(mktemp "${TMPDIR:-/tmp}/fnirs-filter-sublist.XXXXXX")
trap 'rm -f -- "$tmp"' EXIT

awk '
    NR == FNR {
        if (FNR > 1 && $1 != "") available[$1] = 1
        next
    }
    {
        count = split($0, parts, "/")
        split(parts[count], fields, ".")
        subject = fields[2]
        if (subject in available) {
            print
        } else {
            printf "[WARN] excluding subject absent from covariates: %s\n", subject > "/dev/stderr"
        }
    }
' "$covariates" "$sublist" > "$tmp"

mkdir -p -- "$(dirname -- "$output")"
mv -- "$tmp" "$output"
trap - EXIT
printf '[OK] %d subjects -> %s\n' "$(wc -l < "$output")" "$output"
