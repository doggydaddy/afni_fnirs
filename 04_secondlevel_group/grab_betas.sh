#!/usr/bin/env bash
set -euo pipefail

if (( $# < 7 || $# > 8 )); then
    echo "usage: $0 GROUP1_LIST GROUP1_LABEL GROUP2_LIST GROUP2_LABEL VOLUME MASK CHANNEL [OUTPUT]" >&2
    exit 2
fi

group1_list=$1
group1_label=$2
group2_list=$3
group2_label=$4
volume=$5
channel_mask=$6
channel=$7
output=${8:-output.txt}

for path in "$group1_list" "$group2_list" "$channel_mask"; do
    [[ -f $path ]] || { echo "[ERROR] missing input: $path" >&2; exit 1; }
done
for command in 3dcalc 3dROIstats; do
    command -v "$command" >/dev/null || {
        echo "[ERROR] required AFNI command not found on PATH: $command" >&2
        exit 127
    }
done
[[ $volume =~ ^[0-9]+$ ]] || { echo "[ERROR] volume must be an integer" >&2; exit 2; }
[[ $channel =~ ^([1-9]|1[0-6])$ ]] || { echo "[ERROR] channel must be 1-16" >&2; exit 2; }

tmp_dir=$(mktemp -d "${TMPDIR:-/tmp}/fnirs-grab-betas.XXXXXX")
trap 'rm -rf -- "$tmp_dir"' EXIT

mapfile -t group1_files < <(sed "s/$/[$volume]/" "$group1_list")
mapfile -t group2_files < <(sed "s/$/[$volume]/" "$group2_list")
(( ${#group1_files[@]} > 0 )) || { echo "[ERROR] empty subject list: $group1_list" >&2; exit 1; }
(( ${#group2_files[@]} > 0 )) || { echo "[ERROR] empty subject list: $group2_list" >&2; exit 1; }

mask="$tmp_dir/channel_mask.nii"
3dcalc -a "$channel_mask" -expr "within(a,$channel,$channel)" -prefix "$mask"
3dROIstats -mask "$mask" "${group1_files[@]}" > "$tmp_dir/group1.tsv"
3dROIstats -mask "$mask" "${group2_files[@]}" > "$tmp_dir/group2.tsv"

awk -v group="$group1_label" 'BEGIN { OFS="\t" } { print $0, group }' \
    "$tmp_dir/group1.tsv" > "$output"
tail -n +2 "$tmp_dir/group2.tsv" | \
    awk -v group="$group2_label" 'BEGIN { OFS="\t" } { print $0, group }' \
    >> "$output"
