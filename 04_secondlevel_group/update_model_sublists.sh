#!/usr/bin/env bash
set -euo pipefail

if (( $# < 1 || $# > 2 )); then
    echo "usage: $0 {amplitude|na} [OUTPUT_DIR]" >&2
    exit 2
fi

model=$1
output_dir=${2:-.}
script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
repo_root=$(cd -- "$script_dir/.." && pwd)
data_root=${FNIRS_DATA_ROOT:-$repo_root/../data}

case $model in
    amplitude)
        filename_suffix=""
        ;;
    na)
        filename_suffix=".na"
        ;;
    *)
        echo "[ERROR] model must be 'amplitude' or 'na'" >&2
        exit 2
        ;;
esac

[[ -d $data_root ]] || { echo "[ERROR] data root not found: $data_root" >&2; exit 1; }
mkdir -p -- "$output_dir"

group_names=(con pcon lrv)
group_dirs=(friska_kontroller patienter_kontroll patienter_lrv)
for index in "${!group_names[@]}"; do
    group=${group_names[$index]}
    source_dir="$data_root/${group_dirs[$index]}"
    [[ -d $source_dir ]] || { echo "[ERROR] group directory not found: $source_dir" >&2; exit 1; }
    for chromophore in hbo hbr; do
        output="$output_dir/${group}.${chromophore}${filename_suffix}.sublist.txt"
        find "$source_dir" -type f \
            -name "stats.*.${chromophore}.mem${filename_suffix}.nii.gz" \
            | sort > "$output"
        printf '[OK] %4d files -> %s\n' "$(wc -l < "$output")" "$output"
    done
done
