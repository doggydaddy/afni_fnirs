#!/usr/bin/env bash
set -euo pipefail

if (( $# != 6 )); then
    echo "usage: $0 GROUP1_LIST GROUP1_LABEL GROUP2_LIST GROUP2_LABEL VOLUME OUTPUT.csv" >&2
    exit 2
fi

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
repo_root=$(cd -- "$script_dir/.." && pwd)

group1_list=$1
group1_label=$2
group2_list=$3
group2_label=$4
volume=$5
output=$6

template_mask=${TEMPLATE_MASK:-$repo_root/templates/biopac16ch_template_mask.nii}
covariate_file=${COVARIATE_FILE:-$script_dir/select.cov.local.1D}
if [[ -n ${GROUP_ANALYSIS_PYTHON:-} ]]; then
    python=$GROUP_ANALYSIS_PYTHON
elif [[ -x $script_dir/.venv/bin/python ]]; then
    python=$script_dir/.venv/bin/python
else
    python=python3
fi

# Resolve caller-relative inputs before running helpers from a temporary folder.
group1_list=$(realpath "$group1_list")
group2_list=$(realpath "$group2_list")
template_mask=$(realpath "$template_mask")
covariate_file=$(realpath "$covariate_file")
output=$(realpath -m "$output")

tmp_dir=$(mktemp -d "${TMPDIR:-/tmp}/fnirs-dump-channels.XXXXXX")
trap 'rm -rf -- "$tmp_dir"' EXIT

for channel in {1..16}; do
    "$script_dir/grab_betas.sh" \
        "$group1_list" "$group1_label" \
        "$group2_list" "$group2_label" \
        "$volume" "$template_mask" "$channel" "$tmp_dir/output.txt"
    "$python" "$script_dir/add_covariates.py" \
        -input "$tmp_dir/output.txt" \
        -covar "$covariate_file" \
        -output "$tmp_dir/TMP_ch${channel}.csv"
    "$script_dir/relabel_groups.sh" "$tmp_dir/TMP_ch${channel}.csv"
done

"$python" "$script_dir/merge_channel_dumps.py" \
    --input-dir "$tmp_dir" \
    --output "$output"
