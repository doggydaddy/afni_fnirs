#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat >&2 <<'EOF'
usage: run_group_analyses.sh MODEL [--output-prefix PREFIX]
                              [--exclude-symptom-covariates] [--dry-run]

MODEL is "amplitude" or "na". Inputs are read from the current directory;
results are written there. Set GROUP_ANALYSIS_PYTHON to override Python.
EOF
}

(( $# >= 1 )) || { usage; exit 2; }
model=$1
shift
output_prefix=""
exclude_symptom_covariates=false
dry_run=false
while (( $# )); do
    case $1 in
        --output-prefix)
            (( $# >= 2 )) || { usage; exit 2; }
            output_prefix=$2
            shift 2
            ;;
        --dry-run)
            dry_run=true
            shift
            ;;
        --exclude-symptom-covariates)
            exclude_symptom_covariates=true
            shift
            ;;
        *)
            echo "[ERROR] unknown argument: $1" >&2
            usage
            exit 2
            ;;
    esac
done

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
analyzer="$script_dir/analyze_psd_agg_vs_psd_nagg.py"
if [[ -n ${GROUP_ANALYSIS_PYTHON:-} ]]; then
    python=$GROUP_ANALYSIS_PYTHON
elif [[ -x $script_dir/.venv/bin/python ]]; then
    python=$script_dir/.venv/bin/python
else
    python=python3
fi

case $model in
    amplitude)
        input_tokens=(v5 v7)
        output_tokens=(hard_easy corsi)
        contrast_labels=(hard-easy corsi)
        covariates_psd=(sans saps wais_matrix sud dose)
        ;;
    na)
        input_tokens=(na-hard-easy na-corsi)
        output_tokens=(na_hard_easy na_corsi)
        contrast_labels=(na-hard-easy na-corsi)
        covariates_psd=(wais_matrix sud dose)
        ;;
    *)
        echo "[ERROR] model must be 'amplitude' or 'na'" >&2
        usage
        exit 2
        ;;
esac

input_pairs=(lrv-con lrv-pcon pcon-con)
output_pairs=(psd_agg_vs_hc psd_agg_vs_psd_nagg psd_nagg_vs_hc)
positive_groups=(PSD+AGG PSD+AGG PSD-AGG)
negative_groups=(HC PSD-AGG HC)
covariate_sets=(hc psd hc)
covariates_hc=(wais_matrix sud)
if $exclude_symptom_covariates; then
    covariates_psd=(wais_matrix sud dose)
fi

run_one() {
    local pair_index=$1 contrast_index=$2
    local input output label covariate_set
    input="${input_pairs[$pair_index]}_${input_tokens[$contrast_index]}_allchannels.csv"
    output="results_${output_prefix}${output_pairs[$pair_index]}_${output_tokens[$contrast_index]}"
    label="${positive_groups[$pair_index]} vs ${negative_groups[$pair_index]} — ${contrast_labels[$contrast_index]}"
    covariate_set=${covariate_sets[$pair_index]}

    local -a covariates command
    if [[ $covariate_set == hc ]]; then
        covariates=("${covariates_hc[@]}")
    else
        covariates=("${covariates_psd[@]}")
    fi
    command=(
        "$python" "$analyzer"
        --input "$input"
        --output-dir "$output"
        --group-pos "${positive_groups[$pair_index]}"
        --group-neg "${negative_groups[$pair_index]}"
        --covariates "${covariates[@]}"
        --interactions
        --alpha 0.05
    )

    printf '====== %s ======\n' "$label"
    if $dry_run; then
        printf '  %q' "${command[@]}"
        printf '\n'
    else
        "${command[@]}"
    fi
}

for pair_index in "${!input_pairs[@]}"; do
    for contrast_index in "${!input_tokens[@]}"; do
        run_one "$pair_index" "$contrast_index"
    done
done
