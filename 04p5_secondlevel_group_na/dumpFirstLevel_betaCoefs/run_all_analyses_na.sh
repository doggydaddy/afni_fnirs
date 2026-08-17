#!/usr/bin/env zsh
PY=/mnt/highlands/karim_fnirs/afni_fnirs/04_secondlevel_group/.venv/bin/python
SCRIPT=analyze_psd_agg_vs_psd_nagg.py

# na model, group-level comparisons.
#   sans/saps are never included (per instruction).
#   HC comparisons drop dose too (all HCs have dose=0 -> perfect collinearity
#   with group), same as the original amplitude-modulated analysis.
#   PSD+AGG vs PSD-AGG (no HC involved) keeps dose.
COVARS_HC="wais_matrix sud"
COVARS_PSD="wais_matrix sud dose"

echo "====== PSD+AGG vs HC — na-hard-easy ======"
$PY $SCRIPT \
    --input lrv-con_na-hard-easy_allchannels.csv \
    --output-dir results_psd_agg_vs_hc_na_hard_easy \
    --group-pos "PSD+AGG" \
    --group-neg "HC" \
    --covariates ${=COVARS_HC} \
    --interactions \
    --alpha 0.05

echo ""
echo "====== PSD+AGG vs HC — na-corsi ======"
$PY $SCRIPT \
    --input lrv-con_na-corsi_allchannels.csv \
    --output-dir results_psd_agg_vs_hc_na_corsi \
    --group-pos "PSD+AGG" \
    --group-neg "HC" \
    --covariates ${=COVARS_HC} \
    --interactions \
    --alpha 0.05

echo ""
echo "====== PSD+AGG vs PSD-AGG — na-hard-easy ======"
$PY $SCRIPT \
    --input lrv-pcon_na-hard-easy_allchannels.csv \
    --output-dir results_psd_agg_vs_psd_nagg_na_hard_easy \
    --covariates ${=COVARS_PSD} \
    --interactions \
    --alpha 0.05

echo ""
echo "====== PSD+AGG vs PSD-AGG — na-corsi ======"
$PY $SCRIPT \
    --input lrv-pcon_na-corsi_allchannels.csv \
    --output-dir results_psd_agg_vs_psd_nagg_na_corsi \
    --covariates ${=COVARS_PSD} \
    --interactions \
    --alpha 0.05

echo ""
echo "====== PSD-AGG vs HC — na-hard-easy ======"
$PY $SCRIPT \
    --input pcon-con_na-hard-easy_allchannels.csv \
    --output-dir results_psd_nagg_vs_hc_na_hard_easy \
    --group-pos "PSD-AGG" \
    --group-neg "HC" \
    --covariates ${=COVARS_HC} \
    --interactions \
    --alpha 0.05

echo ""
echo "====== PSD-AGG vs HC — na-corsi ======"
$PY $SCRIPT \
    --input pcon-con_na-corsi_allchannels.csv \
    --output-dir results_psd_nagg_vs_hc_na_corsi \
    --group-pos "PSD-AGG" \
    --group-neg "HC" \
    --covariates ${=COVARS_HC} \
    --interactions \
    --alpha 0.05
