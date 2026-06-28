#!/usr/bin/env zsh
PY=/mnt/speyside/karim_fnirs/afni_fnirs/04_secondlevel_group/.venv/bin/python
SCRIPT=analyze_psd_agg_vs_psd_nagg.py

# Covariate sets
#   HC comparisons drop `dose` (all HCs have dose=0 → perfect collinearity with group)
#   PSD+AGG vs PSD-AGG keeps `dose`
COVARS_HC="sans saps wais_matrix sud"
COVARS_PSD="sans saps wais_matrix sud dose"

echo "====== PSD+AGG vs HC — hard-easy ======"
$PY $SCRIPT \
    --input lrv-con_v5_allchannels.csv \
    --output-dir results_psd_agg_vs_hc_hard_easy \
    --group-pos "PSD+AGG" \
    --group-neg "HC" \
    --covariates ${=COVARS_HC} \
    --interactions \
    --alpha 0.05

echo ""
echo "====== PSD+AGG vs HC — corsi ======"
$PY $SCRIPT \
    --input lrv-con_v7_allchannels.csv \
    --output-dir results_psd_agg_vs_hc_corsi \
    --group-pos "PSD+AGG" \
    --group-neg "HC" \
    --covariates ${=COVARS_HC} \
    --interactions \
    --alpha 0.05

echo ""
echo "====== PSD+AGG vs PSD-AGG — hard-easy ======"
$PY $SCRIPT \
    --input lrv-pcon_v5_allchannels.csv \
    --output-dir results_psd_agg_vs_psd_nagg_hard_easy \
    --covariates ${=COVARS_PSD} \
    --interactions \
    --alpha 0.05

echo ""
echo "====== PSD+AGG vs PSD-AGG — corsi ======"
$PY $SCRIPT \
    --input lrv-pcon_v7_allchannels.csv \
    --output-dir results_psd_agg_vs_psd_nagg_corsi \
    --covariates ${=COVARS_PSD} \
    --interactions \
    --alpha 0.05

echo ""
echo "====== PSD-AGG vs HC — hard-easy ======"
$PY $SCRIPT \
    --input pcon-con_v5_allchannels.csv \
    --output-dir results_psd_nagg_vs_hc_hard_easy \
    --group-pos "PSD-AGG" \
    --group-neg "HC" \
    --covariates ${=COVARS_HC} \
    --interactions \
    --alpha 0.05

echo ""
echo "====== PSD-AGG vs HC — corsi ======"
$PY $SCRIPT \
    --input pcon-con_v7_allchannels.csv \
    --output-dir results_psd_nagg_vs_hc_corsi \
    --group-pos "PSD-AGG" \
    --group-neg "HC" \
    --covariates ${=COVARS_HC} \
    --interactions \
    --alpha 0.05
