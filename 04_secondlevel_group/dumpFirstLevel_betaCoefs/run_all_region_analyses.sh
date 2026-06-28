#!/usr/bin/env zsh
# Region-level (4-target) version of run_all_analyses.sh
# Mirrors the 6 comparisons but operates on *_allregions.csv produced by
# pool_to_regions.py, with --targets / --target-labels for the 4 anatomical
# regions: left dlPFC, left mPFC, right mPFC, right dlPFC.
#
# --no-cluster: cluster-based permutation is dropped here because with only 4
# spatially-pre-specified ROIs it is largely redundant with permutation max-stat
# and provides no additional FWER protection. FDR-BH, Hommel, and perm max-stat
# remain the primary multiple-comparison controls.

PY=/mnt/speyside/karim_fnirs/afni_fnirs/04_secondlevel_group/.venv/bin/python
SCRIPT=analyze_psd_agg_vs_psd_nagg.py

# Covariate sets (same logic as the channel runner)
COVARS_HC="sans saps wais_matrix sud"
COVARS_PSD="sans saps wais_matrix sud dose"

# 4-region target spec (column names use underscores for formula safety;
# display labels keep human-readable names)
TARGETS="left_dlPFC left_mPFC right_mPFC right_dlPFC"
LABELS=('Left dlPFC' 'Left mPFC' 'Right mPFC' 'Right dlPFC')

echo "====== PSD+AGG vs HC — hard-easy — regions ======"
$PY $SCRIPT \
    --input lrv-con_v5_allregions.csv \
    --output-dir results_psd_agg_vs_hc_hard_easy_regions \
    --group-pos "PSD+AGG" \
    --group-neg "HC" \
    --covariates ${=COVARS_HC} \
    --targets ${=TARGETS} \
    --target-labels "${LABELS[@]}" \
    --no-cluster \
    --interactions \
    --alpha 0.05

echo ""
echo "====== PSD+AGG vs HC — corsi — regions ======"
$PY $SCRIPT \
    --input lrv-con_v7_allregions.csv \
    --output-dir results_psd_agg_vs_hc_corsi_regions \
    --group-pos "PSD+AGG" \
    --group-neg "HC" \
    --covariates ${=COVARS_HC} \
    --targets ${=TARGETS} \
    --target-labels "${LABELS[@]}" \
    --no-cluster \
    --interactions \
    --alpha 0.05

echo ""
echo "====== PSD+AGG vs PSD-AGG — hard-easy — regions ======"
$PY $SCRIPT \
    --input lrv-pcon_v5_allregions.csv \
    --output-dir results_psd_agg_vs_psd_nagg_hard_easy_regions \
    --covariates ${=COVARS_PSD} \
    --targets ${=TARGETS} \
    --target-labels "${LABELS[@]}" \
    --no-cluster \
    --interactions \
    --alpha 0.05

echo ""
echo "====== PSD+AGG vs PSD-AGG — corsi — regions ======"
$PY $SCRIPT \
    --input lrv-pcon_v7_allregions.csv \
    --output-dir results_psd_agg_vs_psd_nagg_corsi_regions \
    --covariates ${=COVARS_PSD} \
    --targets ${=TARGETS} \
    --target-labels "${LABELS[@]}" \
    --no-cluster \
    --interactions \
    --alpha 0.05

echo ""
echo "====== PSD-AGG vs HC — hard-easy — regions ======"
$PY $SCRIPT \
    --input pcon-con_v5_allregions.csv \
    --output-dir results_psd_nagg_vs_hc_hard_easy_regions \
    --group-pos "PSD-AGG" \
    --group-neg "HC" \
    --covariates ${=COVARS_HC} \
    --targets ${=TARGETS} \
    --target-labels "${LABELS[@]}" \
    --no-cluster \
    --interactions \
    --alpha 0.05

echo ""
echo "====== PSD-AGG vs HC — corsi — regions ======"
$PY $SCRIPT \
    --input pcon-con_v7_allregions.csv \
    --output-dir results_psd_nagg_vs_hc_corsi_regions \
    --group-pos "PSD-AGG" \
    --group-neg "HC" \
    --covariates ${=COVARS_HC} \
    --targets ${=TARGETS} \
    --target-labels "${LABELS[@]}" \
    --no-cluster \
    --interactions \
    --alpha 0.05
