# Automated SUMA Figure Pipeline

This is the SUMA-first replacement for the manual workflow:

1. read one or more analysis result CSV files,
2. copy each channel statistic into the voxels listed in `fNIRS_template.txt`,
3. write a NIfTI overlay using `fNIRS_template.nii` and `fNIRS_template_mask.nii`,
4. optionally project the NIfTI to MNI_N27 surface datasets with `3dVol2Surf`
   and call `drive_suma_render.sh` to make a headless SUMA screenshot.

The older nilearn/matplotlib renderer is still in `render_figure.py`, but this
pipeline is meant to preserve the AFNI/SUMA rendering look.

## Convert CSVs to NIfTI overlays

```bash
./run_auto_suma_figures.sh \
  --out-dir auto_suma_outputs \
  250528_result/lrv-con.hbo.v7.allchannel.csv
```

Outputs:

- `auto_suma_outputs/nii/<input-stem>.nii`
- `auto_suma_outputs/audit/<input-stem>.channels.csv`

The audit CSV records the selected statistic column, selected p-value column,
raw channel value, significance decision, and rendered value for channels 1-16.

## Convert and render with SUMA

```bash
./run_auto_suma_figures.sh \
  --suma-spec /path/to/surface.spec \
  --out-dir auto_suma_outputs \
  250528_result/*.csv
```

Outputs also include:

- `auto_suma_outputs/png/<input-stem>.suma.png`

`drive_suma_render.sh` requires `afni`, `suma`, `DriveSuma`, `plugout_drive`,
`3dVol2Surf`, `Xvfb`, and either `convert` or `pnmtopng`. If AFNI/SUMA is
installed in `$HOME/abin`, the driver adds that directory.

The current SUMA driver assumes the MNI_N27 template layout:

- `MNI_N27_both.spec`
- `MNI_N27_lh.spec`
- `MNI_N27_rh.spec`
- `MNI_N27_SurfVol.nii`

Those files are present in `~/.afni/data/MNI_N27` on this machine.

The SUMA render defaults to the pial surface state and an anterior-facing view
starting from the left cardinal view, then rotating right 18 times. Override
these if needed:

```bash
SUMA_SURFACE_LABELS="lh.pial.gii rh.pial.gii" \
SUMA_VIEW_KEYS="ctrl+left right right right right right right right right right right right right right right right right right right" \
./run_auto_suma_figures.sh --suma-spec ~/.afni/data/MNI_N27/MNI_N27_both.spec ...
```

## Choosing result columns

Two separate column choices control each render:

- `--stat-col`: values used for color/intensity, such as `ancova_t`,
  `cohens_d`, `ancova_beta`, or `int_sans_t`.
- `--significance-p-col` or `--significance-bool-col`: values used to decide
  which channels are colored instead of zeroed.

By default, `auto_suma_figures.py` auto-detects a statistic column such as
`stat`, `ancova_t`, `t`, or an interaction t column, and auto-detects a p-value
column such as `pval`, `ancova_p`, or `ancova_p_fdr`. The audit CSV records
both the statistic column and the significance rule used for each channel.

Color channels by `ancova_t`, but decide significance using FDR-corrected
ANCOVA p-values:

```bash
./run_auto_suma_figures.sh \
  --stat-col ancova_t \
  --significance-p-col ancova_p_fdr \
  --alpha 0.05 \
  results_psd_agg_vs_hc_hard_easy/results_psd_agg_vs_psd_nagg.csv
```

Color channels by Cohen's d, but still decide significance using raw ANCOVA
p-values:

```bash
./run_auto_suma_figures.sh \
  --stat-col cohens_d \
  --significance-p-col ancova_p \
  --alpha 0.05 \
  results_psd_agg_vs_hc_hard_easy/results_psd_agg_vs_psd_nagg.csv
```

For boolean significance columns:

```bash
./run_auto_suma_figures.sh \
  --stat-col ancova_t \
  --significance-bool-col ancova_p_cluster_sig \
  --significance-mode bool \
  results_psd_agg_vs_hc_hard_easy/results_psd_agg_vs_psd_nagg.csv
```

To render every non-zero channel without p-value thresholding:

```bash
./run_auto_suma_figures.sh --no-threshold 250528_result/*.csv
```
