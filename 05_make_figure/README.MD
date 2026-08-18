# Figure-making routines

This directory converts channel-level fNIRS statistics into template-space
NIfTI overlays and brain figures. The recommended workflow is the automated
AFNI/SUMA pipeline:

```text
result CSV -> channel mask/audit -> template-space NIfTI -> SUMA PNG
```

The older notebook and matplotlib/nilearn renderer are retained for manual or
custom figures, but `run_auto_suma_figures.sh` is the main entry point for
reproducible analysis figures.

## Files

- `run_auto_suma_figures.sh`: recommended command-line entry point.
- `auto_figures/auto_suma_figures.py`: reads CSVs, selects statistics and
  significance, writes NIfTI overlays and audit CSVs, and optionally invokes
  SUMA.
- `drive_suma_render.sh`: projects an overlay to the MNI N27 surfaces and
  takes a headless SUMA screenshot.
- `auto_figures/render_figure.py`: optional matplotlib/nilearn renderer.
- `fNIRS_template.txt`: voxel-to-channel mapping for channels 1-16.
- `fNIRS_template.nii`: spatial reference for generated overlays.
- `fNIRS_template_mask.nii`: mask applied to generated overlays.
- `brain.nii`, `results_template.csv`, and `toTemplate.ipynb`: legacy manual
  workflow resources.

## Requirements

NIfTI conversion requires Python with `numpy` and `nibabel`. The repository's
local `.venv` is used automatically when available; otherwise the wrapper
uses `python3` from `PATH`.

SUMA rendering additionally requires:

- AFNI, SUMA, `DriveSuma`, `plugout_drive`, and `3dVol2Surf`
- `Xvfb` for headless rendering
- ImageMagick `convert` or `pnmtopng`
- the AFNI MNI N27 surface package

On the analysis machine the surface specification is normally:

```text
~/.afni/data/MNI_N27/MNI_N27_both.spec
```

The renderer expects the corresponding left/right specifications and
`MNI_N27_SurfVol.nii` in the same directory.

## Input CSV format

Each input CSV must contain one row per channel and a `channel` or `ch`
column. Channel values such as `channel1`, `Ch1`, or `1` are accepted. The
template covers channels 1-16; absent channels are written as zero.

Two independent selections control a figure:

1. The statistic column controls color and intensity, for example
   `ancova_t`, `cohens_d`, or `int_wais_matrix_t`.
2. A p-value or boolean column controls which channels are rendered.

Always choose the significance rule to match the inferential claim in the
report. For primary group effects in the current analyses, the preferred map
uses `ancova_t` masked by the cluster-permutation FWER result
`ancova_p_cluster_sig`. For interaction maps, use the interaction's t-value
and its matching FDR boolean column.

## Recommended usage

Run commands from `05_make_figure`. Write outputs to a private analysis
directory outside the repository because figures and audit files may contain
sensitive results.

### 1. Create NIfTI overlays and audits only

```bash
./run_auto_suma_figures.sh \
  --out-dir /path/to/private-analysis/figures_cluster_fwer \
  --stat-col ancova_t \
  --significance-bool-col ancova_p_cluster_sig \
  --significance-mode bool \
  path/to/results_one.csv path/to/results_two.csv
```

Do not include `*_interactions.csv` in this primary-effect command.

### 2. Create NIfTI overlays, audits, and SUMA PNGs

Add `--suma-spec`:

```bash
./run_auto_suma_figures.sh \
  --out-dir /path/to/private-analysis/figures_cluster_fwer \
  --stat-col ancova_t \
  --significance-bool-col ancova_p_cluster_sig \
  --significance-mode bool \
  --suma-spec ~/.afni/data/MNI_N27/MNI_N27_both.spec \
  path/to/results_one.csv path/to/results_two.csv
```

Multiple CSV paths can be supplied in one invocation. SUMA renders them
sequentially.

### 3. Render an FDR-corrected interaction

Each covariate must be rendered separately so its statistic is paired with
the correct significance column:

```bash
./run_auto_suma_figures.sh \
  --out-dir /path/to/private-analysis/figures_interactions_fdr \
  --stat-col int_wais_matrix_t \
  --significance-bool-col int_wais_matrix_sig \
  --significance-mode bool \
  --suma-spec ~/.afni/data/MNI_N27/MNI_N27_both.spec \
  path/to/results_interactions.csv
```

For SUD or dose interactions, pair columns in the same way, for example
`int_sud_t` with `int_sud_sig`, or `int_dose_t` with `int_dose_sig`.

### 4. Threshold using a p-value column

Boolean corrected-significance columns are preferred when present. A numeric
p-value can instead be thresholded explicitly:

```bash
./run_auto_suma_figures.sh \
  --stat-col ancova_t \
  --significance-p-col ancova_p_fdr \
  --significance-mode p \
  --alpha 0.05 \
  --out-dir output_fdr \
  path/to/results.csv
```

### 5. Exploratory unthresholded map

```bash
./run_auto_suma_figures.sh \
  --stat-col ancova_t \
  --no-threshold \
  --out-dir output_unthresholded \
  path/to/results.csv
```

Label unthresholded output clearly; it does not represent statistically
significant channels.

## Output layout

For `--out-dir OUTPUT`, the pipeline creates:

```text
OUTPUT/
├── audit/
│   └── INPUT_STEM.channels.csv
├── nii/
│   └── INPUT_STEM.nii
└── png/
    └── INPUT_STEM.suma.png
```

The `png` directory is created only when `--suma-spec` is supplied. Audit
files record the source file, statistic and significance columns, raw values,
significance decision, and final rendered value for every channel. Keep them
with published figures as provenance.

A null comparison is valid output: its NIfTI contains zeros and its SUMA PNG
shows an uncolored brain. The command reports the number of rendered channels
for every input.

## Automatic column selection

If columns are not specified, the program attempts to identify common names.
This is convenient for quick inspection but can choose a raw p-value or only
the first interaction family in a multi-interaction CSV. For final figures,
pass `--stat-col` and the significance column explicitly.

Useful options:

- `--stat-col COLUMN`: statistic used for the overlay.
- `--significance-bool-col COLUMN --significance-mode bool`: boolean mask.
- `--significance-p-col COLUMN --significance-mode p`: numeric p-value mask.
- `--alpha VALUE`: numeric p-value threshold; default `0.05`.
- `--no-threshold`: render every nonzero statistic.
- `--template-txt`, `--template-nii`, `--template-mask`: override templates.
- `--suma-driver`: override the SUMA driver script.
- `--suma-spec`: enable SUMA screenshots.

Use `./run_auto_suma_figures.sh --help` for the complete current option list.

## SUMA view customization

The default renderer uses the pial surfaces and an anterior-facing view.
Environment variables can override the surfaces or simulated view keys:

```bash
SUMA_SURFACE_LABELS="lh.pial.gii rh.pial.gii" \
SUMA_VIEW_KEYS="ctrl+left right right right right right" \
./run_auto_suma_figures.sh \
  --suma-spec ~/.afni/data/MNI_N27/MNI_N27_both.spec \
  --out-dir output \
  path/to/results.csv
```

Other advanced variables accepted by the driver are `SUMA_DISPLAY` and
`SUMA_NPB`, which select the X display and AFNI/SUMA port block. They are
mainly useful when another AFNI/SUMA process is already running.

## Optional matplotlib/nilearn renderer

`auto_figures/render_figure.py` renders an existing overlay without launching
SUMA. It may download nilearn's fsaverage data on first use.

```bash
./.venv/bin/python auto_figures/render_figure.py \
  --overlay output/nii/example.nii \
  --out output/example.png \
  --view anterior \
  --surface inflated \
  --surface-mode patch \
  --labels
```

Available views include `anterior`, `posterior`, `dorsal`, `ventral`, `left`,
and `right`; a custom `elevation,azimuth` pair is also accepted. Run
`./.venv/bin/python auto_figures/render_figure.py --help` for surface,
camera, label, patch, and cube controls.

## Legacy manual workflow

`results_template.csv` and `toTemplate.ipynb` implement the older manual
CSV-to-template conversion. They remain for reproducing historical work, but
new analysis figures should use the automated pipeline because it records the
selected statistic, significance rule, and per-channel rendered values.

## Validation checklist

Before using figures in a report or manuscript:

1. Confirm that every input is from the intended analysis date and model.
2. Confirm that the statistic and significance columns describe the same
   effect.
3. Compare the audit's significant-channel count and channel identities with
   the statistical report.
4. Check that all expected `.nii`, audit `.csv`, and `.png` files exist and
   are nonempty.
5. Keep primary effects, interactions, correction methods, and exploratory
   maps in separately named output directories.

More implementation-specific notes are available in
`auto_figures/README.md`.
