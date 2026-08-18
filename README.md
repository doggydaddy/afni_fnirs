# fNIRS processing pipeline using AFNI program package

Corsi-span working-memory task, 16-channel fNIRS (biopac montage), three
groups: healthy controls (`friska_kontroller` / con / HC), psychiatric
patients with aggression history (`patienter_lrv` / lrv / PSD+AGG),
psychiatric patients without (`patienter_kontroll` / pcon / PSD-AGG). Data
lives outside this repo, as a sibling directory: `../data/<group>/<subject>/`.

Pipeline stages are numbered directories, run in order. `04p5_secondlevel_group_na`
is a second, parallel group-analysis branch for an alternate first-level
model (see `03_firstlevel_glm`) — not a later pipeline stage, a sibling to `04`.

```
01_dataConversion       fnirsoft export -> NIfTI + block-timing files
02_preprocessing        scale to %-signal-change, bandpass filter
03_firstlevel_glm       per-subject GLM (two model variants)
04_secondlevel_group    group-level stats, amplitude-modulated model
04p5_secondlevel_group_na   group-level stats, na model
05_make_figure          channel-statistic CSV -> brain-surface figure
templates/              shared 16-channel montage mask + coordinates
```

Each subject's data folder accumulates files with prefixes marking pipeline
stage: raw fnirsoft exports -> `<subj>.<type>.<task>.fnirs.nii` (01) ->
`mean_`/`scl_`/`prp_` (02) -> `stats.<subj>.<type>.<task>[.na].nii.gz` (03).
`<type>` is `hbo` or `hbr` (oxy-/deoxyhemoglobin); `<task>` is `mem` (corsi,
the task analyzed throughout this pipeline) or `gonogo` (converted, not
analyzed downstream).

For the full narrative of a major timing bug found and fixed in this
pipeline on 2026-08-17 (every subject's block onsets were offset from the
actual recording by a per-subject amount), see
`analysis_2026-08-17/analysis_report_2026-08-17.md`. The technical
descriptions below reflect the pipeline as it stands after that fix.

---

## 01_dataConversion

Converts fnirsoft's raw per-channel text exports into AFNI-readable NIfTI
+ timing files.

**`fnirsoft2nii/export2nii.py`** — the core converter.

```
python3 export2nii.py -foldername <subject dir> -tag {hbo,hbr,hbt,oxy} -prefix <subj_id> [-skip_nii]
```

What it does, technically:

1. Reads the subject's `*.Block#.txt` / `*.Time#.txt` / `*Marker#.txt`
   triplet(s) for the requested tag (`find_pairs`). A subject can have
   more than one triplet (repeated/aborted acquisition attempts); each is
   tried in turn.
2. `parse_pair` loads Block (16-channel signal) and Time (sample
   timestamps, seconds since acquisition start) into one array.
3. `parse_marker` loads the marker log and, for each marker type, finds
   the nearest actual sample timestamp (`find_nearest`) — markers are
   logged on their own clock, slightly offset from the signal's sample
   grid.
4. Two task segments are cut out of the continuous recording by marker
   ID: go/no-go is markers 107→108; corsi ("mem") is markers -5→106 (task
   start → task end), with block onsets in between at markers 101–105
   (block 1–5; block 1 is a non-scored practice/instruction block, dropped
   downstream).
5. `marker2onsetdur` builds the per-block onset/duration table from the
   marker timestamps. **Onsets are written relative to marker -5** (task
   start), matching frame 0 of the NIfTI segment `data2nii` builds in step
   6 below — `-5` is subtracted from every marker's timestamp before
   computing onsets. (This is the fix for the 2026-08-17 bug: onsets were
   previously written in the raw recording's absolute clock, which doesn't
   match the exported `.nii`'s own zeroed time axis.)
6. A marker-set validation check (exactly one each of -5, 101–105, 106)
   runs right after the segment is cut; a subject with a malformed set
   (duplicate or missing markers — seen in this dataset) is skipped with a
   clear log line rather than crashing the whole batch or silently writing
   a wrong result.
7. If the subject ID is listed in `../../data/mirror_list.txt`, channels
   are laterally flipped (`mirror_data` — some montages were physically
   mounted mirrored).
8. `data2nii` writes the NIfTI: one 3D volume per timepoint via
   `3dUndump -master biopac16ch_template.nii` (placing each channel's
   value at its montage coordinate from `biopac16ch_template.txt`), then
   concatenated into a 4D dataset with `3dTcat -tr 0.51`. TR is fixed at
   0.51s for the whole dataset.

Outputs per subject/tag: `<subj>.<tag>.mem.fnirs.nii` (+ `.gonogo.fnirs.nii`
if that segment was found), `<subj>.<tag>.mem.timing.txt` (human-readable
marker table, always in absolute clock time — kept as a diagnostic
reference), `<subj>.<tag>.mem.onsetdur.txt` (FSL onset/duration/index
format, feeds `03_firstlevel_glm`).

`-skip_nii` regenerates only the timing/onsetdur files (steps 1–6), skipping
the slow per-timepoint `.nii` reconstruction (step 8, dominated by one
`3dUndump` subprocess call per timepoint — the reason this script is slow
per the original author's note). Use it when only the timing logic changed
and the underlying `.nii` doesn't need to be touched.

Template/mirror-list paths are resolved relative to the script's own file
location (`os.path.dirname` chain up to the `karim_fnirs/` root), so it
runs unchanged regardless of which mount alias (`/mnt/highlands`,
`/mnt/bellevue`, `/mnt/speyside`, ...) the repo happens to be checked out
under.

**`wrapper_export2nii.sh <data_folder>`** — finds every subject directory
two levels under `<data_folder>` and calls `export2nii.py` for each
(currently configured for the `hbt` tag only; edit the `-tag` line to run
`hbo`/`hbr` instead — both need a separate pass).

**`debug_markers.py`** — standalone inspector for one subject's raw marker
file; useful when a subject fails the marker-set validation in step 6.

---

## 02_preprocessing

**`preproc.sh <data_directory>`** — finds every `*.nii` under
`<data_directory>` and, for each:

1. `3dTstat` → `mean_<file>.nii` (per-voxel temporal mean).
2. `3dcalc`, masked to `fNIRS_template_mask.nii` (16-channel montage
   mask): `scl_<file>.nii = mask * min(200, signal/mean*100)` — rescales
   so the run mean is 100, i.e. converts to **percent signal change**,
   with a ceiling at 200 to avoid blow-up from near-zero baseline voxels.
3. `3dBandpass -despike -norm -mask fNIRS_template_mask.nii -band 0.08 0.12`
   → `prp_<file>.nii` — despikes, band-pass filters to 0.08–0.12 Hz, and
   L2-normalizes each channel's timeseries.

`prp_*.nii` is the file every downstream stage (03+) actually analyzes.

---

## 03_firstlevel_glm

Per-subject GLM via AFNI `3dDeconvolve`. **Two model variants**, run in
parallel on the same underlying data — pick one, or run both:

| | Amplitude-modulated (original) | na (non-amplitude-modulated) |
|---|---|---|
| Regressors | 2: `corsi_easy` (blocks 2–3), `corsi_hard` (blocks 4–5), each amplitude-modulated | 4: `b1`–`b4`, one per real block, static amplitude 1 |
| Amplitude source | Corsi performance score (`data_aux/corsi.csv`), block-number fallback if subject missing from CSV, amplitude-1 fallback if CSV score is exactly 0 | none — fixed amplitude 1 |
| Contrasts | `hard-easy` = -corsi_easy+corsi_hard; `corsi` = mean(easy,hard) | `na-hard-easy` = -0.5(b1+b2)+0.5(b3+b4); `na-corsi` = mean(b1..b4) |
| GLM script | `first_level_glm.sh` / `_wrapper.sh` | `first_level_glm_na.sh` / `_na_wrapper.sh` |
| Block prep script | `remodulateBlocks.sh` | `splitTimingIntoStaticBlocks.sh` |
| Output | `stats.<subj>.<type>.mem.nii.gz` | `stats.<subj>.<type>.mem.na.nii.gz` |

Both use `-polort A` (AFNI auto-selects detrending order per subject as
`1+int(run_duration/150)`, rather than a fixed order — a fixed low-enough
order for these ~500–650s runs still needs to clear individual block
durations of up to ~130s without eating into the task response), `-GOFORIT
6 -noFDR`, `-jobs 16`, `dmBLOCK` basis, `-local_times` (onsets relative to
each run's own start — see the 01_dataConversion timing-offset note above),
and `-tout` for per-regressor t-stats.

### Preparation steps (in order)

```
./onsetDur2afniDecon.sh <data dir>                          # onsetdur.txt -> <subj>.<type>.mem.timing.1D (via AFNI's timing_tool.py)
```
Then **either**:
```
./remodulateBlocks.sh <data dir> <data_aux>/corsi.csv        # timing.1D -> model.b2-3.1D / model.b4-5.1D (amp-modulated model)
```
**or**:
```
./splitTimingIntoStaticBlocks.sh <data dir>                  # timing.1D -> na.block1-4.1D (na model)
```

`timing.1D` format is AFNI's `dmBLOCK` AM1 convention:
`onset*amplitude:duration`, one entry per block, space-separated, one line
per subject/type. Block 1 (practice, not scored) is always dropped by both
`remodulateBlocks.sh` and `splitTimingIntoStaticBlocks.sh`; the remaining
4 real blocks are what get modeled.

`splitTimingIntoStaticBlocks.sh` additionally clips each block's modeled
duration to whatever portion of it actually falls within the subject's
scan (`min(onset+duration, scan_duration) - onset`, via `3dinfo -nt -tr`
on the matching `prp_*.nii`), and drops (no file written, causing
`3dDeconvolve` to fail cleanly on that subject rather than silently fit an
unstable regressor) any block where under 10% of its intended duration was
actually captured. Post-2026-08-17-fix this triggers for only 2 of 439
subjects — before the timing fix it looked like it affected nearly all of
them, which was actually the bug, not real truncation.

### Running the GLM

```
./first_level_glm_wrapper.sh <data dir> <mask.nii>       # amp-modulated model — finds prp_*.mem.fnirs.nii, calls first_level_glm.sh per file
./first_level_glm_na_wrapper.sh <data dir> <mask.nii>     # na model
```

Bucket sub-brick layout (`3dinfo -label stats.*.nii.gz`), 0-indexed —
needed for `-stim_times`/`3dROIstats` sub-brick selectors downstream:

- Amp-modulated: `0`=Full_Fstat, `1-2`=corsi_easy coef/tstat,
  `3-4`=corsi_hard coef/tstat, `5-6`=hard-easy_GLT coef/tstat,
  `7-8`=corsi_GLT coef/tstat. **Volume 5 = hard-easy contrast, volume 7 =
  corsi contrast** — these indices are what `04_secondlevel_group`'s dump
  scripts select.
- na: `0`=Full_Fstat, `1-8`=b1–b4 coef/tstat pairs,
  `9-10`=na-hard-easy_GLT coef/tstat, `11-12`=na-corsi_GLT coef/tstat.
  **Volume 9 = na-hard-easy, volume 11 = na-corsi** — selected by
  `04p5_secondlevel_group_na`'s dump scripts.

`splitTimingIntoBlocks.sh` is an older, simpler block-splitter (keeps
block-number as amplitude, no corsi-score lookup) — superseded by
`remodulateBlocks.sh` for the active model, kept for reference.

---

## 04_secondlevel_group — group analysis, amplitude-modulated model

Takes the per-subject `stats.*.nii.gz` buckets from `03_firstlevel_glm`
and produces group-comparison statistics, channel by channel.

### 1. Subject lists

```
cd fix_sublists && ./update_sublists.sh
```
Writes `{con,lrv,pcon}.{hbo,hbr}.sublist.txt` — one path per subject with a
completed `stats.*.nii.gz` (`find ... -name 'stat*<type>*.nii.gz'`; must be
filtered to exclude `*.na.nii.gz` if both models have been run in the same
tree). Copy the `hbo` lists into `dumpFirstLevel_betaCoefs/` for the next
step (that's where the group-comparison scripts actually read them from).

### 2. Dump per-channel beta coefficients to CSV

```
cd dumpFirstLevel_betaCoefs
./dumpBetaCoefs_allChannels.sh <group1 sublist> <group1 label> <group2 sublist> <group2 label> <volume idx> <output.csv>
```

For each of the 16 channels: `grab_betas.sh` builds a single-channel binary
mask from `templates/biopac16ch_template_mask.nii` (`3dcalc ... within(a,ch,ch)`),
runs `3dROIstats` on both groups' files at the given sub-brick index
(5=hard-easy, 7=corsi), and labels rows by group. `add_covariates.py`
left-joins in per-subject clinical covariates from a whitespace-delimited
`.1D` file (`X` column = subject ID; any other columns are covariates,
picked up automatically). `merge_channel_dumps.py` pivots the 16
per-channel temp files into one CSV (`channel1`...`channel16` columns).
`relabel_groups.sh` renames the internal group tags (`lrv`→`PSD+AGG`,
`pcon`→`PSD-AGG`, `con`→`HC`) via `sed`.

Covariate file: `covariate_file` in `dumpBetaCoefs_allChannels.sh` — build
it with `prep_select_cov.py` (merges `data_aux/covars.1D`'s
sans/saps/wais_matrix/dose with a binarized substance-use flag pulled from
the Mindata xlsx) or `prep_select_cov_local.py` (same logic, portable
paths, writes a local copy rather than overwriting the shared
`data_aux/select.cov.1D` — use this if the canonical file is stale/missing
columns on your mount).

### 3. Group-level statistics

```
python3 analyze_psd_agg_vs_psd_nagg.py --input <csv> --output-dir <dir> --group-pos <label> --group-neg <label> --covariates <col ...> [--interactions] [--alpha 0.05]
```

Per channel: Mann-Whitney U (unadjusted) and an ANCOVA
(`channel ~ C(group) + covariates`, OLS, Frisch-Waugh residualization for
the permutation tests below). Multiple-comparison correction, all applied
across the 16 channels: FDR-Benjamini/Hochberg, Simes-Hochberg, Hommel,
permutation max-statistic (5,000 perms — strongest single-stat FWER
control), and cluster-based permutation (5,000 perms, linear Ch1–Ch16
adjacency, cluster-forming threshold = uncorrected p<0.05). `--interactions`
additionally fits `channel ~ C(group) * covariates` and FDR-corrects each
group×covariate interaction term separately. Writes a results CSV, an
interactions CSV, and 3 figures (group means, correction-method heatmap,
volcano plot) per comparison.

`run_all_analyses.sh` (or the dated `run_all_analyses_<date>.sh` variant —
copy and edit the covariate lists per run) drives all 6 standard
comparisons: {PSD+AGG, PSD-AGG} × {HC} and PSD+AGG × PSD-AGG, on both
`hard-easy` and `corsi`. HC comparisons drop `dose` (all HC have dose=0,
perfectly collinear with group).

### 4. Report

```
python3 generate_report.py --base-dir . --output-dir . [--mode channels|regions]
```
Reads all 6 `results_*/*.csv` (+ `*_interactions.csv`), and writes
`report.md`/`report.html` — overview table, auto-generated narrative per
comparison, full per-channel tables, cluster summaries, interaction
tables. Edit the `COMPARISONS` list at the top of the script if result
directory names change.

### Region-pooled variant

`pool_to_regions.py` averages the 16 channels into 4 a-priori regions
(`left_dlPFC`=Ch1-4, `left_mPFC`=Ch5-8, `right_mPFC`=Ch9-12,
`right_dlPFC`=Ch13-16) on an existing `*_allchannels.csv`, producing a
`*_allregions.csv`. `run_all_region_analyses.sh` runs the same
`analyze_psd_agg_vs_psd_nagg.py` on those (with `--no-cluster`, since
cluster-permutation is redundant at k=4). `generate_report.py --mode
regions` reports on those instead.

### One-off utility

`augment_allchannels_csvs.py` — adds `wais_matrix`/`sud` columns to
already-dumped `*_allchannels.csv` files without re-running the (expensive,
one-`3dROIstats`-call-per-channel-per-subject) dump pipeline, by joining on
subject ID. Historical patch tool, not part of the normal per-run pipeline.

`fix_sublists/filter_sublists_cov.sh` / `find_missing_subjects.sh` —
utilities for cross-checking subject lists against the covariate file.

---

## 04p5_secondlevel_group_na — group analysis, na model

Mirrors `04_secondlevel_group` exactly, pointed at the na model's
`stats.*.mem.na.nii.gz` buckets (sub-brick 9=na-hard-easy, 11=na-corsi
instead of 5/7) and its own `.venv`-independent covariate/mask path fixes.

```
./update_sublists_na.sh                                    # con/lrv/pcon.{hbo,hbr}.na.sublist.txt (excludes stats.*.na.nii.gz's non-na counterpart automatically, since it only globs *na.nii.gz)
cd dumpFirstLevel_betaCoefs
./dumpBetaCoefs_allChannels.sh ... 9|11 <output.csv>        # same grab_betas.sh/add_covariates.py/merge_channel_dumps.py/relabel_groups.sh chain as 04
./run_all_analyses_na.sh                                    # same analyze_psd_agg_vs_psd_nagg.py, na model's 6 comparisons
```

A subject with no `stats.*.na.nii.gz` (na model 3dDeconvolve failed or
wasn't run) is automatically absent from the sublists — no special
exclusion handling needed, `find` just won't match it.

`prep_select_cov_na.py` — same portable covariate-file builder as
`04_secondlevel_group/dumpFirstLevel_betaCoefs/prep_select_cov_local.py`.

---

## 05_make_figure

Turns a group-analysis results CSV into a rendered brain figure. Uses its
**own** 16-channel template (`fNIRS_template.nii`/`.txt`/`_mask.nii`) — a
different voxel-coordinate layout than `templates/biopac16ch_template.*`
used by stages 01–04p5; don't mix them up.

```
./run_auto_suma_figures.sh --out-dir <out> <results.csv> [<results2.csv> ...]
```

For each input CSV: auto-detects a statistic column (`ancova_t`, `t`,
`cohens_d`, an interaction-t column, ...) and a significance column
(`ancova_p`, `ancova_p_fdr`, a boolean `*_sig` column, ...) unless
overridden with `--stat-col`/`--significance-p-col`/`--significance-bool-col`/`--alpha`,
writes each channel's (thresholded) statistic into its montage voxel
(`fNIRS_template.txt`) to build a NIfTI overlay (`auto_suma_outputs/nii/`),
and an audit CSV recording exactly which columns/values were used per
channel (`auto_suma_outputs/audit/`). `--no-threshold` skips significance
masking entirely.

Add `--suma-spec <spec>` to also project the overlay onto a cortical
surface with `3dVol2Surf` and render a screenshot via
`drive_suma_render.sh` (needs `afni`, `suma`, `DriveSuma`, `plugout_drive`,
`Xvfb`, `3dVol2Surf`, and a SUMA spec — this machine has `MNI_N27` under
`~/.afni/data/MNI_N27`). Override camera angle via
`SUMA_SURFACE_LABELS`/`SUMA_VIEW_KEYS` env vars (see
`auto_figures/README.md` for the exact rotation-key syntax).

`render_figure.py` is an older nilearn/matplotlib-based renderer (no X
server / SUMA needed), kept as a fallback for a different visual style;
invoke directly with `python3 auto_figures/render_figure.py`.

`toTemplate.ipynb` — notebook used to derive/inspect the montage-to-voxel
coordinate mapping in `fNIRS_template.txt`.

---

## templates/

`biopac16ch_template.nii`/`.txt` — the 16-channel biopac montage: 16
single-voxel "channels" at fixed grid coordinates, used as the `-master`
for `3dUndump` in `01_dataConversion` and as the channel-index reference
for `03_firstlevel_glm`/`04_secondlevel_group`/`04p5_secondlevel_group_na`'s
`biopac16ch_template_mask.nii` (channel-labeled mask, values 1–16, one
voxel-pair per channel — used both as the whole-brain analysis mask in
`3dDeconvolve` and, isolated per channel via `3dcalc within(a,N,N)`, for
`3dROIstats` beta extraction). `biopac16ch_template_compact*.nii` are
denser/alternate-layout variants; `biopac16ch_template_mir.txt` is the
coordinate table for mirrored montages. `brain.nii` is a full-resolution
anatomical reference, used by `05_make_figure` for rendering context.

---

## Project log

`notes.txt` — running, dated log of decisions, data-quality issues found,
and their resolutions; read before assuming any given result is final.
Dated `analysis_<date>/` folders (inside `dumpFirstLevel_betaCoefs/` and at
the repo root) hold point-in-time snapshots of group-analysis output —
treat later dates as superseding earlier ones for the same comparisons
unless a specific reason is given not to (e.g. `notes.txt` 2026-08-17: all
analyses before that date used a since-fixed, systematically misaligned
first-level timing).
