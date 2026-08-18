# Second-level group analysis

This directory contains the complete second-level workflow for both first-level
models. The amplitude-modulated and non-amplitude-modulated (`na`) analyses use
the same channel extraction, covariate, statistics, correction, and reporting
routines; only bucket names, GLT indices, input names, and covariate profiles
differ.

| Model | Subject bucket | Coefficient volumes | Comparison entry point |
|---|---|---|---|
| amplitude | `stats.*.<type>.mem.nii.gz` | 5 hard-easy, 7 corsi | `run_all_analyses.sh` |
| na | `stats.*.<type>.mem.na.nii.gz` | 9 hard-easy, 11 corsi | `run_all_analyses_na.sh` |

Run the following commands from `04_secondlevel_group/` unless noted.

## 1. Generate subject lists

```text
./update_sublists.sh
./update_sublists_na.sh
```

These preserve the former call syntax and write the amplitude and na lists,
respectively, directly into this directory. `update_model_sublists.sh` is the
shared implementation. Set `FNIRS_DATA_ROOT` to override the default sibling
`../data` directory.

## 2. Build the shared local covariate table

```text
./.venv/bin/python prep_select_cov_local.py
# compatibility alias for the na workflow:
./prep_select_cov_na.py
```

Both models consume `select.cov.local.1D`. It is ignored because it contains
subject-level clinical information. `--data-aux` and `--output` can override
the default paths.

## 3. Dump all 16 channels

The positional call syntax is unchanged:

```text
./dumpBetaCoefs_allChannels.sh GROUP1_LIST GROUP1_LABEL GROUP2_LIST GROUP2_LABEL VOLUME OUTPUT.csv
```

Use volume 5 or 7 for the amplitude model and 9 or 11 for the na model. The
same covariate table and helper chain are used for both. `COVARIATE_FILE`,
`TEMPLATE_MASK`, and `GROUP_ANALYSIS_PYTHON` override the defaults.

## 4. Run the six standard comparisons

Invoke the appropriate wrapper from the directory containing the six
`*_allchannels.csv` inputs, normally a private analysis directory outside the
repository:

```text
/path/to/repo/04_secondlevel_group/run_all_analyses.sh [--dry-run]
/path/to/repo/04_secondlevel_group/run_all_analyses_na.sh [--dry-run]
```

Both wrappers call `run_group_analyses.sh`. The legacy amplitude wrapper keeps
SANS/SAPS for the patient-only comparison. The dated amplitude wrapper and na
wrapper use the established symptom-excluded profile.

## Reports and region pooling

- `pool_to_regions.py` pools channels into four predefined regions.
- `run_all_region_analyses.sh` runs the region-level comparisons.
- `generate_report.py` creates Markdown and HTML reports from result folders.
- `augment_allchannels_csvs.py` is a historical utility for adding covariates
  to existing channel dumps without rerunning AFNI.

Subject lists, covariate tables, channel dumps, results, reports, logs, and
dated analysis directories are ignored because they may contain sensitive
information.
