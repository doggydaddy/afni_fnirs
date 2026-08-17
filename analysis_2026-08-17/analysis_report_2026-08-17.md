# fNIRS Corsi Task — Analysis Report, 2026-08-17

**What this document is:** a record of a timing-alignment bug found in the
first-level GLM pipeline, how it was fixed, and the current end-to-end
analysis pipeline as a result. It accompanies `report.md`/`report.html`
in this folder, which holds the actual statistical results (12
comparisons: 2 first-level models x 3 group pairs x 2 contrasts).

---

## 1. Executive summary

While setting up a new non-amplitude-modulated (na) first-level model,
group-level comparisons came back null across the board. Chasing that led
to a genuine bug: **every subject's block-onset timing fed to
`3dDeconvolve` was offset from where those blocks actually occurred in the
exported fNIRS data**, by an amount unique to each subject (12–232s,
median 35.6s). This affected the **entire pipeline** — both the
already-established amplitude-modulated model and the new na model — for
every prior analysis run through it.

The bug was fixed at its source (`01_dataConversion/fnirsoft2nii/export2nii.py`),
all timing metadata was regenerated for all 457 subject-runs, and both
first-level models were rerun end-to-end along with their group-level
comparisons. Two unrelated pre-existing data-quality issues were also
found and fixed along the way (malformed markers in 2 subjects; all-zero
corsi performance scores in 5 subjects). Full sample recovered on both
models (85 HC / 77 PSD+AGG / 67 PSD-AGG, matching raw subject-folder
counts) — see §6 for what changed in the results.

---

## 2. The bug: block-onset timing offset

### 2.1 How it was found

The na model (4 separate per-block regressors, `-polort A`) was showing
apparent block-4 truncation in ~99% of subjects — `3dDeconvolve` either
hard-failed (18 subjects, both blocks 3 and 4 entirely "missing") or, more
insidiously, silently produced numerically unstable beta estimates from a
single subject whose block 4 had a tiny (0.8%) sliver of apparent data
(`1068KP`, beta = -201 on one channel, versus a normal range of
±0.0001–0.01). Investigating that outlier — checking whether it was a
data artifact — led to checking `01_dataConversion` at the user's request.

### 2.2 Root cause

`01_dataConversion/fnirsoft2nii/export2nii.py`:

- `marker2onsetdur()` wrote block onset times using the **raw recording's
  absolute clock** (e.g. block 5 onset = 510.705s for `1068KP`, i.e.
  510.705s after the fNIRS device started recording).
- `data2nii()` builds each subject's exported `.nii` starting at marker
  **-5** (task start) — so frame 0 of the `.nii` is *not* absolute time
  zero. It's whatever the raw clock read when marker -5 fired (96.56s for
  `1068KP`; ranges 12–232s across the dataset, median 35.6s).
- Nothing subtracted that offset anywhere downstream
  (`onsetDur2afniDecon.sh`, `splitTimingIntoBlocks.sh`,
  `remodulateBlocks.sh`, `splitTimingIntoStaticBlocks.sh` all just pass
  the onset values through).
- `first_level_glm.sh` / `first_level_glm_na.sh` call `3dDeconvolve` with
  `-local_times`, which — per AFNI's own documentation — interprets onset
  times as **relative to the start of the analyzed run** (i.e. frame 0 of
  the `.nii`), not the raw recording's absolute clock.

Net effect: every block, for every subject, was scheduled in the design
matrix at a position offset from where it actually happened in the data,
by that subject's own marker -5 time.

### 2.3 Proof

Recomputed, for all 439 then-completed na-model subjects, what fraction of
each block's *intended* duration actually fell within the recording — once
before and once after subtracting the marker -5 offset:

| | before correction | after correction |
|---|---|---|
| block 4 <10% captured | 80 subjects | **0** |
| block 4 <95% captured | 435 subjects | **2** |
| block 3 <95% captured | 66 subjects | **0** |

The "truncated block 4" finding all but disappeared once the clocks were
aligned. Separately, correlating the (buggy) na model's beta estimates
against the (buggy) old model's beta estimates on the same underlying data
showed real, structured correlation at several channels (e.g. na-corsi vs.
old corsi at Ch4: r=0.40, p<0.0001, n=219) — proof that `3dDeconvolve`
itself was working correctly throughout; the two models agreed with each
other on a real signal even while both were misaligned in absolute terms.

---

## 3. Secondary issues found while fixing this

None of these are timing-offset related — found incidentally while
regenerating data end-to-end.

### 3.1 Malformed marker sets (data conversion)

`marker2onsetdur()` assumed exactly one each of markers
-5, 101, 102, 103, 104, 105, 106, with no validation. Two subjects have
marker sets that don't satisfy this and previously crashed the batch
uncaught:

- **`0020P`** — a duplicate `101` marker (recording glitch). Recovered
  automatically via a second acquisition attempt present in the subject's
  raw files.
- **`0028P`** — markers 104 and 105 entirely missing (task administration
  cut short). Already flagged in `notes.txt` (2026-02-02 entry, "Missing
  Markers"). Confirmed this subject has no `prp_0028P.*.mem.fnirs.nii` at
  all — never converted for the corsi task, unrelated to this bug, and
  remains excluded.

Fixed by validating the marker set right after `trim_marker()` and
skipping (not crashing) any subject whose markers don't match, with a
clear log message naming what's missing/unexpected.

### 3.2 All-zero corsi.csv rows (amplitude-modulated model only)

Regenerating the amplitude-modulated model's `model.b2-3.1D`/`model.b4-5.1D`
surfaced 5 subjects (`0053P`, `1013KP`, `1023KP`, `1052KP`, `1072KP`) with a
literal all-zero row in `data_aux/corsi.csv` (`Block_2` through
`Click_tot` all 0). Amplitude-modulating a block by a literal 0 zeroes
that event out of the regressor entirely; with every block zeroed,
`3dDeconvolve` reported `Signal-only matrix condition: UNDEFINED ** VERY
BAD **` for all 10 runs (5 subjects x hbo/hbr). `-GOFORIT 6` let these
through silently — this would not have been caught without inspecting
individual subjects' matrix-condition warnings.

Per instruction: when a block's corsi.csv score is exactly 0,
`remodulateBlocks.sh` now falls back to a static amplitude of 1 for that
block (distinct from its existing "subject not found in csv" fallback,
which still uses block-number amplitude, left unchanged).

### 3.3 `-polort 6` → `-polort A`

Unrelated to the above, raised separately: `-polort 6` was a fixed
detrending order regardless of each subject's actual run length. AFNI's
own guidance (`-polort A`, `pnum = 1 + int(D/150)`) auto-selects per run
duration. At these run lengths (~500-660s), fixed polort 6 set the
drift-removal cutoff at ~140-165s of period — close enough to individual
block durations (up to ~130s) to risk eating into the task response; `A`
roughly doubles that cutoff. Applied to both first-level models.

### 3.4 Hardcoded `/mnt/speyside` paths

Several scripts hardcoded absolute paths to a `/mnt/speyside` mount that
doesn't exist on this machine (this dataset is mirrored across
`/mnt/highlands`, `/mnt/bellevue`, and `/mnt/speyside`, apparently
different machines' views of the same underlying storage).  Fixed in
`export2nii.py` (template/mirror-list paths, now resolved relative to the
script's own location) and in both `dumpBetaCoefs_allChannels.sh` scripts
(template mask path, now relative; covariate file, see below).

### 3.5 Stale covariate file

The canonical `data_aux/select.cov.1D` on this mount is dated 2024-12-21
and only has columns `X sans saps dose` — it predates `wais_matrix`/`sud`
being merged in by `prep_select_cov.py` (from `covars.1D` and the Mindata
xlsx), which a later analysis run (`analysis_2026-07-04`) clearly had
access to on whichever machine produced it. Rather than overwrite the
shared canonical file from this session, built local copies
(`select.cov.na.1D`, `select.cov.local.1D`) via a portable-path variant of
the same merge logic (`prep_select_cov_na.py`, `prep_select_cov_local.py`).

---

## 4. Files changed

| File | Change |
|---|---|
| `01_dataConversion/fnirsoft2nii/export2nii.py` | Fixed timing offset (`marker2onsetdur`); added `-skip_nii` flag; added marker-set validation; fixed hardcoded paths |
| `01_dataConversion/fnirsoft2nii/README.md` | Documented `-skip_nii` and the timing-offset fix |
| `03_firstlevel_glm/first_level_glm.sh` | `-polort 6` → `-polort A` |
| `03_firstlevel_glm/first_level_glm_na.sh` | New (na model GLM); `-polort A` |
| `03_firstlevel_glm/first_level_glm_na_wrapper.sh` | New (na model batch wrapper) |
| `03_firstlevel_glm/splitTimingIntoStaticBlocks.sh` | New (na model block prep); duration-clipping safety net (now rarely triggers post-fix) |
| `03_firstlevel_glm/remodulateBlocks.sh` | Zero-score fallback to amplitude 1 (§3.2) |
| `03_firstlevel_glm/README.md` | Documented na model scripts |
| `04_secondlevel_group/dumpFirstLevel_betaCoefs/dumpBetaCoefs_allChannels.sh` | Fixed hardcoded paths |
| `04_secondlevel_group/dumpFirstLevel_betaCoefs/add_covariates.py` | pandas 3.0 compatibility (`delim_whitespace` removed) |
| `04_secondlevel_group/dumpFirstLevel_betaCoefs/prep_select_cov_local.py` | New (local covariate file builder) |
| `04_secondlevel_group/dumpFirstLevel_betaCoefs/run_all_analyses_2026-08-17.sh` | New (this run's group-analysis driver, sans/saps excluded) |
| `04p5_secondlevel_group_na/` | New (entire na model group-analysis setup — sublists, beta dumps, group comparisons) |
| `analysis_2026-08-17/generate_combined_report.py` | New (this report's generator) |
| `notes.txt` | Full narrative log, dated entries |

Pre-fix timing metadata (`onsetdur.txt`, `timing.txt`, `timing.1D`) for all
subjects backed up to
`data_aux/backup_pre_offset_fix_2026-08-17/timing_metadata_backup.tar.gz`
before any regeneration.

---

## 5. Current end-to-end pipeline

Steps to reproduce either model's results from scratch, in order:

### 5.1 Data conversion (only needed once, already done)

```
cd 01_dataConversion/fnirsoft2nii
python3 export2nii.py -foldername <subject dir> -tag hbo -prefix <subj> [-skip_nii]
```
(`-skip_nii` regenerates only timing metadata — use this for any future
timing-logic fix; omit it for a genuinely new subject's first conversion.)

### 5.2 First-level GLM — amplitude-modulated model

```
cd 03_firstlevel_glm
./onsetDur2afniDecon.sh <data dir>                        # onsetdur.txt -> timing.1D
./remodulateBlocks.sh <data dir> <data_aux>/corsi.csv      # timing.1D -> model.b2-3/b4-5.1D
./first_level_glm_wrapper.sh <data dir> <mask>             # -> stats.<subj>.<type>.mem.nii.gz
```

### 5.3 First-level GLM — na model

```
cd 03_firstlevel_glm
./onsetDur2afniDecon.sh <data dir>                         # (shared with 5.2, if not already run)
./splitTimingIntoStaticBlocks.sh <data dir>                # timing.1D -> na.block1-4.1D
./first_level_glm_na_wrapper.sh <data dir> <mask>          # -> stats.<subj>.<type>.mem.na.nii.gz
```

### 5.4 Group-level: amplitude-modulated model

```
cd 04_secondlevel_group/fix_sublists
./update_sublists.sh                                       # (or build directly with find, excluding *.na.nii.gz)
cd ../dumpFirstLevel_betaCoefs
./dumpBetaCoefs_allChannels.sh <group1 sublist> <label1> <group2 sublist> <label2> <vol idx> <out.csv>
  # vol idx 5 = hard-easy_GLT_Coef, 7 = corsi_GLT_Coef
./run_all_analyses_2026-08-17.sh                            # 6 comparisons, sans/saps excluded
```

### 5.5 Group-level: na model

```
cd 04p5_secondlevel_group_na
./update_sublists_na.sh
cd dumpFirstLevel_betaCoefs
./dumpBetaCoefs_allChannels.sh <group1 sublist> <label1> <group2 sublist> <label2> <vol idx> <out.csv>
  # vol idx 9 = na-hard-easy_GLT_Coef, 11 = na-corsi_GLT_Coef
./run_all_analyses_na.sh                                    # 6 comparisons, sans/saps excluded
```

### 5.6 Combined report (this document's companion)

```
cd analysis_2026-08-17
python3 generate_combined_report.py                         # -> report.md, report.html
```

---

## 6. Impact on results

Both models now use correctly time-aligned data. Full sample recovered on
both (85 HC / 77 PSD+AGG / 67 PSD-AGG on hbo, matching raw subject-folder
counts — vs. 82/74/64 mid-debugging before the offset fix was found).

**na model** (see `report.md` §2, 6, 9, 11): PSD+AGG vs HC shows 3
significant clusters on hard-easy (Ch3, Ch6-7, Ch9-16); PSD-AGG vs HC
shows a significant cluster on hard-easy (Ch9-12); PSD+AGG vs PSD-AGG
shows a wide group x covariate interaction effect on corsi.

**Amplitude-modulated model** (see `report.md` §1, 3, 5, 7): much
narrower — only PSD+AGG vs PSD-AGG on corsi shows a significant cluster
(Ch15-16).

Both models are now correctly time-aligned, so this difference reflects
the two models' different sensitivity/power (na's 4 separate
static-amplitude regressors vs. this model's 2 merged, corsi-score
amplitude-modulated regressors), not one being right and the other wrong.
Not yet reconciled or further interpreted — see `report.md` for full
per-channel tables.

Prior to this fix, both models' group comparisons were **null across all
6 comparisons** — the timing bug was actively suppressing a real signal,
not just adding noise.

---

## 7. Caveats / not yet done

- Prior group-level analyses under `04_secondlevel_group/dumpFirstLevel_
  betaCoefs/analysis_2026-06-26`, `2026-06-27`, `2026-07-04` were built on
  the same misaligned timing and should be treated as superseded by this
  report, not as independent confirmation of anything.
- `gonogo` task data was not touched by this investigation — same
  `export2nii.py` machinery is used for it, and the same class of offset
  bug plausibly applies there too, but it hasn't been checked.
- The two models' differing results (§6) haven't been reconciled or
  interpreted beyond noting the discrepancy.
- Region-pooled analyses (`*_regions` variants elsewhere in the repo)
  were not rerun as part of this fix.
