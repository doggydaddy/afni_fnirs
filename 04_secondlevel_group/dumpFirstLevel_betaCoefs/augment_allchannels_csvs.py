#!/usr/bin/env python3
"""
augment_allchannels_csvs.py
===========================
One-shot utility to add `wais_matrix` and `sud` covariate columns to existing
`*_allchannels.csv` files **without** re-running the AFNI dump pipeline.

Why
---
The dump pipeline (`dumpBetaCoefs_allChannels.sh`) is expensive because it
calls 3dROIstats once per (channel × subject). Beta coefficients in the
existing CSVs are byte-for-byte correct; the only thing missing is the two new
covariate columns. So we just join them in by subject ID.

Sources
-------
* wais_matrix:  /mnt/speyside/karim_fnirs/data_aux/covars.1D
* sud (raw):    /mnt/speyside/karim_fnirs/data_aux/Mindata1 (kopia).xlsx
                  column "substance abuse "  →  binarize 0→0, 1-6→1
                  missing IDs → impute 0 (population mode; HCs default)

Behaviour
---------
For each <input_dir>/<base>_allchannels.csv:
    read   →  insert wais_matrix and sud after "dose"  →  write <out_dir>/<base>_allchannels.csv

Usage
-----
    python augment_allchannels_csvs.py \
        --input-dir  analysis_2026-06-26 \
        --output-dir .
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd


COVARS_1D = Path("/mnt/speyside/karim_fnirs/data_aux/covars.1D")
XLSX      = Path("/mnt/speyside/karim_fnirs/data_aux/Mindata1 (kopia).xlsx")

# CSV basenames to augment (the .csv suffix is added automatically)
CSV_BASES = [
    "lrv-con_v5_allchannels",
    "lrv-con_v7_allchannels",
    "lrv-pcon_v5_allchannels",
    "lrv-pcon_v7_allchannels",
    "pcon-con_v5_allchannels",
    "pcon-con_v7_allchannels",
]


def load_lookups() -> pd.DataFrame:
    """Return DataFrame indexed by subject id with columns wais_matrix, sud."""
    # ── wais_matrix from covars.1D ───────────────────────────────────────────
    covars = pd.read_csv(COVARS_1D, sep=r"\s+")
    if "wais_matrix" not in covars.columns:
        sys.exit(f"[ERROR] 'wais_matrix' not in {COVARS_1D}; columns: {list(covars.columns)}")
    covars["X"] = covars["X"].astype(str).str.strip()
    wais = covars[["X", "wais_matrix"]].rename(columns={"X": "id"})

    # ── sud from xlsx (binarize) ─────────────────────────────────────────────
    xlsx   = pd.read_excel(XLSX)
    sa_col = "substance abuse "
    if sa_col not in xlsx.columns:
        candidates = [c for c in xlsx.columns if "substance" in c.lower()]
        sys.exit(f"[ERROR] '{sa_col}' not in {XLSX}; candidates: {candidates}")
    sud = xlsx[["Unnamed: 0", sa_col]].copy()
    sud.columns = ["id", "sud_raw"]
    sud["id"]   = sud["id"].astype(str).str.strip()
    # Binarize: 0→0, 1-6→1, NaN→NaN
    sud["sud"] = sud["sud_raw"].where(sud["sud_raw"].isna(),
                                       (sud["sud_raw"] >= 1).astype(int))

    # ── merge into a single lookup ───────────────────────────────────────────
    lookup = wais.merge(sud[["id", "sud"]], on="id", how="left")
    return lookup


def augment_one(in_path: Path, out_path: Path, lookup: pd.DataFrame) -> None:
    """Read in_path, add wais_matrix + sud columns after 'dose', write out_path."""
    df = pd.read_csv(in_path)
    if "id" not in df.columns:
        sys.exit(f"[ERROR] {in_path} has no 'id' column")

    # Drop any pre-existing wais_matrix / sud columns to keep the merge idempotent
    df = df.drop(columns=[c for c in ("wais_matrix", "sud") if c in df.columns],
                 errors="ignore")

    merged = df.merge(lookup, on="id", how="left")

    # Report (and impute) missing values per file
    miss_wais = merged.loc[merged["wais_matrix"].isna(), "id"].tolist()
    miss_sud  = merged.loc[merged["sud"].isna(), "id"].tolist()
    if miss_wais:
        print(f"  [WARN] {len(miss_wais)} subj missing wais_matrix "
              f"(left as NaN — will be dropped by ANCOVA): {miss_wais}", file=sys.stderr)
    if miss_sud:
        print(f"  [WARN] {len(miss_sud)} subj missing sud → imputing 0 "
              f"(population mode; HCs default): {miss_sud}", file=sys.stderr)
        merged["sud"] = merged["sud"].fillna(0)

    merged["sud"] = merged["sud"].astype("Int64")  # nullable int (NaN-safe)

    # Reorder: insert wais_matrix, sud right after "dose"
    cols = list(merged.columns)
    for c in ("wais_matrix", "sud"):
        cols.remove(c)
    dose_idx = cols.index("dose") + 1
    new_order = cols[:dose_idx] + ["wais_matrix", "sud"] + cols[dose_idx:]
    merged = merged[new_order]

    out_path.parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(out_path, index=False)
    print(f"  [✓] {in_path.name}  →  {out_path}   "
          f"({len(merged)} rows × {len(merged.columns)} cols)")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input-dir",  default="analysis_2026-06-26",
                    help="Directory containing the original *_allchannels.csv files")
    ap.add_argument("--output-dir", default=".",
                    help="Where to write the augmented *_allchannels.csv files")
    args = ap.parse_args()

    in_dir  = Path(args.input_dir).resolve()
    out_dir = Path(args.output_dir).resolve()

    print(f"[INFO] Input  dir : {in_dir}")
    print(f"[INFO] Output dir : {out_dir}")
    print(f"[INFO] Loading lookups from {COVARS_1D.name} and {XLSX.name} ...")
    lookup = load_lookups()
    print(f"[INFO] Lookup table: {len(lookup)} subjects "
          f"(non-null sud: {lookup['sud'].notna().sum()}, "
          f"non-null wais_matrix: {lookup['wais_matrix'].notna().sum()})\n")

    for base in CSV_BASES:
        in_path  = in_dir  / f"{base}.csv"
        out_path = out_dir / f"{base}.csv"
        if not in_path.exists():
            print(f"  [SKIP] missing: {in_path}", file=sys.stderr)
            continue
        augment_one(in_path, out_path, lookup)

    print("\n[DONE]")
    return 0


if __name__ == "__main__":
    sys.exit(main())
