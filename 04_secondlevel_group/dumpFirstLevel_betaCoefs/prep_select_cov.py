#!/usr/bin/env python3
"""
prep_select_cov.py
==================
Build the master covariate file consumed by the channel-dump pipeline
(`dumpBetaCoefs_allChannels.sh` → `add_covariates.py`).

Sources
-------
1. /mnt/speyside/karim_fnirs/data_aux/covars.1D
     columns: X sans saps wais_matrix wais_info dose
2. /mnt/speyside/karim_fnirs/data_aux/Mindata1 (kopia).xlsx
     column "substance abuse " (0-6, NaN allowed) → binarized to 0/1

Output
------
/mnt/speyside/karim_fnirs/data_aux/select.cov.1D
    Whitespace-delimited table with columns:
        X  sans  saps  wais_matrix  sud  dose

Missing SUD values are imputed to 0 (population mode; also the safe default
for HCs which are typically screened to exclude substance abuse). The imputed
subject IDs are printed to stderr.

Run
---
    python prep_select_cov.py [--dry-run]
"""

from __future__ import annotations

import argparse
import shutil
import sys
from datetime import datetime
from pathlib import Path

import pandas as pd


DATA_AUX  = Path("/mnt/speyside/karim_fnirs/data_aux")
COVARS_1D = DATA_AUX / "covars.1D"
XLSX      = DATA_AUX / "Mindata1 (kopia).xlsx"
OUT_1D    = DATA_AUX / "select.cov.1D"

# Final column order (must match the runner's --covariates flags)
COLUMNS = ["X", "sans", "saps", "wais_matrix", "sud", "dose"]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dry-run", action="store_true",
                    help="Show what would be written but don't overwrite select.cov.1D")
    args = ap.parse_args()

    # ── 1. base covariates ────────────────────────────────────────────────────
    covars = pd.read_csv(COVARS_1D, sep=r"\s+")
    needed = {"X", "sans", "saps", "wais_matrix", "dose"}
    missing = needed - set(covars.columns)
    if missing:
        sys.exit(f"[ERROR] {COVARS_1D} is missing columns: {sorted(missing)}")
    covars["X"] = covars["X"].astype(str).str.strip()

    # ── 2. SUD from xlsx ──────────────────────────────────────────────────────
    xlsx = pd.read_excel(XLSX)
    sa_col = "substance abuse "
    if sa_col not in xlsx.columns:
        candidates = [c for c in xlsx.columns if "substance" in c.lower()]
        sys.exit(f"[ERROR] Column '{sa_col}' not found in {XLSX}. "
                 f"Candidates: {candidates}")
    sud = xlsx[["Unnamed: 0", sa_col]].copy()
    sud.columns = ["X", "sud_raw"]
    sud["X"] = sud["X"].astype(str).str.strip()

    # Binarize: 0 → 0, any of 1-6 → 1, NaN → NaN (handled later)
    sud["sud"] = sud["sud_raw"].apply(
        lambda v: 0 if pd.isna(v) is False and v == 0
                  else (1 if pd.notna(v) and v >= 1 else pd.NA)
    )
    print("[INFO] SUD raw distribution (Mindata xlsx):", file=sys.stderr)
    print(sud["sud_raw"].value_counts(dropna=False).sort_index().to_string(), file=sys.stderr)
    print("[INFO] SUD binarized counts (0 / 1 / NaN):", file=sys.stderr)
    print(sud["sud"].value_counts(dropna=False).sort_index().to_string(), file=sys.stderr)

    # ── 3. merge ──────────────────────────────────────────────────────────────
    merged = covars.merge(sud[["X", "sud"]], on="X", how="left")

    # ── 4. impute missing SUD with 0 (mode; safe default for HCs) ─────────────
    missing_sud = merged.loc[merged["sud"].isna(), "X"].tolist()
    if missing_sud:
        print(f"\n[WARN] {len(missing_sud)} subject(s) have no SUD entry in "
              f"{XLSX.name}; imputing SUD=0 (population mode):", file=sys.stderr)
        for sid in missing_sud:
            print(f"         {sid}", file=sys.stderr)
        merged["sud"] = merged["sud"].fillna(0).infer_objects(copy=False)

    # Cast SUD to integer for clean output
    merged["sud"] = merged["sud"].astype(int)

    # ── 5. select & order columns ─────────────────────────────────────────────
    out = merged[COLUMNS].copy()
    print(f"\n[INFO] Final table: {len(out)} subjects × {len(COLUMNS)} columns")
    print(out.head().to_string(index=False))

    # ── 6. write ──────────────────────────────────────────────────────────────
    if args.dry_run:
        print("\n[DRY-RUN] Not writing. Re-run without --dry-run to overwrite "
              f"{OUT_1D}", file=sys.stderr)
        return 0

    # Backup existing
    if OUT_1D.exists():
        stamp = datetime.now().strftime("%Y%m%d-%H%M%S")
        backup = OUT_1D.with_suffix(f".1D.bak-{stamp}")
        shutil.copy2(OUT_1D, backup)
        print(f"\n[INFO] Backed up existing {OUT_1D.name} → {backup.name}", file=sys.stderr)

    # Whitespace-delimited, no quoting, header on
    out.to_csv(OUT_1D, sep=" ", index=False, header=True)
    print(f"[✓] wrote {OUT_1D}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
