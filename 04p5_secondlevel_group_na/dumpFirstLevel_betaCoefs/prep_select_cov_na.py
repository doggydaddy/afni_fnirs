#!/usr/bin/env python3
"""
prep_select_cov_na.py
======================
Local variant of ../../04_secondlevel_group/dumpFirstLevel_betaCoefs/prep_select_cov.py.

The canonical select.cov.1D on this mount (data_aux/select.cov.1D, dated
2024-12-21) predates the wais_matrix/sud columns and only has
X sans saps dose - it was never regenerated here after prep_select_cov.py
started merging in wais_matrix (from covars.1D) and sud (from the Mindata
xlsx). Rather than overwrite the shared canonical file from this analysis,
build a local copy (select.cov.na.1D) for the na group-level comparisons.

Sources (local paths, both present under data_aux/):
    covars.1D                 columns: X sans saps wais_matrix wais_info dose
    Mindata1 (kopia).xlsx     column "substance abuse " (0-6, NaN allowed)

Output: select.cov.na.1D, columns X sans saps wais_matrix sud dose
"""

import sys
from pathlib import Path

import pandas as pd

DATA_AUX  = Path("/mnt/highlands/karim_fnirs/data_aux")
COVARS_1D = DATA_AUX / "covars.1D"
XLSX      = DATA_AUX / "Mindata1 (kopia).xlsx"
OUT_1D    = Path(__file__).parent / "select.cov.na.1D"

COLUMNS = ["X", "sans", "saps", "wais_matrix", "sud", "dose"]


def main() -> int:
    covars = pd.read_csv(COVARS_1D, sep=r"\s+")
    needed = {"X", "sans", "saps", "wais_matrix", "dose"}
    missing = needed - set(covars.columns)
    if missing:
        sys.exit(f"[ERROR] {COVARS_1D} is missing columns: {sorted(missing)}")
    covars["X"] = covars["X"].astype(str).str.strip()

    xlsx = pd.read_excel(XLSX)
    sa_col = "substance abuse "
    if sa_col not in xlsx.columns:
        candidates = [c for c in xlsx.columns if "substance" in c.lower()]
        sys.exit(f"[ERROR] Column '{sa_col}' not found in {XLSX}. Candidates: {candidates}")
    sud = xlsx[["Unnamed: 0", sa_col]].copy()
    sud.columns = ["X", "sud_raw"]
    sud["X"] = sud["X"].astype(str).str.strip()

    sud["sud"] = sud["sud_raw"].apply(
        lambda v: 0 if pd.isna(v) is False and v == 0
                  else (1 if pd.notna(v) and v >= 1 else pd.NA)
    )
    print("[INFO] SUD raw distribution (Mindata xlsx):", file=sys.stderr)
    print(sud["sud_raw"].value_counts(dropna=False).sort_index().to_string(), file=sys.stderr)

    merged = covars.merge(sud[["X", "sud"]], on="X", how="left")

    missing_sud = merged.loc[merged["sud"].isna(), "X"].tolist()
    if missing_sud:
        print(f"\n[WARN] {len(missing_sud)} subject(s) have no SUD entry in "
              f"{XLSX.name}; imputing SUD=0 (population mode):", file=sys.stderr)
        for sid in missing_sud:
            print(f"         {sid}", file=sys.stderr)
        merged["sud"] = merged["sud"].fillna(0).infer_objects(copy=False)

    merged["sud"] = merged["sud"].astype(int)

    out = merged[COLUMNS].copy()
    print(f"\n[INFO] Final table: {len(out)} subjects x {len(COLUMNS)} columns")
    print(out.head().to_string(index=False))

    out.to_csv(OUT_1D, sep=" ", index=False, header=True)
    print(f"[OK] wrote {OUT_1D}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
