#!/usr/bin/env python3
"""Build a local covariate table for either second-level model.

The input directory must contain ``covars.1D`` and
``Mindata1 (kopia).xlsx``. The output is deliberately local and ignored by
Git because it contains subject-level clinical information.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd


COLUMNS = ["X", "sans", "saps", "wais_matrix", "sud", "dose"]


def main() -> int:
    script_dir = Path(__file__).resolve().parent
    repo_root = script_dir.parent

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--data-aux",
        type=Path,
        default=repo_root.parent / "data_aux",
        help="Directory containing covars.1D and the Mindata workbook",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=script_dir / "select.cov.local.1D",
        help="Destination .1D file",
    )
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    covars_path = args.data_aux / "covars.1D"
    workbook_path = args.data_aux / "Mindata1 (kopia).xlsx"
    covars = pd.read_csv(covars_path, sep=r"\s+")
    needed = {"X", "sans", "saps", "wais_matrix", "dose"}
    missing = needed - set(covars.columns)
    if missing:
        sys.exit(f"[ERROR] {covars_path} is missing columns: {sorted(missing)}")
    covars["X"] = covars["X"].astype(str).str.strip()

    workbook = pd.read_excel(workbook_path)
    substance_column = "substance abuse "
    if substance_column not in workbook.columns:
        candidates = [column for column in workbook.columns
                      if "substance" in column.lower()]
        sys.exit(
            f"[ERROR] column {substance_column!r} not found in {workbook_path}; "
            f"candidates: {candidates}"
        )

    sud = workbook[["Unnamed: 0", substance_column]].copy()
    sud.columns = ["X", "sud_raw"]
    sud["X"] = sud["X"].astype(str).str.strip()
    sud["sud"] = sud["sud_raw"].apply(
        lambda value: 0 if pd.notna(value) and value == 0
        else (1 if pd.notna(value) and value >= 1 else pd.NA)
    )

    print("[INFO] SUD raw distribution:", file=sys.stderr)
    print(sud["sud_raw"].value_counts(dropna=False).sort_index().to_string(),
          file=sys.stderr)

    merged = covars.merge(sud[["X", "sud"]], on="X", how="left")
    missing_ids = merged.loc[merged["sud"].isna(), "X"].tolist()
    if missing_ids:
        print(
            f"[WARN] imputing SUD=0 for {len(missing_ids)} subject(s) with no entry: "
            f"{missing_ids}",
            file=sys.stderr,
        )
        merged["sud"] = merged["sud"].fillna(0).infer_objects(copy=False)

    output = merged[COLUMNS].copy()
    output["sud"] = output["sud"].astype(int)
    print(f"[INFO] final table: {len(output)} subjects x {len(COLUMNS)} columns")

    if args.dry_run:
        print(f"[DRY-RUN] would write {args.output}")
        return 0

    args.output.parent.mkdir(parents=True, exist_ok=True)
    output.to_csv(args.output, sep=" ", index=False)
    print(f"[OK] wrote {args.output}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
