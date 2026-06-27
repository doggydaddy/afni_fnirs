#!/usr/bin/python
"""
add_covariates.py
=================
Merge a per-subject covariate table (.1D, whitespace-delimited) onto a beta-values
channel dump produced by `grab_betas.sh`.

Generalized: the output column names are derived from the covariate file header
(the first column, conventionally "X", is renamed to "id"). Adding new covariates
no longer requires editing this script — just add a column to the .1D file.

Input table format (tab-separated, header line skipped):
    path   volume   value   group

Covariate file format (whitespace-separated, with header):
    X   sans   saps   wais_matrix   sud   dose
    0001P 14 53 13 1 25
    ...

Output CSV columns:
    path, volume, value, group, id, <all covariate columns except X>
"""

import csv
import sys
import pandas as pd
import argparse

desc = "add covariates to a beta values channel dump by grab_betas.sh"
parser = argparse.ArgumentParser(prog='add_covariates.py', description=desc)
parser.add_argument('-input',  required=True, help="input channel dump table, output from grab_betas.sh")
parser.add_argument('-covar',  required=True, help="covariates (.1D) file")
parser.add_argument('-output', required=True, help="output (csv) file")
args = parser.parse_args()

data_table       = args.input
covars_file      = args.covar
output_filename  = args.output

# ── Load covariates ────────────────────────────────────────────────────────────
# Tolerate one-or-more spaces between columns; AFNI's .1D files vary.
covars = pd.read_csv(covars_file, delim_whitespace=True)
if "X" not in covars.columns:
    sys.exit(f"[ERROR] Covariate file '{covars_file}' must have an 'X' column (subject IDs).")

# Strip whitespace from the subject ID column (just in case)
covars["X"] = covars["X"].astype(str).str.strip()

# Derive output header dynamically: input columns + covariate columns (X → id)
input_cols    = ["path", "volume", "value", "group"]
covar_cols    = list(covars.columns)                              # e.g. ["X","sans","saps","wais_matrix","sud","dose"]
covar_renamed = ["id" if c == "X" else c for c in covar_cols]     # e.g. ["id","sans","saps","wais_matrix","sud","dose"]
headers       = input_cols + covar_renamed

output_df = pd.DataFrame(columns=headers)
missing_ids = []

with open(data_table, mode="r", newline="", encoding="utf-8") as file:
    reader = csv.reader(file, delimiter="\t")
    next(reader)  # skip header

    for row in reader:
        if not row:
            continue
        full_path = row[0]
        bn        = full_path.split("/")[-1]
        subj_id   = bn.split(".")[1]

        matching_rows = covars[covars["X"] == subj_id]
        if matching_rows.size > 0:
            new_row = row + matching_rows.values.flatten().tolist()
            output_df.loc[len(output_df)] = new_row
        else:
            missing_ids.append(subj_id)

if missing_ids:
    print(f"[WARN] {len(missing_ids)} subject(s) not found in covariate file: "
          f"{sorted(set(missing_ids))}", file=sys.stderr)

output_df.to_csv(output_filename, index=False)
print(f"[OK] wrote {len(output_df)} rows × {len(headers)} cols → {output_filename}")