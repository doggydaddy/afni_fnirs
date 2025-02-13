#!/usr/bin/python

import csv 
import pandas as pd
import numpy as np
import argparse

desc="add covariates to a beta values channel dump by grab_betas.sh"
parser = argparse.ArgumentParser(prog='add_covariates.py', description=desc)
parser.add_argument('-input', help="input channel dump table, output from grab_betas.sh")
parser.add_argument('-covar', help="covariates (.1D) files")
parser.add_argument('-output', help="output (csv) files")
args = parser.parse_args()


data_table = args.input
covars_file = args.covar
output_filename = args.output

covars = pd.read_csv(covars_file, delimiter=" ")

headers=["path", "volume", "value", "group", "id", "sans", "saps", "dose"]

output_df = pd.DataFrame(columns=headers)
with open(data_table, mode="r", newline="", encoding="utf-8") as file:
    reader = csv.reader(file, delimiter="\t")
    next(reader) 

    for row in reader:
        # Check if the row is not empty
        if row:
            full_path = row[0]  # Read the first column

            split_res_1 = full_path.split("/")
            bn = split_res_1[-1]

            split_res_2 = bn.split(".")
            subj_id = split_res_2[1]

            # Find rows where the name "X" matches subject ID
            matching_rows = covars[covars["X"] == subj_id]

            # Print the result
            if matching_rows.size > 0:
                new_row = row + matching_rows.values.flatten().tolist()
                output_df.loc[len(output_df)] = new_row
            else:
                print("row not found")

output_df.to_csv(output_filename, index=False)