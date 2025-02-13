#!/usr/bin/python
import pandas as pd
import glob
import os
from functools import reduce

# Example: files named TMP_ch1.csv, TMP_ch2.csv, ..., TMP_ch16.csv
file_pattern = "TMP_ch*.csv"

# Gather all filenames
csv_files = glob.glob(file_pattern)
csv_files.sort()

# Read the first CSV to identify the common columns
temp_df = pd.read_csv(csv_files[0])
common_columns = [col for col in temp_df.columns if col != 'value']

df_list = []
for file_name in csv_files:
    base = os.path.basename(file_name)
    channel_number = base.replace("TMP_ch", "").replace(".csv", "")
    
    df = pd.read_csv(file_name)
    df = df.rename(columns={'value': f'channel{channel_number}'})
    df_list.append(df)

merged_df = reduce(lambda left, right: pd.merge(left, right, on=common_columns), df_list)

# --- Reorder columns so that channel columns are last and in ascending channel order ---
all_cols = list(merged_df.columns)

# Identify channel columns
channel_cols = [c for c in all_cols if c.startswith('channel')]

# Identify the non-channel (common) columns
non_channel_cols = [c for c in all_cols if not c.startswith('channel')]

# Sort the channel columns by their numeric suffix (channel1, channel2, etc.)
# If the channel suffix can be multi-digit, parse it as integer:
channel_cols_sorted = sorted(channel_cols, key=lambda x: int(x.replace('channel', '')))

# Combine them: non-channel columns first, then sorted channel columns
final_column_order = non_channel_cols + channel_cols_sorted
merged_df = merged_df[final_column_order]

# Write out the merged dataset
merged_df.to_csv("merged.csv", index=False)