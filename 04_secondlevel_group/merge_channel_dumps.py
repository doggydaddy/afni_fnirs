#!/usr/bin/env python3
"""Merge one CSV per fNIRS channel into a single wide subject table."""

from __future__ import annotations

import argparse
import re
import sys
from functools import reduce
from pathlib import Path

import pandas as pd


CHANNEL_RE = re.compile(r"TMP_ch(\d+)\.csv$")


def channel_number(path: Path) -> int:
    match = CHANNEL_RE.search(path.name)
    if not match:
        raise ValueError(f"cannot extract channel number from {path}")
    return int(match.group(1))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", default=".",
                        help="Directory containing TMP_ch<N>.csv files")
    parser.add_argument("--output", default="merged.csv",
                        help="Merged CSV path (default: merged.csv)")
    args = parser.parse_args()

    input_dir = Path(args.input_dir)
    csv_files = sorted(input_dir.glob("TMP_ch*.csv"), key=channel_number)
    if not csv_files:
        sys.exit(f"[ERROR] no TMP_ch*.csv files found in {input_dir}")

    frames = []
    common_columns: list[str] | None = None
    for path in csv_files:
        frame = pd.read_csv(path)
        if "value" not in frame.columns:
            sys.exit(f"[ERROR] {path} has no 'value' column")
        current_common = [column for column in frame.columns if column != "value"]
        if common_columns is None:
            common_columns = current_common
        elif current_common != common_columns:
            sys.exit(f"[ERROR] metadata columns differ in {path}")
        frame = frame.rename(columns={"value": f"channel{channel_number(path)}"})
        frames.append(frame)

    assert common_columns is not None
    merged = reduce(
        lambda left, right: pd.merge(left, right, on=common_columns, validate="one_to_one"),
        frames,
    )
    channel_columns = sorted(
        (column for column in merged.columns if column.startswith("channel")),
        key=lambda column: int(column.removeprefix("channel")),
    )
    merged = merged[common_columns + channel_columns]

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(output, index=False)
    print(f"[OK] merged {len(csv_files)} channels into {output}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
