#!/usr/bin/env python3
"""
pool_to_regions.py
==================
Pool the 16 channel columns of each `*_allchannels.csv` file into 4 anatomical
region columns by averaging:

    left_dlPFC  = mean(channel1,  channel2,  channel3,  channel4 )
    left_mPFC   = mean(channel5,  channel6,  channel7,  channel8 )
    right_mPFC  = mean(channel9,  channel10, channel11, channel12)
    right_dlPFC = mean(channel13, channel14, channel15, channel16)

Output files have the same name with `allchannels` replaced by `allregions`.

Note: column names use underscores (not hyphens) so they are safe to use as
identifiers in statsmodels formulas downstream. Display labels remain
"Left dlPFC", "Left mPFC", etc.

Usage
-----
    python pool_to_regions.py [--input-dir .] [--output-dir .]
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd


REGION_MAP = {
    "left_dlPFC":  [f"channel{i}" for i in range(1,  5)],   # 1-4
    "left_mPFC":   [f"channel{i}" for i in range(5,  9)],   # 5-8
    "right_mPFC":  [f"channel{i}" for i in range(9, 13)],   # 9-12
    "right_dlPFC": [f"channel{i}" for i in range(13, 17)],  # 13-16
}

CSV_BASES = [
    "lrv-con_v5",
    "lrv-con_v7",
    "lrv-pcon_v5",
    "lrv-pcon_v7",
    "pcon-con_v5",
    "pcon-con_v7",
]


def pool_one(in_path: Path, out_path: Path) -> None:
    df = pd.read_csv(in_path)
    all_chs = sum(REGION_MAP.values(), [])
    missing = [c for c in all_chs if c not in df.columns]
    if missing:
        sys.exit(f"[ERROR] {in_path.name}: missing channels {missing}")

    # Compute region means (subject-wise, i.e. row-wise)
    regions = pd.DataFrame({
        region: df[chs].mean(axis=1) for region, chs in REGION_MAP.items()
    })

    # Drop the 16 channel columns, append the 4 region columns
    meta_cols = [c for c in df.columns if c not in all_chs]
    out = pd.concat([df[meta_cols].reset_index(drop=True),
                     regions.reset_index(drop=True)], axis=1)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(out_path, index=False)
    print(f"  [✓] {in_path.name}  →  {out_path.name}   "
          f"({len(out)} rows × {len(out.columns)} cols)")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input-dir",  default=".",
                    help="Directory containing *_allchannels.csv files")
    ap.add_argument("--output-dir", default=".",
                    help="Where to write *_allregions.csv files")
    args = ap.parse_args()

    in_dir  = Path(args.input_dir).resolve()
    out_dir = Path(args.output_dir).resolve()

    print(f"[INFO] Input  dir : {in_dir}")
    print(f"[INFO] Output dir : {out_dir}")
    print(f"[INFO] Pooling map:")
    for r, chs in REGION_MAP.items():
        print(f"           {r:11s} ← {', '.join(chs)}")
    print()

    for base in CSV_BASES:
        in_path  = in_dir  / f"{base}_allchannels.csv"
        out_path = out_dir / f"{base}_allregions.csv"
        if not in_path.exists():
            print(f"  [SKIP] missing: {in_path}", file=sys.stderr)
            continue
        pool_one(in_path, out_path)

    print("\n[DONE]")
    return 0


if __name__ == "__main__":
    sys.exit(main())
