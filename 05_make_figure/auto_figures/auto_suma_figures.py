#!/usr/bin/env python3
"""
Batch fNIRS result CSVs into template-space NIfTIs and optional SUMA renders.

This replaces the manual results_template.csv -> toTemplate.ipynb step:
each input CSV is read directly, channel statistics are copied into the
template voxels listed in fNIRS_template.txt, and a .nii is written beside
a small audit CSV. If --suma-spec is supplied, the generated overlay is then
sent to drive_suma_render.sh for a headless SUMA screenshot.
"""
from __future__ import annotations

import argparse
import csv
import math
import re
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path

import nibabel as nib
import numpy as np


HERE = Path(__file__).resolve().parent
REPO = HERE.parent

DEFAULT_TEMPLATE_TXT = REPO / "fNIRS_template.txt"
DEFAULT_TEMPLATE_NII = REPO / "fNIRS_template.nii"
DEFAULT_TEMPLATE_MASK = REPO / "fNIRS_template_mask.nii"
DEFAULT_SUMA_DRIVER = REPO / "drive_suma_render.sh"

AUTO_VALUE_COLUMNS = (
    "stat",
    "ancova_t",
    "t",
    "int_sans_t",
    "int_saps_t",
    "int_wais_matrix_t",
    "int_sud_t",
    "cohens_d",
    "ancova_beta",
    "beta",
)
AUTO_P_COLUMNS = (
    "pval",
    "p",
    "ancova_p",
    "ancova_p_fdr",
    "ancova_p_hoch",
    "ancova_p_homm",
    "ancova_p_perm_max",
    "ancova_p_cluster",
    "mw_p",
    "mw_p_fdr",
    "mw_p_hoch",
    "mw_p_homm",
    "int_sans_p",
    "int_sans_p_fdr",
    "int_saps_p",
    "int_saps_p_fdr",
    "int_wais_matrix_p",
    "int_wais_matrix_p_fdr",
    "int_sud_p",
    "int_sud_p_fdr",
)


@dataclass(frozen=True)
class ChannelValue:
    channel: int
    statistic_value: float
    rendered_value: float
    significance_value: float | str | None
    significant: bool


def parse_channel(value: str) -> int:
    match = re.search(r"(\d+)", str(value))
    if not match:
        raise ValueError(f"Cannot parse channel from {value!r}")
    channel = int(match.group(1))
    if channel < 1:
        raise ValueError(f"Channel numbers must be positive, got {channel}")
    return channel


def parse_float(value: str | None, default: float | None = None) -> float | None:
    if value is None:
        return default
    text = str(value).strip()
    if text == "" or text.upper() == "NA":
        return default
    try:
        out = float(text)
    except ValueError:
        return default
    return out if math.isfinite(out) else default


def parse_bool(value: str | None) -> bool | None:
    if value is None:
        return None
    text = str(value).strip().lower()
    if text in {"true", "t", "1", "yes", "y"}:
        return True
    if text in {"false", "f", "0", "no", "n"}:
        return False
    return None


def choose_column(fieldnames: list[str], requested: str, candidates: tuple[str, ...],
                  label: str) -> str | None:
    if requested.lower() != "auto":
        if requested not in fieldnames:
            raise ValueError(
                f"Requested {label} column {requested!r} is not in CSV. "
                f"Available columns: {', '.join(fieldnames)}"
            )
        return requested
    lower_to_name = {name.lower(): name for name in fieldnames}
    for candidate in candidates:
        if candidate.lower() in lower_to_name:
            return lower_to_name[candidate.lower()]
    return None


def read_result_csv(path: Path, statistic_col: str, significance_p_col: str,
                    significance_bool_col: str | None, significance_mode: str,
                    alpha: float) -> tuple[list[ChannelValue], str, str | None, str]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames:
            raise ValueError(f"{path} has no header row")
        fieldnames = list(reader.fieldnames)
        channel_col = choose_column(fieldnames, "auto", ("channel", "ch"), "channel")
        if channel_col is None:
            raise ValueError(f"{path} has no channel column")
        selected_statistic_col = choose_column(fieldnames, statistic_col, AUTO_VALUE_COLUMNS, "statistic")
        if selected_statistic_col is None:
            raise ValueError(
                f"Could not auto-detect a statistic column in {path}. "
                "Pass --stat-col explicitly."
            )
        selected_p_col = choose_column(
            fieldnames, significance_p_col, AUTO_P_COLUMNS, "significance p-value")

        if significance_bool_col is not None and significance_bool_col not in fieldnames:
            raise ValueError(
                f"Requested significance column {significance_bool_col!r} is not in CSV")

        selected_mode = significance_mode
        if selected_mode == "auto":
            if significance_bool_col is not None:
                selected_mode = "bool"
            elif selected_p_col is not None:
                selected_mode = "p"
            else:
                selected_mode = "nonzero"
        if selected_mode == "p" and selected_p_col is None:
            raise ValueError(
                f"Could not auto-detect a p-value column in {path}. "
                "Pass --significance-p-col or choose another --significance-mode."
            )
        if selected_mode == "bool" and significance_bool_col is None:
            raise ValueError("--significance-mode bool requires --significance-bool-col")

        by_channel: dict[int, ChannelValue] = {}
        for row in reader:
            channel = parse_channel(row[channel_col])
            statistic_value = parse_float(row.get(selected_statistic_col), default=0.0)
            assert statistic_value is not None

            significance_value: float | str | None = None
            if selected_mode == "none":
                significant = abs(statistic_value) > 0
            elif selected_mode == "nonzero":
                significant = abs(statistic_value) > 0
            elif selected_mode == "bool":
                raw_sig = row.get(significance_bool_col)
                significance_value = raw_sig
                sig_from_column = parse_bool(raw_sig)
                if sig_from_column is None:
                    raise ValueError(
                        f"Could not parse boolean significance value {raw_sig!r} "
                        f"for channel {channel} in {path}"
                    )
                significant = sig_from_column
            else:
                significance_value = parse_float(row.get(selected_p_col), default=None)
                if significance_value is None:
                    significant = False
                else:
                    significant = significance_value < alpha

            rendered_value = statistic_value if significant else 0.0
            by_channel[channel] = ChannelValue(
                channel=channel,
                statistic_value=statistic_value,
                rendered_value=rendered_value,
                significance_value=significance_value,
                significant=significant,
            )

    values = [by_channel.get(ch, ChannelValue(ch, 0.0, 0.0, None, False))
              for ch in range(1, 17)]
    significance_col = (
        significance_bool_col if selected_mode == "bool" else selected_p_col
    )
    return values, selected_statistic_col, significance_col, selected_mode


def load_template_rows(path: Path) -> np.ndarray:
    rows = np.loadtxt(path)
    if rows.ndim == 1:
        rows = rows[None, :]
    if rows.shape[1] != 4:
        raise ValueError(f"{path} must have four columns: i j k channel")
    return rows


def write_overlay(values: list[ChannelValue], out_nii: Path, template_txt: Path,
                  template_nii: Path, template_mask: Path | None) -> None:
    template_img = nib.load(str(template_nii))
    data = np.zeros(template_img.shape, dtype=np.float32)
    channel_to_value = {item.channel: item.rendered_value for item in values}

    rows = load_template_rows(template_txt)
    for i_float, j_float, k_float, channel_float in rows:
        i, j, k = int(i_float), int(j_float), int(k_float)
        channel = int(channel_float)
        if not (0 <= i < data.shape[0] and 0 <= j < data.shape[1] and 0 <= k < data.shape[2]):
            raise ValueError(f"Template voxel {(i, j, k)} is outside image shape {data.shape}")
        data[i, j, k] = float(channel_to_value.get(channel, 0.0))

    if template_mask is not None:
        mask = np.asarray(nib.load(str(template_mask)).dataobj) != 0
        if mask.shape != data.shape:
            raise ValueError(
                f"Mask shape {mask.shape} does not match template shape {data.shape}"
            )
        data[~mask] = 0.0

    out_nii.parent.mkdir(parents=True, exist_ok=True)
    nib.save(nib.Nifti1Image(data, template_img.affine, template_img.header), str(out_nii))


def write_audit_csv(values: list[ChannelValue], out_csv: Path, source: Path,
                    statistic_col: str, significance_col: str | None,
                    significance_mode: str, alpha: float) -> None:
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow([
            "source",
            "channel",
            "statistic_col",
            "statistic_value",
            "significance_mode",
            "significance_col",
            "significance_value",
            "alpha",
            "significant",
            "rendered_value",
        ])
        for item in values:
            writer.writerow([
                source,
                item.channel,
                statistic_col,
                item.statistic_value,
                significance_mode,
                significance_col or "",
                "" if item.significance_value is None else item.significance_value,
                alpha,
                item.significant,
                item.rendered_value,
            ])


def run_suma(driver: Path, spec: Path, overlay: Path, out_png: Path) -> None:
    out_png.parent.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        [str(driver), str(spec), str(overlay), str(out_png)],
        check=True,
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Convert fNIRS result CSVs to NIfTI overlays and optional SUMA renders."
    )
    parser.add_argument("csvs", nargs="+", type=Path, help="Result CSV file(s).")
    parser.add_argument("--out-dir", type=Path, default=REPO / "auto_suma_outputs")
    parser.add_argument("--template-txt", type=Path, default=DEFAULT_TEMPLATE_TXT)
    parser.add_argument("--template-nii", type=Path, default=DEFAULT_TEMPLATE_NII)
    parser.add_argument("--template-mask", type=Path, default=DEFAULT_TEMPLATE_MASK)
    parser.add_argument(
        "--stat-col", "--value-col",
        dest="statistic_col",
        default="auto",
        help="Column whose values set overlay color/intensity, e.g. ancova_t or cohens_d. Default: auto.",
    )
    parser.add_argument(
        "--significance-p-col", "--p-col",
        dest="significance_p_col",
        default="auto",
        help="P-value column deciding which channels are colored. Default: auto.",
    )
    parser.add_argument(
        "--significance-bool-col", "--sig-col",
        dest="significance_bool_col",
        default=None,
        help="Boolean column deciding which channels are colored, e.g. ancova_p_cluster_sig.",
    )
    parser.add_argument(
        "--significance-mode",
        choices=["auto", "p", "bool", "nonzero", "none"],
        default="auto",
        help=(
            "How to decide whether a channel is colored. auto uses a boolean "
            "column if supplied, else p < alpha if a p column exists, else nonzero. "
            "none/nonzero color every non-zero statistic."
        ),
    )
    parser.add_argument("--alpha", type=float, default=0.05)
    parser.add_argument(
        "--no-threshold",
        action="store_true",
        help="Alias for --significance-mode none: render all non-zero statistics.",
    )
    parser.add_argument(
        "--suma-spec",
        type=Path,
        default=None,
        help="SUMA .spec file. If omitted, only NIfTI overlays and audit CSVs are written.",
    )
    parser.add_argument("--suma-driver", type=Path, default=DEFAULT_SUMA_DRIVER)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    significance_mode = "none" if args.no_threshold else args.significance_mode
    rendered = 0

    for csv_path in args.csvs:
        values, statistic_col, significance_col, selected_significance_mode = read_result_csv(
            csv_path,
            statistic_col=args.statistic_col,
            significance_p_col=args.significance_p_col,
            significance_bool_col=args.significance_bool_col,
            significance_mode=significance_mode,
            alpha=args.alpha,
        )
        stem = csv_path.stem
        out_nii = args.out_dir / "nii" / f"{stem}.nii"
        out_audit = args.out_dir / "audit" / f"{stem}.channels.csv"
        write_overlay(
            values,
            out_nii=out_nii,
            template_txt=args.template_txt,
            template_nii=args.template_nii,
            template_mask=args.template_mask,
        )
        write_audit_csv(
            values, out_audit, csv_path, statistic_col, significance_col,
            selected_significance_mode, args.alpha)

        nonzero = sum(1 for item in values if abs(item.rendered_value) > 0)
        print(f"[auto_suma] wrote {out_nii} ({nonzero} rendered channels)")
        print(f"[auto_suma] wrote {out_audit}")

        if args.suma_spec is not None:
            out_png = args.out_dir / "png" / f"{stem}.suma.png"
            run_suma(args.suma_driver, args.suma_spec, out_nii, out_png)
            print(f"[auto_suma] wrote {out_png}")
            rendered += 1

    if args.suma_spec is None:
        print("[auto_suma] no --suma-spec supplied, skipped SUMA screenshots")
    else:
        print(f"[auto_suma] rendered {rendered} SUMA screenshot(s)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
