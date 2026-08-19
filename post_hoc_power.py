#!/usr/bin/env python3
"""
Post-hoc power analysis for the model V.1 (performance-weighted, corsi contrast)
ANCOVA group comparisons.

Reproduces the power numbers reported in the manuscript's Limitations section
(Frontiers v.1.5 KCL.docx) directly from the group-level result CSVs, rather
than from numbers typed by hand.

For each of the three V.1 group comparisons this script:
  1. Reads the per-channel ANCOVA group-effect t-statistic from the result CSV.
  2. Derives the residual degrees of freedom from the achieved sample size and
     the number of model parameters (group term + covariates + intercept),
     matching the covariate sets used in run_group_analyses.sh:
       - PSD+AGG vs PSD-AGG: wais_matrix, sud, dose  (3 covariates)
       - PSD vs HC          : wais_matrix, sud        (2 covariates)
  3. Computes observed ("post-hoc") power for each channel using the
     noncentral t-distribution, treating the observed |t| as the
     noncentrality parameter for a two-sided test at alpha=0.05.
  4. Computes the minimum |t| required for 80% power at each comparison's df.
  5. Estimates the total sample size that would be required for 80% power at
     an effect the size of a given channel (default: channel 15, the weaker
     of the two channels surviving cluster-based permutation correction),
     under the approximation that the coefficient standard error scales as
     1/sqrt(n) with the covariate structure and residual variance held fixed.

Run with:
    04_secondlevel_group/.venv/bin/python3 post_hoc_power.py

Requires: pandas, scipy (both already present in 04_secondlevel_group/.venv).
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

import pandas as pd
from scipy import optimize, stats

REPO_ROOT = Path(__file__).resolve().parent
RESULTS_DIR = REPO_ROOT / "04_secondlevel_group" / "analysis_2026-08-17"

ALPHA = 0.05
TARGET_POWER = 0.80


@dataclass
class Comparison:
    label: str
    csv_path: Path
    n_covariates: int  # excludes group term and intercept


COMPARISONS = [
    Comparison(
        label="PSD+AGG vs PSD-AGG",
        csv_path=RESULTS_DIR
        / "results_2026-08-17_psd_agg_vs_psd_nagg_corsi"
        / "results_2026-08-17_psd_agg_vs_psd_nagg_corsi.csv",
        n_covariates=3,  # wais_matrix, sud, dose
    ),
    Comparison(
        label="PSD+AGG vs HC",
        csv_path=RESULTS_DIR
        / "results_2026-08-17_psd_agg_vs_hc_corsi"
        / "results_2026-08-17_psd_agg_vs_hc_corsi.csv",
        n_covariates=2,  # wais_matrix, sud
    ),
    Comparison(
        label="PSD-AGG vs HC",
        csv_path=RESULTS_DIR
        / "results_2026-08-17_psd_nagg_vs_hc_corsi"
        / "results_2026-08-17_psd_nagg_vs_hc_corsi.csv",
        n_covariates=2,  # wais_matrix, sud
    ),
]

N_MODEL_TERMS_EXCL_COVARIATES = 2  # intercept + group term


def observed_power(t_obs: float, df: float, alpha: float = ALPHA) -> float:
    """Two-sided post-hoc power for the observed |t|, given residual df."""
    t_crit = stats.t.ppf(1 - alpha / 2, df)
    nc = abs(t_obs)
    return 1 - stats.nct.cdf(t_crit, df, nc) + stats.nct.cdf(-t_crit, df, nc)


def t_required_for_power(df: float, target_power: float = TARGET_POWER,
                          alpha: float = ALPHA) -> float:
    """Minimum |t| (noncentrality parameter) needed to reach target_power."""
    t_crit = stats.t.ppf(1 - alpha / 2, df)

    def f(nc: float) -> float:
        return (1 - stats.nct.cdf(t_crit, df, nc) + stats.nct.cdf(-t_crit, df, nc)) - target_power

    return optimize.brentq(f, 1e-6, 20)


def required_n_for_power(t_obs: float, n_current: int, df: float,
                          target_power: float = TARGET_POWER,
                          alpha: float = ALPHA) -> float:
    """
    Total sample size needed for target_power at an effect the size of
    t_obs, assuming se ~ 1/sqrt(n) (covariate structure and residual
    variance held fixed as n grows).
    """
    nc_target = t_required_for_power(df, target_power, alpha)
    return n_current * (nc_target / abs(t_obs)) ** 2


def load_comparison(comp: Comparison) -> pd.DataFrame:
    df = pd.read_csv(comp.csv_path)
    n_total = int(df["n_pos"].iloc[0] + df["n_neg"].iloc[0])
    n_params = N_MODEL_TERMS_EXCL_COVARIATES + comp.n_covariates
    resid_df = n_total - n_params
    df["comparison"] = comp.label
    df["n_total"] = n_total
    df["resid_df"] = resid_df
    df["observed_power"] = df["ancova_t"].apply(lambda t: observed_power(t, resid_df))
    return df


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--reference-channel", default="channel15",
        help="Channel used for the required-sample-size illustration "
             "(default: channel15, the weaker of the two cluster-corrected "
             "channels in PSD+AGG vs PSD-AGG).",
    )
    parser.add_argument(
        "--reference-comparison", default="PSD+AGG vs PSD-AGG",
        help="Comparison label the reference channel is drawn from.",
    )
    args = parser.parse_args()

    print(f"alpha = {ALPHA}, target power = {TARGET_POWER}\n")

    all_rows = []
    for comp in COMPARISONS:
        d = load_comparison(comp)
        all_rows.append(d)

        t_req = t_required_for_power(d["resid_df"].iloc[0])
        print(f"=== {comp.label} ===")
        print(f"  n = {d['n_total'].iloc[0]}  (n_pos={d['n_pos'].iloc[0]}, "
              f"n_neg={d['n_neg'].iloc[0]}), covariates={comp.n_covariates}, "
              f"resid_df={d['resid_df'].iloc[0]}")
        print(f"  |t| required for {int(TARGET_POWER*100)}% power: {t_req:.2f}")

        top = d.reindex(d["ancova_t"].abs().sort_values(ascending=False).index)
        print("  channels ranked by |t| (top 5):")
        for _, row in top.head(5).iterrows():
            sig = "SURVIVED cluster correction" if bool(row["ancova_p_cluster_sig"]) else "did not survive correction"
            print(f"    {row['channel']:>10s}: t={row['ancova_t']:+.3f}  "
                  f"observed power={row['observed_power']:.3f}  ({sig})")
        print()

    combined = pd.concat(all_rows, ignore_index=True)
    out_path = REPO_ROOT / "post_hoc_power_results.csv"
    combined.to_csv(out_path, index=False)
    print(f"Full per-channel table written to {out_path.relative_to(REPO_ROOT)}\n")

    ref_comp = next(c for c in COMPARISONS if c.label == args.reference_comparison)
    ref_df = combined[combined["comparison"] == args.reference_comparison]
    ref_row = ref_df[ref_df["channel"] == args.reference_channel].iloc[0]
    n_req = required_n_for_power(ref_row["ancova_t"], int(ref_row["n_total"]), ref_row["resid_df"])
    print(
        f"Sample size needed for {int(TARGET_POWER*100)}% power at an effect the size of "
        f"{args.reference_channel} in '{args.reference_comparison}' "
        f"(t={ref_row['ancova_t']:.3f}, current n={int(ref_row['n_total'])}):\n"
        f"  estimated total n ~= {n_req:.0f} "
        f"(assuming se scales as 1/sqrt(n), covariate structure and residual "
        f"variance held fixed)."
    )


if __name__ == "__main__":
    main()
