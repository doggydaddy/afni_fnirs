#!/usr/bin/env python3
"""
analyze_psd_agg_vs_psd_nagg.py
================================
Second-level group analysis: PSD+AGG vs PSD-AGG
All 16 fNIRS channels, covariate-corrected (SANS, SAPS, DOSE)

Analysis pipeline
-----------------
1. Load data (default: lrv-pcon_v5_allchannels.csv)
2. Descriptive statistics per group per channel
3. For each channel:
   a. Mann-Whitney U-test (non-parametric, unadjusted)
   b. ANCOVA via OLS: channel ~ C(group) + sans + saps + dose
4. FDR correction (Benjamini-Hochberg) over 16 channels
   — applied separately to Mann-Whitney and ANCOVA p-values
5. Effect-size: rank-biserial r (Mann-Whitney) and Cohen's d (adjusted means)
6. Output:
   • Console table
   • results_psd_agg_vs_psd_nagg.csv
   • figures/  (group means bar chart, volcano plot, FDR-corrected p-value heatmap)

Usage
-----
python analyze_psd_agg_vs_psd_nagg.py [--input CSV] [--output-dir DIR] [--alpha 0.05]
"""

import argparse
import os
import sys
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, t as t_dist
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests
import matplotlib
matplotlib.use("Agg")   # headless-safe; swap to "TkAgg" / "Qt5Agg" for interactive
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns

# ──────────────────────────────────────────────
# CLI
# ──────────────────────────────────────────────
def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input",
                   default="lrv-pcon_v5_allchannels.csv",
                   help="Input CSV (default: lrv-pcon_v5_allchannels.csv)")
    p.add_argument("--output-dir",
                   default=".",
                   help="Directory for CSV results and figures (default: .)")
    p.add_argument("--alpha",
                   type=float,
                   default=0.05,
                   help="Significance threshold after FDR correction (default: 0.05)")
    p.add_argument("--group-pos",
                   default="PSD+AGG",
                   help="Label for the positive / reference group (default: PSD+AGG)")
    p.add_argument("--group-neg",
                   default="PSD-AGG",
                   help="Label for the negative / comparison group (default: PSD-AGG)")
    p.add_argument("--covariates",
                   nargs="+",
                   default=["sans", "saps", "dose"],
                   help="Covariates to include in ANCOVA (default: sans saps dose). "
                        "Tip: drop 'dose' when comparing against HC (dose=0 for all HC).")
    p.add_argument("--n-perms",
                   type=int,
                   default=5000,
                   help="Number of permutations for max-stat and cluster tests (default: 5000)")
    p.add_argument("--cluster-thresh",
                   type=float,
                   default=0.05,
                   help="Uncorrected p-threshold for cluster-forming (default: 0.05)")
    p.add_argument("--interactions",
                   action="store_true",
                   default=False,
                   help="Run secondary interaction analysis: "
                        "channel ~ C(group) × (covariates). "
                        "Saves interaction_results.csv alongside main results.")
    p.add_argument("--targets",
                   nargs="+",
                   default=None,
                   help="Names of the dependent-variable columns (channels/regions). "
                        "If omitted, auto-detects columns matching r'^channel\\d+$' from the input CSV.")
    p.add_argument("--target-labels",
                   nargs="+",
                   default=None,
                   help="Display labels for --targets (same length). "
                        "Defaults to the targets with 'channel' → 'Ch' substitution.")
    p.add_argument("--no-cluster",
                   action="store_true",
                   default=False,
                   help="Skip cluster-based permutation correction. Recommended when the "
                        "number of targets is small (e.g. a handful of ROIs) since cluster-perm "
                        "is largely redundant with permutation max-stat in that regime. "
                        "Cluster columns in the output CSV will be NaN.")
    return p.parse_args()


# ──────────────────────────────────────────────
# Helpers
# ──────────────────────────────────────────────
GROUP_COL = "group"


def derive_targets(df, cli_targets):
    """Return (targets, labels) given the dataframe and optional CLI override."""
    import re
    if cli_targets:
        targets = list(cli_targets)
    else:
        targets = sorted(
            [c for c in df.columns if re.fullmatch(r"channel\d+", c)],
            key=lambda c: int(c[len("channel"):]),
        )
    if not targets:
        sys.exit("[ERROR] Could not auto-detect any target columns "
                 "(expected 'channelN'). Use --targets to specify explicitly.")
    return targets


def derive_labels(targets, cli_labels):
    """Short labels for plotting axes; falls back to a sensible auto-mapping."""
    if cli_labels:
        if len(cli_labels) != len(targets):
            sys.exit(f"[ERROR] --target-labels has {len(cli_labels)} entries "
                     f"but --targets has {len(targets)}.")
        return list(cli_labels)
    return [t.replace("channel", "Ch") for t in targets]


def linear_adjacency(n):
    """Linear-chain adjacency: target_i ~ target_{i-1}, target_{i+1}."""
    return {i: [j for j in (i - 1, i + 1) if 0 <= j < n] for i in range(n)}


def rank_biserial_r(group1_vals, group2_vals):
    """Rank-biserial correlation as effect size for Mann-Whitney U."""
    n1, n2 = len(group1_vals), len(group2_vals)
    if n1 == 0 or n2 == 0:
        return np.nan
    u_stat, _ = mannwhitneyu(group1_vals, group2_vals, alternative="two-sided")
    r = 1 - (2 * u_stat) / (n1 * n2)
    return r


def cohens_d(group1_vals, group2_vals):
    """Cohen's d: (mean1 - mean2) / pooled SD."""
    m1, m2 = np.mean(group1_vals), np.mean(group2_vals)
    s1, s2 = np.std(group1_vals, ddof=1), np.std(group2_vals, ddof=1)
    n1, n2 = len(group1_vals), len(group2_vals)
    pooled_sd = np.sqrt(((n1 - 1) * s1**2 + (n2 - 1) * s2**2) / (n1 + n2 - 2))
    return (m1 - m2) / pooled_sd if pooled_sd > 0 else np.nan


def run_ancova(df, channel, group_pos, group_neg, covariates):
    """
    OLS ANCOVA: channel ~ C(group, Treatment(ref)) + covariate1 + covariate2 ...
    Returns (beta_group, se_group, t_group, p_ancova, r_squared)
    """
    df_sub = df.dropna(subset=[GROUP_COL, channel] + covariates)
    # keep only the two groups of interest
    df_sub = df_sub[df_sub[GROUP_COL].isin([group_pos, group_neg])].copy()
    if len(df_sub) < 6:
        return (np.nan,) * 5

    cov_terms = " + ".join(covariates)
    formula = (
        f"Q('{channel}') ~ C({GROUP_COL}, Treatment(reference='{group_neg}'))"
        f" + {cov_terms}"
    )
    try:
        model = smf.ols(formula, data=df_sub).fit()
        term = f"C({GROUP_COL}, Treatment(reference='{group_neg}'))[T.{group_pos}]"
        beta   = model.params.get(term, np.nan)
        se     = model.bse.get(term, np.nan)
        t_val  = model.tvalues.get(term, np.nan)
        p_val  = model.pvalues.get(term, np.nan)
        r2     = model.rsquared
        return beta, se, t_val, p_val, r2
    except Exception as e:
        print(f"  [WARN] ANCOVA failed for {channel}: {e}", file=sys.stderr)
        return (np.nan,) * 5


def run_interaction_model(df, channel, group_pos, group_neg, covariates):
    """
    Fits the full interaction model:
      channel ~ C(group, Treatment(ref)) * (cov1 + cov2 + ...)
    Returns dict {cov: (beta, t, p)} for each group × covariate interaction term.
    """
    df_sub = df.dropna(subset=[GROUP_COL, channel] + covariates)
    df_sub = df_sub[df_sub[GROUP_COL].isin([group_pos, group_neg])].copy()
    if len(df_sub) < 8:
        return {c: (np.nan, np.nan, np.nan) for c in covariates}

    cov_terms = " + ".join(covariates)
    formula = (
        f"Q('{channel}') ~ C({GROUP_COL}, Treatment(reference='{group_neg}'))"
        f" * ({cov_terms})"
    )
    out = {}
    try:
        model = smf.ols(formula, data=df_sub).fit()
        for cov in covariates:
            term = (f"C({GROUP_COL}, Treatment(reference='{group_neg}'))"
                    f"[T.{group_pos}]:{cov}")
            out[cov] = (
                model.params.get(term,   np.nan),
                model.tvalues.get(term,  np.nan),
                model.pvalues.get(term,  np.nan),
            )
    except Exception as e:
        print(f"  [WARN] Interaction model failed for {channel}: {e}", file=sys.stderr)
        out = {c: (np.nan, np.nan, np.nan) for c in covariates}
    return out


# ──────────────────────────────────────────────
# Permutation test helpers  (Frisch-Waugh based, vectorised over channels)
# ──────────────────────────────────────────────

# Default adjacency builder: linear chain  target_1 – target_2 – ... – target_n
# Override by passing a custom dict (0-based indices) to run_perm_cluster, e.g.
# for a 4×4 grid:
#   {i: [j for j in
#       [i-1 if i%4 != 0 else -1,  i+1 if (i+1)%4 != 0 else -1, i-4, i+4]
#       if 0 <= j < 16] for i in range(16)}


def _fw_precompute(df_sub, covariates, targets):
    """
    Frisch-Waugh-Lovell setup.
    Returns (Z, ZtZinv, MY, g, n_params) where:
      Z     – (n, q) nuisance matrix [1 | covariates]
      ZtZinv– precomputed (Z'Z)^-1
      MY    – (n, n_targets) covariate-residualised target matrix
      g     – (n,)  binary group indicator (1 = group_pos)
      n_params – total params in the full model (used for df_resid)
    Uses complete-case rows only.
    """
    Y = df_sub[targets].values.astype(float)
    g = df_sub["_group_bin"].values.astype(float)
    Z = np.column_stack(
        [np.ones(len(df_sub))] + [df_sub[c].values.astype(float) for c in covariates]
    )
    ZtZinv = np.linalg.pinv(Z.T @ Z)
    MY = Y - Z @ (ZtZinv @ (Z.T @ Y))
    n_params = Z.shape[1] + 1          # intercept + n_covariates + group
    return Z, ZtZinv, MY, g, n_params


def _t_stats(MY, Z, ZtZinv, g, n_params):
    """
    Vectorised group t-stats for all channels using Frisch-Waugh.
    Only the group column is permuted between calls; Z / ZtZinv are fixed.
    """
    n = len(g)
    Mg     = g - Z @ (ZtZinv @ (Z.T @ g))   # residualise group   (n,)
    Mg_sq  = Mg @ Mg
    beta_g = (Mg @ MY) / Mg_sq              # group betas          (n_ch,)
    e      = MY - np.outer(Mg, beta_g)      # full-model residuals (n, n_ch)
    sigma2 = np.sum(e ** 2, axis=0) / (n - n_params)
    return beta_g / np.sqrt(sigma2 / Mg_sq)  # t-stats              (n_ch,)


def run_perm_max_stat(MY, Z, ZtZinv, g, n_params, n_perms, seed=42):
    """
    Max-statistic permutation test (strong FWER control).
    For each permutation, shuffles the group label, computes t-stats for all
    channels, records max|t|.  Channel p = proportion(max|t|_null >= |t|_obs).
    Returns (t_obs, p_perm)  both shape (n_ch,).
    """
    rng    = np.random.default_rng(seed)
    t_obs  = _t_stats(MY, Z, ZtZinv, g, n_params)
    null   = np.empty(n_perms)
    for i in range(n_perms):
        null[i] = np.max(np.abs(_t_stats(MY, Z, ZtZinv,
                                         rng.permutation(g), n_params)))
    p_perm = np.array([(null >= abs(t)).mean() for t in t_obs])
    return t_obs, p_perm


def _find_clusters(t_stats, adj, t_thresh):
    """BFS over channels whose |t| > t_thresh.  Returns [(indices, mass), ...]."""
    active  = set(np.where(np.abs(t_stats) > t_thresh)[0])
    visited = set()
    clusters = []
    for seed in list(active):
        if seed in visited:
            continue
        cluster, queue = [], [seed]
        while queue:
            node = queue.pop(0)
            if node in visited:
                continue
            visited.add(node)
            if node in active:
                cluster.append(node)
                queue.extend(nb for nb in adj[node] if nb not in visited)
        if cluster:
            clusters.append((cluster, float(np.sum(t_stats[cluster]))))
    return clusters


def run_perm_cluster(MY, Z, ZtZinv, g, n_params, adj,
                     n_perms, cluster_thresh=0.05, seed=42):
    """
    Cluster-based permutation test (strong FWER control).
    Cluster-forming threshold = t_dist.ppf(1 - cluster_thresh/2, df=n-n_params).
    Returns (ch_p, obs_clusters, t_thresh).
    """
    n       = len(g)
    df_res  = n - n_params
    t_thr   = t_dist.ppf(1 - cluster_thresh / 2, df=df_res)
    rng     = np.random.default_rng(seed)

    t_obs        = _t_stats(MY, Z, ZtZinv, g, n_params)
    obs_clusters = _find_clusters(t_obs, adj, t_thr)

    null_max = np.zeros(n_perms)
    for i in range(n_perms):
        t_p = _t_stats(MY, Z, ZtZinv, rng.permutation(g), n_params)
        pc  = _find_clusters(t_p, adj, t_thr)
        null_max[i] = max((abs(m) for _, m in pc), default=0.0)

    ch_p = np.ones(len(t_obs))
    for idx_list, mass in obs_clusters:
        pv = (null_max >= abs(mass)).mean()
        for idx in idx_list:
            ch_p[idx] = pv

    return ch_p, obs_clusters, t_thr


# ──────────────────────────────────────────────
# Main analysis
# ──────────────────────────────────────────────
def main():
    args = parse_args()

    # ── 1. Load data ──────────────────────────
    input_path = os.path.join(os.path.dirname(__file__), args.input) \
                 if not os.path.isabs(args.input) else args.input
    if not os.path.exists(input_path):
        sys.exit(f"[ERROR] Input file not found: {input_path}")

    df = pd.read_csv(input_path)
    df.columns = df.columns.str.strip().str.lower()   # normalise column names

    # Confirm groups exist
    present_groups = df[GROUP_COL].unique()
    for g in [args.group_pos, args.group_neg]:
        if g not in present_groups:
            sys.exit(f"[ERROR] Group '{g}' not found in column '{GROUP_COL}'. "
                     f"Present: {list(present_groups)}")

    df_analysis = df[df[GROUP_COL].isin([args.group_pos, args.group_neg])].copy()

    # Use covariates from CLI (allows dropping dose for HC comparisons)
    covariates = args.covariates

    # Resolve target columns (channels or regions). Header was lowercased above,
    # so apply the same to CLI overrides for a clean match.
    cli_targets = [t.lower() for t in args.targets] if args.targets else None
    targets = derive_targets(df, cli_targets)
    labels  = derive_labels(targets, args.target_labels)
    missing = [t for t in targets if t not in df.columns]
    if missing:
        sys.exit(f"[ERROR] Target columns not in {args.input}: {missing}\n"
                 f"        Available: {[c for c in df.columns if c not in ('path','volume','group','id') and c not in covariates]}")

    n_pos = (df_analysis[GROUP_COL] == args.group_pos).sum()
    n_neg = (df_analysis[GROUP_COL] == args.group_neg).sum()
    print(f"\n{'='*60}")
    print(f"  {args.group_pos} vs {args.group_neg}  |  fNIRS beta-coefficient analysis")
    print(f"{'='*60}")
    print(f"  Input : {args.input}")
    print(f"  Groups: {args.group_pos} (n={n_pos})  vs  {args.group_neg} (n={n_neg})")
    print(f"  Covariates : {', '.join(covariates)}")
    print(f"  Targets    : {len(targets)}  ({', '.join(labels)})")
    print(f"  Alpha (FDR): {args.alpha}")
    print(f"{'='*60}\n")

    # ── 2. Descriptive statistics ─────────────
    print("── Descriptive statistics (covariates) ──")
    desc = df_analysis.groupby(GROUP_COL)[covariates].agg(["mean", "std"]).round(2)
    print(desc.to_string())
    print()

    # ── 3. Per-channel tests ──────────────────
    rows = []
    for ch in targets:
        g_pos = df_analysis[df_analysis[GROUP_COL] == args.group_pos][ch].dropna()
        g_neg = df_analysis[df_analysis[GROUP_COL] == args.group_neg][ch].dropna()

        # Descriptives
        mean_pos, sd_pos = g_pos.mean(), g_pos.std()
        mean_neg, sd_neg = g_neg.mean(), g_neg.std()

        # Mann-Whitney U
        if len(g_pos) > 0 and len(g_neg) > 0:
            mw_stat, mw_p = mannwhitneyu(g_pos, g_neg, alternative="two-sided")
            rb_r           = rank_biserial_r(g_pos.values, g_neg.values)
            cd             = cohens_d(g_pos.values, g_neg.values)
        else:
            mw_stat = mw_p = rb_r = cd = np.nan

        # ANCOVA
        beta, se, t_val, anc_p, r2 = run_ancova(df_analysis, ch,
                                                  args.group_pos, args.group_neg,
                                                  covariates)

        rows.append({
            "channel"      : ch,
            "n_pos"        : len(g_pos),
            "n_neg"        : len(g_neg),
            "mean_pos"     : mean_pos,
            "sd_pos"       : sd_pos,
            "mean_neg"     : mean_neg,
            "sd_neg"       : sd_neg,
            "mw_stat"      : mw_stat,
            "mw_p"         : mw_p,
            "rb_r"         : rb_r,           # rank-biserial r
            "cohens_d"     : cd,
            "ancova_beta"  : beta,           # PSD+AGG – PSD-AGG, adjusted
            "ancova_se"    : se,
            "ancova_t"     : t_val,
            "ancova_p"     : anc_p,
            "ancova_r2"    : r2,
        })

    results = pd.DataFrame(rows)

    # ── 4. Multiple-comparison corrections ───────────────────────────────
    # Three methods, all over all 16 channels:
    #   fdr_bh   – Benjamini/Hochberg (FDR control, least conservative)
    #   hochberg – Hochberg step-up   (FWER control)
    #   hommel   – Hommel             (FWER control, most powerful of the three)
    METHODS = [
        ("fdr_bh",        "fdr"),
        ("simes-hochberg","hoch"),
        ("hommel",        "homm"),
    ]
    for test_col in ("ancova_p", "mw_p"):
        valid = results[test_col].notna()
        for method, tag in METHODS:
            adj = np.full(len(results), np.nan)
            if valid.sum() > 0:
                _, padj, _, _ = multipletests(
                    results.loc[valid, test_col], alpha=args.alpha, method=method)
                adj[valid] = padj
            results[f"{test_col}_{tag}"]     = adj
            results[f"{test_col}_{tag}_sig"] = adj < args.alpha

    # ── 4c. Permutation tests ─────────────────────────────────────────────
    # Use complete cases across ALL targets for a single consistent design matrix
    cols_needed = [GROUP_COL] + covariates + targets
    df_perm = (df_analysis
               .dropna(subset=cols_needed)
               .query(f"{GROUP_COL} in ['{args.group_pos}', '{args.group_neg}']")
               .copy())
    df_perm["_group_bin"] = (df_perm[GROUP_COL] == args.group_pos).astype(float)
    n_cc = len(df_perm)

    Z, ZtZinv, MY, g, n_params = _fw_precompute(df_perm, covariates, targets)

    print(f"\n── Permutation tests  ({args.n_perms:,} perms, n={n_cc} complete cases) ──")

    t_obs, p_max = run_perm_max_stat(MY, Z, ZtZinv, g, n_params,
                                     n_perms=args.n_perms)

    if args.no_cluster:
        print("  [info] Cluster-based permutation SKIPPED (--no-cluster). "
              "Cluster columns will be NaN in the output.")
        p_clust       = np.full(len(t_obs), np.nan)
        obs_clusters  = []
        t_thr         = np.nan
    else:
        p_clust, obs_clusters, t_thr = run_perm_cluster(
            MY, Z, ZtZinv, g, n_params, linear_adjacency(len(targets)),
            n_perms=args.n_perms, cluster_thresh=args.cluster_thresh)

    results["ancova_t"]               = t_obs
    results["ancova_p_perm_max"]      = p_max
    results["ancova_p_perm_max_sig"]  = p_max < args.alpha
    results["ancova_p_cluster"]       = p_clust
    results["ancova_p_cluster_sig"]   = np.where(np.isnan(p_clust), False, p_clust < args.alpha)

    # Report observed clusters
    if args.no_cluster:
        pass  # already announced above
    elif obs_clusters:
        print(f"  Cluster-forming |t| threshold: {t_thr:.3f}")
        for k, (idx_list, mass) in enumerate(obs_clusters, 1):
            chs = ", ".join(labels[i] for i in sorted(idx_list))
            pv  = p_clust[idx_list[0]]
            print(f"  Cluster {k}: [{chs}]  mass={mass:.3f}  p={pv:.4f}"
                  f"  {'✓' if pv < args.alpha else '—'}")
    else:
        print(f"  No clusters above |t|={t_thr:.3f} threshold.")

    # ── 5. Print summary table ────────────────
    print("\n── Per-channel results (ANCOVA — all corrections) ──")
    tbl = results[[
        "channel", "n_pos", "n_neg",
        "ancova_beta", "ancova_p",
        "ancova_p_fdr",  "ancova_p_fdr_sig",
        "ancova_p_hoch", "ancova_p_hoch_sig",
        "ancova_p_homm", "ancova_p_homm_sig",
        "ancova_p_perm_max", "ancova_p_perm_max_sig",
        "ancova_p_cluster",  "ancova_p_cluster_sig",
        "cohens_d",
    ]].copy()
    for col in ["ancova_beta"]:
        tbl[col] = tbl[col].map(lambda v: f"{v:.3e}" if pd.notna(v) else "NA")
    for col in ["ancova_p", "ancova_p_fdr", "ancova_p_hoch", "ancova_p_homm",
                "ancova_p_perm_max", "ancova_p_cluster"]:
        tbl[col] = tbl[col].map(lambda v: f"{v:.4f}" if pd.notna(v) else "NA")
    tbl["cohens_d"] = tbl["cohens_d"].map(lambda v: f"{v:.3f}" if pd.notna(v) else "NA")
    for col in ["ancova_p_fdr_sig", "ancova_p_hoch_sig", "ancova_p_homm_sig",
                "ancova_p_perm_max_sig", "ancova_p_cluster_sig"]:
        tbl[col] = tbl[col].map({True: "YES", False: "---"})
    print(tbl.rename(columns={
        "ancova_p":              "p_raw",
        "ancova_p_fdr":          "p_FDR",    "ancova_p_fdr_sig":      "sig_FDR",
        "ancova_p_hoch":         "p_Hoch",   "ancova_p_hoch_sig":     "sig_Hoch",
        "ancova_p_homm":         "p_Homm",   "ancova_p_homm_sig":     "sig_Homm",
        "ancova_p_perm_max":     "p_Pmax",   "ancova_p_perm_max_sig": "sig_Pmax",
        "ancova_p_cluster":      "p_Clust",  "ancova_p_cluster_sig":  "sig_Clust",
        "ancova_beta":           "beta",
    }).to_string(index=False))

    print("\n── Per-channel results (Mann-Whitney — all three corrections) ──")
    tbl2 = results[[
        "channel", "mw_p",
        "mw_p_fdr",  "mw_p_fdr_sig",
        "mw_p_hoch", "mw_p_hoch_sig",
        "mw_p_homm", "mw_p_homm_sig",
        "rb_r",
    ]].copy()
    for col in ["mw_p", "mw_p_fdr", "mw_p_hoch", "mw_p_homm"]:
        tbl2[col] = tbl2[col].map(lambda v: f"{v:.4f}" if pd.notna(v) else "NA")
    for col in ["rb_r"]:
        tbl2[col] = tbl2[col].map(lambda v: f"{v:.3f}" if pd.notna(v) else "NA")
    for col in ["mw_p_fdr_sig", "mw_p_hoch_sig", "mw_p_homm_sig"]:
        tbl2[col] = tbl2[col].map({True: "YES", False: "---"})
    print(tbl2.rename(columns={
        "mw_p":          "p_raw",
        "mw_p_fdr":      "p_FDR-BH",  "mw_p_fdr_sig":  "sig_FDR",
        "mw_p_hoch":     "p_Hochberg", "mw_p_hoch_sig": "sig_Hoch",
        "mw_p_homm":     "p_Hommel",   "mw_p_homm_sig": "sig_Homm",
    }).to_string(index=False))

    # ── 6. Save CSV ───────────────────────────
    os.makedirs(args.output_dir, exist_ok=True)
    csv_out = os.path.join(args.output_dir, "results_psd_agg_vs_psd_nagg.csv")
    results.to_csv(csv_out, index=False, float_format="%.6e")
    print(f"\n[✓] Results saved → {csv_out}")

    # ── 7. Figures — use Hommel as the primary highlight ──────────────────
    fig_dir = os.path.join(args.output_dir, "figures")
    os.makedirs(fig_dir, exist_ok=True)
    _plot_group_means(results, args.group_pos, args.group_neg, args.alpha, fig_dir, labels)
    _plot_correction_heatmap(results, args.alpha, fig_dir, labels)
    _plot_volcano(results, args.alpha, fig_dir, labels)
    print(f"[✓] Figures saved → {fig_dir}/")

    # ── 8. Interaction analysis (secondary, opt-in) ───────────────────────────
    if args.interactions:
        print(f"\n── Interaction analysis: {args.group_pos} vs {args.group_neg} ──")
        print(f"   Model : channel ~ C(group) × ({' + '.join(covariates)})")
        print(f"   FDR-BH correction across {len(targets)} targets per interaction term\n")

        int_rows = []
        for ch in targets:
            row = {"channel": ch}
            for cov, (b, t, p) in run_interaction_model(
                    df_analysis, ch, args.group_pos, args.group_neg, covariates
            ).items():
                row[f"int_{cov}_beta"] = b
                row[f"int_{cov}_t"]    = t
                row[f"int_{cov}_p"]    = p
            int_rows.append(row)

        int_df = pd.DataFrame(int_rows)

        # FDR-BH correction: one family per covariate interaction term (16 channels each)
        for cov in covariates:
            p_col = f"int_{cov}_p"
            valid = int_df[p_col].notna()
            adj   = np.full(len(int_df), np.nan)
            if valid.sum() > 0:
                _, padj, _, _ = multipletests(int_df.loc[valid, p_col],
                                              alpha=args.alpha, method="fdr_bh")
                adj[valid] = padj
            int_df[f"int_{cov}_p_fdr"] = adj
            int_df[f"int_{cov}_sig"]   = adj < args.alpha

        # Print one compact table per covariate
        for cov in covariates:
            any_sig = int_df[f"int_{cov}_sig"].any()
            print(f"  group × {cov.upper()}  {'← significant channels found' if any_sig else ''}")
            tbl = int_df[["channel",
                           f"int_{cov}_beta", f"int_{cov}_t",
                           f"int_{cov}_p",    f"int_{cov}_p_fdr",
                           f"int_{cov}_sig"]].copy()
            tbl[f"int_{cov}_beta"]  = tbl[f"int_{cov}_beta"].map(
                lambda v: f"{v:.3e}" if pd.notna(v) else "NA")
            tbl[f"int_{cov}_t"]     = tbl[f"int_{cov}_t"].map(
                lambda v: f"{v:.2f}"  if pd.notna(v) else "NA")
            for pc in [f"int_{cov}_p", f"int_{cov}_p_fdr"]:
                tbl[pc] = tbl[pc].map(lambda v: f"{v:.4f}" if pd.notna(v) else "NA")
            tbl[f"int_{cov}_sig"] = tbl[f"int_{cov}_sig"].map({True: "YES ✓", False: "---"})
            print(tbl.rename(columns={
                f"int_{cov}_beta":  "beta",
                f"int_{cov}_t":     "t",
                f"int_{cov}_p":     "p_raw",
                f"int_{cov}_p_fdr": "p_FDR",
                f"int_{cov}_sig":   "sig",
            }).to_string(index=False))
            print()

        int_csv = os.path.join(args.output_dir, "interaction_results.csv")
        int_df.to_csv(int_csv, index=False, float_format="%.6e")
        print(f"[✓] Interaction results saved → {int_csv}")

    print()


# ──────────────────────────────────────────────
# Plotting helpers
# ──────────────────────────────────────────────

def _plot_group_means(results, group_pos, group_neg, alpha, out_dir, labels):
    """Bar chart: mean betas per target. Stars = Hommel-significant (most powerful FWER)."""
    n = len(labels)
    fig, ax = plt.subplots(figsize=(max(6, 0.85 * n + 2), 5))
    x = np.arange(n)
    w = 0.35

    sig_homm = results["ancova_p_homm_sig"].values
    colors_pos = ["#d62728" if s else "#ff9896" for s in sig_homm]
    colors_neg = ["#1f77b4" if s else "#aec7e8" for s in sig_homm]

    ax.bar(x - w/2, results["mean_pos"],  width=w, color=colors_pos,
           label=group_pos,  edgecolor="black", linewidth=0.5)
    ax.bar(x + w/2, results["mean_neg"],  width=w, color=colors_neg,
           label=group_neg,  edgecolor="black", linewidth=0.5)

    # Mark Hommel-significant targets with *
    for i, sig in enumerate(sig_homm):
        if sig:
            ymax = max(results.loc[i, "mean_pos"], results.loc[i, "mean_neg"])
            ax.text(i, ymax * 1.1 if ymax > 0 else ymax - abs(ymax)*0.1,
                    "*", ha="center", va="bottom", fontsize=14, color="black")

    ax.axhline(0, color="black", linewidth=0.8, linestyle="--")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=9, rotation=30 if n <= 8 else 0, ha="right" if n <= 8 else "center")
    ax.set_ylabel("Mean beta coefficient (ΔHbO, a.u.)")
    ax.set_title(f"{group_pos} vs {group_neg} — mean beta per target\n"
                 f"(* = Hommel FWER < {alpha}, ANCOVA-adjusted)")

    # Legend patches
    sig_patch   = mpatches.Patch(color="#d62728",  label=f"{group_pos} (sig.)")
    nsig_patch  = mpatches.Patch(color="#ff9896",  label=f"{group_pos} (n.s.)")
    sig_patch2  = mpatches.Patch(color="#1f77b4",  label=f"{group_neg} (sig.)")
    nsig_patch2 = mpatches.Patch(color="#aec7e8",  label=f"{group_neg} (n.s.)")
    ax.legend(handles=[sig_patch, nsig_patch, sig_patch2, nsig_patch2],
              fontsize=8, loc="upper right")

    plt.tight_layout()
    out = os.path.join(out_dir, "group_means_per_channel.png")
    fig.savefig(out, dpi=150)
    plt.close(fig)


def _plot_correction_heatmap(results, alpha, out_dir, labels):
    """Heatmap: −log10(corrected p) for all correction methods (ANCOVA only).
    Any panel whose column is entirely NaN (e.g. cluster-perm under --no-cluster)
    is omitted so it doesn't appear as spuriously significant."""
    thresh = -np.log10(alpha)
    candidate_panels = [
        ("ancova_p_fdr",      "ANCOVA — FDR-BH"),
        ("ancova_p_hoch",     "ANCOVA — Hochberg"),
        ("ancova_p_homm",     "ANCOVA — Hommel"),
        ("ancova_p_perm_max", "ANCOVA — Perm max-stat"),
        ("ancova_p_cluster",  "ANCOVA — Cluster perm"),
        ("mw_p_fdr",          "Mann-Whitney — FDR-BH"),
        ("mw_p_homm",         "Mann-Whitney — Hommel"),
    ]
    panels = [(c, lab) for c, lab in candidate_panels
              if c in results.columns and results[c].notna().any()]
    if not panels:
        return  # nothing to plot

    mat = np.vstack([
        np.clip(-np.log10(np.where(results[col].values > 0, results[col].values, 1e-10)), 0, 5)
        for col, _ in panels
    ])
    row_labels = [label for _, label in panels]
    n = len(labels)

    fig, ax = plt.subplots(figsize=(max(6, 0.85 * n + 2), 3.5))
    im = ax.imshow(mat, cmap="YlOrRd", vmin=0, vmax=5, aspect="auto")
    ax.set_xticks(range(n))
    ax.set_xticklabels(labels, fontsize=9, rotation=30 if n <= 8 else 0, ha="right" if n <= 8 else "center")
    ax.set_yticks(range(len(row_labels)))
    ax.set_yticklabels(row_labels, fontsize=8)
    plt.colorbar(im, ax=ax, label="−log₁₀(corrected p)", shrink=0.7)
    ax.set_title(f"Correction comparison  |  blue border = significant at α={alpha}")

    for row in range(len(row_labels)):
        for col in range(n):
            if mat[row, col] >= thresh:
                rect = plt.Rectangle((col - 0.5, row - 0.5), 1, 1,
                                     fill=False, edgecolor="blue", linewidth=2)
                ax.add_patch(rect)

    plt.tight_layout()
    fig.savefig(os.path.join(out_dir, "correction_heatmap.png"), dpi=150)
    plt.close(fig)


def _plot_volcano(results, alpha, out_dir, labels, sig_col="ancova_p_homm_sig"):
    """Volcano plot: ANCOVA beta (x) vs −log10(Hommel p) (y)."""
    betas = results["ancova_beta"].values.astype(float)
    pvals = results["ancova_p_homm"].values.astype(float)
    nl_p  = -np.log10(np.where(pvals > 0, pvals, 1e-10))
    sigs  = results[sig_col].values

    fig, ax = plt.subplots(figsize=(7, 6))
    colors = np.where(sigs, "#d62728", "#888888")

    ax.scatter(betas, nl_p, c=colors, s=80, edgecolors="black", linewidths=0.5, zorder=3)
    ax.axhline(-np.log10(alpha), color="blue", linestyle="--", linewidth=1)
    ax.axvline(0, color="black", linestyle="--", linewidth=0.8)

    for i, (b, p, ch) in enumerate(zip(betas, nl_p, labels)):
        if not np.isnan(b):
            ax.text(b, p + 0.05, ch, fontsize=7, ha="center", va="bottom")

    ax.set_xlabel("ANCOVA β  (adjusted)\n[ΔHbO, a.u.]")
    ax.set_ylabel("−log₁₀(Hommel corrected p)")
    ax.set_title("Volcano plot — Hommel FWER correction")

    sig_patch = mpatches.Patch(color="#d62728", label="Hommel significant")
    ns_patch  = mpatches.Patch(color="#888888", label="Not significant")
    ax.legend(handles=[sig_patch, ns_patch,
                       plt.Line2D([0], [0], color="blue", linestyle="--",
                                  label=f"α = {alpha}")],
              fontsize=8)

    plt.tight_layout()
    fig.savefig(os.path.join(out_dir, "volcano_plot.png"), dpi=150)
    plt.close(fig)


# ──────────────────────────────────────────────
if __name__ == "__main__":
    main()
