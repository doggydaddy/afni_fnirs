#!/usr/bin/env python3
"""
generate_report.py
==================
Reads all 6 comparison result CSVs produced by analyze_psd_agg_vs_psd_nagg.py
and generates two report files:

  • report.md   — Markdown (readable in VS Code / GitHub)
  • report.html — Styled HTML (open in browser)

Each report contains:
  - Overview table of significant findings across all 6 comparisons
  - Per-comparison section: narrative summary, cluster details,
    full ANCOVA channel table (all corrections), Mann-Whitney table
  - Cross-comparison discussion paragraph

Usage (from dumpFirstLevel_betaCoefs/):
  python generate_report.py [--base-dir .] [--output-dir .] [--alpha 0.05]
"""

import argparse
import datetime
import os
import numpy as np
import pandas as pd


# ──────────────────────────────────────────────────────────────────────────────
# Comparison manifest — edit result_dir paths here if you move folders
# ──────────────────────────────────────────────────────────────────────────────
COMPARISONS = [
    dict(label="PSD+AGG vs HC",      contrast="v5 (LRV hard-easy)",
         group_pos="PSD+AGG",  group_neg="HC",      covariates="SANS, SAPS, WAIS-matrix, SUD",
         result_dir="results_psd_agg_vs_hc_lrv_v5"),
    dict(label="PSD+AGG vs HC",      contrast="v7 (LRV corsi_GLT)",
         group_pos="PSD+AGG",  group_neg="HC",      covariates="SANS, SAPS, WAIS-matrix, SUD",
         result_dir="results_psd_agg_vs_hc_lrv_v7"),
    dict(label="PSD+AGG vs PSD-AGG", contrast="v5 (LRV hard-easy)",
         group_pos="PSD+AGG",  group_neg="PSD-AGG", covariates="SANS, SAPS, WAIS-matrix, SUD, DOSE",
         result_dir="results_psd_agg_lrv_v5"),
    dict(label="PSD+AGG vs PSD-AGG", contrast="v7 (LRV corsi_GLT)",
         group_pos="PSD+AGG",  group_neg="PSD-AGG", covariates="SANS, SAPS, WAIS-matrix, SUD, DOSE",
         result_dir="results_psd_agg_lrv_v7"),
    dict(label="PSD-AGG vs HC",      contrast="v5 (pcon hard-easy)",
         group_pos="PSD-AGG",  group_neg="HC",      covariates="SANS, SAPS, WAIS-matrix, SUD",
         result_dir="results_psd_nagg_vs_hc_pcon_v5"),
    dict(label="PSD-AGG vs HC",      contrast="v7 (pcon corsi_GLT)",
         group_pos="PSD-AGG",  group_neg="HC",      covariates="SANS, SAPS, WAIS-matrix, SUD",
         result_dir="results_psd_nagg_vs_hc_pcon_v7"),
]

# Mode-specific configuration. The region runner writes to result_dirs with the
# `_regions` suffix; the analyze script writes a `channel` column whose values
# match the original target column names (e.g. "left_dlPFC").
TARGET_LABEL_MAP_CHANNELS = {f"channel{i}": f"Ch{i}" for i in range(1, 17)}
TARGET_LABEL_MAP_REGIONS  = {
    "left_dlpfc":  "Left dlPFC",
    "left_mpfc":   "Left mPFC",
    "right_mpfc":  "Right mPFC",
    "right_dlpfc": "Right dlPFC",
}


def get_mode_config(mode):
    """Return (comparisons, label_map, report_filename_stem)."""
    if mode == "channels":
        return COMPARISONS, TARGET_LABEL_MAP_CHANNELS, "report"
    if mode == "regions":
        comps_r = [dict(c, result_dir=c["result_dir"] + "_regions") for c in COMPARISONS]
        return comps_r, TARGET_LABEL_MAP_REGIONS, "report_regions"
    raise ValueError(f"unknown mode: {mode}")

CSV_NAME   = "results_psd_agg_vs_psd_nagg.csv"
INT_CSV    = "interaction_results.csv"
ALPHA      = 0.05


# ──────────────────────────────────────────────────────────────────────────────
# Data helpers
# ──────────────────────────────────────────────────────────────────────────────
def _relabel(series, label_map):
    """Map raw target names → display labels; passthrough if unknown."""
    return series.map(lambda v: label_map.get(str(v).lower(), str(v)))


def load(comp, base_dir, label_map):
    path = os.path.join(base_dir, comp["result_dir"], CSV_NAME)
    if not os.path.exists(path):
        return None
    df = pd.read_csv(path)
    df["ch_label"] = _relabel(df["channel"], label_map)
    return df


def load_interactions(comp, base_dir, label_map):
    path = os.path.join(base_dir, comp["result_dir"], INT_CSV)
    if not os.path.exists(path):
        return None
    df = pd.read_csv(path)
    df["ch_label"] = _relabel(df["channel"], label_map)
    # Detect which covariates are present from column names
    cols = df.columns.tolist()
    covs = sorted(set(
        c.replace("int_", "").replace("_beta", "").replace("_t", "")
         .replace("_p_fdr", "").replace("_p", "").replace("_sig", "")
        for c in cols if c.startswith("int_")
    ))
    # Preferred display order; any not in this list go at the end alphabetically
    preferred = ["sans", "saps", "wais_matrix", "sud", "dose"]
    ordered = [c for c in preferred if c in covs] + sorted(c for c in covs if c not in preferred)
    return df, ordered


def extract_clusters(df):
    """
    Reconstruct clusters from the saved CSV.
    Channels sharing the same ancova_p_cluster value and ancova_p_cluster_sig=True
    belong to the same cluster.
    """
    sig = df[df["ancova_p_cluster_sig"] == True].copy()
    if sig.empty:
        return []
    clusters = []
    for pval, grp in sig.groupby(sig["ancova_p_cluster"].round(6)):
        mass = grp["ancova_t"].sum()
        clusters.append(dict(
            channels  = grp["ch_label"].tolist(),
            p         = float(pval),
            mass      = float(mass),
            direction = "↓" if mass < 0 else "↑",
        ))
    return sorted(clusters, key=lambda c: c["p"])


def sig_channels(df, col_sig):
    return df[df[col_sig] == True]["ch_label"].tolist()


# ──────────────────────────────────────────────────────────────────────────────
# Narrative generator
# ──────────────────────────────────────────────────────────────────────────────
def narrative(comp, df, clusters, alpha=ALPHA):
    gpos     = comp["group_pos"]
    gneg     = comp["group_neg"]
    n_pos    = int(df["n_pos"].iloc[0])
    n_neg    = int(df["n_neg"].iloc[0])
    covs     = comp["covariates"]
    contrast = comp["contrast"]

    best     = df.loc[df["ancova_p"].idxmin()]
    best_ch  = best["ch_label"]
    best_p   = best["ancova_p"]
    best_t   = best["ancova_t"]
    best_d   = best["cohens_d"]

    fdr_chs  = sig_channels(df, "ancova_p_fdr_sig")
    pmax_chs = sig_channels(df, "ancova_p_perm_max_sig")

    parts = [
        f"Comparison of {gpos} (n={n_pos}) vs {gneg} (n={n_neg}) on the "
        f"{contrast} contrast, with {covs} as covariates."
    ]

    if clusters:
        for cl in clusters:
            ch_str  = ", ".join(cl["channels"])
            dir_str = "reduced" if cl["direction"] == "↓" else "elevated"
            parts.append(
                f"A spatially contiguous cluster comprising {ch_str} reached "
                f"significance under cluster-based permutation FWER correction "
                f"(cluster mass = {cl['mass']:.2f}, p = {cl['p']:.4f}), "
                f"reflecting {dir_str} ΔHbO in {gpos} relative to {gneg}."
            )
    else:
        parts.append(
            f"No channel cluster survived FWER correction. "
            f"The strongest single-channel effect was {best_ch} "
            f"(t = {best_t:.2f}, p_raw = {best_p:.4f}, d = {best_d:.2f}), uncorrected."
        )

    if pmax_chs:
        parts.append(
            f"{', '.join(pmax_chs)} also survived permutation max-statistic correction (FWER, p < {alpha})."
        )

    if fdr_chs:
        parts.append(
            f"{', '.join(fdr_chs)} survived FDR-BH correction (p < {alpha})."
        )

    # Effect direction summary
    neg_betas = (df["ancova_beta"] < 0).sum()
    if neg_betas > 12:
        parts.append(
            f"The majority of channels ({neg_betas}/16) showed negative β values, "
            f"indicating a broad pattern of lower ΔHbO in {gpos} than {gneg}."
        )

    return " ".join(parts)


def interaction_narrative(int_df, covariates, gpos, gneg, alpha=ALPHA):
    """One-paragraph summary of interaction effects."""
    sig_findings = []
    for cov in covariates:
        sig_chs = int_df[int_df[f"int_{cov}_sig"] == True]["ch_label"].tolist()
        if sig_chs:
            t_vals = int_df[int_df[f"int_{cov}_sig"] == True][f"int_{cov}_t"].tolist()
            p_vals = int_df[int_df[f"int_{cov}_sig"] == True][f"int_{cov}_p_fdr"].tolist()
            details = "; ".join(
                f"{ch} (t={t:.2f}, FDR p={p:.3f})"
                for ch, t, p in zip(sig_chs, t_vals, p_vals)
            )
            sig_findings.append(
                f"group × {cov.upper()}: {details}"
            )
    if not sig_findings:
        return (
            f"No significant group × covariate interaction effects were found after "
            f"FDR-BH correction (α={alpha}), suggesting that the relationships between "
            f"{', '.join(c.upper() for c in covariates)} and ΔHbO do not differ "
            f"significantly between {gpos} and {gneg} in any channel."
        )
    return (
        f"Secondary interaction analysis revealed the following significant "
        f"group × covariate effects (FDR-BH corrected, α={alpha}): "
        + " | ".join(sig_findings) + ". "
        "This indicates that the modulating effect of these clinical variables on "
        "cortical haemodynamic response differs between groups."
    )


def cross_comparison_narrative(all_data, alpha=ALPHA):
    """Auto-generate a brief cross-comparison discussion (purely data-driven)."""
    lines = []
    for comp, df, clusters in all_data:
        if clusters:
            for cl in clusters:
                chs = "+".join(cl["channels"])
                lines.append(
                    f"- **{comp['label']}** ({comp['contrast']}): cluster [{chs}], "
                    f"p={cl['p']:.3f}, mass={cl['mass']:.2f} {cl['direction']}"
                )
    if not lines:
        return "No cluster-corrected significant findings were observed across any of the six comparisons."

    intro = (
        "Across the six comparisons, cluster-based permutation correction (FWER) identified "
        "spatially contiguous target clusters in the following contrasts:"
    )
    return "\n\n".join([intro, "\n".join(lines)])


# ──────────────────────────────────────────────────────────────────────────────
# Mode-dependent boilerplate for the report header / methods note
# ──────────────────────────────────────────────────────────────────────────────
def _mode_methods_note_md(mode):
    if mode == "regions":
        return (
            "> **Methods note (region analysis).** Channels are pooled into four "
            "anatomically pre-specified regions of interest (ROIs) by averaging: "
            "**Left dlPFC** = mean(Ch1–4), **Left mPFC** = mean(Ch5–8), "
            "**Right mPFC** = mean(Ch9–12), **Right dlPFC** = mean(Ch13–16). "
            "Because these four regions are *a priori* hypotheses and the number of "
            "tests is small (k=4), cluster-based permutation correction is **not** "
            "applied here — it is largely redundant with permutation max-statistic "
            "at this resolution and contributes no additional FWER protection. "
            "Multiple-comparison control is provided by FDR-BH, Simes-Hochberg, "
            "Hommel, and permutation max-statistic. Raw p-values are reported "
            "alongside corrected p-values; under an ROI-confirmatory framing they "
            "may be interpreted as the primary statistic for each region."
        )
    return ""


def _mode_corrections_line(mode):
    if mode == "regions":
        return ("**Corrections:** FDR-BH · Simes-Hochberg · Hommel · "
                "Permutation max-statistic (5,000 perms). "
                "*Cluster-based permutation omitted for ROI analysis — see Methods note below.*  ")
    return ("**Corrections:** FDR-BH · Simes-Hochberg · Hommel · "
            "Permutation max-statistic · Cluster-based permutation "
            "(5,000 perms, linear adjacency)  ")


# ──────────────────────────────────────────────────────────────────────────────
# Markdown builder
# ──────────────────────────────────────────────────────────────────────────────
def build_markdown(all_data, alpha, base_dir, label_map, mode="channels"):
    today = datetime.date.today().isoformat()
    L = []

    L += [
        "# fNIRS Second-Level Group Analysis — Report",
        "",
        f"**Generated:** {today}  ",
        f"**Alpha:** {alpha}  ",
        _mode_corrections_line(mode),
        "**Effect sizes:** Cohen's *d* (ANCOVA) · rank-biserial *r* (Mann-Whitney)  ",
        "",
        "---",
        "",
    ]

    methods_note = _mode_methods_note_md(mode)
    if methods_note:
        L += [methods_note, "", "---", ""]

    # Overview table — drop the "Sig. cluster" column in regions mode
    if mode == "regions":
        L += [
            "## Overview of Significant Findings",
            "",
            "| # | Comparison | Contrast | Perm max-stat ROIs | FDR-BH ROIs |",
            "|---|---|---|---|---|",
        ]
        for i, (comp, df, clusters) in enumerate(all_data, 1):
            pmax = sig_channels(df, "ancova_p_perm_max_sig")
            fdr  = sig_channels(df, "ancova_p_fdr_sig")
            L.append(
                f"| {i} | **{comp['label']}** | {comp['contrast']} "
                f"| {', '.join(pmax) or '—'} | {', '.join(fdr) or '—'} |"
            )
    else:
        L += [
            "## Overview of Significant Findings",
            "",
            "| # | Comparison | Contrast | Sig. cluster (p) | Perm max-stat chs | FDR-BH chs |",
            "|---|---|---|---|---|---|",
        ]
        for i, (comp, df, clusters) in enumerate(all_data, 1):
            clust_str = " ; ".join(
                f"{'+'.join(c['channels'])} (p={c['p']:.3f})" for c in clusters
            ) or "—"
            pmax = sig_channels(df, "ancova_p_perm_max_sig")
            fdr  = sig_channels(df, "ancova_p_fdr_sig")
            L.append(
                f"| {i} | **{comp['label']}** | {comp['contrast']} "
                f"| {clust_str} | {', '.join(pmax) or '—'} | {', '.join(fdr) or '—'} |"
            )

    L += ["", "---", "", "## Cross-Comparison Discussion", ""]
    L.append(cross_comparison_narrative(all_data, alpha))
    L += ["", "---", ""]

    for i, (comp, df, clusters) in enumerate(all_data, 1):
        L += [
            f"## {i}. {comp['label']} — {comp['contrast']}",
            "",
            f"> {narrative(comp, df, clusters, alpha)}",
            "",
        ]

        # Cluster box (only meaningful when cluster-perm was actually run)
        if clusters and mode != "regions":
            L += ["### Cluster Summary", "",
                  "| Channels | Mass | Direction | p (cluster-perm) |",
                  "|---|---|---|---|"]
            for cl in clusters:
                L.append(
                    f"| {', '.join(cl['channels'])} | {cl['mass']:.3f} "
                    f"| {cl['direction']} | **{cl['p']:.4f}** |"
                )
            L.append("")

        # ANCOVA table — narrower in regions mode (no cluster column)
        if mode == "regions":
            L += [
                "### ANCOVA Region Results",
                "",
                "| Region | β (ΔHbO) | t | d | p raw | FDR-BH | Hochberg | Hommel | Perm-max |",
                "|---|---|---|---|---|---|---|---|---|",
            ]
        else:
            L += [
                "### ANCOVA Channel Results",
                "",
                "| Ch | β (ΔHbO) | t | d | p raw | FDR-BH | Hochberg | Hommel | Perm-max | Cluster |",
                "|---|---|---|---|---|---|---|---|---|---|",
            ]
        for _, r in df.iterrows():
            def fp(col, sig=None):
                v = r[col]
                if pd.isna(v): return "—"
                s = f"{v:.4f}"
                if sig and r.get(sig, False): s = f"**{s}\\***"
                return s
            row_md = (
                f"| {r['ch_label']} "
                f"| {r['ancova_beta']:.2e} "
                f"| {r['ancova_t']:.2f} "
                f"| {r['cohens_d']:.2f} "
                f"| {fp('ancova_p')} "
                f"| {fp('ancova_p_fdr',      'ancova_p_fdr_sig')} "
                f"| {fp('ancova_p_hoch',     'ancova_p_hoch_sig')} "
                f"| {fp('ancova_p_homm',     'ancova_p_homm_sig')} "
                f"| {fp('ancova_p_perm_max', 'ancova_p_perm_max_sig')} "
            )
            if mode == "regions":
                row_md += "|"
            else:
                row_md += f"| {fp('ancova_p_cluster',  'ancova_p_cluster_sig')} |"
            L.append(row_md)

        # Mann-Whitney table
        L += [
            "",
            "### Mann-Whitney (non-parametric, unadjusted for covariates)",
            "",
            "| " + ("Region" if mode == "regions" else "Ch") + " | p raw | FDR-BH | Hommel | rb *r* |",
            "|---|---|---|---|---|",
        ]
        for _, r in df.iterrows():
            def fp(col, sig=None):
                v = r[col]
                if pd.isna(v): return "—"
                s = f"{v:.4f}"
                if sig and r.get(sig, False): s = f"**{s}\\***"
                return s
            L.append(
                f"| {r['ch_label']} "
                f"| {fp('mw_p')} "
                f"| {fp('mw_p_fdr',  'mw_p_fdr_sig')} "
                f"| {fp('mw_p_homm', 'mw_p_homm_sig')} "
                f"| {r['rb_r']:.3f} |"
            )

        # Interaction section (if CSV exists)
        int_result = load_interactions(comp, base_dir, label_map)
        if int_result is not None:
            int_df, int_covs = int_result
            int_df["ch_label"] = _relabel(int_df["channel"], label_map)
            L += [
                "",
                "### Secondary Analysis: Group × Covariate Interaction Effects",
                "",
                f"> {interaction_narrative(int_df, int_covs, comp['group_pos'], comp['group_neg'], alpha)}",
                "",
            ]
            for cov in int_covs:
                L += [f"#### group × {cov.upper()}", "",
                      f"| Ch | β (interaction) | t | p raw | p FDR-BH | sig |",
                      "|---|---|---|---|---|---|"]
                for _, r in int_df.iterrows():
                    def fp(col):
                        v = r[col]; return f"{v:.4f}" if pd.notna(v) else "—"
                    sig = r.get(f"int_{cov}_sig", False)
                    beta_str = f"{r[f'int_{cov}_beta']:.3e}" if pd.notna(r[f"int_{cov}_beta"]) else "—"
                    t_str    = f"{r[f'int_{cov}_t']:.2f}"   if pd.notna(r[f"int_{cov}_t"])    else "—"
                    p_fdr    = fp(f"int_{cov}_p_fdr")
                    if sig: p_fdr = f"**{p_fdr}\\***"
                    L.append(
                        f"| {r['ch_label']} | {beta_str} | {t_str} "
                        f"| {fp(f'int_{cov}_p')} | {p_fdr} | {'YES ✓' if sig else '---'} |"
                    )
                L.append("")

        L += ["", "---", ""]

    return "\n".join(L)


# ──────────────────────────────────────────────────────────────────────────────
# HTML builder
# ──────────────────────────────────────────────────────────────────────────────
CSS = """
* { box-sizing: border-box; }
body { font-family: 'Segoe UI', Arial, sans-serif; max-width: 1300px;
       margin: 32px auto; padding: 0 24px; color: #2c2c2c; line-height: 1.6; }
h1 { font-size: 1.8em; border-bottom: 3px solid #c0392b; padding-bottom: 8px; color: #1a1a2e; }
h2 { font-size: 1.3em; margin-top: 44px; padding: 6px 14px;
     border-left: 5px solid #c0392b; background: #fafafa; color: #16213e; }
h3 { font-size: 1.05em; color: #0f3460; margin-top: 20px; }
table { border-collapse: collapse; width: 100%; font-size: 12.5px; margin: 12px 0; }
th { background: #16213e; color: #fff; padding: 7px 10px; text-align: left; white-space: nowrap; }
td { padding: 5px 10px; border-bottom: 1px solid #e0e0e0; white-space: nowrap; }
tr:nth-child(even) td { background: #f7f8fb; }
tr:hover td { background: #edf0ff; }
.neg { color: #c0392b; font-weight: 600; }
.pos { color: #1e8449; font-weight: 600; }
.sig-cluster td { background: #fce8e8 !important; font-weight: bold; }
.sig-cluster td.p-val { color: #922b21; }
.sig-pmax   td { background: #fef5e7 !important; }
.sig-pmax   td.p-val { color: #7d6608; font-weight: bold; }
.sig-fdr    td { background: #eaf6fb !important; }
.sig-fdr    td.p-val { color: #1a5276; font-weight: bold; }
.p-sig  { font-weight: bold; }
.p-cluster { color: #922b21; font-weight: bold; }
.p-pmax    { color: #7d6608; font-weight: bold; }
.p-fdr     { color: #1a5276; font-weight: bold; }
.narrative { background: #f0f5ff; border-left: 4px solid #2980b9;
             padding: 12px 18px; border-radius: 4px; margin: 12px 0;
             font-size: 13.5px; }
.discussion { background: #f9fbe7; border-left: 4px solid #7dbb00;
              padding: 12px 18px; border-radius: 4px; margin: 12px 0; }
.cluster-badge { display: inline-block; background: #c0392b; color: #fff;
                 border-radius: 3px; padding: 1px 6px; font-size: 11px;
                 font-weight: bold; margin-right: 4px; }
.none { color: #aaa; font-style: italic; }
footer { margin-top: 40px; padding-top: 10px; border-top: 1px solid #ddd;
         font-size: 11px; color: #999; }
"""

def _td_p(val, sig, cls_name):
    if pd.isna(val) or val > 1.5:
        return "<td>—</td>"
    s = f"{val:.4f}"
    if sig:
        return f'<td class="p-val {cls_name}">{s} ✓</td>'
    return f"<td>{s}</td>"


def build_html(all_data, alpha, base_dir, label_map, mode="channels"):
    today = datetime.date.today().isoformat()
    B = []

    target_word = "ROIs" if mode == "regions" else "chs"
    target_word_singular = "Region" if mode == "regions" else "Ch"

    if mode == "regions":
        corr_html = ("FDR-BH · Simes-Hochberg · Hommel (FWER) · "
                     "Permutation max-statistic (FWER). "
                     "<em>Cluster-perm omitted for ROI analysis.</em>")
    else:
        corr_html = ("FDR-BH · Simes-Hochberg · Hommel (FWER) · "
                     "Permutation max-statistic (FWER) · Cluster-based permutation (FWER)")

    B += [
        f"<h1>fNIRS Second-Level Group Analysis — Report</h1>",
        f"<p><b>Generated:</b> {today} &nbsp;|&nbsp; <b>Alpha:</b> {alpha} &nbsp;|&nbsp; "
        f"<b>Permutations:</b> 5,000</p>",
        f"<p><b>Corrections applied:</b> {corr_html}</p>",
        "<hr>",
    ]

    if mode == "regions":
        B += [
            '<div class="discussion"><b>Methods note (region analysis).</b> '
            "Channels are pooled into four anatomically pre-specified regions of "
            "interest (ROIs) by averaging: "
            "<b>Left dlPFC</b> = mean(Ch1–4), <b>Left mPFC</b> = mean(Ch5–8), "
            "<b>Right mPFC</b> = mean(Ch9–12), <b>Right dlPFC</b> = mean(Ch13–16). "
            "Because these four regions are <em>a priori</em> hypotheses and the "
            "number of tests is small (k=4), cluster-based permutation correction "
            "is <b>not</b> applied here — it is largely redundant with permutation "
            "max-statistic at this resolution and contributes no additional FWER "
            "protection. Raw p-values are reported alongside corrected p-values; "
            "under an ROI-confirmatory framing they may be interpreted as the "
            "primary statistic for each region.</div>",
            "<hr>",
        ]

    if mode == "regions":
        B += [
            "<h2>Overview of Significant Findings</h2>",
            "<table><tr><th>#</th><th>Comparison</th><th>Contrast</th>"
            f"<th>Perm max-stat {target_word}</th><th>FDR-BH {target_word}</th></tr>",
        ]
        for i, (comp, df, clusters) in enumerate(all_data, 1):
            pmax = sig_channels(df, "ancova_p_perm_max_sig")
            fdr  = sig_channels(df, "ancova_p_fdr_sig")
            B.append(
                f"<tr><td>{i}</td><td><b>{comp['label']}</b></td><td>{comp['contrast']}</td>"
                f"<td>{'<b>' + ', '.join(pmax) + '</b>' if pmax else '<span class=\"none\">—</span>'}</td>"
                f"<td>{', '.join(fdr) if fdr else '<span class=\"none\">—</span>'}</td></tr>"
            )
    else:
        B += [
            "<h2>Overview of Significant Findings</h2>",
            "<table><tr><th>#</th><th>Comparison</th><th>Contrast</th>"
            f"<th>Significant cluster</th><th>Perm max-stat {target_word}</th><th>FDR-BH {target_word}</th></tr>",
        ]
        for i, (comp, df, clusters) in enumerate(all_data, 1):
            clust_td = " ".join(
                f'<span class="cluster-badge">{"+".join(c["channels"])}: p={c["p"]:.3f}</span>'
                for c in clusters
            ) or '<span class="none">—</span>'
            pmax = sig_channels(df, "ancova_p_perm_max_sig")
            fdr  = sig_channels(df, "ancova_p_fdr_sig")
            B.append(
                f"<tr><td>{i}</td><td><b>{comp['label']}</b></td><td>{comp['contrast']}</td>"
                f"<td>{clust_td}</td>"
                f"<td>{'<b>' + ', '.join(pmax) + '</b>' if pmax else '<span class=\"none\">—</span>'}</td>"
                f"<td>{', '.join(fdr) if fdr else '<span class=\"none\">—</span>'}</td></tr>"
            )
    B += ["</table>", "<hr>",
          "<h2>Cross-Comparison Discussion</h2>",
          f'<div class="discussion">{cross_comparison_narrative(all_data, alpha).replace(chr(10), "<br>")}</div>',
          "<hr>"]

    for i, (comp, df, clusters) in enumerate(all_data, 1):
        B += [
            f"<h2>{i}. {comp['label']} &mdash; {comp['contrast']}</h2>",
            f'<div class="narrative">{narrative(comp, df, clusters, alpha)}</div>',
        ]

        if clusters and mode != "regions":
            B += ["<h3>Cluster Summary</h3>",
                  "<table><tr><th>Channels</th><th>Cluster mass</th><th>Direction</th><th>p (cluster-perm)</th></tr>"]
            for cl in clusters:
                B.append(
                    f"<tr class='sig-cluster'>"
                    f"<td>{', '.join(cl['channels'])}</td>"
                    f"<td>{cl['mass']:.3f}</td>"
                    f"<td>{cl['direction']}</td>"
                    f"<td class='p-cluster'>{cl['p']:.4f}</td></tr>"
                )
            B.append("</table>")

        # ANCOVA table — narrower in regions mode (no cluster column)
        ancova_title = "ANCOVA Region Results" if mode == "regions" else "ANCOVA Channel Results"
        ancova_hdr = (
            f"<tr><th>{target_word_singular}</th><th>β (ΔHbO)</th><th>t</th><th>d</th>"
            f"<th>p raw</th><th>FDR-BH</th><th>Hochberg</th><th>Hommel</th>"
            f"<th>Perm max-stat</th>"
        )
        if mode != "regions":
            ancova_hdr += "<th>Cluster-perm</th>"
        ancova_hdr += "</tr>"
        B += [f"<h3>{ancova_title}</h3>", "<table>", ancova_hdr]
        for _, r in df.iterrows():
            row_cls = ""
            if mode != "regions" and r.get("ancova_p_cluster_sig", False):
                row_cls = "sig-cluster"
            elif r["ancova_p_perm_max_sig"]:   row_cls = "sig-pmax"
            elif r["ancova_p_fdr_sig"]:        row_cls = "sig-fdr"
            bc = "neg" if r["ancova_beta"] < 0 else "pos"
            tc = "neg" if r["ancova_t"]    < 0 else "pos"
            dc = "neg" if r["cohens_d"]    < 0 else "pos"
            row_html = (
                f"<tr class='{row_cls}'>"
                f"<td><b>{r['ch_label']}</b></td>"
                f"<td class='{bc}'>{r['ancova_beta']:.2e}</td>"
                f"<td class='{tc}'>{r['ancova_t']:.2f}</td>"
                f"<td class='{dc}'>{r['cohens_d']:.2f}</td>"
                + _td_p(r["ancova_p"],          False,                       "")
                + _td_p(r["ancova_p_fdr"],       r["ancova_p_fdr_sig"],       "p-fdr")
                + _td_p(r["ancova_p_hoch"],      r["ancova_p_hoch_sig"],      "p-pmax")
                + _td_p(r["ancova_p_homm"],      r["ancova_p_homm_sig"],      "p-pmax")
                + _td_p(r["ancova_p_perm_max"],  r["ancova_p_perm_max_sig"],  "p-pmax")
            )
            if mode != "regions":
                row_html += _td_p(r["ancova_p_cluster"], r["ancova_p_cluster_sig"], "p-cluster")
            row_html += "</tr>"
            B.append(row_html)
        B.append("</table>")

        # Mann-Whitney table
        B += [f"<h3>Mann-Whitney (non-parametric, unadjusted for covariates)</h3>",
              f"<table><tr><th>{target_word_singular}</th><th>p raw</th><th>FDR-BH</th><th>Hommel</th><th>rb <em>r</em></th></tr>"]
        for _, r in df.iterrows():
            B.append(
                f"<tr><td><b>{r['ch_label']}</b></td>"
                + _td_p(r["mw_p"],       False,               "")
                + _td_p(r["mw_p_fdr"],   r["mw_p_fdr_sig"],   "p-fdr")
                + _td_p(r["mw_p_homm"],  r["mw_p_homm_sig"],  "p-pmax")
                + f"<td>{r['rb_r']:.3f}</td></tr>"
            )
        B += ["</table>"]

        # Interaction section (if CSV exists)
        int_result = load_interactions(comp, base_dir, label_map)
        if int_result is not None:
            int_df, int_covs = int_result
            int_df["ch_label"] = _relabel(int_df["channel"], label_map)
            narr = interaction_narrative(int_df, int_covs, comp["group_pos"], comp["group_neg"], alpha)
            B += [
                "<h3>Secondary Analysis: Group &times; Covariate Interaction Effects</h3>",
                f'<div class="narrative">{narr}</div>',
            ]
            for cov in int_covs:
                any_sig = int_df[f"int_{cov}_sig"].any()
                B += [
                    f"<h4>group &times; {cov.upper()}"
                    + (" <span style='color:#c0392b'>← significant</span>" if any_sig else "")
                    + "</h4>",
                    "<table>",
                    f"<tr><th>Ch</th><th>β (interaction)</th><th>t</th>"
                    f"<th>p raw</th><th>p FDR-BH</th><th>sig</th></tr>",
                ]
                for _, r in int_df.iterrows():
                    sig = bool(r.get(f"int_{cov}_sig", False))
                    rc  = "sig-cluster" if sig else ""
                    def fv(c): return f"{r[c]:.4f}" if pd.notna(r[c]) else "—"
                    beta_str = f"{r[f'int_{cov}_beta']:.3e}" if pd.notna(r[f"int_{cov}_beta"]) else "—"
                    t_str    = f"{r[f'int_{cov}_t']:.2f}"   if pd.notna(r[f"int_{cov}_t"])    else "—"
                    p_fdr    = fv(f"int_{cov}_p_fdr")
                    p_raw    = fv(f"int_{cov}_p")
                    sig_str  = "✓ YES" if sig else "—"
                    B.append(
                        f"<tr class='{rc}'><td><b>{r['ch_label']}</b></td>"
                        f"<td>{beta_str}</td><td>{t_str}</td>"
                        f"<td>{p_raw}</td>"
                        f"<td class=\"{'p-cluster' if sig else ''}\">{p_fdr}</td>"
                        f"<td>{sig_str}</td></tr>"
                    )
                B += ["</table>"]

        B += ["<hr>"]

    B.append(
        f"<footer>Report generated by <code>generate_report.py</code> on {today}. "
        "Permutation tests: 5,000 permutations. Channel adjacency: linear chain Ch1–Ch16. "
        "ANCOVA: covariate-adjusted OLS (Frisch-Waugh). "
        "Cluster forming threshold: uncorrected p &lt; 0.05.</footer>"
    )

    return (
        f'<!DOCTYPE html><html lang="en"><head><meta charset="UTF-8">'
        f'<title>fNIRS Report {today}</title>'
        f"<style>{CSS}</style></head><body>"
        + "".join(B)
        + "</body></html>"
    )


# ──────────────────────────────────────────────────────────────────────────────
# Entry point
# ──────────────────────────────────────────────────────────────────────────────
def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--base-dir",    default=".",
                   help="Folder containing all results_* directories (default: .)")
    p.add_argument("--output-dir",  default=".",
                   help="Where to write report.md and report.html (default: .)")
    p.add_argument("--alpha",       type=float, default=ALPHA)
    p.add_argument("--mode",        choices=["channels", "regions"], default="channels",
                   help="Which comparison flavor to report on. "
                        "'channels' reads results_*/, 'regions' reads results_*_regions/. "
                        "Report filenames are auto-suffixed in regions mode.")
    return p.parse_args()


def main():
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    comparisons, label_map, stem = get_mode_config(args.mode)

    # Load all comparisons
    all_data = []
    for comp in comparisons:
        df = load(comp, args.base_dir, label_map)
        if df is None:
            print(f"  [WARN] Missing results for: {comp['result_dir']}")
            continue
        clusters = extract_clusters(df)
        all_data.append((comp, df, clusters))

    print(f"Loaded {len(all_data)}/{len(comparisons)} comparisons  (mode={args.mode}).")

    md   = build_markdown(all_data, args.alpha, args.base_dir, label_map, mode=args.mode)
    html = build_html(all_data, args.alpha, args.base_dir, label_map, mode=args.mode)

    md_path   = os.path.join(args.output_dir, f"{stem}.md")
    html_path = os.path.join(args.output_dir, f"{stem}.html")

    with open(md_path,   "w", encoding="utf-8") as f:
        f.write(md)
    with open(html_path, "w", encoding="utf-8") as f:
        f.write(html)

    print(f"[✓] Markdown → {md_path}")
    print(f"[✓] HTML     → {html_path}")


if __name__ == "__main__":
    main()
