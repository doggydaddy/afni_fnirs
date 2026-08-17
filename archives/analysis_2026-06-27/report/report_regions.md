# fNIRS Second-Level Group Analysis — Report

**Generated:** 2026-06-27  
**Alpha:** 0.05  
**Corrections:** FDR-BH · Simes-Hochberg · Hommel · Permutation max-statistic (5,000 perms). *Cluster-based permutation omitted for ROI analysis — see Methods note below.*  
**Effect sizes:** Cohen's *d* (ANCOVA) · rank-biserial *r* (Mann-Whitney)  

---

> **Methods note (region analysis).** Channels are pooled into four anatomically pre-specified regions of interest (ROIs) by averaging: **Left dlPFC** = mean(Ch1–4), **Left mPFC** = mean(Ch5–8), **Right mPFC** = mean(Ch9–12), **Right dlPFC** = mean(Ch13–16). Because these four regions are *a priori* hypotheses and the number of tests is small (k=4), cluster-based permutation correction is **not** applied here — it is largely redundant with permutation max-statistic at this resolution and contributes no additional FWER protection. Multiple-comparison control is provided by FDR-BH, Simes-Hochberg, Hommel, and permutation max-statistic. Raw p-values are reported alongside corrected p-values; under an ROI-confirmatory framing they may be interpreted as the primary statistic for each region.

---

## Overview of Significant Findings

| # | Comparison | Contrast | Perm max-stat ROIs | FDR-BH ROIs |
|---|---|---|---|---|
| 1 | **PSD+AGG vs HC** | v5 (LRV hard-easy) | — | — |
| 2 | **PSD+AGG vs HC** | v7 (LRV corsi_GLT) | — | — |
| 3 | **PSD+AGG vs PSD-AGG** | v5 (LRV hard-easy) | — | — |
| 4 | **PSD+AGG vs PSD-AGG** | v7 (LRV corsi_GLT) | — | — |
| 5 | **PSD-AGG vs HC** | v5 (pcon hard-easy) | — | — |
| 6 | **PSD-AGG vs HC** | v7 (pcon corsi_GLT) | — | — |

---

## Cross-Comparison Discussion

No cluster-corrected significant findings were observed across any of the six comparisons.

---

## 1. PSD+AGG vs HC — v5 (LRV hard-easy)

> Comparison of PSD+AGG (n=77) vs HC (n=85) on the v5 (LRV hard-easy) contrast, with SANS, SAPS, WAIS-matrix, SUD as covariates. No channel cluster survived FWER correction. The strongest single-channel effect was Left dlPFC (t = -2.32, p_raw = 0.0214, d = -0.29), uncorrected.

### ANCOVA Region Results

| Region | β (ΔHbO) | t | d | p raw | FDR-BH | Hochberg | Hommel | Perm-max |
|---|---|---|---|---|---|---|---|---|
| Left dlPFC | -3.78e-05 | -2.32 | -0.29 | 0.0214 | 0.0856 | 0.0856 | 0.0856 | 0.0646 |
| Left mPFC | -1.70e-05 | -1.17 | -0.18 | 0.2424 | 0.4847 | 0.7271 | 0.7271 | 0.5714 |
| Right mPFC | -2.81e-06 | -0.25 | -0.10 | 0.8057 | 0.8057 | 0.8057 | 0.8057 | 0.9974 |
| Right dlPFC | -7.98e-06 | -0.69 | -0.12 | 0.4912 | 0.6550 | 0.8057 | 0.8057 | 0.8954 |

### Mann-Whitney (non-parametric, unadjusted for covariates)

| Region | p raw | FDR-BH | Hommel | rb *r* |
|---|---|---|---|---|
| Left dlPFC | 0.2961 | 0.4750 | 0.4750 | 0.095 |
| Left mPFC | 0.2854 | 0.4750 | 0.4750 | 0.097 |
| Right mPFC | 0.4257 | 0.4750 | 0.4750 | 0.073 |
| Right dlPFC | 0.4750 | 0.4750 | 0.4750 | 0.065 |

### Secondary Analysis: Group × Covariate Interaction Effects

> No significant group × covariate interaction effects were found after FDR-BH correction (α=0.05), suggesting that the relationships between SANS, SAPS, WAIS_MATRIX, SUD and ΔHbO do not differ significantly between PSD+AGG and HC in any channel.

#### group × SANS

| Ch | β (interaction) | t | p raw | p FDR-BH | sig |
|---|---|---|---|---|---|
| Left dlPFC | 1.214e-06 | 0.76 | 0.4499 | 0.8999 | --- |
| Left mPFC | -1.102e-07 | -0.08 | 0.9387 | 0.9892 | --- |
| Right mPFC | -1.526e-08 | -0.01 | 0.9892 | 0.9892 | --- |
| Right dlPFC | -1.770e-06 | -1.57 | 0.1195 | 0.4779 | --- |

#### group × SAPS

| Ch | β (interaction) | t | p raw | p FDR-BH | sig |
|---|---|---|---|---|---|
| Left dlPFC | 4.975e-07 | 0.20 | 0.8403 | 0.8403 | --- |
| Left mPFC | 1.347e-06 | 0.61 | 0.5413 | 0.8403 | --- |
| Right mPFC | 4.538e-07 | 0.26 | 0.7940 | 0.8403 | --- |
| Right dlPFC | 6.796e-07 | 0.39 | 0.6963 | 0.8403 | --- |

#### group × WAIS_MATRIX

| Ch | β (interaction) | t | p raw | p FDR-BH | sig |
|---|---|---|---|---|---|
| Left dlPFC | 7.267e-07 | 0.27 | 0.7858 | 0.9229 | --- |
| Left mPFC | -2.310e-07 | -0.10 | 0.9229 | 0.9229 | --- |
| Right mPFC | 4.301e-07 | 0.23 | 0.8193 | 0.9229 | --- |
| Right dlPFC | -2.202e-06 | -1.17 | 0.2440 | 0.9229 | --- |

#### group × SUD

| Ch | β (interaction) | t | p raw | p FDR-BH | sig |
|---|---|---|---|---|---|
| Left dlPFC | 2.557e-05 | 0.93 | 0.3528 | 0.7663 | --- |
| Left mPFC | 1.450e-05 | 0.59 | 0.5547 | 0.7663 | --- |
| Right mPFC | -6.253e-06 | -0.32 | 0.7466 | 0.7663 | --- |
| Right dlPFC | 5.760e-06 | 0.30 | 0.7663 | 0.7663 | --- |


---

## 2. PSD+AGG vs HC — v7 (LRV corsi_GLT)

> Comparison of PSD+AGG (n=77) vs HC (n=85) on the v7 (LRV corsi_GLT) contrast, with SANS, SAPS, WAIS-matrix, SUD as covariates. No channel cluster survived FWER correction. The strongest single-channel effect was Left dlPFC (t = -2.08, p_raw = 0.0391, d = -0.21), uncorrected.

### ANCOVA Region Results

| Region | β (ΔHbO) | t | d | p raw | FDR-BH | Hochberg | Hommel | Perm-max |
|---|---|---|---|---|---|---|---|---|
| Left dlPFC | -2.19e-05 | -2.08 | -0.21 | 0.0391 | 0.1563 | 0.1563 | 0.1563 | 0.1180 |
| Left mPFC | -1.43e-05 | -1.49 | -0.25 | 0.1373 | 0.1831 | 0.2746 | 0.2746 | 0.3696 |
| Right mPFC | -7.25e-06 | -0.89 | -0.21 | 0.3763 | 0.3763 | 0.3763 | 0.3763 | 0.7864 |
| Right dlPFC | -1.38e-05 | -1.72 | -0.03 | 0.0869 | 0.1738 | 0.2607 | 0.2059 | 0.2470 |

### Mann-Whitney (non-parametric, unadjusted for covariates)

| Region | p raw | FDR-BH | Hommel | rb *r* |
|---|---|---|---|---|
| Left dlPFC | 0.3383 | 0.4510 | 0.6766 | 0.087 |
| Left mPFC | 0.0629 | 0.2517 | 0.2517 | 0.170 |
| Right mPFC | 0.2534 | 0.4510 | 0.5074 | 0.104 |
| Right dlPFC | 0.9639 | 0.9639 | 0.9639 | -0.004 |

### Secondary Analysis: Group × Covariate Interaction Effects

> Secondary interaction analysis revealed the following significant group × covariate effects (FDR-BH corrected, α=0.05): group × SANS: Left mPFC (t=2.40, FDR p=0.039); Right mPFC (t=2.36, FDR p=0.039). This indicates that the modulating effect of these clinical variables on cortical haemodynamic response differs between groups.

#### group × SANS

| Ch | β (interaction) | t | p raw | p FDR-BH | sig |
|---|---|---|---|---|---|
| Left dlPFC | 1.274e-06 | 1.24 | 0.2179 | 0.2905 | --- |
| Left mPFC | 2.218e-06 | 2.40 | 0.0174 | **0.0388\*** | YES ✓ |
| Right mPFC | 1.864e-06 | 2.36 | 0.0194 | **0.0388\*** | YES ✓ |
| Right dlPFC | 2.190e-07 | 0.28 | 0.7800 | 0.7800 | --- |

#### group × SAPS

| Ch | β (interaction) | t | p raw | p FDR-BH | sig |
|---|---|---|---|---|---|
| Left dlPFC | 4.914e-07 | 0.31 | 0.7567 | 0.7567 | --- |
| Left mPFC | 1.824e-06 | 1.29 | 0.2005 | 0.2673 | --- |
| Right mPFC | 1.627e-06 | 1.34 | 0.1820 | 0.2673 | --- |
| Right dlPFC | 1.963e-06 | 1.63 | 0.1049 | 0.2673 | --- |

#### group × WAIS_MATRIX

| Ch | β (interaction) | t | p raw | p FDR-BH | sig |
|---|---|---|---|---|---|
| Left dlPFC | -1.075e-06 | -0.63 | 0.5315 | 0.9649 | --- |
| Left mPFC | 4.600e-07 | 0.30 | 0.7651 | 0.9649 | --- |
| Right mPFC | -5.799e-08 | -0.04 | 0.9649 | 0.9649 | --- |
| Right dlPFC | -3.163e-07 | -0.24 | 0.8086 | 0.9649 | --- |

#### group × SUD

| Ch | β (interaction) | t | p raw | p FDR-BH | sig |
|---|---|---|---|---|---|
| Left dlPFC | 1.680e-05 | 0.95 | 0.3421 | 0.7700 | --- |
| Left mPFC | 8.397e-06 | 0.53 | 0.5957 | 0.7943 | --- |
| Right mPFC | -1.837e-06 | -0.14 | 0.8920 | 0.8920 | --- |
| Right dlPFC | 1.167e-05 | 0.87 | 0.3850 | 0.7700 | --- |


---

## 3. PSD+AGG vs PSD-AGG — v5 (LRV hard-easy)

> Comparison of PSD+AGG (n=77) vs PSD-AGG (n=66) on the v5 (LRV hard-easy) contrast, with SANS, SAPS, WAIS-matrix, SUD, DOSE as covariates. No channel cluster survived FWER correction. The strongest single-channel effect was Left dlPFC (t = -1.24, p_raw = 0.2177, d = -0.12), uncorrected.

### ANCOVA Region Results

| Region | β (ΔHbO) | t | d | p raw | FDR-BH | Hochberg | Hommel | Perm-max |
|---|---|---|---|---|---|---|---|---|
| Left dlPFC | -1.36e-05 | -1.24 | -0.12 | 0.2177 | 0.8710 | 0.8710 | 0.8710 | 0.5638 |
| Left mPFC | -2.19e-06 | -0.23 | -0.07 | 0.8148 | 0.9678 | 0.9678 | 0.9678 | 0.9984 |
| Right mPFC | 2.17e-06 | 0.26 | 0.06 | 0.7991 | 0.9678 | 0.9678 | 0.9678 | 0.9978 |
| Right dlPFC | 3.57e-07 | 0.04 | 0.04 | 0.9678 | 0.9678 | 0.9678 | 0.9678 | 1.0000 |

### Mann-Whitney (non-parametric, unadjusted for covariates)

| Region | p raw | FDR-BH | Hommel | rb *r* |
|---|---|---|---|---|
| Left dlPFC | 0.7444 | 0.8300 | 0.8300 | 0.032 |
| Left mPFC | 0.2921 | 0.8300 | 0.8300 | 0.103 |
| Right mPFC | 0.5915 | 0.8300 | 0.8300 | 0.052 |
| Right dlPFC | 0.8300 | 0.8300 | 0.8300 | 0.021 |

### Secondary Analysis: Group × Covariate Interaction Effects

> No significant group × covariate interaction effects were found after FDR-BH correction (α=0.05), suggesting that the relationships between SANS, SAPS, WAIS_MATRIX, SUD, DOSE and ΔHbO do not differ significantly between PSD+AGG and PSD-AGG in any channel.

#### group × SANS

| Ch | β (interaction) | t | p raw | p FDR-BH | sig |
|---|---|---|---|---|---|
| Left dlPFC | 2.321e-08 | 0.03 | 0.9769 | 0.9769 | --- |
| Left mPFC | -1.998e-07 | -0.30 | 0.7679 | 0.9769 | --- |
| Right mPFC | -3.987e-07 | -0.64 | 0.5221 | 0.9769 | --- |
| Right dlPFC | -1.517e-06 | -2.42 | 0.0167 | 0.0668 | --- |

#### group × SAPS

| Ch | β (interaction) | t | p raw | p FDR-BH | sig |
|---|---|---|---|---|---|
| Left dlPFC | 7.669e-07 | 0.86 | 0.3910 | 0.7764 | --- |
| Left mPFC | 1.172e-06 | 1.56 | 0.1216 | 0.4863 | --- |
| Right mPFC | -1.968e-07 | -0.28 | 0.7764 | 0.7764 | --- |
| Right dlPFC | 3.371e-07 | 0.48 | 0.6292 | 0.7764 | --- |

#### group × WAIS_MATRIX

| Ch | β (interaction) | t | p raw | p FDR-BH | sig |
|---|---|---|---|---|---|
| Left dlPFC | 2.274e-06 | 0.82 | 0.4125 | 0.9292 | --- |
| Left mPFC | 2.078e-07 | 0.09 | 0.9292 | 0.9292 | --- |
| Right mPFC | -7.453e-07 | -0.35 | 0.7289 | 0.9292 | --- |
| Right dlPFC | -9.845e-07 | -0.46 | 0.6497 | 0.9292 | --- |

#### group × SUD

| Ch | β (interaction) | t | p raw | p FDR-BH | sig |
|---|---|---|---|---|---|
| Left dlPFC | 2.213e-05 | 1.01 | 0.3127 | 0.8283 | --- |
| Left mPFC | -1.801e-06 | -0.10 | 0.9223 | 0.9223 | --- |
| Right mPFC | -1.332e-05 | -0.79 | 0.4333 | 0.8283 | --- |
| Right dlPFC | -8.454e-06 | -0.50 | 0.6213 | 0.8283 | --- |

#### group × DOSE

| Ch | β (interaction) | t | p raw | p FDR-BH | sig |
|---|---|---|---|---|---|
| Left dlPFC | -5.531e-07 | -0.49 | 0.6284 | 0.6284 | --- |
| Left mPFC | 5.495e-07 | 0.57 | 0.5688 | 0.6284 | --- |
| Right mPFC | 6.172e-07 | 0.70 | 0.4866 | 0.6284 | --- |
| Right dlPFC | 1.879e-06 | 2.11 | 0.0370 | 0.1478 | --- |


---

## 4. PSD+AGG vs PSD-AGG — v7 (LRV corsi_GLT)

> Comparison of PSD+AGG (n=77) vs PSD-AGG (n=66) on the v7 (LRV corsi_GLT) contrast, with SANS, SAPS, WAIS-matrix, SUD, DOSE as covariates. No channel cluster survived FWER correction. The strongest single-channel effect was Left dlPFC (t = -1.15, p_raw = 0.2531, d = -0.15), uncorrected.

### ANCOVA Region Results

| Region | β (ΔHbO) | t | d | p raw | FDR-BH | Hochberg | Hommel | Perm-max |
|---|---|---|---|---|---|---|---|---|
| Left dlPFC | -8.15e-06 | -1.15 | -0.15 | 0.2531 | 0.6173 | 0.7117 | 0.7036 | 0.6364 |
| Left mPFC | -4.30e-06 | -0.73 | -0.11 | 0.4691 | 0.6254 | 0.7117 | 0.7117 | 0.8948 |
| Right mPFC | 1.83e-06 | 0.37 | 0.06 | 0.7117 | 0.7117 | 0.7117 | 0.7117 | 0.9924 |
| Right dlPFC | -5.30e-06 | -1.02 | -0.13 | 0.3087 | 0.6173 | 0.7117 | 0.7036 | 0.7294 |

### Mann-Whitney (non-parametric, unadjusted for covariates)

| Region | p raw | FDR-BH | Hommel | rb *r* |
|---|---|---|---|---|
| Left dlPFC | 0.5735 | 0.9060 | 0.9371 | 0.055 |
| Left mPFC | 0.2352 | 0.9060 | 0.9060 | 0.116 |
| Right mPFC | 0.9371 | 0.9371 | 0.9371 | -0.008 |
| Right dlPFC | 0.6795 | 0.9060 | 0.9371 | 0.040 |

### Secondary Analysis: Group × Covariate Interaction Effects

> No significant group × covariate interaction effects were found after FDR-BH correction (α=0.05), suggesting that the relationships between SANS, SAPS, WAIS_MATRIX, SUD, DOSE and ΔHbO do not differ significantly between PSD+AGG and PSD-AGG in any channel.

#### group × SANS

| Ch | β (interaction) | t | p raw | p FDR-BH | sig |
|---|---|---|---|---|---|
| Left dlPFC | 1.031e-06 | 2.04 | 0.0434 | 0.1735 | --- |
| Left mPFC | 1.398e-07 | 0.33 | 0.7443 | 0.7443 | --- |
| Right mPFC | 3.089e-07 | 0.86 | 0.3925 | 0.5233 | --- |
| Right dlPFC | 4.319e-07 | 1.16 | 0.2501 | 0.5002 | --- |

#### group × SAPS

| Ch | β (interaction) | t | p raw | p FDR-BH | sig |
|---|---|---|---|---|---|
| Left dlPFC | 1.877e-07 | 0.33 | 0.7392 | 0.7392 | --- |
| Left mPFC | 7.353e-07 | 1.54 | 0.1250 | 0.5000 | --- |
| Right mPFC | 1.694e-07 | 0.42 | 0.6733 | 0.7392 | --- |
| Right dlPFC | 4.376e-07 | 1.05 | 0.2950 | 0.5901 | --- |

#### group × WAIS_MATRIX

| Ch | β (interaction) | t | p raw | p FDR-BH | sig |
|---|---|---|---|---|---|
| Left dlPFC | -9.490e-07 | -0.54 | 0.5877 | 0.6111 | --- |
| Left mPFC | -1.327e-06 | -0.90 | 0.3711 | 0.6111 | --- |
| Right mPFC | -6.340e-07 | -0.51 | 0.6111 | 0.6111 | --- |
| Right dlPFC | 1.051e-06 | 0.81 | 0.4174 | 0.6111 | --- |

#### group × SUD

| Ch | β (interaction) | t | p raw | p FDR-BH | sig |
|---|---|---|---|---|---|
| Left dlPFC | 1.704e-05 | 1.24 | 0.2187 | 0.6124 | --- |
| Left mPFC | 5.112e-07 | 0.04 | 0.9651 | 0.9651 | --- |
| Right mPFC | -6.650e-06 | -0.68 | 0.4994 | 0.6659 | --- |
| Right dlPFC | 1.048e-05 | 1.03 | 0.3062 | 0.6124 | --- |

#### group × DOSE

| Ch | β (interaction) | t | p raw | p FDR-BH | sig |
|---|---|---|---|---|---|
| Left dlPFC | -1.503e-06 | -2.09 | 0.0387 | 0.1546 | --- |
| Left mPFC | -5.002e-07 | -0.82 | 0.4131 | 0.7333 | --- |
| Right mPFC | 1.751e-07 | 0.34 | 0.7333 | 0.7333 | --- |
| Right dlPFC | 2.474e-07 | 0.46 | 0.6430 | 0.7333 | --- |


---

## 5. PSD-AGG vs HC — v5 (pcon hard-easy)

> Comparison of PSD-AGG (n=66) vs HC (n=85) on the v5 (pcon hard-easy) contrast, with SANS, SAPS, WAIS-matrix, SUD as covariates. No channel cluster survived FWER correction. The strongest single-channel effect was Right dlPFC (t = -2.36, p_raw = 0.0197, d = -0.14), uncorrected.

### ANCOVA Region Results

| Region | β (ΔHbO) | t | d | p raw | FDR-BH | Hochberg | Hommel | Perm-max |
|---|---|---|---|---|---|---|---|---|
| Left dlPFC | -2.31e-05 | -2.19 | -0.18 | 0.0298 | 0.0596 | 0.0894 | 0.0894 | 0.0872 |
| Left mPFC | -9.35e-06 | -0.99 | -0.12 | 0.3214 | 0.3214 | 0.3214 | 0.3214 | 0.6860 |
| Right mPFC | -1.21e-05 | -1.16 | -0.14 | 0.2476 | 0.3214 | 0.3214 | 0.3214 | 0.5640 |
| Right dlPFC | -2.55e-05 | -2.36 | -0.14 | 0.0197 | 0.0596 | 0.0788 | 0.0596 | 0.0594 |

### Mann-Whitney (non-parametric, unadjusted for covariates)

| Region | p raw | FDR-BH | Hommel | rb *r* |
|---|---|---|---|---|
| Left dlPFC | 0.4841 | 0.9357 | 0.9357 | 0.067 |
| Left mPFC | 0.9357 | 0.9357 | 0.9357 | 0.008 |
| Right mPFC | 0.7684 | 0.9357 | 0.9357 | 0.028 |
| Right dlPFC | 0.6364 | 0.9357 | 0.9357 | 0.045 |

### Secondary Analysis: Group × Covariate Interaction Effects

> No significant group × covariate interaction effects were found after FDR-BH correction (α=0.05), suggesting that the relationships between SANS, SAPS, WAIS_MATRIX, SUD and ΔHbO do not differ significantly between PSD-AGG and HC in any channel.

#### group × SANS

| Ch | β (interaction) | t | p raw | p FDR-BH | sig |
|---|---|---|---|---|---|
| Left dlPFC | 1.284e-06 | 0.93 | 0.3561 | 0.9892 | --- |
| Left mPFC | 1.684e-08 | 0.01 | 0.9892 | 0.9892 | --- |
| Right mPFC | 3.403e-07 | 0.25 | 0.8057 | 0.9892 | --- |
| Right dlPFC | -4.994e-07 | -0.35 | 0.7276 | 0.9892 | --- |

#### group × SAPS

| Ch | β (interaction) | t | p raw | p FDR-BH | sig |
|---|---|---|---|---|---|
| Left dlPFC | -2.960e-07 | -0.14 | 0.8884 | 0.9453 | --- |
| Left mPFC | 1.296e-07 | 0.07 | 0.9453 | 0.9453 | --- |
| Right mPFC | 4.585e-07 | 0.22 | 0.8272 | 0.9453 | --- |
| Right dlPFC | 1.818e-07 | 0.08 | 0.9335 | 0.9453 | --- |

#### group × WAIS_MATRIX

| Ch | β (interaction) | t | p raw | p FDR-BH | sig |
|---|---|---|---|---|---|
| Left dlPFC | -1.798e-06 | -0.80 | 0.4231 | 0.9031 | --- |
| Left mPFC | -2.444e-07 | -0.12 | 0.9031 | 0.9031 | --- |
| Right mPFC | 1.286e-06 | 0.58 | 0.5649 | 0.9031 | --- |
| Right dlPFC | -5.564e-07 | -0.24 | 0.8100 | 0.9031 | --- |

#### group × SUD

| Ch | β (interaction) | t | p raw | p FDR-BH | sig |
|---|---|---|---|---|---|
| Left dlPFC | 5.951e-06 | 0.25 | 0.8034 | 0.8530 | --- |
| Left mPFC | 1.385e-05 | 0.65 | 0.5178 | 0.8530 | --- |
| Right mPFC | 4.410e-06 | 0.19 | 0.8530 | 0.8530 | --- |
| Right dlPFC | 5.845e-06 | 0.24 | 0.8127 | 0.8530 | --- |


---

## 6. PSD-AGG vs HC — v7 (pcon corsi_GLT)

> Comparison of PSD-AGG (n=66) vs HC (n=85) on the v7 (pcon corsi_GLT) contrast, with SANS, SAPS, WAIS-matrix, SUD as covariates. No channel cluster survived FWER correction. The strongest single-channel effect was Right dlPFC (t = 0.46, p_raw = 0.6434, d = 0.10), uncorrected.

### ANCOVA Region Results

| Region | β (ΔHbO) | t | d | p raw | FDR-BH | Hochberg | Hommel | Perm-max |
|---|---|---|---|---|---|---|---|---|
| Left dlPFC | 2.05e-06 | 0.28 | -0.04 | 0.7804 | 0.7804 | 0.7804 | 0.7804 | 0.9928 |
| Left mPFC | -2.33e-06 | -0.37 | -0.16 | 0.7122 | 0.7804 | 0.7804 | 0.7804 | 0.9822 |
| Right mPFC | -2.79e-06 | -0.41 | -0.24 | 0.6794 | 0.7804 | 0.7804 | 0.7804 | 0.9770 |
| Right dlPFC | 2.99e-06 | 0.46 | 0.10 | 0.6434 | 0.7804 | 0.7804 | 0.7804 | 0.9668 |

### Mann-Whitney (non-parametric, unadjusted for covariates)

| Region | p raw | FDR-BH | Hommel | rb *r* |
|---|---|---|---|---|
| Left dlPFC | 0.7684 | 0.7684 | 0.7684 | 0.028 |
| Left mPFC | 0.5150 | 0.7682 | 0.7684 | 0.062 |
| Right mPFC | 0.2456 | 0.7682 | 0.7682 | 0.111 |
| Right dlPFC | 0.5762 | 0.7682 | 0.7684 | -0.053 |

### Secondary Analysis: Group × Covariate Interaction Effects

> Secondary interaction analysis revealed the following significant group × covariate effects (FDR-BH corrected, α=0.05): group × SANS: Left mPFC (t=2.64, FDR p=0.037). This indicates that the modulating effect of these clinical variables on cortical haemodynamic response differs between groups.

#### group × SANS

| Ch | β (interaction) | t | p raw | p FDR-BH | sig |
|---|---|---|---|---|---|
| Left dlPFC | 4.169e-07 | 0.43 | 0.6679 | 0.7724 | --- |
| Left mPFC | 2.143e-06 | 2.64 | 0.0091 | **0.0365\*** | YES ✓ |
| Right mPFC | 1.530e-06 | 1.74 | 0.0834 | 0.1669 | --- |
| Right dlPFC | -2.450e-07 | -0.29 | 0.7724 | 0.7724 | --- |

#### group × SAPS

| Ch | β (interaction) | t | p raw | p FDR-BH | sig |
|---|---|---|---|---|---|
| Left dlPFC | 5.199e-07 | 0.35 | 0.7246 | 0.7246 | --- |
| Left mPFC | 1.137e-06 | 0.92 | 0.3572 | 0.4763 | --- |
| Right mPFC | 1.451e-06 | 1.09 | 0.2780 | 0.4763 | --- |
| Right dlPFC | 1.502e-06 | 1.17 | 0.2439 | 0.4763 | --- |

#### group × WAIS_MATRIX

| Ch | β (interaction) | t | p raw | p FDR-BH | sig |
|---|---|---|---|---|---|
| Left dlPFC | -5.881e-07 | -0.38 | 0.7078 | 0.7078 | --- |
| Left mPFC | 1.615e-06 | 1.23 | 0.2192 | 0.6984 | --- |
| Right mPFC | 6.444e-07 | 0.45 | 0.6500 | 0.7078 | --- |
| Right dlPFC | -1.282e-06 | -0.94 | 0.3492 | 0.6984 | --- |

#### group × SUD

| Ch | β (interaction) | t | p raw | p FDR-BH | sig |
|---|---|---|---|---|---|
| Left dlPFC | 6.395e-06 | 0.38 | 0.7022 | 0.9948 | --- |
| Left mPFC | 1.011e-05 | 0.72 | 0.4698 | 0.9948 | --- |
| Right mPFC | 4.027e-06 | 0.27 | 0.7901 | 0.9948 | --- |
| Right dlPFC | 9.534e-08 | 0.01 | 0.9948 | 0.9948 | --- |


---
