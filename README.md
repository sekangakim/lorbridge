# lorbridge

<!-- badges: start -->
<!-- [![CRAN status](https://www.r-pkg.org/badges/version/lorbridge)](https://CRAN.R-project.org/package=lorbridge) -->
[![R-CMD-check](https://github.com/sekangakim/lorbridge/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/sekangakim/lorbridge/actions/workflows/R-CMD-check.yaml)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
<!-- badges: end -->

## Overview

**lorbridge** provides a unified analytical framework that connects two 
worlds that clinical and medical researchers often treat separately: 
conventional logistic regression and geometric correspondence analysis.

### The Problem

When a logistic regression model produces a log-odds ratio (LOR) of 1.47 
for a predictor, most clinical faculty find it difficult to judge whether 
this represents a small, moderate, or large association — because the 
log-odds ratio scale is unbounded and asymmetric. The exponentiated form 
(odds ratio = 4.35) is easier to communicate but still not directly 
comparable across predictors or across studies.

### The Solution

Building on Kim and Grochowalski (2019), **lorbridge** re-expresses every 
LOR from binary and multinomial logistic regression as a **cosine theta** 
estimate — a measure bounded between −1 and +1, interpretable exactly like 
a Pearson correlation. A cosine theta of +0.68 is immediately recognisable 
as a strong positive association to any researcher familiar with 
correlations, without requiring statistical training in odds ratios or 
log-transformations.

In addition to cosine theta, **lorbridge** computes three further 
**closeness-of-concordance measures (CCMs)** — Yule's Q, Yule's Y, and 
r_meta (a meta-analytic correlation) — all on the same −1 to +1 scale, 
with bootstrap confidence intervals throughout.

### The Bridge to Correspondence Analysis

The same cosine theta metric connects logistic regression output to the 
geometry of **Singly-Ordered Nonsymmetric Correspondence Analysis 
(SONSCA)** and **Doubly-Ordered Nonsymmetric Correspondence Analysis 
(DONSCA)**. These methods extend conventional correspondence analysis to 
tables where one or both variable orderings carry substantive meaning — 
for example, a table of racial groups (nominal rows) by discretised IQ 
or Vocabulary Meaning score bins (ordered columns).

SONSCA and DONSCA produce geometric biplots that make complex categorical 
relationships immediately visible, while the cosine theta values anchored 
at user-specified reference categories provide a single, interpretable 
effect-size summary for each pairwise contrast — directly comparable to 
the LOR-derived CCMs from the logistic regression in the same analysis.

### Who Is This For?

**lorbridge** is designed for:

- Clinical and medical researchers who routinely report odds ratios and 
  want an intuitive, correlation-like effect-size companion measure
- Methodologists working with ordered categorical data who want to connect 
  logistic regression and correspondence analysis in a single coherent 
  workflow
- Graduate students and faculty learning to interpret association strength 
  beyond p-values and raw odds ratios

### Methodological Lineage

This package implements and extends methods from:

- Kim, S.-K., & Grochowalski, J. H. (2019) — the original cosine theta 
  bridge between LORs and CA geometry
- Kim, S.-K. (2020, 2021, 2022) — clinical applications of matched and 
  nonsymmetric CA in psychiatric and pediatric research
- Kim, S.-K. (2024) — factor analytic profile analysis and interpretable 
  machine learning for biomedical data

---

## Installation

```r
# From CRAN (once published):
install.packages("lorbridge")

# Development version from GitHub:
# install.packages("remotes")
remotes::install_github("sekangakim/lorbridge")
```

---

## Quick Start

```r
library(lorbridge)
data(lorbridge_data)

## --- Binary logistic regression: continuous predictor ---
res_cont <- blr_continuous(
  outcome   = lorbridge_data$minority,
  predictor = lorbridge_data$VM
)
print(res_cont$summary_table[, c("LOR","OR","OR_lo","OR_hi","p","r_meta")])

## --- Binary logistic regression: categorical predictor ---
res_cat <- blr_categorical(
  outcome   = lorbridge_data$minority,
  predictor = lorbridge_data$VMbin,
  ref_level = "VM4"
)
print(res_cat$results[, c("Category","LOR","OR","p","YuleQ","cos_theta")])

## --- SONSCA pairwise CCMs ---
data(tab_IQ)
sonsca_ccm(tab_IQ,
           row_k = "Race1", bin_j = "IQ1",
           row_anchor = "Race2", col_anchor = "IQ4")

## --- DONSCA cosine theta ---
data(tab_IQ_VM)
fit <- donsca_fit(tab_IQ_VM)
donsca_cosines(fit, col_anchor_idx = 4, row_anchor_idx = 4)
```

---

## Package Structure

| Function | Description |
|---|---|
| `lor_ci_2x2()` | Haldane–Anscombe-corrected 2×2 LOR with Wald CI |
| `ccm_row()` | Full CCM row: OR, LOR, YuleQ, YuleY, r_meta with CIs |
| `blr_continuous()` | Binary logistic regression, continuous predictor (per 1 SD) |
| `blr_categorical()` | Binary logistic regression, categorical predictor |
| `cosine_theta_2row()` | Cosine theta from 2-row CA via SVD (LOR bridge) |
| `sonsca_coords()` | SONSCA column-isometric coordinates |
| `sonsca_cosines()` | Doubly-anchored cosine theta matrix (SONSCA) |
| `sonsca_bootstrap()` | Bootstrap CIs for SONSCA cosine theta |
| `sonsca_ccm()` | Pairwise CCMs for a single SONSCA contrast |
| `inertia_pct()` | Percent inertia per CA dimension |
| `donsca_fit()` | Fit a DONSCA model |
| `donsca_cosines()` | Doubly-anchored cosine theta (DONSCA) |
| `mlr_ccm()` | Multinomial logistic regression with CCMs |

---

## Datasets

| Dataset | Description |
|---|---|
| `lorbridge_data` | Individual-level data (N = 900): VM scores, bins, minority/majority |
| `tab_IQ` | 4 races × 6 IQ bins (SONSCA) |
| `tab_VM` | 4 races × 6 VM bins (SONSCA) |
| `tab_IQ_VM` | 6 IQ bins × 6 VM bins (DONSCA) |

---

## References

Kim, S.-K., & Grochowalski, J. H. (2019). Gaining from discretization of
continuous data: The correspondence analysis biplot approach.
*Behavior Research Methods*, 51(2), 589–601.
https://doi.org/10.3758/s13428-018-1161-1

Kim, S.-K. (2020). Test treatment effect differences in repeatedly measured
symptoms with binary values: The matched correspondence analysis approach.
*Behavior Research Methods*, 52, 1480–1490.
https://doi.org/10.3758/s13428-019-01328-9

Kim, S.-K., McKay, D., Murphy, T. K., Bussing, R., McNamara, J. P.,
Goodman, W. K., & Storch, E. C. (2021). Age moderated–anxiety mediation
for multimodal treatment outcome among children with obsessive-compulsive
disorder: An evaluation with correspondence analysis.
*Journal of Affective Disorders*, 282, 766–775.
https://doi.org/10.1016/j.jad.2020.12.198

Kim, S.-K. (2022). Assessment of improvement in anxiety severity for children
with autism spectrum disorders: The matched correspondence analysis approach.
*Journal of Psychiatric Research*, 145, 175–181.
https://doi.org/10.1016/j.jpsychires.2021.12.004

Kim, S.-K. (2024). Factorization of person response profiles to identify
summative profiles carrying central response patterns.
*Psychological Methods*, 29(4), 723–730.
https://doi.org/10.1037/met0000568

---

## Author

**Se-Kang Kim, Ph.D.**
Professor (Tenured), Division of Pediatric Psychology
Baylor College of Medicine, Houston, TX
se-kang.kim@bcm.edu
ORCID: [0000-0003-0928-3396](https://orcid.org/0000-0003-0928-3396)

---

## License

GPL-3 © Se-Kang Kim
