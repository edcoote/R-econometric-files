# Econometrics Toolkit — R

*Built for applied health economics research. Adapt by swapping in your own data frame and variable names.*

Companion to the Python toolkit. Restructured to follow how R's package ecosystem actually divides the work, rather than mirroring the Python file split — which means data preparation, evidence synthesis and economic evaluation get their own files, since R is where that work usually happens.

---

## Files

| File | Methods |
|---|---|
| `00_setup.R` | Packages, session hygiene, import conventions, data integrity checks, numeric-coercion diagnostics |
| `01_data_prep.R` | Type fixing, factors and reference levels, derived variables, plausibility ranges, missing data (MICE), Table 1 |
| `02_linear_regression.R` | OLS, functional form, interactions and marginal effects, full diagnostic suite, robust/clustered SEs (HC0–HC3, CR, HAC), WLS, survey weights, quantile regression |
| `03_glm_outcomes.R` | Logit/probit/cloglog, risk ratios and risk differences, ordered logit (+ Brant test), multinomial, Poisson, negative binomial, hurdle, zero-inflated, cost models (Gamma-log, two-part, modified Park test) |
| `04_survival.R` | Kaplan–Meier, log-rank, RMST, Cox PH (+ Schoenfeld, time-varying), parametric models and extrapolation (TSD 14), AFT, competing risks (CIF, Fine–Gray), transition probabilities for Markov models |
| `05_causal_inference.R` | PSM (MatchIt), IPTW and doubly robust (WeightIt), balance diagnostics (SMD, Love plots), IV/2SLS (+ weak-instrument tests), DiD and event studies (+ staggered adoption), RDD (+ McCrary density), panel FE/RE, E-values |
| `06_evidence_synthesis.R` | Pairwise meta-analysis (binary, continuous, generic IV, proportions), heterogeneity and prediction intervals, publication bias, meta-regression, multilevel/RVE, network meta-analysis, diagnostic accuracy, feeding pooled estimates into PSA |
| `07_economic_evaluation.R` | Markov cohort models, discounting and half-cycle correction, ICERs and NMB, PSA with distribution selection, CE plane, CEAC, tornado/OWSA, EVPI and EVPPI, budget impact, CHEERS checklist |
| `08_reporting.R` | Publication tables (gtsummary, modelsummary), export to Word/Excel, figure themes and journal sizing, inline reporting helpers, reproducibility |

---

## Placeholder convention

The scripts are templates, not runnable demos. Substitute your own names for:

| Placeholder | Meaning |
|---|---|
| `d` | data frame / tibble |
| `y` | continuous outcome |
| `ybin` | binary outcome (factor, `"No"` / `"Yes"`) |
| `yord` | ordered categorical outcome |
| `ycount` | count outcome |
| `time`, `status` | follow-up time; event indicator (1 = event, 0 = censored) |
| `treat` | exposure or treatment (factor, reference = control) |
| `x1`, `x2` | continuous covariates |
| `grp` | categorical covariate (factor) |
| `id`, `clust` | individual identifier; clustering variable |

---

## Suggested project structure

```
project/
  project.Rproj
  data/raw/          read-only, never written to by a script
  data/derived/      analysis-ready, produced by 01_data_prep.R
  R/                 these scripts
  output/tables/
  output/figures/
```

Use an RStudio Project and `here()` rather than `setwd()`. Paths stay relative and the work is portable.

---

## Conventions

**Tidyverse for wrangling, base and specialist packages for modelling.** `dplyr` and `ggplot2` throughout for data handling and plots; `lm`, `glm`, `survival`, `flexsurv`, `metafor` for the models themselves, since that's where the methods actually live.

**`MASS` masks `dplyr::select`.** Load `MASS` before `tidyverse`, or write `dplyr::select()` explicitly.

**Fix the analysis sample once.** `lm()` and `glm()` drop rows with `NA` on any model variable, so models with different covariate sets have different `n` — which invalidates `anova()`, `AIC()` and likelihood ratio comparisons between them. `01_data_prep.R` sets `d_analysis` for this reason.

**Data prep runs once, from a fresh import.** `factor()` applied to an existing factor maps level *indices*, not original codes. Re-running a prep script on already-converted objects destroys data silently.

**Restart R and run top to bottom before believing a result.** Anything that only works because of what happens to be loaded isn't a result yet.

---

## Not covered here

Deep learning, and most general-purpose machine learning. For ML in R use `tidymodels` or `mlr3`, which are ecosystems rather than snippets and don't compress usefully into a reference file.

---

## Standards referenced

STROBE · CONSORT · PRISMA 2020 · CHEERS 2022 · TRIPOD+AI · RECORD · NICE DSU TSD 14 (survival extrapolation) · NICE reference case
