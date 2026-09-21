# =============================================================================
# 06_evidence_synthesis.R  —  Meta-analysis, meta-regression, NMA
# =============================================================================
# R's metafor/meta ecosystem is more complete than anything in Python, which
# is why this gets its own file rather than sitting under "advanced methods".
# Pairs with PRISMA 2020 for the review itself.
# =============================================================================

library(tidyverse)
library(meta); library(metafor); library(dmetar)


# --- 1. Input data -----------------------------------------------------------
# One row per study. Columns depend on the effect measure:
#   binary   : event_t, n_t, event_c, n_c
#   continuous: mean_t, sd_t, n_t, mean_c, sd_c, n_c
#   generic  : effect estimate (on the analysis scale) + its standard error

dat <- readr::read_csv(here::here("data", "raw", "extraction.csv"))


# --- 2. Pairwise meta-analysis -----------------------------------------------

## Binary outcomes
m_bin <- metabin(
  event.e = event_t, n.e = n_t,
  event.c = event_c, n.c = n_c,
  studlab = study, data = dat,
  sm      = "RR",            # "RR", "OR", "RD"
  method  = "MH",            # Mantel-Haenszel; "Inverse" for generic
  random  = TRUE, common = FALSE,
  method.tau = "REML",       # between-study variance estimator
  method.random.ci = "HK"    # Hartung-Knapp: wider, better coverage with
)                            # few studies. Use it by default.
summary(m_bin)

## Continuous outcomes
m_cont <- metacont(
  n.e = n_t, mean.e = mean_t, sd.e = sd_t,
  n.c = n_c, mean.c = mean_c, sd.c = sd_c,
  studlab = study, data = dat,
  sm = "SMD",                # "MD" if all studies share a scale
  method.smd = "Hedges",     # small-sample corrected
  random = TRUE, common = FALSE, method.random.ci = "HK"
)

## Generic inverse variance — for adjusted estimates, HRs, anything pre-computed
## NOTE: ratio measures must be entered on the LOG scale.
m_gen <- metagen(
  TE = log(hr), seTE = se_log_hr,
  studlab = study, data = dat,
  sm = "HR", random = TRUE, common = FALSE, method.random.ci = "HK"
)

## Single-arm proportions (prevalence, utility values, complication rates)
m_prop <- metaprop(
  event = n_event, n = n_total, studlab = study, data = dat,
  sm = "PLOGIT",             # logit transform keeps CIs inside [0,1]
  random = TRUE, method.random.ci = "HK"
)

## Means (e.g. pooling OHIP-14 or EQ-5D scores across studies)
m_mean <- metamean(n = n, mean = mean_score, sd = sd_score,
                   studlab = study, data = dat, random = TRUE)


# --- 3. Heterogeneity --------------------------------------------------------

# I-squared  % of variability due to heterogeneity rather than chance
# tau-squared between-study variance, on the analysis scale
# Q          chi-squared test; underpowered with few studies
#
# The PREDICTION INTERVAL is the most useful single statistic — it gives the
# range a new study's true effect would fall in, and is what tells you whether
# the pooled estimate means anything.

m_bin$I2; m_bin$tau2; m_bin$lower.predict; m_bin$upper.predict

# Do not treat I2 thresholds mechanically. High I2 with a narrow prediction
# interval is unimportant; low I2 across three studies is uninformative.


# --- 4. Forest and funnel plots ----------------------------------------------

forest(
  m_bin,
  sortvar    = TE,
  prediction = TRUE,
  print.tau2 = TRUE,
  leftcols   = c("studlab", "event.e", "n.e", "event.c", "n.c"),
  leftlabs   = c("Study", "Events", "Total", "Events", "Total"),
  label.e    = "Intervention", label.c = "Control",
  xlab       = "Risk ratio"
)

# Publication bias. Needs >= 10 studies; below that the tests have no power.
funnel(m_bin, studlab = TRUE)
metabias(m_bin, method.bias = "Egger")       # continuous outcomes
metabias(m_bin, method.bias = "peters")      # binary outcomes; Egger is
                                             # miscalibrated for log RR/OR
trimfill(m_bin)                              # sensitivity, not a correction

# Funnel asymmetry has causes other than publication bias: small-study effects,
# heterogeneity, poor methodology in small trials. Do not report it as proof.


# --- 5. Subgroup analysis and meta-regression --------------------------------

update(m_bin, subgroup = risk_of_bias, tau.common = FALSE)

# Meta-regression needs roughly 10 studies per covariate. Aggregation bias
# means study-level covariates cannot be read as individual-level effects.
mr <- metareg(m_bin, ~ mean_age + follow_up_months)
summary(mr)
bubble(mr, studlab = TRUE)


# --- 6. metafor — when you need more control ---------------------------------

es <- escalc(measure = "RR",
             ai = event_t, n1i = n_t,
             ci = event_c, n2i = n_c, data = dat)

res <- rma(yi, vi, data = es, method = "REML", test = "knha")
summary(res)
predict(res, transf = exp, digits = 2)

## Influence and leave-one-out
inf <- influence(res); plot(inf)
leave1out(res, transf = exp)
baujat(res)                          # heterogeneity vs influence

## Multilevel / multiple outcomes per study — accounts for the dependence
## you create when one trial contributes several estimates
rma.mv(yi, vi, random = ~ 1 | study_id / effect_id, data = es)

## Robust variance estimation for dependent effect sizes
library(clubSandwich)
coef_test(res_mv, vcov = "CR2")


# --- 7. Network meta-analysis ------------------------------------------------
# For >2 comparators connected through a network of trials.

nma <- netmeta(
  TE = TE, seTE = seTE,
  treat1 = treat1, treat2 = treat2, studlab = study,
  data = dat_nma, sm = "RR",
  reference.group = "Usual care",
  common = FALSE, random = TRUE
)
summary(nma)

netgraph(nma, plastic = FALSE, thickness = "number.of.studies")
forest(nma, reference.group = "Usual care")

## Transitivity is the core assumption: studies must be similar enough that
## indirect comparison is valid. Check the distribution of effect modifiers
## across comparisons before fitting.

## Inconsistency — direct vs indirect evidence must agree
netsplit(nma)                        # node splitting
decomp.design(nma)

## Ranking. Report SUCRA with its uncertainty; a rank order alone is
## overconfident and routinely over-interpreted.
netrank(nma, small.values = "good")
rankogram(nma)


# --- 8. Diagnostic test accuracy ---------------------------------------------
# Bivariate model: sensitivity and specificity are negatively correlated
# across studies (threshold effect), so they must be pooled jointly.

library(mada)
madad(dat_dta)                                   # per-study estimates
fit_dta <- reitsma(dat_dta)                      # bivariate / HSROC
summary(fit_dta)
plot(fit_dta, sroclwd = 2)


# --- 9. Feeding results into an economic model -------------------------------
# A cost-effectiveness model needs a distribution, not a point estimate.
# Carry the full uncertainty through — and prefer the predictive distribution
# when the model represents a new setting rather than the average trial.

pooled_lor <- res$b[1]
pooled_se  <- res$se

set.seed(1234)
psa_or <- exp(rnorm(10000, pooled_lor, pooled_se))       # parameter uncertainty
quantile(psa_or, c(0.025, 0.5, 0.975))

# Predictive distribution — includes between-study heterogeneity
psa_or_pred <- exp(rnorm(10000, pooled_lor, sqrt(pooled_se^2 + res$tau2)))
quantile(psa_or_pred, c(0.025, 0.5, 0.975))

# Beta distribution for a pooled probability (utilities, event risks):
#   method of moments from pooled mean and variance
mu <- 0.75; v <- 0.01
alpha <- mu * (mu * (1 - mu) / v - 1)
beta  <- alpha * (1 - mu) / mu
rbeta(10000, alpha, beta)
