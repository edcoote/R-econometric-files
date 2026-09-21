# =============================================================================
# 05_causal_inference.R  —  Matching, weighting, IV, DiD, RDD, panel
# =============================================================================
# Every method here rests on an identifying assumption that the data cannot
# verify. State the assumption before the estimate, and test what is testable.
# =============================================================================

library(tidyverse); library(broom)
library(MatchIt); library(WeightIt); library(cobalt)
library(ivreg); library(fixest); library(rdrobust); library(plm)
library(sandwich); library(lmtest); library(marginaleffects)


# =============================================================================
# 1. PROPENSITY SCORE MATCHING
# =============================================================================
# Assumption: no unmeasured confounding (conditional ignorability). Matching
# balances what you measured; it does nothing about what you did not.

## Estimate and match in one step
m.out <- matchit(
  treat ~ age + sex + bmi + n_comorbid + grp,
  data     = d,
  method   = "nearest",
  distance = "glm",          # logistic propensity score
  link     = "logit",
  caliper  = 0.2,            # in SD of the linear propensity score
  std.caliper = TRUE,
  ratio    = 1,
  replace  = FALSE
)

summary(m.out, un = TRUE)

## Balance assessment — the actual output of a matching exercise
love.plot(m.out, thresholds = c(m = 0.1), abs = TRUE,
          var.order = "unadjusted")
bal.tab(m.out, un = TRUE, m.threshold = 0.1, v.threshold = 2)
plot(m.out, type = "density", interactive = FALSE)

# Judge on standardised mean differences (<0.1) and variance ratios (0.5-2).
# NEVER use p-values for balance: they conflate imbalance with sample size,
# and matching reduces n, so balance appears to improve automatically.

## Overlap / common support
plot(m.out, type = "jitter", interactive = FALSE)
# If treated units have no comparable controls, you are extrapolating.
# Restrict to the region of common support and say so.

## Outcome model on the matched sample — cluster on subclass
d_m <- match.data(m.out)

fit_m <- lm(y ~ treat + age + sex + bmi, data = d_m, weights = weights)
coeftest(fit_m, vcov = vcovCL(fit_m, cluster = ~ subclass))
avg_comparisons(fit_m, variables = "treat",
                vcov = ~ subclass, wts = "weights")

## Other matching approaches
matchit(treat ~ age + bmi, data = d, method = "full")       # full matching
matchit(treat ~ age + bmi, data = d, method = "optimal")    # optimal pairs
matchit(treat ~ age + bmi, data = d, method = "cem")        # coarsened exact
matchit(treat ~ age + bmi, data = d, method = "nearest",
        distance = "mahalanobis")


# =============================================================================
# 2. INVERSE PROBABILITY WEIGHTING
# =============================================================================
# Keeps the whole sample rather than discarding unmatched units.

w.out <- weightit(
  treat ~ age + sex + bmi + n_comorbid,
  data      = d,
  method    = "glm",
  estimand  = "ATE"        # ATT if the question is about the treated only
)

summary(w.out)             # check the weight distribution
bal.tab(w.out, m.threshold = 0.1)
love.plot(w.out, thresholds = c(m = 0.1))

## Extreme weights destabilise the estimate. Trim or stabilise:
w.out_t <- trim(w.out, at = 0.99)
# Stabilised weights: numerator = marginal P(treat), denominator = PS.

d$ipw <- w.out$weights
fit_w <- glm(y ~ treat, data = d, weights = ipw)
coeftest(fit_w, vcov = vcovHC(fit_w, "HC3"))    # robust SEs REQUIRED

## Doubly robust — consistent if EITHER the PS or the outcome model is correct
fit_dr <- glm(y ~ treat + age + sex + bmi, data = d, weights = ipw)
avg_comparisons(fit_dr, variables = "treat", wts = "ipw",
                vcov = "HC3")

## Targeted maximum likelihood, for a machine-learning nuisance model
# library(tmle); tmle(Y = d$y, A = d$treat_num, W = W_matrix)


# =============================================================================
# 3. INSTRUMENTAL VARIABLES
# =============================================================================
# Assumptions: (1) relevance — Z predicts treatment; (2) exclusion — Z affects
# Y only through treatment; (3) independence; (4) monotonicity.
# (2) is untestable and is where IV analyses usually fail. Argue it explicitly.
# Estimates the LATE — the effect among compliers, not the population.

iv <- ivreg(y ~ treat + age + sex | z + age + sex, data = d)
summary(iv, diagnostics = TRUE)
coeftest(iv, vcov = vcovHC(iv, "HC3"))

## Diagnostics printed by the above:
##   Weak instruments — first-stage F. Rule of thumb F > 10; below that the
##                      2SLS estimate is biased towards OLS and CIs undercover.
##   Wu-Hausman       — is OLS actually inconsistent? If not, prefer OLS.
##   Sargan           — overidentifying restrictions; only with >1 instrument.

## Inspect the first stage directly
fs <- lm(treat ~ z + age + sex, data = d)
summary(fs)
linearHypothesis(fs, "z = 0")          # F statistic on the excluded instrument

## Weak-instrument-robust inference
# library(ivmodel); AR.test(...)        Anderson-Rubin confidence set

## Common health economics instruments: distance to specialist centre,
## regional prescribing variation, calendar time of a guideline change,
## physician preference. Each needs the exclusion restriction defended.


# =============================================================================
# 4. DIFFERENCE-IN-DIFFERENCES
# =============================================================================
# Assumption: parallel trends. Untestable for the post period; test it in the
# pre period and show the plot.

## Canonical 2x2
did <- lm(y ~ treated * post, data = d)
coeftest(did, vcov = vcovCL(did, cluster = ~ clust))
# The interaction coefficient is the DiD estimate.

## Two-way fixed effects
did_fe <- feols(y ~ treated_post | unit + period, data = d, cluster = ~ unit)
summary(did_fe)

## Event study — the standard way to show parallel trends
es <- feols(y ~ i(period_rel, treated, ref = -1) | unit + period,
            data = d, cluster = ~ unit)
iplot(es, main = "Event study", xlab = "Periods relative to treatment")
# Pre-period coefficients should be flat and near zero. Trends before
# treatment mean parallel trends is implausible.

## STAGGERED ADOPTION — two-way FE is BIASED when units are treated at
## different times and effects vary over time (already-treated units serve as
## controls). Use a heterogeneity-robust estimator instead.
library(did)
cs <- att_gt(yname = "y", tname = "period", idname = "unit",
             gname = "first_treat_period", data = d,
             control_group = "notyettreated")
aggte(cs, type = "dynamic")            # event-study aggregation
ggdid(aggte(cs, type = "dynamic"))

# Alternatives: fixest::sunab() (Sun & Abraham), didimputation, DIDmultiplegt

## Inference: cluster at the level of treatment assignment. With few clusters
## (<40) use wild cluster bootstrap.
fwildclusterboot::boottest(did, clustid = "clust", param = "treated:post",
                           B = 9999)


# =============================================================================
# 5. REGRESSION DISCONTINUITY
# =============================================================================
# Assumption: units cannot precisely manipulate the running variable at the
# cutoff. Estimates a LATE at the threshold only.
# Health economics uses: age-based eligibility, risk-score treatment thresholds,
# QOF payment thresholds, NICE cost-effectiveness cutoffs.

## Sharp RDD
rd <- rdrobust(y = d$y, x = d$running, c = 0)
summary(rd)                            # use the robust bias-corrected estimate
rdplot(y = d$y, x = d$running, c = 0,
       x.label = "Running variable", y.label = "Outcome")

## Optimal bandwidth
rdbwselect(y = d$y, x = d$running, c = 0, bwselect = "mserd")

## Manipulation test — density of the running variable at the cutoff
library(rddensity)
dt <- rddensity(X = d$running, c = 0)
summary(dt)
rdplotdensity(dt, X = d$running)
# A discontinuity in density means sorting; the design fails.

## Covariate balance at the cutoff — should show NO jump
rdrobust(y = d$age, x = d$running, c = 0)

## Placebo cutoffs — should show no effect
rdrobust(y = d$y, x = d$running, c = -10)
rdrobust(y = d$y, x = d$running, c =  10)

## Fuzzy RDD — crossing the threshold changes the probability of treatment
rdrobust(y = d$y, x = d$running, c = 0, fuzzy = d$treat)


# =============================================================================
# 6. PANEL DATA
# =============================================================================

pd <- pdata.frame(d, index = c("id", "period"))

fe <- plm(y ~ treat + age, data = pd, model = "within")   # fixed effects
re <- plm(y ~ treat + age, data = pd, model = "random")   # random effects
pool <- plm(y ~ treat + age, data = pd, model = "pooling")

## Hausman test: is RE consistent?
phtest(fe, re)
# p < 0.05 -> unit effects correlate with regressors -> use FE.
# FE removes all time-invariant confounding but also cannot estimate any
# time-invariant covariate (sex, baseline genotype).

## Clustered SEs are essential in panel data
coeftest(fe, vcov = vcovHC(fe, type = "HC1", cluster = "group"))

## fixest — faster, handles high-dimensional fixed effects
feols(y ~ treat + age | id + period, data = d, cluster = ~ id)

## Serial correlation and cross-sectional dependence
pbgtest(fe)
pcdtest(fe, test = "cd")


# =============================================================================
# 7. REPORTING CAUSAL ESTIMATES
# =============================================================================
# State, in this order:
#   - the estimand (ATE, ATT, LATE, RMST difference) and for whom
#   - the identifying assumption and why it is plausible HERE
#   - what was tested (balance, pre-trends, first-stage F, density)
#   - sensitivity to unmeasured confounding

## E-value: how strong would unmeasured confounding have to be?
library(EValue)
evalues.RR(est = 1.8, lo = 1.3, hi = 2.5)
evalues.OLS(est = 0.3, se = 0.1, sd = sd(d$y, na.rm = TRUE))

## Rosenbaum bounds for matched designs
# library(rbounds); psens(...)

modelsummary::modelsummary(
  list("OLS" = m, "PSM" = fit_m, "IPW" = fit_w, "IV" = iv),
  vcov = list("HC3", ~subclass, "HC3", "HC3"),
  stars = TRUE, statistic = "conf.int"
)
