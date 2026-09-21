# =============================================================================
# 03_glm_outcomes.R  —  Binary, ordinal, multinomial, count, cost outcomes
# =============================================================================
# Pick the family from the outcome, then report on a scale people can read:
# odds ratios, risk ratios, incidence rate ratios, or marginal effects.
# =============================================================================

library(tidyverse); library(broom)
library(MASS)        # glm.nb, polr — load BEFORE tidyverse or use dplyr::select
library(nnet)        # multinom
library(pscl)        # hurdle, zeroinfl, pR2
library(marginaleffects)
library(sandwich); library(lmtest)


# =============================================================================
# BINARY OUTCOMES
# =============================================================================
# glm() models the probability of the SECOND factor level. Check it:
levels(d$ybin)       # want c("No", "Yes") so that "Yes" is the event

mb <- glm(ybin ~ treat + age + sex + bmi, data = d, family = binomial)
summary(mb)

## Odds ratios
exp(cbind(OR = coef(mb), confint(mb)))                # profile likelihood CI
tidy(mb, exponentiate = TRUE, conf.int = TRUE)

## Risk ratios — usually what a clinical audience wants, and OR overstates RR
## badly when the outcome is common (>10%).
mrr <- glm(ybin ~ treat + age, data = d, family = poisson(link = "log"))
coeftest(mrr, vcov = vcovHC(mrr, "HC3"))              # robust SEs are REQUIRED
exp(cbind(RR = coef(mrr), coefci(mrr, vcov = vcovHC(mrr, "HC3"))))

## Risk difference — the most directly interpretable for policy
avg_comparisons(mb, variables = "treat")              # average marginal effect
avg_comparisons(mb, variables = "treat", type = "response")

## Predicted probabilities
predictions(mb, newdata = datagrid(age = 40:80, treat = levels(d$treat)))
plot_predictions(mb, condition = c("age", "treat"))

## Other links
glm(ybin ~ treat, data = d, family = binomial(link = "probit"))
glm(ybin ~ treat, data = d, family = binomial(link = "cloglog"))

## Fit and discrimination
pscl::pR2(mb)                                          # pseudo R2
performance::performance_hosmer(mb, n_bins = 10)       # calibration
library(pROC)
roc_obj <- roc(d$ybin, fitted(mb)); auc(roc_obj); plot(roc_obj)
performance::check_collinearity(mb)

## Separation: "fitted probabilities numerically 0 or 1 occurred"
## A predictor splits the outcome perfectly. Check table(d$ybin, d$x), then:
library(logistf)
logistf(ybin ~ treat + age, data = d)                  # Firth penalised


# =============================================================================
# ORDERED CATEGORICAL  (severity grades, Likert, EQ-5D dimensions)
# =============================================================================
d$yord <- factor(d$yord, levels = c("None","Mild","Moderate","Severe"),
                 ordered = TRUE)

mo <- MASS::polr(yord ~ treat + age + sex, data = d,
                 method = "logistic", Hess = TRUE)   # Hess = TRUE or no SEs
summary(mo)
exp(cbind(OR = coef(mo), confint(mo)))               # proportional odds ratios

## Test the proportional odds assumption — the whole model rests on it
library(brant)
brant::brant(mo)                                     # p < 0.05 -> violated

## If violated: partial proportional odds, or a multinomial model
library(VGAM)
vglm(yord ~ treat + age, family = cumulative(parallel = FALSE), data = d)


# =============================================================================
# NOMINAL CATEGORICAL  (unordered: treatment chosen, site of care)
# =============================================================================
d$ynom <- relevel(factor(d$ynom), ref = "Reference category")

mm <- nnet::multinom(ynom ~ treat + age + sex, data = d, trace = FALSE)
summary(mm)

z <- summary(mm)$coefficients / summary(mm)$standard.errors
(1 - pnorm(abs(z))) * 2                               # Wald p-values
exp(coef(mm))                                         # relative risk ratios
tidy(mm, exponentiate = TRUE, conf.int = TRUE)


# =============================================================================
# COUNT OUTCOMES  (admissions, GP visits, events per person)
# =============================================================================

## Poisson — assumes mean = variance (equidispersion)
mp <- glm(ycount ~ treat + age + sex, data = d, family = poisson)
exp(cbind(IRR = coef(mp), confint(mp)))

## Rate model — when follow-up time differs between people. The offset enters
## with a fixed coefficient of 1, making the outcome events per person-time.
mp_rate <- glm(ycount ~ treat + age + offset(log(pt_years)),
               data = d, family = poisson)

## Is Poisson adequate?
AER::dispersiontest(mp)                                # H0: equidispersion
summary(mp)$deviance / summary(mp)$df.residual         # >> 1 -> overdispersed
performance::check_overdispersion(mp)

## Overdispersion does NOT bias the Poisson coefficients, but it understates
## their standard errors — everything looks more significant than it is.

## Negative binomial — the usual fix
mnb <- MASS::glm.nb(ycount ~ treat + age + sex, data = d)
exp(cbind(IRR = coef(mnb), confint(mnb)))
mnb$theta                                              # dispersion; small = more

## Quasi-Poisson — rescales SEs only. No likelihood, so no AIC.
mqp <- glm(ycount ~ treat + age, data = d, family = quasipoisson)

AIC(mp, mnb)                                           # NB vs Poisson
lmtest::lrtest(mp, mnb)


# =============================================================================
# EXCESS ZEROS
# =============================================================================
# Distinguish the two mechanisms before choosing:
#   Hurdle       — one process decides any/none, a second decides how many
#                  (all zeros come from the first stage)
#   Zero-inflated— a structurally never-at-risk group mixed with a count process
#                  (zeros arise from both)
# For healthcare utilisation, hurdle usually matches the story better: the
# decision to seek care differs from the intensity of care received.

mh <- pscl::hurdle(ycount ~ treat + age | treat + age,
                   data = d, dist = "negbin")
summary(mh)

mzi <- pscl::zeroinfl(ycount ~ treat + age | treat + age,
                      data = d, dist = "negbin")
summary(mzi)
# Left of | = count model; right of | = zero model.

AIC(mnb, mh, mzi)
pscl::vuong(mnb, mzi)          # non-nested comparison

## Predicted vs observed zero counts — the honest check
sum(d$ycount == 0, na.rm = TRUE)
sum(dpois(0, fitted(mp)))      # Poisson badly under-predicts if zeros matter
sum(predict(mh, type = "prob")[, 1])


# =============================================================================
# COST AND UTILISATION OUTCOMES
# =============================================================================
# Right-skewed, non-negative, often with a zero spike. Log-OLS is biased on
# retransformation to the natural scale (the smearing problem) — GLM with a
# log link avoids it by modelling E[y|x] directly.

## GLM with log link. Gamma for costs; check the variance function.
mg <- glm(cost ~ treat + age + sex, data = d,
          family = Gamma(link = "log"))
exp(coef(mg))                                          # cost ratios

## Choosing the family: modified Park test regresses squared raw-scale
## residuals on the fitted mean. Slope ~1 Poisson, ~2 gamma, ~3 inverse Gaussian.
r2 <- residuals(mg, type = "response")^2
coef(glm(r2 ~ log(fitted(mg)), family = Gamma(link = "log")))[2]

## Two-part model: P(any cost) x E[cost | cost > 0]
p1 <- glm(I(cost > 0) ~ treat + age, data = d, family = binomial)
p2 <- glm(cost ~ treat + age, data = filter(d, cost > 0),
          family = Gamma(link = "log"))
# Expected cost = predicted probability x predicted conditional mean
e_cost <- predict(p1, type = "response") *
          predict(p2, newdata = d, type = "response")

## Report incremental cost on the natural scale with a bootstrap CI —
## coefficients from a log-link model are ratios, not differences.
avg_comparisons(mg, variables = "treat", type = "response")


# =============================================================================
# REPORTING
# =============================================================================
library(gtsummary)

tbl_regression(mb, exponentiate = TRUE) %>%
  add_global_p() %>%                    # factors tested as a block
  add_n() %>%
  bold_p(t = 0.05) %>%
  modify_caption("**Adjusted odds ratios**")

# Side by side
tbl_merge(
  list(tbl_regression(mb, exponentiate = TRUE),
       tbl_regression(mnb, exponentiate = TRUE)),
  tab_spanner = c("**Binary outcome (OR)**", "**Count outcome (IRR)**")
)

modelsummary::modelsummary(
  list("Poisson" = mp, "Negative binomial" = mnb, "Hurdle" = mh),
  exponentiate = TRUE, statistic = "conf.int", stars = TRUE
)
