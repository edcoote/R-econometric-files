# =============================================================================
# 02_linear_regression.R  —  OLS, robust/clustered SEs, diagnostics
# =============================================================================
# Order that saves time: fit, plot residuals, test, then re-estimate.
# Rule out misspecification BEFORE "correcting" heteroskedasticity or
# autocorrelation — a wrong functional form produces both.
# =============================================================================

library(tidyverse); library(broom)
library(lmtest); library(sandwich); library(car)
library(marginaleffects); library(emmeans)


# --- 1. Fit ------------------------------------------------------------------

m <- lm(y ~ treat + age + sex + bmi + grp, data = d)

summary(m)
tidy(m, conf.int = TRUE)      # coefficients as a data frame
glance(m)                     # R2, adj R2, AIC, BIC, df, n
augment(m)                    # .fitted .resid .hat .cooksd .std.resid

nobs(m)                       # ALWAYS check: did lm() drop rows to NA?

# Formula syntax
#   y ~ x1 + x2        main effects
#   y ~ x1 * x2        x1 + x2 + x1:x2
#   y ~ x1:x2          interaction only
#   y ~ . - id         everything except id
#   y ~ x - 1          no intercept
#   y ~ I(x^2)         literal square (^ means interaction inside a formula)
#   y ~ poly(x, 2)     orthogonal polynomial
#   y ~ x + offset(log(pt))   fixed coefficient of 1, for rate models


# --- 2. Functional form ------------------------------------------------------

m_quad  <- lm(y ~ age + I(age^2) + treat, data = d)
m_log   <- lm(log(y) ~ treat + age, data = d)        # β ≈ 100β % change in y
m_loglog<- lm(log(y) ~ log(x1), data = d)            # elasticity
m_spline<- lm(y ~ splines::ns(age, df = 4) + treat, data = d)

# Interpretation
#   level-level  y ~ x        1 unit x  -> β units y
#   log-level    log(y) ~ x   1 unit x  -> 100·β % in y
#   level-log    y ~ log(x)   1% in x   -> β/100 units y
#   log-log      log(y)~log(x) 1% in x  -> β% in y  (elasticity)
#   quadratic                 turning point at -β1 / (2·β2)

# Inverse hyperbolic sine — log-like but defined at zero (costs, utilisation)
ihs <- function(x) log(x + sqrt(x^2 + 1))

# Plot the fitted non-linearity rather than reading coefficients
plot_predictions(m_spline, condition = "age")


# --- 3. Interactions ---------------------------------------------------------

m_int <- lm(y ~ treat * age + sex, data = d)

# The interaction coefficient is not the effect of treatment. Report the
# marginal effect across the moderator instead.
plot_slopes(m_int, variables = "treat", condition = "age")
avg_slopes(m_int, variables = "treat")                       # average ME
slopes(m_int, variables = "treat",
       newdata = datagrid(age = c(50, 60, 70)))              # at set values

emmeans(m_int, ~ treat | age, at = list(age = c(50, 60, 70)))


# --- 4. Categorical predictors -----------------------------------------------

levels(d$grp)                 # [1] is the reference
contrasts(d$grp)
d$grp <- relevel(d$grp, ref = "Chosen reference")

# Test a multi-level factor as a single block, not level by level
car::Anova(m, type = "II")    # type III only if you have interactions AND
                              # sum-to-zero contrasts set

emmeans(m, ~ grp)                                   # adjusted means
pairs(emmeans(m, ~ grp), adjust = "tukey")          # all pairwise, corrected


# --- 5. Diagnostics ----------------------------------------------------------

par(mfrow = c(2, 2)); plot(m); par(mfrow = c(1, 1))
#   Residuals vs Fitted  flat band  | curve -> functional form; funnel -> hetero
#   Q-Q residuals        on the line| fat tails / skew
#   Scale-Location       flat       | rising -> heteroskedasticity
#   Residuals vs Leverage           | points past Cook's contours

## Multicollinearity
car::vif(m)                   # >5 watch, >10 serious
                              # factors: use GVIF^(1/(2*Df)) column
# Inflates SEs; does NOT bias β. If prediction is the goal it may not matter.

## Heteroskedasticity
bptest(m)                                              # Breusch-Pagan
bptest(m, ~ age*bmi + I(age^2) + I(bmi^2), data = d)   # White
skedastic::white(m)
skedastic::goldfeld_quandt(m)

## Autocorrelation (time series, repeated measures)
dwtest(m)                     # AR(1) only; invalid with a lagged DV
bgtest(m, order = 2)          # Breusch-Godfrey, more general
acf(resid(m))

## Normality of residuals — least urgent; CLT covers t and F at decent n
tseries::jarque.bera.test(resid(m))
shapiro.test(resid(m))        # n < 5000

## Influence
car::influencePlot(m)
which(cooks.distance(m) > 4 / nobs(m))
which(hatvalues(m) > 2 * length(coef(m)) / nobs(m))
car::outlierTest(m)
# High leverage is not the same as influential. Report fits with and without
# rather than deleting quietly.

## Specification
resettest(m, power = 2:3, type = "regressor")   # rejection -> something omitted
car::crPlots(m)                                 # component + residual
car::avPlots(m)                                 # added variable


# --- 6. Robust and clustered standard errors ---------------------------------
# Keep the OLS coefficients; correct the variance matrix.

coeftest(m, vcov = vcovHC(m, type = "HC3"))
#   HC0  large n        HC1  = Stata's robust
#   HC2                 HC3  best below n ~ 250; safe default

coefci(m, vcov = vcovHC(m, type = "HC3"))

# Clustered (practice, site, region). Needs enough clusters — below ~40,
# use wild cluster bootstrap instead.
coeftest(m, vcov = vcovCL(m, cluster = ~ clust))
fwildclusterboot::boottest(m, clustid = "clust", param = "treat", B = 9999)

# Heteroskedasticity- and autocorrelation-consistent
coeftest(m, vcov = NeweyWest(m))

# Robust SEs into a tidy frame
tidy(coeftest(m, vcov = vcovHC(m, "HC3")), conf.int = TRUE)


# --- 7. Weighted least squares -----------------------------------------------
# If you know the source of the non-constant variance, WLS is more efficient
# than robust SEs. If you do not, use robust SEs.

m_wls <- lm(y ~ treat + age, data = d, weights = 1 / x1^2)

# Survey weights are a different thing — use the survey package
library(survey)
des <- svydesign(ids = ~ psu, strata = ~ stratum, weights = ~ wt,
                 data = d, nest = TRUE)
svyglm(y ~ treat + age, design = des)


# --- 8. Model comparison -----------------------------------------------------

anova(m_reduced, m_full)                             # nested F test
linearHypothesis(m, c("treat = 0", "age = 0"))       # Wald
lrtest(m_reduced, m_full)                            # likelihood ratio
AIC(m1, m2); BIC(m1, m2)                             # lower is better

# Nested tests require IDENTICAL rows. Fit on d_analysis (fixed in 01) or the
# comparison is meaningless.

# Stepwise: acceptable for exploration, but the reported p-values are not
# adjusted for having searched the model space. Do not present them as inference.
step(m_full, direction = "backward")


# --- 9. Prediction -----------------------------------------------------------

nd <- data.frame(
  treat = factor("Treatment", levels = levels(d$treat)),
  age   = 60,
  sex   = factor("Female",    levels = levels(d$sex)),
  bmi   = 27,
  grp   = factor("Normal",    levels = levels(d$grp))
)

predict(m, nd, interval = "confidence")   # mean response — narrow
predict(m, nd, interval = "prediction")   # single new case — wide

# Always build newdata factors with levels = levels(d$f), or you get
# "factor has new levels".

predictions(m, newdata = datagrid(age = 40:80, treat = levels(d$treat)))


# --- 10. Quantile regression -------------------------------------------------
# For skewed outcomes (cost, length of stay) where the mean is not the target,
# and when you want effects across the distribution rather than at the mean.

library(quantreg)
mq <- rq(cost ~ treat + age, data = d, tau = c(0.25, 0.50, 0.75, 0.90))
summary(mq, se = "boot")
plot(summary(mq))
