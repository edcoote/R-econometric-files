# =============================================================================
# 04_survival.R  —  KM, Cox, parametric models, extrapolation, competing risks
# =============================================================================
# For economic evaluation the parametric extrapolation matters more than the
# in-sample fit: models that agree over observed follow-up routinely diverge
# two- or three-fold in the tail, and the tail drives lifetime results.
# =============================================================================

library(tidyverse); library(broom)
library(survival); library(ggsurvfit); library(survminer)
library(flexsurv); library(tidycmprsk); library(gtsummary)


# --- 1. The Surv object ------------------------------------------------------
# Survival data needs two pieces per person: time observed, and whether the
# event occurred or they were censored.
#
#   status: 1 = event, 0 = censored      <- standard
#   Some datasets (survival::lung) use 1 = censored, 2 = event. Recode ONCE:
#   d <- d %>% mutate(status = recode(status, `1` = 0, `2` = 1))
#   Re-running that maps 1 -> 0 again and silently censors all your events.

Surv(d$time, d$status)[1:10]        # "+" marks censored observations

table(d$status, useNA = "ifany")    # confirm event count before modelling

# Interval censored, or delayed entry (left truncation)
# Surv(time1, time2, type = "interval2")
# Surv(entry_time, exit_time, status)     # counting-process form


# --- 2. Kaplan-Meier ---------------------------------------------------------

km <- survfit2(Surv(time, status) ~ treat, data = d)

km                                        # median survival with 95% CI
summary(km, times = c(365, 730, 1825))    # survival at fixed timepoints

km %>%
  ggsurvfit() +
  add_confidence_interval() +
  add_risktable() +
  add_quantile(y_value = 0.5, linetype = "dashed") +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Days since index", y = "Overall survival probability") +
  theme_classic()

## Log-rank test
survdiff(Surv(time, status) ~ treat, data = d)
# Assumes proportional hazards. If the curves cross it loses power badly and
# can return a null result despite a real difference. Plot before you test.

## Weighted alternatives when hazards are non-proportional
# survMisc::comp()      Gehan-Wilcoxon (early differences), Fleming-Harrington

## Restricted mean survival time — no PH assumption, and directly interpretable
## as life-years gained within the horizon. Often the better estimand for HTA.
library(survRM2)
rmst2(d$time, d$status, arm = as.numeric(d$treat) - 1, tau = 1825)


# --- 3. Cox proportional hazards ---------------------------------------------

cox <- coxph(Surv(time, status) ~ treat + age + sex + bmi, data = d)
summary(cox)
tidy(cox, exponentiate = TRUE, conf.int = TRUE)      # hazard ratios

tbl_regression(cox, exponentiate = TRUE) %>% add_global_p()

## Test the proportional hazards assumption — this is not optional
ph <- cox.zph(cox)
ph                                    # p < 0.05 -> violated for that covariate
plot(ph)                              # Schoenfeld residuals; want a flat line
survminer::ggcoxzph(ph)

## If PH is violated, options in rough order of preference:
# (a) stratify on the offending variable — removes the assumption for it,
#     but you lose its hazard ratio
cox_s <- coxph(Surv(time, status) ~ treat + age + strata(sex), data = d)

# (b) time-varying coefficient
cox_tv <- coxph(Surv(time, status) ~ treat + age + tt(age), data = d,
                tt = function(x, t, ...) x * log(t))

# (c) report RMST or a fully parametric AFT model instead
# (d) fit separately within time periods (episode splitting via survSplit)

## Other diagnostics
ggcoxdiagnostics(cox, type = "martingale")   # functional form of covariates
ggcoxdiagnostics(cox, type = "dfbeta")       # influential observations
concordance(cox)                             # Harrell's C

## Time-varying covariates (exposure that changes during follow-up)
d_long <- tmerge(d, d, id = id, endpt = event(time, status))
d_long <- tmerge(d_long, exposure_data, id = id, on_treat = tdc(switch_time))
coxph(Surv(tstart, tstop, endpt) ~ on_treat + age, data = d_long)

## Clustered / recurrent events
coxph(Surv(time, status) ~ treat + cluster(clust), data = d)
coxph(Surv(time, status) ~ treat + frailty(id), data = d)


# --- 4. Parametric survival models -------------------------------------------
# Required for extrapolation beyond trial follow-up. Fit the standard set and
# compare on fit AND plausibility of the extrapolated hazard.

fits <- list(
  exp      = flexsurvreg(Surv(time, status) ~ treat, data = d, dist = "exp"),
  weibull  = flexsurvreg(Surv(time, status) ~ treat, data = d, dist = "weibull"),
  gompertz = flexsurvreg(Surv(time, status) ~ treat, data = d, dist = "gompertz"),
  lnorm    = flexsurvreg(Surv(time, status) ~ treat, data = d, dist = "lnorm"),
  llogis   = flexsurvreg(Surv(time, status) ~ treat, data = d, dist = "llogis"),
  gengamma = flexsurvreg(Surv(time, status) ~ treat, data = d, dist = "gengamma")
)

map_dfr(fits, ~ tibble(AIC = AIC(.x), BIC = BIC(.x)), .id = "dist") %>%
  arrange(AIC)

## Hazard shapes — pick on clinical plausibility, not AIC alone
#   exponential   constant hazard
#   Weibull       monotonic increasing or decreasing (shape >1 / <1)
#   Gompertz      exponentially increasing or decreasing; common for mortality
#   log-normal    hazard rises then falls
#   log-logistic  hazard rises then falls; heavier tail than log-normal
#   gen. gamma    flexible, nests Weibull / log-normal / gamma

## Visual fit against the KM
plot(fits$weibull, ci = FALSE, col = "red", xlab = "Days", ylab = "Survival")
lines(fits$gompertz, ci = FALSE, col = "blue")
legend("topright", c("Kaplan-Meier","Weibull","Gompertz"),
       col = c("black","red","blue"), lty = 1, bty = "n")

## Formal shape diagnostics
plot(survfit(Surv(time, status) ~ 1, data = d), fun = "cloglog")
# straight line          -> Weibull
# straight with slope 1  -> exponential
plot(survfit(Surv(time, status) ~ treat, data = d), fun = "cloglog")
# parallel lines         -> proportional hazards holds

## EXTRAPOLATION — where the money is
t_horizon <- seq(0, 365.25 * 40, by = 30)

extrap <- map_dfr(fits, function(f) {
  s <- summary(f, t = t_horizon, ci = FALSE, tidy = TRUE)
  as_tibble(s)
}, .id = "dist")

ggplot(extrap, aes(time / 365.25, est, colour = dist)) +
  geom_line() +
  geom_vline(xintercept = max(d$time) / 365.25, linetype = "dashed") +
  labs(x = "Years", y = "Survival", colour = "Distribution",
       caption = "Dashed line = end of observed follow-up") +
  theme_classic()

## Checks before adopting an extrapolation (NICE DSU TSD 14):
##  - Does the implied long-term hazard make clinical sense?
##  - Does modelled survival ever exceed general population survival?
##    Blend with a lifetable if so.
##  - How much do lifetime QALYs change across distributions? Report it.
##  - Is the treatment effect assumed to persist forever? Justify or waive it.

## General population mortality constraint
# library(demography); lifetable <- ...
# s_model <- pmin(s_model, s_genpop)


# --- 5. Accelerated failure time ---------------------------------------------
# Models time directly rather than the hazard. Coefficients are time ratios,
# which some audiences find more intuitive than hazard ratios.

aft <- survreg(Surv(time, status) ~ treat + age, data = d, dist = "weibull")
summary(aft)
exp(coef(aft))                       # time ratio: >1 = longer survival


# --- 6. Competing risks ------------------------------------------------------
# When a competing event prevents the event of interest (death from other
# causes blocking cancer death). Standard KM OVERESTIMATES incidence here
# because it treats competing events as censored, i.e. as if those people
# remained at risk.

d$status_cr <- factor(d$status_cr, levels = c(0, 1, 2),
                      labels = c("Censored", "Event", "Competing"))

## Cumulative incidence function — use instead of 1 - KM
cif <- tidycmprsk::cuminc(Surv(time, status_cr) ~ treat, data = d)
cif
ggcuminc(cif, outcome = "Event") +
  add_confidence_interval() +
  add_risktable() +
  labs(x = "Days", y = "Cumulative incidence")

## Fine-Gray subdistribution hazard — for predicting absolute risk
fg <- tidycmprsk::crr(Surv(time, status_cr) ~ treat + age, data = d)
tbl_regression(fg, exponentiate = TRUE)

## Cause-specific hazard — for aetiology. Fit a Cox model per cause,
## censoring the other. Use this when asking "does treatment affect this
## cause?"; use Fine-Gray when asking "what is this patient's absolute risk?"
coxph(Surv(time, status_cr == "Event") ~ treat + age, data = d)


# --- 7. Transition probabilities for a Markov model --------------------------
# Converting survival output into cycle transition probabilities.
# p = 1 - exp(-rate * cycle_length) for a constant rate.
# For a time-varying hazard, derive per cycle from the fitted survivor function:

cycle_length <- 1/12                 # monthly cycles, time in years
cycles <- 0:(40 * 12)
tt <- cycles * cycle_length

s <- summary(fits$weibull, t = tt, ci = FALSE, tidy = TRUE)$est
tp <- 1 - (s[-1] / s[-length(s)])    # conditional probability of transition
                                     # within each cycle
head(tp)

# Sanity: probabilities must lie in [0,1] and be non-NaN at t = 0.
stopifnot(all(tp >= 0 & tp <= 1, na.rm = TRUE))
