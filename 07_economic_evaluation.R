# =============================================================================
# 07_economic_evaluation.R  —  Markov models, CUA, PSA, CEAC, VOI, BIA
# =============================================================================
# Written to match the NICE reference case; adapt discount rates, perspective
# and threshold for other jurisdictions.
# =============================================================================

library(tidyverse)
library(heemod)      # Markov cohort models
library(dampack)     # CEA analysis, CEAC, EVPI, tornado
library(BCEA)        # Bayesian-flavoured CEA output


# --- 1. Parameters -----------------------------------------------------------
# Hold everything in one place. Every parameter used in the model should be
# named here, with a source and a distribution for PSA.

params <- list(
  # horizon and discounting (NICE reference case: 3.5% both)
  cycle_length   = 1,          # years
  n_cycles       = 40,
  disc_c         = 0.035,
  disc_e         = 0.035,
  wtp            = 20000,      # £/QALY

  # transition probabilities (annual)
  p_prog         = 0.12,
  p_death_well   = 0.008,
  p_death_prog   = 0.15,
  rr_treat       = 0.72,       # from 06_evidence_synthesis.R

  # costs (£, current price year — state the year)
  c_treat        = 4200,
  c_well         = 320,
  c_prog         = 8900,
  c_death        = 2100,

  # utilities
  u_well         = 0.82,
  u_prog         = 0.54,
  u_dead         = 0
)


# --- 2. Discounting and half-cycle correction --------------------------------

disc_factor <- function(rate, cycles, cycle_length = 1) {
  1 / (1 + rate)^((cycles - 1) * cycle_length)
}

# Half-cycle correction: transitions occur throughout a cycle, not at its
# boundary. Apply via the trapezoid rule, or use shorter cycles and skip it.
hcc <- function(x) (head(x, -1) + tail(x, -1)) / 2


# --- 3. Markov cohort model, written directly --------------------------------
# Transparent and debuggable. Use heemod (section 4) for anything larger.

run_markov <- function(p, treat = FALSE) {

  states <- c("Well", "Progressed", "Dead")
  n_s <- length(states)

  # treatment effect applied to progression only
  p_prog <- if (treat) p$p_prog * p$rr_treat else p$p_prog

  # transition matrix — rows MUST sum to 1
  tm <- matrix(0, n_s, n_s, dimnames = list(states, states))
  tm["Well", "Progressed"] <- p_prog
  tm["Well", "Dead"]       <- p$p_death_well
  tm["Well", "Well"]       <- 1 - p_prog - p$p_death_well
  tm["Progressed", "Dead"] <- p$p_death_prog
  tm["Progressed", "Progressed"] <- 1 - p$p_death_prog
  tm["Dead", "Dead"]       <- 1

  stopifnot(all(abs(rowSums(tm) - 1) < 1e-10), all(tm >= 0))

  # trace
  trace <- matrix(0, p$n_cycles, n_s, dimnames = list(NULL, states))
  trace[1, ] <- c(1, 0, 0)
  for (i in 2:p$n_cycles) trace[i, ] <- trace[i - 1, ] %*% tm

  # payoffs per cycle
  c_state <- c(p$c_well + if (treat) p$c_treat else 0, p$c_prog, 0)
  u_state <- c(p$u_well, p$u_prog, p$u_dead)

  cost_cyc <- as.vector(trace %*% c_state)
  qaly_cyc <- as.vector(trace %*% u_state) * p$cycle_length

  df <- disc_factor(p$disc_c, 1:p$n_cycles, p$cycle_length)
  dfe <- disc_factor(p$disc_e, 1:p$n_cycles, p$cycle_length)

  list(
    trace = trace,
    cost  = sum(hcc(cost_cyc) * head(df, -1)),
    qaly  = sum(hcc(qaly_cyc) * head(dfe, -1))
  )
}

res_soc <- run_markov(params, treat = FALSE)
res_trt <- run_markov(params, treat = TRUE)

## Markov trace — always plot it. Most model errors are visible here.
as_tibble(res_trt$trace) %>%
  mutate(cycle = row_number()) %>%
  pivot_longer(-cycle) %>%
  ggplot(aes(cycle, value, fill = name)) +
  geom_area() +
  labs(x = "Cycle", y = "Proportion", fill = "State") +
  theme_classic()


# --- 4. Incremental analysis -------------------------------------------------

inc_cost <- res_trt$cost - res_soc$cost
inc_qaly <- res_trt$qaly - res_soc$qaly
icer     <- inc_cost / inc_qaly
nmb      <- inc_qaly * params$wtp - inc_cost     # net monetary benefit

tibble(
  Strategy = c("Standard care", "Treatment"),
  Cost     = c(res_soc$cost, res_trt$cost),
  QALYs    = c(res_soc$qaly, res_trt$qaly)
) %>%
  mutate(IncCost = Cost - lag(Cost),
         IncQALY = QALYs - lag(QALYs),
         ICER    = IncCost / IncQALY)

# With >2 strategies, rank by cost, remove dominated and extendedly dominated
# options, then compute sequential ICERs on the efficiency frontier.
# An ICER against a non-adjacent comparator is meaningless.
# dampack::calculate_icers(cost, effect, strategies)

# Report NMB alongside the ICER. ICERs behave badly when the increment is
# negative or near zero, and cannot be averaged across PSA runs.


# --- 5. Probabilistic sensitivity analysis -----------------------------------
# Distribution choice follows the parameter's support:
#   probability / utility (0-1)   beta
#   relative risk, HR, OR         lognormal (normal on the log scale)
#   cost                          gamma
#   multi-state transitions       Dirichlet
#   utility that can be < 0       normal, or beta shifted onto [-0.594, 1]

beta_mom <- function(mu, se) {
  v <- se^2
  a <- mu * (mu * (1 - mu) / v - 1)
  list(shape1 = a, shape2 = a * (1 - mu) / mu)
}

gamma_mom <- function(mu, se) {
  list(shape = (mu / se)^2, scale = se^2 / mu)
}

n_sim <- 10000
set.seed(1234)

psa_params <- tibble(
  p_prog  = do.call(rbeta,  c(list(n_sim), beta_mom(0.12, 0.02))),
  u_well  = do.call(rbeta,  c(list(n_sim), beta_mom(0.82, 0.03))),
  u_prog  = do.call(rbeta,  c(list(n_sim), beta_mom(0.54, 0.05))),
  c_prog  = do.call(rgamma, c(list(n_sim), gamma_mom(8900, 1200))),
  c_treat = do.call(rgamma, c(list(n_sim), gamma_mom(4200, 400))),
  rr_treat = exp(rnorm(n_sim, log(0.72), 0.12))
)

psa <- psa_params %>%
  mutate(row = row_number()) %>%
  rowwise() %>%
  mutate({
    p <- modifyList(params, list(
      p_prog = p_prog, u_well = u_well, u_prog = u_prog,
      c_prog = c_prog, c_treat = c_treat, rr_treat = rr_treat))
    s <- run_markov(p, FALSE); t <- run_markov(p, TRUE)
    tibble(c_soc = s$cost, q_soc = s$qaly,
           c_trt = t$cost, q_trt = t$qaly)
  }) %>%
  ungroup() %>%
  mutate(inc_c = c_trt - c_soc,
         inc_q = q_trt - q_soc,
         nmb   = inc_q * params$wtp - inc_c)

## Cost-effectiveness plane
ggplot(psa, aes(inc_q, inc_c)) +
  geom_point(alpha = 0.15) +
  geom_abline(slope = params$wtp, linetype = "dashed", colour = "red") +
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
  labs(x = "Incremental QALYs", y = "Incremental cost (£)",
       caption = "Dashed line = £20,000/QALY threshold") +
  theme_classic()

## CEAC — probability cost-effective across thresholds
thresholds <- seq(0, 100000, by = 1000)
ceac <- map_dfr(thresholds, ~ tibble(
  wtp = .x,
  p_ce = mean(psa$inc_q * .x - psa$inc_c > 0)
))

ggplot(ceac, aes(wtp, p_ce)) +
  geom_line() +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Willingness to pay (£/QALY)",
       y = "Probability cost-effective") +
  theme_classic()

## Report the mean of the PSA results as the base case, not the deterministic
## run — non-linearity means they differ, and the mean is the decision-relevant
## quantity.
psa %>% summarise(across(c(inc_c, inc_q, nmb), mean))


# --- 6. Deterministic sensitivity analysis -----------------------------------

owsa <- function(param, lo, hi) {
  map_dfr(c(lo, hi), function(v) {
    p <- modifyList(params, setNames(list(v), param))
    s <- run_markov(p, FALSE); t <- run_markov(p, TRUE)
    tibble(param = param, value = v,
           nmb = (t$qaly - s$qaly) * params$wtp - (t$cost - s$cost))
  })
}

tornado <- bind_rows(
  owsa("p_prog",   0.08,  0.16),
  owsa("u_prog",   0.44,  0.64),
  owsa("c_prog",   6000,  12000),
  owsa("rr_treat", 0.60,  0.88)
)

tornado %>%
  group_by(param) %>%
  summarise(lo = min(nmb), hi = max(nmb), range = hi - lo) %>%
  mutate(param = fct_reorder(param, range)) %>%
  ggplot(aes(y = param)) +
  geom_segment(aes(x = lo, xend = hi, yend = param), linewidth = 6,
               colour = "steelblue") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x = "Incremental net monetary benefit (£)", y = NULL) +
  theme_classic()


# --- 7. Value of information -------------------------------------------------
# EVPI: the value of eliminating ALL parameter uncertainty. An upper bound on
# what further research could be worth.

evpi_calc <- function(nmb_soc, nmb_trt) {
  mean(pmax(nmb_soc, nmb_trt)) - max(mean(nmb_soc), mean(nmb_trt))
}

psa <- psa %>%
  mutate(nmb_soc = q_soc * params$wtp - c_soc,
         nmb_trt = q_trt * params$wtp - c_trt)

evpi_pp <- evpi_calc(psa$nmb_soc, psa$nmb_trt)     # per patient
n_eligible  <- 25000
years_tech  <- 10
evpi_pop <- evpi_pp * n_eligible *
  sum(1 / (1 + params$disc_e)^(0:(years_tech - 1)))

# EVPI across thresholds — EVPI peaks where the decision is closest
map_dfr(thresholds, function(k) {
  s <- psa$q_soc * k - psa$c_soc
  t <- psa$q_trt * k - psa$c_trt
  tibble(wtp = k, evpi = evpi_calc(s, t))
}) %>%
  ggplot(aes(wtp, evpi)) + geom_line() +
  labs(x = "Willingness to pay (£/QALY)", y = "EVPI per patient (£)") +
  theme_classic()

## EVPPI — value of resolving uncertainty in a SUBSET of parameters, which is
## what tells you where to spend research money. Use a regression-based method.
library(voi)
evppi(outputs = list(e = cbind(psa$q_soc, psa$q_trt),
                     c = cbind(psa$c_soc, psa$c_trt),
                     k = params$wtp),
      inputs = psa_params,
      pars = list("u_prog", "rr_treat"))

# EVSI (value of a specific proposed study) via voi::evsi()


# --- 8. heemod, for larger models --------------------------------------------
# Handles time-dependent transitions, tunnel states, state residence time,
# and PSA natively.

# tm <- define_transition(
#   state_names = c("Well","Progressed","Dead"),
#   C, p_prog, p_death_well,
#   0, C,      p_death_prog,
#   0, 0,      1
# )
# s_well <- define_state(cost = discount(c_well, .035),
#                        qaly = discount(u_well, .035))
# strat  <- define_strategy(transition = tm, Well = s_well, ...)
# res    <- run_model(soc = strat_soc, trt = strat_trt,
#                     cycles = 40, cost = cost, effect = qaly)
# summary(res); plot(res, type = "counts")


# --- 9. Budget impact --------------------------------------------------------
# Distinct from cost-effectiveness: short horizon (typically 5 years),
# NO discounting of costs in most guidance, actual population not a cohort,
# and uptake modelled explicitly.

bia <- tibble(
  year    = 1:5,
  pop     = 25000 * 1.02^(0:4),
  uptake  = c(0.05, 0.15, 0.30, 0.45, 0.55)
) %>%
  mutate(
    n_treated   = pop * uptake,
    cost_new    = n_treated * params$c_treat,
    cost_offset = n_treated * params$p_prog * params$rr_treat * params$c_prog,
    net_impact  = cost_new - cost_offset
  )

bia
sum(bia$net_impact)


# --- 10. Reporting checklist -------------------------------------------------
# CHEERS 2022. State explicitly:
#   perspective (NHS/PSS vs societal) and price year
#   time horizon and why
#   discount rates for costs and outcomes
#   source of every transition probability, cost and utility
#   how uncertainty was characterised (PSA + DSA + scenarios)
#   structural assumptions and structural sensitivity analyses
#   whether treatment effect was assumed to persist, and for how long
