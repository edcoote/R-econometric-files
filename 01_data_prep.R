# =============================================================================
# 01_data_prep.R  —  Types, factors, missing data, derived variables, Table 1
# =============================================================================
# The stage where most analyses go wrong silently. Every operation here is
# written to be run ONCE, top to bottom, from a fresh import.
#
# CARDINAL RULE: factor() applied to an existing factor maps LEVEL INDICES
# (1, 2, 3...), not your original codes. Re-running this script on an already
# converted object turns data into NA without warning. Always start from a
# fresh read of the raw file.
# =============================================================================

library(tidyverse)
library(gtsummary)


# --- 1. Fix types ------------------------------------------------------------

# Character column that should be numeric. Check blockers() from 00_setup first.
d$x1 <- as.numeric(d$x1)

# Factor that should be numeric — as.character() step is ESSENTIAL.
# as.numeric(f) alone returns level indices and produces plausible wrong numbers.
d$x2 <- as.numeric(as.character(d$x2))

# Text that needs parsing rather than coercion
d <- d %>%
  mutate(
    cost = readr::parse_number(cost_raw),              # strips £ and commas
    sbp  = as.numeric(sub("/.*", "", bp_raw)),         # "140/90" -> 140
    dbp  = as.numeric(sub(".*/", "", bp_raw))
  )

# Dates
d <- d %>%
  mutate(
    date_index = as.Date(date_index, format = "%d/%m/%Y"),
    date_event = as.Date(date_event, format = "%d/%m/%Y"),
    fu_days    = as.numeric(date_event - date_index),
    fu_years   = fu_days / 365.25
  )


# --- 2. Factors --------------------------------------------------------------

# ALWAYS specify levels as well as labels. With levels given, any value not
# listed becomes NA — which is what you want for "don't know" codes. Without
# levels, labels are applied in sorted order and a stray code throws:
#   "invalid 'labels'; length 2 should be 1 or 3"

d$treat <- factor(d$treat, levels = c(0, 1), labels = c("Control", "Treatment"))

# Several variables sharing one coding scheme
yn_vars <- c("cvd", "t2dm", "cancer", "currsmoker")

stopifnot(all(yn_vars %in% names(d)))                    # fail early on typos
stopifnot(all(sapply(d[yn_vars], is.numeric)))           # guard double-conversion

d[yn_vars] <- lapply(d[yn_vars], function(x)
  factor(x, levels = c(0, 1), labels = c("No", "Yes")))

# Ordered categorical
d$bmi_grp <- factor(
  d$bmi_grp,
  levels = 1:4,
  labels = c("Underweight", "Normal", "Overweight", "Obese")
)

# Reference level decides what every coefficient is compared against.
# Choose it deliberately — the default is whichever level sorts first.
d$bmi_grp <- relevel(d$bmi_grp, ref = "Normal")

# Ordered factor, where the ordering is substantive (for polr, trend tests)
d$yord <- factor(d$yord, levels = c("Mild", "Moderate", "Severe"), ordered = TRUE)

# forcats helpers
d <- d %>%
  mutate(
    grp = fct_relevel(grp, "Reference level"),
    grp = fct_lump_min(grp, min = 20, other_level = "Other"),  # rare categories
    grp = fct_recode(grp, "New name" = "old_name"),
    grp = fct_drop(grp)                                        # unused levels
  )

# Verify
sapply(d[c(yn_vars, "treat", "bmi_grp")], levels)
map(d[yn_vars], ~ table(.x, useNA = "ifany"))


# --- 3. Derived variables ----------------------------------------------------

d <- d %>%
  mutate(
    bmi      = weight_kg / (height_m^2),
    age_c    = age - mean(age, na.rm = TRUE),              # centre before
    age_grp  = cut(age, breaks = c(0, 45, 55, 65, 75, Inf),
                   labels = c("<45","45-54","55-64","65-74","75+"),
                   right = FALSE),
    log_cost = log(cost + 1)                               # +1 handles zeros
  )

# Comorbidity count. Counting factors requires comparison, not addition:
# "Yes" + "No" is not defined and returns NA with a warning.
d <- d %>%
  mutate(
    n_comorbid = (cvd == "Yes") + (t2dm == "Yes") + (cancer == "Yes")
  )
# NA on any component makes the whole count NA. That is usually correct —
# rowSums(..., na.rm = TRUE) would silently treat unknown as absent.

# Conditional recoding
d <- d %>%
  mutate(
    risk = case_when(
      n_comorbid == 0              ~ "Low",
      n_comorbid %in% 1:2          ~ "Moderate",
      n_comorbid >= 3              ~ "High",
      TRUE                         ~ NA_character_
    ),
    risk = factor(risk, levels = c("Low", "Moderate", "High"))
  )


# --- 4. Plausibility ---------------------------------------------------------
# Coercion produces numbers, not errors. Range-check every clinical variable.

d %>%
  select(age, bmi, sbp, chol) %>%
  summary()

d %>% summarise(
  age_bad = sum(age  < 18  | age  > 110, na.rm = TRUE),
  bmi_bad = sum(bmi  < 12  | bmi  > 70,  na.rm = TRUE),
  sbp_bad = sum(sbp  < 60  | sbp  > 260, na.rm = TRUE)
)

# Flag rather than delete, so the decision stays visible and reversible
d <- d %>% mutate(bmi = if_else(bmi < 12 | bmi > 70, NA_real_, bmi))


# --- 5. Missing data ---------------------------------------------------------

colSums(is.na(d))
naniar::vis_miss(d)                    # visual pattern
naniar::gg_miss_upset(d)               # co-occurrence of missingness

# Is missingness related to the outcome? If so, complete-case analysis is biased.
d %>%
  mutate(miss_x1 = is.na(x1)) %>%
  group_by(miss_x1) %>%
  summarise(across(c(age, y), ~ mean(.x, na.rm = TRUE)), n = n())

# Multiple imputation by chained equations
library(mice)
imp <- mice(d, m = 20, method = "pmm", seed = 1234, printFlag = FALSE)
plot(imp)                              # convergence
fit <- with(imp, lm(y ~ treat + age + sex))
pooled <- pool(fit)                    # Rubin's rules
summary(pooled, conf.int = TRUE)

# Include the outcome in the imputation model but not in downstream analysis
# decisions; exclude identifiers and perfectly collinear derived variables.


# --- 6. Analysis set ---------------------------------------------------------
# Fix the analysis sample ONCE. lm() and glm() drop rows with NA on any model
# variable, so models with different covariate sets have different n — which
# invalidates anova(), AIC() and likelihood ratio comparisons between them.

model_vars <- c("y", "treat", "age", "sex", "bmi", "grp")

d_analysis <- d %>%
  filter(complete.cases(select(., all_of(model_vars))))

nrow(d); nrow(d_analysis)

saveRDS(d_analysis, here::here("data", "derived", "analysis.rds"))


# --- 7. Table 1 --------------------------------------------------------------

tbl1 <- d_analysis %>%
  select(age, sex, bmi, bmi_grp, n_comorbid, treat) %>%
  tbl_summary(
    by = treat,
    statistic = list(all_continuous()  ~ "{mean} ({sd})",
                     all_categorical() ~ "{n} ({p}%)"),
    digits    = all_continuous() ~ 1,
    missing   = "ifany",
    missing_text = "Missing"
  ) %>%
  add_n() %>%
  add_overall() %>%
  modify_caption("**Table 1. Baseline characteristics**") %>%
  bold_labels()

tbl1

# For a randomised comparison report standardised mean differences, not
# p-values — a significance test on baseline imbalance tests randomisation,
# which you already know occurred.
# tbl1 %>% add_difference(everything() ~ "smd")

tbl1 %>% as_gt() %>% gt::gtsave(here::here("output","tables","table1.docx"))
