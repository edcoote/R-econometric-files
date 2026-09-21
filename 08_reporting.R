# =============================================================================
# 08_reporting.R  —  Tables, figures, export, reproducibility
# =============================================================================
# Getting results out of R and into a manuscript, HTA submission or slide deck
# without manual retyping. Manual transcription is where numbers go wrong.
# =============================================================================

library(tidyverse)
library(gtsummary); library(gt); library(modelsummary); library(flextable)
library(patchwork)


# --- 1. Global table settings ------------------------------------------------

theme_gtsummary_journal("jama")       # or "lancet", "nejm", "qjecon"
theme_gtsummary_compact()

# Reset with: reset_gtsummary_theme()


# --- 2. Regression tables ----------------------------------------------------

tbl_regression(m, exponentiate = TRUE) %>%
  add_global_p() %>%                  # factors tested as a block, not level-wise
  add_n(location = "level") %>%
  add_nevent() %>%
  bold_labels() %>%
  italicize_levels() %>%
  modify_header(estimate ~ "**aOR (95% CI)**") %>%
  modify_caption("**Table 2. Adjusted associations**") %>%
  modify_footnote(estimate ~ "Adjusted for age, sex and comorbidity count.")

## Unadjusted and adjusted side by side — what most journals expect
tbl_merge(
  list(
    tbl_uvregression(d, y = ybin, method = glm,
                     method.args = list(family = binomial),
                     exponentiate = TRUE,
                     include = c(treat, age, sex, bmi)),
    tbl_regression(mb, exponentiate = TRUE)
  ),
  tab_spanner = c("**Unadjusted**", "**Adjusted**")
)

## Several models across columns
modelsummary(
  list("OLS" = m1, "Robust SE" = m2, "IPW" = m3),
  vcov      = list(NULL, "HC3", "HC3"),
  statistic = "conf.int",
  stars     = c("*" = .05, "**" = .01, "***" = .001),
  gof_map   = c("nobs", "r.squared", "adj.r.squared", "aic"),
  coef_rename = c("treatTreatment" = "Treatment",
                  "age" = "Age (years)"),
  output    = here::here("output", "tables", "models.docx")
)

# Report confidence intervals rather than stars where the journal allows it.
# Effect size and precision carry more information than a significance marker.


# --- 3. Export ---------------------------------------------------------------

tbl %>% as_gt()        %>% gt::gtsave(here::here("output","tables","t1.docx"))
tbl %>% as_gt()        %>% gt::gtsave(here::here("output","tables","t1.html"))
tbl %>% as_flex_table() %>%
  flextable::save_as_docx(path = here::here("output","tables","t1.docx"))

# Excel, for onward analysis or a model appendix
writexl::write_xlsx(
  list(Baseline = tbl1_df, Results = results_df, PSA = psa_summary),
  here::here("output", "tables", "results.xlsx")
)

# Always also write the underlying numbers as CSV. A formatted table is not
# a data source, and reviewers ask for the raw values.
readr::write_csv(results_df, here::here("output","tables","results_raw.csv"))


# --- 4. Figure theme ---------------------------------------------------------

theme_pub <- function(base_size = 10) {
  theme_classic(base_size = base_size) +
    theme(
      panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3),
      axis.title  = element_text(face = "plain"),
      plot.title  = element_text(face = "bold", size = rel(1.1)),
      plot.caption = element_text(colour = "grey40", hjust = 0),
      legend.position = "bottom",
      legend.title = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold")
    )
}
theme_set(theme_pub())

# Colourblind-safe palettes. Never let colour alone carry meaning —
# vary linetype or shape as well, so the figure survives greyscale printing.
scale_colour_viridis_d(option = "D", end = 0.9)
# c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00")   Okabe-Ito


# --- 5. Common figures -------------------------------------------------------

## Forest plot from any tidy model
tidy(m, conf.int = TRUE, exponentiate = TRUE) %>%
  filter(term != "(Intercept)") %>%
  mutate(term = fct_reorder(term, estimate)) %>%
  ggplot(aes(estimate, term)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_pointrange(aes(xmin = conf.low, xmax = conf.high)) +
  scale_x_log10() +
  labs(x = "Odds ratio (95% CI)", y = NULL)

## Predicted values across a covariate
marginaleffects::plot_predictions(m, condition = c("age", "treat")) +
  labs(x = "Age (years)", y = "Predicted outcome")

## Panel figure
(p1 | p2) / p3 +
  plot_annotation(tag_levels = "A") +
  plot_layout(guides = "collect")


# --- 6. Saving figures -------------------------------------------------------
# Specify dimensions in the units the journal asks for. Resizing in Word
# rescales the text and ruins the typography.

ggsave(here::here("output","figures","fig1.png"), p1,
       width = 180, height = 120, units = "mm", dpi = 600)

ggsave(here::here("output","figures","fig1.pdf"), p1,
       width = 180, height = 120, units = "mm", device = cairo_pdf)

# Vector (PDF/EPS/SVG) for line art, 600 dpi PNG/TIFF for anything rasterised.
# Typical journal widths: 90 mm single column, 180 mm double.


# --- 7. Numbers in text ------------------------------------------------------
# Inline R code in Quarto/R Markdown, so the manuscript text updates when the
# analysis does. This removes the single largest source of reporting errors.

fmt_est <- function(est, lo, hi, digits = 2) {
  sprintf("%.*f (95%% CI %.*f to %.*f)", digits, est, digits, lo, digits, hi)
}

fmt_p <- function(p) {
  ifelse(p < 0.001, "p < 0.001", sprintf("p = %.3f", p))
}

# In the .qmd:
#   The adjusted odds ratio was `r fmt_est(or, lo, hi)`, `r fmt_p(pval)`.


# --- 8. Reproducibility ------------------------------------------------------

sessionInfo()
# renv::snapshot()          pin package versions
# renv::restore()           rebuild the environment elsewhere

## Before you believe any result:
##   1. Session > Restart R
##   2. Run the whole pipeline top to bottom from the raw import
##   3. Confirm the numbers are unchanged
##
## Anything that only works because of what happens to be in your environment
## is not a result yet.

# Full pipeline
# source(here::here("R","00_setup.R"))
# source(here::here("R","01_data_prep.R"))
# ...
# Or use targets:: for a dependency-aware pipeline that reruns only what changed.


# --- 9. Reporting guidelines -------------------------------------------------
#   STROBE    observational studies
#   CONSORT   randomised trials
#   PRISMA    systematic reviews and meta-analyses (2020)
#   CHEERS    economic evaluations (2022)
#   TRIPOD+AI prediction models
#   RECORD    routinely collected health data
#   SQUIRE    quality improvement
#
# Fill the checklist while writing, not afterwards. It reliably surfaces
# analyses you meant to run and forgot.
