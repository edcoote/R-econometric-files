# =============================================================================
# 00_setup.R  —  Packages, import conventions, data integrity checks
# =============================================================================
# Run this first in every session. Everything downstream assumes it.
#
# PLACEHOLDER CONVENTION used across this repo:
#   d        data frame / tibble
#   y        continuous outcome
#   ybin     binary outcome (factor, "No" / "Yes")
#   yord     ordered categorical outcome (ordered factor)
#   ycount   count outcome (integer)
#   time     follow-up time
#   status   event indicator, 1 = event, 0 = censored
#   treat    exposure or treatment (factor, reference level = control)
#   x1, x2   continuous covariates
#   grp      categorical covariate (factor)
#   id       individual identifier
#   clust    clustering variable (practice, site, region)
# =============================================================================


# --- Packages ----------------------------------------------------------------
# install.packages(c(
#   "tidyverse", "broom", "janitor", "here",
#   "lmtest", "sandwich", "car", "marginaleffects", "emmeans",
#   "MASS", "nnet", "pscl", "quantreg",
#   "survival", "survminer", "ggsurvfit", "flexsurv", "cmprsk", "tidycmprsk",
#   "MatchIt", "WeightIt", "cobalt", "ivreg", "fixest", "rdrobust", "plm",
#   "meta", "metafor",
#   "heemod", "dampack", "BCEA",
#   "gtsummary", "gt", "modelsummary"
# ))

library(tidyverse)     # dplyr, ggplot2, tidyr, readr, forcats
library(broom)         # tidy() / glance() / augment()
library(janitor)       # clean_names(), tabyl()
library(here)          # project-relative paths

# NOTE: MASS::select() masks dplyr::select(). If you load MASS, either load it
# BEFORE tidyverse, or always write dplyr::select() explicitly.


# --- Session hygiene ---------------------------------------------------------
# In RStudio: Tools > Global Options > General
#   "Restore .RData into workspace at startup"  -> UNCHECK
#   "Save workspace to .RData on exit"          -> Never
# Otherwise objects silently persist between sessions and your script stops
# being the single source of truth for your results.

set.seed(1234)         # any analysis involving randomness
options(scipen = 999)  # suppress scientific notation in output


# --- Import ------------------------------------------------------------------
# Declare missing-value codes AT IMPORT. This is the single highest-value line
# in the whole pipeline: it prevents the two most common downstream failures
# (a numeric column read as character; a factor gaining phantom levels).

d <- readr::read_csv(
  here("data", "raw", "FILENAME.csv"),
  na = c("", "NA", "N/A", ".", "-", "99", "999", "-99")
)

# Other formats
# d <- readxl::read_excel(here("data","raw","FILE.xlsx"), sheet = "Data",
#                         na = c("", "NA", "99"))
# d <- haven::read_dta(here("data","raw","FILE.dta"))   # Stata, keeps labels
# d <- haven::read_sav(here("data","raw","FILE.sav"))   # SPSS

d <- janitor::clean_names(d)   # snake_case, no spaces or punctuation


# --- Integrity check ---------------------------------------------------------
# Run on every new dataset before modelling. Most later errors are visible here.

check_data <- function(d) {
  cat("rows:", nrow(d), "  cols:", ncol(d), "\n\n")

  cat("--- column types ---\n")
  print(sapply(d, function(x) class(x)[1]))

  cat("\n--- missing per column ---\n")
  miss <- colSums(is.na(d))
  print(miss[miss > 0])
  cat("complete rows:", sum(complete.cases(d)),
      sprintf("(%.1f%%)\n", 100 * mean(complete.cases(d))))

  # Numeric-looking columns stored as character = a disguised missing code
  suspect <- names(d)[sapply(d, function(x) {
    if (!is.character(x)) return(FALSE)
    v <- x[!is.na(x)]
    length(v) > 0 && mean(!is.na(suppressWarnings(as.numeric(v)))) > 0.8
  })]
  if (length(suspect))
    cat("\nCHARACTER BUT MOSTLY NUMERIC (check for stray codes):\n  ",
        paste(suspect, collapse = ", "), "\n")

  # Constant columns break contrasts in model.matrix()
  const <- names(d)[sapply(d, function(x) length(unique(na.omit(x))) < 2)]
  if (length(const))
    cat("\nCONSTANT COLUMNS:", paste(const, collapse = ", "), "\n")

  # Distinct values of every non-numeric column
  cat("\n--- categorical levels ---\n")
  for (v in names(d)[!sapply(d, is.numeric)]) {
    cat("\n", v, ":\n", sep = "")
    print(table(d[[v]], useNA = "ifany"))
  }
  invisible(NULL)
}

check_data(d)


# --- Find what blocks a numeric conversion -----------------------------------
# Run BEFORE any as.numeric() on a character column. as.numeric() converts
# whatever it can and silently turns the rest into NA — this tells you exactly
# what you would be destroying.

blockers <- function(x) unique(x[is.na(suppressWarnings(as.numeric(x))) & !is.na(x)])

blockers(d$x1)
# "NA", "", "."           -> safe to coerce
# "<90", "140/90", "1,240" -> real data; parse properly, do not coerce
# factor                   -> use as.numeric(as.character(x)), never as.numeric(x)


# --- Project structure -------------------------------------------------------
# Use an RStudio Project (.Rproj) and here() rather than setwd(). Paths then
# stay relative and the work is portable.
#
#   project/
#     project.Rproj
#     data/raw/         read-only, never edited by script
#     data/derived/     analysis-ready, written by 01_data_prep.R
#     R/                these scripts
#     output/tables/
#     output/figures/
#
# renv::init()      # pin package versions per project
# renv::snapshot()  # after adding packages
# sessionInfo()     # record versions in any output you share
