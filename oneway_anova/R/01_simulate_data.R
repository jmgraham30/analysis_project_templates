# =============================================================================
# 01_simulate_data.R
# Project:  One-way ANOVA (drug dose and locomotor activity)
# Purpose:  Document HOW the teaching data set was created.
#
# WHY SIMULATED DATA?
#   These data are computer-generated (simulated); no animals were tested. When
#   we simulate, we know the "truth" (the real group means and spread), so we
#   can check how well our statistics recover it and demonstrate power and error
#   rates honestly. Real raw data are never created by a script like this; they
#   come from your instrument or notebook and are saved untouched.
#
# ABOUT THE RANDOM SEED
#   A random sample can, by chance, look very different from the population it
#   came from. For a clear teaching example we wanted a REPRESENTATIVE sample, so
#   we used this rule, decided before looking at any p-values: use the first seed
#   (counting 1, 2, 3, ...) whose sample effect size (eta-squared) is within 0.04
#   of the TRUE effect size (about 0.24). That was seed 1. (Our first attempt, with
#   seed 2026, happened to draw a sample that nearly missed the real effect; the
#   notebook shows what such unlucky samples look like through simulation.)
#   Choosing a seed is a form of selection, so we disclose it here; a real
#   experiment gets one sample and no chance to choose.
#
# YOU DO NOT NEED TO RUN THIS FILE to do the analysis. The data already exist in
# Data/locomotor_dose.csv. Running it again re-creates an identical file because
# the seed is set (same seed = same "random" numbers).
# =============================================================================

library(tidyverse)
library(here)

set.seed(1)   # makes the random numbers reproducible (see the note above)

# ---- 1. The scenario --------------------------------------------------------
# A stimulant drug is expected to increase locomotor (movement) activity. 48
# adult mice are randomly assigned to one of FOUR independent groups of 12:
#   Vehicle (saline injection, 0 mg/kg), Low (2 mg/kg), Medium (5 mg/kg),
#   High (10 mg/kg).
# Each mouse is injected once and placed in an open-field arena for 30 minutes.
# Outcome: total distance traveled, in meters.

# ---- 2. The "true" population values used to generate the data --------------
doses <- tibble(
  dose_group = c("Vehicle", "Low", "Medium", "High"),
  dose_mgkg  = c(0, 2, 5, 10),
  true_mean  = c(40, 45, 50, 55)    # true mean distance (m): rises with dose
)
true_sd <- 10        # the same spread in every group (equal variances)
n_group <- 12        # mice per group (a "balanced" design)
# The true overall effect is about f = 0.56 (Cohen's f), a large effect. We chose
# a large effect on purpose so the example gives a clear result that can be
# followed up with post hoc comparisons. (Small effects are explored in the
# notebook's power and simulation sections.)

# ---- 3. Generate each mouse's data -----------------------------------------
locomotor <- doses |>
  slice(rep(row_number(), each = n_group)) |>         # 12 rows per group
  mutate(
    mouse_id   = sprintf("M%02d", row_number()),
    sex        = unlist(map(1:4, ~ sample(rep(c("F", "M"), length.out = n_group)))),
    distance_m = round(pmax(rnorm(n(), mean = true_mean, sd = true_sd), 0), 1)
  ) |>
  select(mouse_id, dose_group, dose_mgkg, sex, distance_m)

write_csv(locomotor, here("Data", "locomotor_dose.csv"))
message("Saved Data/locomotor_dose.csv (", nrow(locomotor), " rows)")
