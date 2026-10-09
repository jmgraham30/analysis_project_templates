# =============================================================================
# 01_simulate_data.R
# Project:  Two-way ANOVA (housing, stress, and hippocampal BDNF)
# Purpose:  Document HOW the teaching data set was created.
#
# WHY SIMULATED DATA?
#   These data are computer-generated (simulated); no animals were tested. When
#   we simulate, we know the "truth" (the real cell means and spread), so we can
#   check how well our statistics recover it and demonstrate power and error
#   rates honestly. Real raw data are never created by a script like this; they
#   come from your instrument or notebook and are saved untouched.
#
# ABOUT THE RANDOM SEED
#   A random sample can, by chance, look very different from the population it
#   came from. For a clear teaching example we wanted a REPRESENTATIVE sample, so
#   we used this rule, decided before looking at any p-values: use the first seed
#   (counting 1, 2, 3, ...) for which the sample partial eta-squared of EACH of
#   the three effects (housing, stress, interaction) is within 0.05 of its TRUE
#   value (.31, .26, and .14). That was seed 14. Choosing a seed is a form of
#   selection, so we disclose it here; a real experiment gets one sample and no
#   chance to choose.
#
# YOU DO NOT NEED TO RUN THIS FILE to do the analysis. The data already exist in
# Data/bdnf_housing_stress.csv. Running it again re-creates an identical file
# because the seed is set (same seed = same "random" numbers).
# =============================================================================

library(tidyverse)
library(here)

set.seed(14)   # makes the random numbers reproducible (see the note above)

# ---- 1. The scenario --------------------------------------------------------
# Chronic stress lowers brain-derived neurotrophic factor (BDNF) in the
# hippocampus. Does living in an enriched environment (running wheel, tunnels,
# nesting material) protect against this? 40 adult mice are randomly assigned to
# one of FOUR groups of 10, formed by crossing two factors:
#   Housing: Standard cage or Enriched cage
#   Stress : Control (handled only) or Chronic mild stress for 3 weeks
# Outcome: hippocampal BDNF protein, in pg per mg of total protein.
# This is a 2 x 2 "factorial" design: every combination of the two factors is
# present, and each mouse is in exactly one of the four cells.

# ---- 2. The "true" population values used to generate the data --------------
cells <- tribble(
  ~housing,   ~stress,   ~true_mean,
  "Standard", "Control", 100,
  "Standard", "Stress",   70,   # stress lowers BDNF a lot in standard cages...
  "Enriched", "Control", 108,
  "Enriched", "Stress",  102    # ...but only a little in enriched cages
)
true_sd <- 15     # the same spread in every cell (equal variances)
n_cell  <- 10     # mice per cell (a "balanced" design)
# The drop caused by stress is 30 units in standard cages but only 6 in enriched
# cages. That difference in differences is the INTERACTION: the effect of one
# factor depends on the level of the other.

# ---- 3. Generate each mouse's data -----------------------------------------
bdnf <- cells |>
  slice(rep(row_number(), each = n_cell)) |>             # 10 rows per cell
  mutate(
    mouse_id  = sprintf("M%02d", row_number()),
    sex       = unlist(map(1:4, ~ sample(rep(c("F", "M"), length.out = n_cell)))),
    bdnf_pgmg = round(rnorm(n(), mean = true_mean, sd = true_sd), 1)
  ) |>
  select(mouse_id, housing, stress, sex, bdnf_pgmg)

write_csv(bdnf, here("Data", "bdnf_housing_stress.csv"))
message("Saved Data/bdnf_housing_stress.csv (", nrow(bdnf), " rows)")
