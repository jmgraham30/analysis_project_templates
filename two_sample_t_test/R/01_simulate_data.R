# =============================================================================
# 01_simulate_data.R
# Project:  Two-sample t-test (elevated plus maze: stress vs. control)
# Purpose:  Document HOW the teaching data set was created.
#
# WHY SIMULATED DATA?
#   These data are computer-generated (simulated); no animals were tested. When
#   we simulate, we know the "truth" (the real group means and spreads), so we
#   can check how well our statistics recover it and demonstrate power and error
#   rates honestly. Real raw data are never created by a script like this; they
#   come from your instrument or notebook and are saved untouched.
#
# ABOUT THE RANDOM SEED
#   The seed below (2026) was the FIRST one tried. We did not search for a seed
#   that gives a "nice" result, so the example data are simply one realistic
#   random sample. (Choosing a seed because of how its results look would be a
#   subtle form of cherry-picking; see the notebook for more on this idea.)
#
# YOU DO NOT NEED TO RUN THIS FILE to do the analysis. The data already exist in
# Data/elevated_plus_maze.csv. Running it again re-creates an identical file
# because the seed is set (same seed = same "random" numbers).
# =============================================================================

library(tidyverse)
library(here)

set.seed(2026)   # makes the random numbers reproducible

# ---- 1. The scenario --------------------------------------------------------
# The elevated plus maze (EPM) is a plus-shaped maze raised off the floor, with
# two OPEN arms and two CLOSED (walled) arms. Mice naturally prefer the safe,
# closed arms; the more anxious a mouse is, the LESS time it spends in the open
# arms. We compare two independent groups of adult mice:
#   - Control : handled normally
#   - Stress  : exposed to two weeks of chronic mild stress
# Outcome: percent of the 5-minute test spent in the open arms.

# ---- 2. The "true" population values used to generate the data --------------
true_mean <- c(Control = 34, Stress = 25)   # % time in open arms
true_sd   <- c(Control =  9, Stress = 11)   # note: slightly unequal spreads
n_group   <- c(Control = 15, Stress = 14)   # mice tested per group (unequal n)
# The true standardized difference is about d = 0.9 (a large effect).

# ---- 3. Generate each mouse's data -----------------------------------------
make_group <- function(grp, prefix) {
  n <- n_group[[grp]]
  tibble(
    mouse_id     = sprintf("%s%02d", prefix, seq_len(n)),
    group        = grp,
    sex          = sample(rep(c("F", "M"), length.out = n)),   # shuffled, balanced
    open_arm_pct = rnorm(n, mean = true_mean[[grp]], sd = true_sd[[grp]]),
    # Total arm entries: a rough measure of general activity
    total_entries = rpois(n, lambda = 24)
  )
}

epm <- bind_rows(
  make_group("Control", "C"),
  make_group("Stress",  "S")
) |>
  mutate(
    # A percentage must stay between 0 and 100; record to 0.1 percentage point
    open_arm_pct = round(pmin(pmax(open_arm_pct, 0), 100), 1)
  )

# ---- 4. Add a realistic data problem: one missing value ---------------------
# Real data are rarely perfect. Suppose mouse S09 jumped off the maze and the
# test could not be completed. We record this as a missing value (NA), NOT as 0.
epm <- epm |>
  mutate(across(c(open_arm_pct, total_entries),
                ~ replace(.x, mouse_id == "S09", NA)))

write_csv(epm, here("Data", "elevated_plus_maze.csv"))
message("Saved Data/elevated_plus_maze.csv (", nrow(epm), " rows)")
