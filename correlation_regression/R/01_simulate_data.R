# =============================================================================
# 01_simulate_data.R
# Project:  Correlation and regression (running and adult neurogenesis)
# Purpose:  Document HOW the teaching data set was created.
#
# WHY SIMULATED DATA?
#   These data are computer-generated (simulated); no animals were tested. When
#   we simulate, we know the "truth" (the real relationship), so we can check
#   how well our statistics recover it and demonstrate power and error rates
#   honestly. Real raw data come from your instrument or notebook and are saved
#   untouched; they are never created by a script like this one.
#
# ABOUT THE RANDOM SEED
#   A random sample can, by chance, look quite different from the population it
#   came from. For a clear teaching example we wanted a REPRESENTATIVE sample, so
#   we used this rule, decided before looking at any p-values from the analysis:
#   use the first seed (counting 1, 2, 3, ...) for which
#     (1) the sample correlation between running distance and new neurons is
#         within 0.03 of the TRUE population correlation (r = 0.62, found by
#         simulating one million mice), AND
#     (2) the random "noise" added to each mouse is not noticeably non-normal
#         (Shapiro-Wilk p > .20 on the true errors).
#   Our first attempt used only rule (1) and picked seed 1, whose noise happened
#   to be strongly right-skewed (p = .003), so we added rule (2) and re-ran the
#   search; the first seed meeting both rules is seed 4. The notebook shows what
#   assumption violations look like using separate demonstrations. Choosing a
#   seed is a form of selection, so we disclose it here; a real experiment gets
#   one sample and no chance to choose.
#
# YOU DO NOT NEED TO RUN THIS FILE to do the analysis. The data already exist in
# Data/neurogenesis_running.csv. Running it again re-creates an identical file
# because the seed is set (same seed = same "random" numbers).
# =============================================================================

library(tidyverse)
library(here)

set.seed(4)   # makes the random numbers reproducible (see the note above)

# ---- 1. The scenario --------------------------------------------------------
# Voluntary wheel running is known to increase adult neurogenesis in the
# hippocampal dentate gyrus. 60 adult mice (8 to 20 weeks old) lived with a
# running wheel for two weeks. For each mouse we recorded:
#   wheel_km  - average distance run per night (km)      [predictor, X]
#   dcx_cells - doublecortin-positive (DCX+) cells counted in the dentate gyrus,
#               a marker of newly born neurons             [outcome, Y]
#   age_weeks - age of the mouse; neurogenesis declines with age (a covariate)
#   sex       - F or M

# ---- 2. The "true" relationship used to generate the data -------------------
#   dcx_cells = 110 + 12 * wheel_km - 3 * (age_weeks - 14) + random noise (SD 25)
# So, in the population, each extra km/night adds about 12 DCX+ cells and each
# extra week of age removes about 3. The population correlation between running
# and cell count is r = 0.62 (R-squared about 0.38).
n <- 60

# ---- 3. Generate each mouse's data -----------------------------------------
neurogenesis <- tibble(
  mouse_id  = sprintf("M%02d", seq_len(n)),
  sex       = sample(rep(c("F", "M"), length.out = n)),        # shuffled 30 F / 30 M
  age_weeks = round(runif(n, 8, 20)),
  wheel_km  = round(pmax(rnorm(n, mean = 4, sd = 1.8), 0.2), 1) # distances can't be negative
) |>
  mutate(
    dcx_cells = round(pmax(rnorm(n, mean = 110 + 12 * wheel_km - 3 * (age_weeks - 14),
                                 sd = 25), 0))                 # counts can't be negative
  )

write_csv(neurogenesis, here("Data", "neurogenesis_running.csv"))
message("Saved Data/neurogenesis_running.csv (", nrow(neurogenesis), " rows)")
