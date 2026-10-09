# =============================================================================
# 01_simulate_data.R
# Project:  One-sample t-test (novel object recognition)
# Purpose:  Document HOW the teaching data set was created.
#
# WHY SIMULATED DATA?
#   These data are computer-generated (simulated); no animals were tested. When
#   we simulate, we know the "truth" (the real population mean and spread), so
#   we can check how well our statistics recover it and demonstrate power and
#   error rates honestly. Real raw data are never created by a script like this;
#   they come from your instrument or notebook and are saved untouched.
#
# YOU DO NOT NEED TO RUN THIS FILE to do the analysis. The data already exist in
# Data/novel_object_recognition.csv. Running it again re-creates an identical
# file because a random seed is set below (same seed = same "random" numbers).
# =============================================================================

library(tidyverse)
library(here)

# The seed makes the random numbers reproducible. (Seed 5 was chosen, after
# trying a few, because it gives a sample whose effect size resembles the true
# one. With a different seed you would get a different sample, and in a small
# study that sample could easily look more or less impressive. That is exactly
# the sampling variability the Type I / Type II error demonstration in the lab
# notebook explores.)
set.seed(5)

# ---- 1. The scenario --------------------------------------------------------
# 26 adult mice complete a novel object recognition (NOR) test. During the test
# phase each mouse explores one FAMILIAR and one NOVEL object. A mouse that
# remembers the familiar object should spend MORE time with the novel one.
#
# Discrimination index (DI) = (time_novel - time_familiar) /
#                             (time_novel + time_familiar)
#   DI =  0  -> no preference (no memory, or chance)
#   DI >  0  -> prefers the novel object (memory for the familiar object)
#   DI <  0  -> prefers the familiar object

n_mice <- 26

# ---- 2. The "true" population values used to generate the data --------------
true_mean_DI <- 0.15   # a modest real preference for the novel object
true_sd_DI   <- 0.28   # mouse-to-mouse variability in the index
# (true_mean / true_sd gives a "true" Cohen's d of about 0.54: a medium effect)

# ---- 3. Generate each mouse's data -----------------------------------------
mice <- tibble(
  mouse_id  = sprintf("M%02d", 1:n_mice),
  sex       = rep(c("F", "M"), each = n_mice / 2),
  age_weeks = sample(10:14, n_mice, replace = TRUE),
  # Total exploration time (s): positively skewed, typically ~40 s
  total_s   = rlnorm(n_mice, meanlog = log(40), sdlog = 0.35),
  di_true   = rnorm(n_mice, mean = true_mean_DI, sd = true_sd_DI)
)

# Two mice barely explored the objects. We deliberately make them low explorers
# so the analysis can demonstrate a pre-specified exclusion rule.
mice$total_s[mice$mouse_id %in% c("M07", "M19")] <- c(11.4, 14.9)

# Keep the index inside its possible range (-1 to 1) and convert to raw times.
mice <- mice |>
  mutate(
    di_true         = pmin(pmax(di_true, -0.9), 0.9),
    time_novel_s    = round(total_s * (1 + di_true) / 2, 1),
    time_familiar_s = round(total_s * (1 - di_true) / 2, 1)
  )

# ---- 4. Keep only what a real experimenter would have recorded -------------
# Notice we save the RAW observations (the two stopwatch times), NOT the
# discrimination index. The index is calculated later in the analysis script.
nor_data <- mice |>
  select(mouse_id, sex, age_weeks, time_novel_s, time_familiar_s)

write_csv(nor_data, here("Data", "novel_object_recognition.csv"))
message("Saved Data/novel_object_recognition.csv (", nrow(nor_data), " rows)")
