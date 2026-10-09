# =============================================================================
# 01_simulate_data.R
# Project:  One-sample t-test (novel object recognition)
# Purpose:  Document HOW the teaching data sets were created.
#
# WHY SIMULATED DATA?
#   These data are computer-generated (simulated); no animals were tested. When
#   we simulate, we know the "truth" (the real population mean and spread), so
#   we can check how well our statistics recover it and demonstrate power and
#   error rates honestly. Real raw data are never created by a script like this;
#   they come from your instrument or notebook and are saved untouched.
#
# TWO COHORTS, ONE POPULATION
#   Both cohorts below are drawn from the SAME population with the SAME true
#   effect. They differ only in the random draw (the "seed"). This lets the
#   lab notebook show how much results can differ from one sample to the next:
#     - Cohort 1 (seed 5):    Data/novel_object_recognition.csv
#     - Cohort 2 (seed 2026): Data/novel_object_recognition_cohort2.csv
#   Seed 5 was chosen, after trying a few seeds, because it gives a sample whose
#   effect resembles the true one. Seed 2026 was the first seed we tried, and it
#   gives a sample that misses the real effect. We keep BOTH on purpose: choosing
#   only the "nicest" sample would be misleading, and showing both illustrates
#   sampling variability and Type II errors.
#
# YOU DO NOT NEED TO RUN THIS FILE to do the analysis. The data already exist in
# the Data folder. Running it again re-creates identical files because random
# seeds are set (same seed = same "random" numbers).
# =============================================================================

library(tidyverse)
library(here)

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

# ---- 3. A function that simulates one cohort of mice ------------------------
simulate_cohort <- function(seed, id_prefix) {
  set.seed(seed)   # makes the random numbers reproducible

  mice <- tibble(
    mouse_id  = sprintf("%s%02d", id_prefix, 1:n_mice),
    sex       = rep(c("F", "M"), each = n_mice / 2),
    age_weeks = sample(10:14, n_mice, replace = TRUE),
    # Total exploration time (s): positively skewed, typically ~40 s
    total_s   = rlnorm(n_mice, meanlog = log(40), sdlog = 0.35),
    di_true   = rnorm(n_mice, mean = true_mean_DI, sd = true_sd_DI)
  )

  # Two mice barely explored the objects. We deliberately make them low
  # explorers so the analysis can demonstrate a pre-specified exclusion rule.
  low <- mice$mouse_id %in% sprintf(c("%s07", "%s19"), id_prefix)
  mice$total_s[low] <- c(11.4, 14.9)

  # Keep the index inside its possible range (-1 to 1) and convert to raw times.
  mice |>
    mutate(
      di_true         = pmin(pmax(di_true, -0.9), 0.9),
      time_novel_s    = round(total_s * (1 + di_true) / 2, 1),
      time_familiar_s = round(total_s * (1 - di_true) / 2, 1)
    ) |>
    # Keep only what a real experimenter would have recorded. Notice we save the
    # RAW observations (the two stopwatch times), NOT the discrimination index.
    # The index is calculated later, in the analysis script.
    select(mouse_id, sex, age_weeks, time_novel_s, time_familiar_s)
}

# ---- 4. Create and save the two cohorts -------------------------------------
cohort1 <- simulate_cohort(seed = 5,    id_prefix = "M")
cohort2 <- simulate_cohort(seed = 2026, id_prefix = "R")   # "R" = second run

write_csv(cohort1, here("Data", "novel_object_recognition.csv"))
write_csv(cohort2, here("Data", "novel_object_recognition_cohort2.csv"))
message("Saved Data/novel_object_recognition.csv (", nrow(cohort1), " rows)")
message("Saved Data/novel_object_recognition_cohort2.csv (", nrow(cohort2), " rows)")
