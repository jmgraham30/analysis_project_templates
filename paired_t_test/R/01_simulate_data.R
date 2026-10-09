# =============================================================================
# 01_simulate_data.R
# Project:  Paired-samples t-test (open-field crossover: drug vs. vehicle)
# Purpose:  Document HOW the two teaching data sets were created.
#
# WHY SIMULATED DATA?
#   These data are computer-generated (simulated); no animals were tested. When
#   we simulate, we know the "truth" (the real effects and spreads), so we can
#   check how well our statistics recover it and demonstrate error rates
#   honestly. Real raw data come from your instrument or notebook and are saved
#   untouched; they are never created by a script like this.
#
# ABOUT THE RANDOM SEED (please read; it is a little different from the other
# projects)
#   The seed below (5) was chosen by a rule fixed BEFORE looking at any
#   results: "use the first seed, trying 1, 2, 3, ..., for which (a) the mean of the
#   analyzed paired differences (open field) and the mean neuron-level change
#   (slices) are each within 10% of the true effect built into the simulation,
#   and (b) a Shapiro-Wilk test does not flag the analyzed
#   differences as non-normal (p > .20)." This keeps the examples typical of
#   the true situation, with no unusual draw. Even so, a random sample will
#   not match the truth exactly, and the notebook shows how much that matters.
#
# YOU DO NOT NEED TO RUN THIS FILE to do the analysis. The data already exist
# in Data/. Running it again re-creates identical files because the seed is set.
# =============================================================================

library(tidyverse)
library(here)

set.seed(5)   # makes the random numbers reproducible

# =============================================================================
# DATA SET 1: open-field crossover (the main analysis)
# =============================================================================
# Scenario: 14 adult mice are each tested in an open-field arena twice, one week
# apart. Before one session a mouse gets a low dose of a stimulant (Drug); before
# the other it gets a saline injection (Vehicle). Half of the mice get the drug
# first and half get the vehicle first (a counterbalanced "crossover" design), so
# the order of the sessions is not confounded with the treatment.
# Outcome: distance traveled (meters) in a 10-minute session.

n_mice          <- 14
true_effect     <- 5      # true mean increase in distance caused by the drug (m)
sd_mouse        <- 6      # SD of mouse-to-mouse differences in baseline activity
sd_drug_resp    <- 4      # SD of mouse-to-mouse differences in drug response
sd_session      <- 3      # session-to-session noise within a mouse (m)
order_effect    <- -2     # session 2 is 2 m lower than session 1 (habituation)

mice <- tibble(
  mouse_id  = sprintf("M%02d", seq_len(n_mice)),
  sex       = sample(rep(c("F", "M"), length.out = n_mice)),
  sequence  = sample(rep(c("Drug first", "Vehicle first"), length.out = n_mice)),
  baseline  = rnorm(n_mice, mean = 30, sd = sd_mouse),    # this mouse's usual level
  drug_resp = rnorm(n_mice, mean = true_effect, sd = sd_drug_resp)
)

open_field <- mice |>
  # One row per mouse per session (tidy "long" format)
  mutate(session_1 = if_else(sequence == "Drug first", "Drug", "Vehicle"),
         session_2 = if_else(sequence == "Drug first", "Vehicle", "Drug")) |>
  pivot_longer(c(session_1, session_2), names_to = "session",
               names_prefix = "session_", values_to = "treatment") |>
  mutate(
    session  = as.integer(session),
    distance = baseline +
      if_else(treatment == "Drug", drug_resp, 0) +
      if_else(session == 2, order_effect, 0) +
      rnorm(n(), 0, sd_session),
    distance_m = round(pmax(distance, 1), 1)
  ) |>
  select(mouse_id, sex, sequence, session, treatment, distance_m) |>
  arrange(mouse_id, session)

# A realistic data problem: for ONE mouse the video tracking failed in one
# session, so that distance is missing (NA, never 0). A paired analysis needs
# BOTH sessions, so this mouse will be left out of the paired test.
lost_row <- open_field |> filter(mouse_id == sample(mouse_id, 1), session == 2) |>
  slice(1)
open_field <- open_field |>
  mutate(distance_m = replace(distance_m,
                              mouse_id == lost_row$mouse_id & session == 2, NA))

write_csv(open_field, here("Data", "open_field_crossover.csv"))
message("Saved Data/open_field_crossover.csv (", nrow(open_field), " rows)")

# =============================================================================
# DATA SET 2: neurons recorded in brain slices (the "many measurements per
# animal" lesson)
# =============================================================================
# Scenario: 8 mice provide brain slices. In each slice several neurons are
# recorded (3 to 6 per mouse). Each neuron's firing rate (Hz) is measured at
# baseline and again after bath application of the drug, so the data are paired
# at the level of the NEURON. But neurons from the same mouse are not
# independent of each other: the mouse is the real unit of analysis.

neurons_per_mouse <- c(4, 3, 5, 4, 6, 3, 4, 5)       # 34 neurons in all
n_slice_mice      <- length(neurons_per_mouse)
true_change       <- 1.3   # true mean increase in firing rate after the drug (Hz)
sd_mouse_base     <- 1.0   # mouse-to-mouse differences in baseline firing
sd_mouse_change   <- 1.2   # mouse-to-mouse differences in the drug's effect
sd_neuron_base    <- 1.5   # neuron-to-neuron differences in baseline firing
sd_neuron_change  <- 1.5   # neuron-to-neuron noise in the change

slice_mice <- tibble(
  mouse_id    = sprintf("S%02d", seq_len(n_slice_mice)),
  mouse_base  = rnorm(n_slice_mice, 0, sd_mouse_base),
  mouse_shift = rnorm(n_slice_mice, 0, sd_mouse_change)   # this mouse's own shift
)

slice_recordings <- slice_mice |>
  mutate(n_neurons = neurons_per_mouse) |>
  uncount(n_neurons) |>
  group_by(mouse_id) |>
  mutate(neuron_id = sprintf("%s-N%d", mouse_id, row_number())) |>
  ungroup() |>
  mutate(
    Baseline = pmax(5 + mouse_base + rnorm(n(), 0, sd_neuron_base), 0.5),
    Drug     = pmax(Baseline + true_change + mouse_shift +
                      rnorm(n(), 0, sd_neuron_change), 0.1)
  ) |>
  pivot_longer(c(Baseline, Drug), names_to = "condition", values_to = "firing_hz") |>
  mutate(firing_hz = round(firing_hz, 1)) |>
  select(mouse_id, neuron_id, condition, firing_hz) |>
  arrange(mouse_id, neuron_id, factor(condition, levels = c("Baseline", "Drug")))

write_csv(slice_recordings, here("Data", "slice_recordings.csv"))
message("Saved Data/slice_recordings.csv (", nrow(slice_recordings), " rows)")
