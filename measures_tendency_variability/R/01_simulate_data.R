# =============================================================================
# 01_simulate_data.R
# Project:  Measures of central tendency and variability (neuron firing rates)
# Purpose:  Document HOW the teaching data set was created.
#
# WHY SIMULATED DATA?
#   These data are computer-generated (simulated); no animals were tested. When
#   we simulate, we know the "truth" (the real distribution the numbers came from),
#   so we can check how well our summary statistics recover it and demonstrate
#   ideas such as sampling variability honestly. Real raw data are never created
#   by a script like this; they come from your instrument and are saved untouched.
#
# ABOUT THE RANDOM SEED
#   The seed (2026) was the first one tried; it was NOT chosen to give any
#   particular pattern. This project runs no hypothesis tests, so there is no
#   result to "select" for.
#
# YOU DO NOT NEED TO RUN THIS FILE to do the analysis. The data already exist in
# Data/neuron_firing.csv. Running it again re-creates an identical file because
# the seed is set (same seed = same "random" numbers).
# =============================================================================

library(tidyverse)
library(here)

set.seed(2026)

# ---- 1. The scenario --------------------------------------------------------
# Extracellular recordings from 120 single neurons: 40 in each of three brain
# regions (somatosensory cortex, hippocampus, and striatum), from 12 mice. For each neuron we
# have:
#   firing_rate_hz   average firing rate (spikes per second); a RATIO variable
#   spike_width_ms   width of the average spike waveform (milliseconds); a RATIO
#                    variable that is much more symmetric than firing rate
#   cell_type        Regular spiking (usually excitatory) or Fast spiking
#                    (usually inhibitory): a NOMINAL variable
#   tuning_rating    how clearly the neuron responded to the stimulus, rated by
#                    an observer from 1 (none) to 5 (very clear): an ORDINAL variable

# ---- 2. The "true" population values used to generate the data --------------
regions <- tribble(
  ~region,         ~median_rs_hz,   # typical firing rate of a regular-spiking neuron
  "Cortex",         3.0,
  "Hippocampus",    1.8,
  "Striatum",       4.5
)
n_per_region <- 40
p_fast       <- 0.20     # 20% of neurons are fast spiking
sigma_log    <- 0.60     # spread of log(firing rate): firing rates are SKEWED
fast_factor  <- 3.5      # fast-spiking neurons fire about 3.5 times faster

# ---- 3. Generate each neuron's data -----------------------------------------
neurons <- regions |>
  slice(rep(row_number(), each = n_per_region)) |>
  mutate(
    neuron_id = sprintf("N%03d", row_number()),
    cell_type = if_else(runif(n()) < p_fast, "Fast spiking", "Regular spiking"),
    # Firing rates are always positive and right-skewed: a LOG-NORMAL distribution
    # (the logarithm of the rate is normal) is a common, realistic model.
    firing_rate_hz = rlnorm(n(), meanlog = log(median_rs_hz) + log(fast_factor) * (cell_type == "Fast spiking"),
                            sdlog = sigma_log),
    # Spike width is roughly symmetric (normal); fast-spiking neurons have narrower spikes
    spike_width_ms = rnorm(n(), mean = if_else(cell_type == "Fast spiking", 0.22, 0.42), sd = 0.06),
    # An ordinal rating from 1 to 5; clearer responses are more likely in fast-spiking cells
    tuning_rating  = map_int(cell_type, ~ sample(1:5, 1, prob = if (.x == "Fast spiking") c(.05, .10, .20, .30, .35)
                                                  else c(.20, .30, .25, .15, .10)))
  ) |>
  mutate(across(c(firing_rate_hz, spike_width_ms), ~ round(.x, 2))) |>
  mutate(spike_width_ms = pmax(spike_width_ms, 0.10)) |>
  # Neurons were recorded from 12 mice (10 neurons per mouse, 4 mice per region).
  # In this simulation mice do not differ from each other; in real data, neurons
  # from the same mouse tend to be more alike than neurons from different mice.
  mutate(mouse_id = sprintf("M%02d", ceiling(row_number() / 10))) |>
  select(neuron_id, mouse_id, region, cell_type, firing_rate_hz, spike_width_ms, tuning_rating)

write_csv(neurons, here("Data", "neuron_firing.csv"))
message("Saved Data/neuron_firing.csv (", nrow(neurons), " rows)")
