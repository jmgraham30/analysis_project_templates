# =============================================================================
# 01_simulate_data.R
# Project:  Organization and display of data (messy lab data to tidy data and figures)
# Purpose:  Document HOW the two teaching data files were created, including the
#           problems that were added ON PURPOSE.
#
# WHY SIMULATED DATA?
#   These data are computer-generated (simulated); no animals were tested. We
#   first create a perfectly clean data set, so we KNOW the right answer, and
#   then deliberately add the kinds of problems that show up in real lab
#   spreadsheets (inconsistent labels, missing-value codes, typos, a duplicated
#   row, and so on). The notebook then shows how to find and fix each one in code,
#   and we can check that the cleaned data match the original truth.
#
# ABOUT THE RANDOM SEED
#   The seed (2026) was the first one tried; it was NOT chosen to give any
#   particular pattern. This project does no hypothesis tests, so there is no
#   result to "select" for.
#
# YOU DO NOT NEED TO RUN THIS FILE to do the analysis. The data already exist in
# Data/. Running it again re-creates identical files because the seed is set.
# =============================================================================

library(tidyverse)
library(here)

set.seed(2026)

# ---- 1. The scenario --------------------------------------------------------
# 40 adult mice from two genotypes (wild type, WT, and a knockout, KO) are given
# either vehicle (a control injection) or a drug daily for 3 weeks. Body weight is
# recorded at the start and then weekly (weeks 0 to 3). At the end, each mouse
# is tested for 10 minutes in an open field. Outcomes:
#   - body weight (g), 4 times
#   - distance traveled in the open field (m)
#   - time spent in the center of the open field (s)
# Mice are housed in cages of 5, all of the same genotype and sex.

# ---- 2. The clean "truth" ----------------------------------------------------
cages <- tibble(
  cage     = sprintf("C%02d", 1:8),
  genotype = rep(c("WT", "KO"), each = 4),
  sex      = rep(c("F", "F", "M", "M"), times = 2),
  birth    = as.Date("2025-11-03") + c(0, 2, 5, 1, 3, 6, 0, 4)   # one litter per cage
)

mice <- cages |>
  slice(rep(row_number(), each = 5)) |>
  mutate(mouse_id = sprintf("M%02d", row_number())) |>
  group_by(cage) |>
  mutate(
    # Odd cages: 3 vehicle + 2 drug; even cages: 2 vehicle + 3 drug (balanced overall)
    treatment = sample(if (as.integer(str_sub(first(cage), 2)) %% 2 == 1)
                         c(rep("Vehicle", 3), rep("Drug", 2))
                       else c(rep("Vehicle", 2), rep("Drug", 3)))
  ) |>
  ungroup() |>
  mutate(
    ko   = genotype == "KO",
    male = sex == "M",
    drug = treatment == "Drug",
    # Weight at week 0, then weekly growth (the drug slows weight gain)
    w0 = 20.5 + 5 * male - 1.0 * ko + rnorm(n(), 0, 1.3),
    growth = 0.45 + 0.10 * male - 0.70 * drug,
    weight_wk0 = w0,
    weight_wk1 = w0 + 1 * growth + rnorm(n(), 0, 0.5),
    weight_wk2 = w0 + 2 * growth + rnorm(n(), 0, 0.5),
    weight_wk3 = w0 + 3 * growth + rnorm(n(), 0, 0.5),
    distance_m    = 40 + 6 * ko + 12 * drug + 5 * male + rnorm(n(), 0, 10),
    center_time_s = rgamma(n(), shape = 4, rate = 4 / (45 * if_else(drug, 0.8, 1))),
    test_date = as.Date("2026-02-09") + (as.integer(str_sub(cage, 2)) - 1) %% 4
  ) |>
  mutate(across(c(weight_wk0:weight_wk3, distance_m, center_time_s), ~ round(.x, 1)))

# ---- 3. File 1: the animal records (a clean-ish table kept by the colony manager)
animals <- mice |>
  transmute(
    animal_id = mouse_id, genotype, cage,
    sex = map_chr(sex, ~ sample(if (.x == "F") c("F", "f", "Female") else c("M", "m", "Male"),
                                1, prob = c(0.6, 0.2, 0.2))),   # inconsistent labels
    birth_date = format(birth, "%Y-%m-%d")
  )

# ---- 4. File 2: the "raw" spreadsheet from the lab (messy on purpose) -------
# Everything is stored as TEXT first so that we can put in problems that a number
# column could not hold.
raw <- mice |>
  transmute(
    `Animal ID`        = mouse_id,
    `Test Date`        = format(test_date, "%Y-%m-%d"),
    Treatment          = map_chr(treatment, ~ sample(
                           if (.x == "Vehicle") c("Vehicle", "vehicle", "VEH", "Vehicle ")
                           else                 c("Drug", "drug", "DRUG", "Drug "),
                           1, prob = c(0.5, 0.2, 0.15, 0.15))),
    `Weight Wk0 (g)`   = as.character(weight_wk0),
    `Weight Wk1 (g)`   = as.character(weight_wk1),
    `Weight Wk2 (g)`   = as.character(weight_wk2),
    `Weight Wk3 (g)`   = as.character(weight_wk3),
    `Distance (m)`     = as.character(distance_m),
    `Center Time (s)`  = as.character(center_time_s),
    Notes              = NA_character_
  )

# Problems added ON PURPOSE (each one is found and fixed in the notebook):
raw$`Animal ID`[raw$`Animal ID` == "M14"] <- "m14"                      # lowercase ID
raw$`Weight Wk0 (g)`[raw$`Animal ID` == "M09"] <- paste(raw$`Weight Wk0 (g)`[raw$`Animal ID` == "M09"], "g")  # unit typed in the cell
raw$`Weight Wk2 (g)`[raw$`Animal ID` == "M05"] <- "-999"                # missing-value code
raw$`Weight Wk1 (g)`[raw$`Animal ID` == "M21"] <- "NA"                  # text "NA"
raw$`Weight Wk3 (g)`[raw$`Animal ID` == "M26"] <- NA                    # blank cell (scale error)
raw$`Weight Wk2 (g)`[raw$`Animal ID` == "M18"] <-                       # decimal point typo: 26.0 -> 260
  as.character(10 * as.numeric(raw$`Weight Wk2 (g)`[raw$`Animal ID` == "M18"]))
raw$`Distance (m)`[raw$`Animal ID` == "M12"] <- NA                      # tracking failed
raw$`Center Time (s)`[raw$`Animal ID` == "M12"] <- NA
raw$Notes[raw$`Animal ID` == "M12"] <- "Tracking software failed"
# One mouse died: no weights after week 1 and no open-field test
died <- raw$`Animal ID` == "M33"
raw$`Weight Wk2 (g)`[died] <- NA
raw$`Weight Wk3 (g)`[died] <- NA
raw$`Distance (m)`[died] <- NA
raw$`Center Time (s)`[died] <- NA
raw$Notes[died] <- "Died day 12 (not tested)"
# One row entered twice
raw <- bind_rows(raw, raw[raw$`Animal ID` == "M27", ]) |> arrange(`Animal ID`)

# ---- 5. Save ------------------------------------------------------------------
write_csv(animals, here("Data", "animals.csv"))
write_csv(raw,     here("Data", "open_field_raw.csv"), na = "")
# (The clean truth is NOT saved in the project: the notebook rebuilds it from the
#  messy files. It is saved to a scratch file only so we can compare.)
message("Saved Data/animals.csv (", nrow(animals), " rows) and ",
        "Data/open_field_raw.csv (", nrow(raw), " rows)")
