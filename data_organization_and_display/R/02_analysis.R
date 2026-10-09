# =============================================================================
# 02_analysis.R
# Project:  Organization and display of data (messy lab data to tidy data and figures)
# Purpose:  Complete, step-by-step script that takes two MESSY lab files to a
#           clean, documented, TIDY data set and then displays it. Run it from top
#           to bottom (or line by line with Ctrl/Cmd + Enter) to reproduce every
#           table and figure in the lab notebook.
#
# WHAT THIS PROJECT IS (AND IS NOT)
#   This project is about getting data READY and LOOKING AT THEM HONESTLY. It
#   uses only descriptive statistics and graphs; no hypothesis tests are run.
#   Every later project in this course starts with the steps you practice here.
#
# BEFORE YOU RUN: open data_organization_and_display.Rproj so the working
#                 directory is the project folder.
# INPUT : Data/open_field_raw.csv, Data/animals.csv   (raw, NEVER edited by hand)
# OUTPUT: cleaned data in output/data/ and figures in output/figures/
#         (re-created each time you run this)
# =============================================================================


# ---- 1. Set up ---------------------------------------------------------------
source("R/00_setup.R")   # loads packages, plot theme, and helper functions

dir.create(here("output", "data"),    recursive = TRUE, showWarnings = FALSE)
dir.create(here("output", "figures"), recursive = TRUE, showWarnings = FALSE)


# ---- 2. Cleaning rules decided BEFORE touching the data ------------------------
# Writing the rules first (and sticking to them) keeps cleaning honest: you never
# decide what to delete because of what it does to a result.
#  (1) Never edit the raw files. Every change is made in code, and logged.
#  (2) Missing-value codes (-999, "NA", blank) all become R's NA. Missing is NOT 0.
#  (3) Labels are standardized (one spelling per category).
#  (4) A weight outside 10 to 50 g is implausible. If dividing it by 10 gives a
#      value that agrees with the same mouse's other weights (within 3 g), treat
#      it as a decimal-point typo and correct it; otherwise set it to NA. In a
#      real study, check the original lab notebook first.
#  (5) Exact duplicate rows are removed (keep one).
#  (6) Animals that did not complete a measurement stay in the data with NA for
#      that measurement. Nobody is dropped, and the reason is recorded.
weight_min <- 10    # grams
weight_max <- 50    # grams


# ---- 3. Import the raw behavior file and LOOK at it ---------------------------
raw <- read_csv(here("Data", "open_field_raw.csv"), show_col_types = FALSE)

glimpse(raw)    # column names and types: anything strange?
# Notice: the weight column for week 0 was read as TEXT (<chr>), not numbers.
# Something in it is not a number. Also notice the awkward column names.

problems(raw)   # (empty here: readr guessed a type for every column)
summary(raw)    # ranges: the minimum of -999 and a maximum of 260 are red flags
raw |> filter(!str_detect(`Weight Wk0 (g)`, "^[0-9.]+$"))   # which cells are not numbers?

# Look at the categories: every distinct spelling is a different "group" to R.
count(raw, Treatment)
count(raw, `Animal ID`) |> filter(n > 1)    # duplicated IDs?
raw |> filter(`Animal ID` %in% c("M27")) |> distinct()      # are the duplicates identical?


# ---- 4. Clean the column names -----------------------------------------------
# Good names: lowercase, no spaces or symbols, and units in the name when helpful.
# We do it with a rule (lowercase; everything that is not a letter or digit becomes
# an underscore), so it works the same for every column.
clean_names <- function(x) {
  x |>
    str_to_lower() |>
    str_replace_all("[^a-z0-9]+", "_") |>
    str_remove("^_|_$")
}
dat <- raw |> rename_with(clean_names)
names(dat)
# We want shorter, clearer names for a few columns:
dat <- dat |>
  rename(mouse_id = animal_id,
         weight_g_wk0 = weight_wk0_g, weight_g_wk1 = weight_wk1_g,
         weight_g_wk2 = weight_wk2_g, weight_g_wk3 = weight_wk3_g)


# ---- 5. Fix the problems, one at a time, and log each one --------------------
# We keep a "cleaning log": what we found, how many rows it touched, what we did.
cleaning_log <- tibble(problem = character(), n_affected = integer(), action = character())
log_step <- function(problem, n, action) {
  cleaning_log <<- bind_rows(cleaning_log,
                             tibble(problem = problem, n_affected = as.integer(n), action = action))
}

# 5a. IDs: whitespace and capitalization
n_id <- sum(dat$mouse_id != str_to_upper(str_trim(dat$mouse_id)))
dat <- dat |> mutate(mouse_id = str_to_upper(str_trim(mouse_id)))
log_step("ID in lowercase or with spaces", n_id, "Converted to uppercase, trimmed spaces")

# 5b. Treatment labels: several spellings of the same two categories
n_trt <- sum(!dat$treatment %in% c("Vehicle", "Drug"))
dat <- dat |>
  mutate(
    treatment = case_when(
      str_detect(str_to_lower(treatment), "^veh")  ~ "Vehicle",
      str_detect(str_to_lower(treatment), "^drug") ~ "Drug",
      .default = NA_character_
    )
  )
count(dat, treatment)                        # exactly two categories now
stopifnot(!anyNA(dat$treatment))             # an unmatched label would show up as NA
log_step("Inconsistent treatment labels (vehicle, VEH, 'Drug ', ...)", n_trt,
         "Mapped every spelling to 'Vehicle' or 'Drug'")

# 5c. Weights: text in a number column, a missing-value code, and an impossible value
# parse_number() pulls the number out of text such as "20.7 g".
n_text <- sum(!is.na(dat$weight_g_wk0) & !str_detect(dat$weight_g_wk0, "^[0-9.]+$"))
dat <- dat |>
  mutate(across(starts_with("weight_g"), ~ parse_number(as.character(.x))))   # "20.7 g" -> 20.7
log_step("Unit typed inside a weight cell ('20.7 g')", n_text, "Extracted the number")

n_code <- sum(dat |> select(starts_with("weight_g"), distance_m, center_time_s) == -999, na.rm = TRUE)
dat <- dat |>
  mutate(across(c(starts_with("weight_g"), distance_m, center_time_s),
                ~ na_if(.x, -999)))                        # -999 -> NA
log_step("Missing-value code -999", n_code, "Recoded as NA")

# Implausible weights: flag them, look at them, apply rule (4)
implausible <- dat |>
  pivot_longer(starts_with("weight_g"), names_to = "week", values_to = "weight_g") |>
  filter(!is.na(weight_g), weight_g < weight_min | weight_g > weight_max)
implausible    # one value: mouse M18, week 2

# Compare with that mouse's other weights
dat |> filter(mouse_id == "M18") |> select(mouse_id, starts_with("weight_g"))

# 260 / 10 = 26.0, which agrees with 24, 24, 26.3 -> a decimal-point typo.
dat <- dat |>
  mutate(
    weight_g_wk2 = if_else(mouse_id == "M18" & weight_g_wk2 == 260, 26.0, weight_g_wk2)
  )
log_step("Implausible weight (260 g, mouse M18, week 2)", nrow(implausible),
         "Decimal-point typo (26.0 g); corrected (rule 4) and flagged in the notebook")

# Make sure nothing implausible is left
dat |>
  pivot_longer(starts_with("weight_g"), values_to = "weight_g") |>
  summarize(min = min(weight_g, na.rm = TRUE), max = max(weight_g, na.rm = TRUE))

# 5d. Exact duplicate rows
n_dup <- nrow(dat) - nrow(distinct(dat))
dat <- distinct(dat)
log_step("Duplicated row (mouse M27 entered twice)", n_dup, "Removed the extra copy")
stopifnot(!anyDuplicated(dat$mouse_id))      # now one row per mouse

# 5e. Dates: readr already recognized test_date as a real Date (see glimpse above)
# because it is written in the international standard YYYY-MM-DD. That is one
# reason to ALWAYS enter dates this way: "03/04/2026" could mean March 4 or April 3.
class(dat$test_date)
range(dat$test_date)

# 5f. Remaining missing values: count them, find out WHY, keep a record
missing_by_var <- dat |>
  summarize(across(everything(), ~ sum(is.na(.x)))) |>
  pivot_longer(everything(), names_to = "variable", values_to = "n_missing") |>
  filter(n_missing > 0)
missing_by_var
dat |> filter(!is.na(notes)) |> select(mouse_id, notes)    # the notes explain some
dat |> filter(if_any(c(starts_with("weight_g"), distance_m, center_time_s), is.na)) |>
  select(mouse_id, treatment, starts_with("weight_g"), distance_m, notes)
# M26's week-3 weight is missing with no note. Ask the person who ran the animal
# (a scale error?). We cannot invent it. We keep the mouse and leave the NA.
log_step("Missing measurements (died, tracking failed, scale error, -999)",
         nrow(filter(dat, if_any(c(starts_with("weight_g"), distance_m, center_time_s), is.na))),
         "Kept the animals; NA for the missing values (never 0); reasons recorded")


# ---- 6. The animal records: another messy table -------------------------------
animals_raw <- read_csv(here("Data", "animals.csv"), show_col_types = FALSE)
count(animals_raw, sex)    # F, f, Female, M, m, Male

animals <- animals_raw |>
  rename(mouse_id = animal_id) |>
  mutate(
    mouse_id   = str_to_upper(str_trim(mouse_id)),
    sex        = case_when(str_to_upper(str_sub(sex, 1, 1)) == "F" ~ "Female",
                           str_to_upper(str_sub(sex, 1, 1)) == "M" ~ "Male"),
    genotype   = factor(genotype, levels = c("WT", "KO"))   # birth_date is already a Date
  )
count(animals, sex)
stopifnot(!anyNA(animals$sex), !anyDuplicated(animals$mouse_id))
log_step("Inconsistent sex labels in the animal records (F, f, Female, ...)",
         sum(!animals_raw$sex %in% c("Female", "Male")),
         "Mapped to 'Female' or 'Male'")


# ---- 7. Combine the two tables with a JOIN ------------------------------------
# Both tables share the key mouse_id. Check that the keys match BEFORE joining.
anti_join(dat, animals, by = "mouse_id")      # in the behavior file, not in the records
anti_join(animals, dat, by = "mouse_id")      # in the records, not in the behavior file
# (Both are empty: every mouse is in both files. If M14 had stayed 'm14', it would
#  have shown up here, which is exactly why we checked.)

mice <- dat |>
  left_join(animals, by = "mouse_id") |>      # keep every row of dat, add the animal info
  mutate(
    treatment = factor(treatment, levels = c("Vehicle", "Drug")),
    sex       = factor(sex, levels = c("Female", "Male")),
    age_days  = as.numeric(test_date - birth_date)
  ) |>
  select(mouse_id, genotype, sex, treatment, cage, birth_date, test_date, age_days,
         starts_with("weight_g"), distance_m, center_time_s, notes)

# ---- 8. Validation checks: does the clean data look the way it should? ---------
stopifnot(
  nrow(mice) == 40,                                             # 40 mice
  !anyDuplicated(mice$mouse_id),                                # one row per mouse
  all(count(mice, genotype, treatment)$n == 10),                # balanced design
  all(mice$age_days > 60, mice$age_days < 200),                 # plausible ages
  all(between(unlist(select(mice, starts_with("weight_g"))), weight_min, weight_max), na.rm = TRUE)
)
count(mice, genotype, treatment)
count(mice, cage, genotype, sex)   # cages hold one genotype and one sex each

glimpse(mice)


# ---- 9. Reshape: WIDE to LONG ("tidy") ----------------------------------------
# WIDE: one row per mouse, one column per week. Good for data entry and for a
# per-mouse summary. LONG (tidy): one row per mouse per week, with ONE column for
# the weight and ONE for the week. ggplot2 and most analyses want long data.
weights_long <- mice |>
  pivot_longer(
    cols         = starts_with("weight_g_wk"),
    names_to     = "week",
    names_prefix = "weight_g_wk",
    values_to    = "weight_g"
  ) |>
  mutate(week = as.integer(week)) |>
  select(mouse_id, genotype, sex, treatment, cage, week, weight_g)

head(weights_long, 8)
dim(mice); dim(weights_long)     # 40 rows -> 160 rows (40 mice x 4 weeks)

# Going back the other way (long to wide) is pivot_wider():
weights_long |>
  pivot_wider(names_from = week, names_prefix = "wk", values_from = weight_g) |>
  head(3)


# ---- 10. Save the CLEAN data (never over-write the raw files) --------------------
write_csv(mice,         here("output", "data", "mice_clean.csv"))
write_csv(weights_long, here("output", "data", "weights_long.csv"))
write_csv(cleaning_log, here("output", "data", "cleaning_log.csv"))
cleaning_log


# ---- 11. Descriptive tables ---------------------------------------------------
# Summary statistics by group. na.rm = TRUE skips missing values; n_valid counts
# the values actually used, so the reader can see how much data each mean rests on.
desc <- mice |>
  group_by(genotype, treatment) |>
  summarize(
    n          = n(),
    across(c(weight_g_wk3, distance_m, center_time_s),
           list(n_valid = ~ sum(!is.na(.x)), mean = ~ mean(.x, na.rm = TRUE),
                sd = ~ sd(.x, na.rm = TRUE), median = ~ median(.x, na.rm = TRUE)),
           .names = "{.col}__{.fn}"),
    .groups = "drop"
  )
desc

# Frequency table of a categorical variable
count(mice, genotype, sex) |>
  group_by(genotype) |>
  mutate(percent = 100 * n / sum(n))


# ---- 12. Displaying the data --------------------------------------------------
# GOLDEN RULES: show the data (not just summaries); pick the display to match the
# question and the variable type; label axes with units; use colorblind-friendly
# colors AND a second cue (shape or line type); do not truncate bar charts.

# 12a. Categorical data: a bar chart of counts (bars start at zero)
fig_counts <- mice |>
  count(genotype, treatment) |>
  ggplot(aes(treatment, n, fill = treatment)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = n), vjust = -0.4) +
  facet_wrap(~ genotype) +
  scale_fill_manual(values = unname(pal[c("blue", "vermillion")]), guide = "none") +
  scale_y_continuous(limits = c(0, 12), expand = c(0, 0)) +
  labs(x = "Treatment", y = "Number of mice")
print(fig_counts)

# 12b. One numeric variable: its DISTRIBUTION (shape, center, spread, outliers)
fig_hist <- ggplot(mice, aes(distance_m)) +
  geom_histogram(binwidth = 5, fill = pal["blue"], color = "white", na.rm = TRUE) +
  labs(x = "Distance traveled (m)", y = "Number of mice")

fig_hist_center <- ggplot(mice, aes(center_time_s)) +
  geom_histogram(binwidth = 10, fill = pal["orange"], color = "white", na.rm = TRUE) +
  labs(x = "Time in center (s)", y = "Number of mice")

fig_dist <- fig_hist + fig_hist_center + plot_annotation(tag_levels = "A")
print(fig_dist)
# Distance is roughly symmetric. Center time is skewed to the right (a few mice
# spend much longer): the mean is pulled up, so report the median too.

# 12c. A numeric variable by GROUP: points + box, never bars alone
fig_groups <- ggplot(mice, aes(treatment, distance_m, color = treatment, shape = treatment)) +
  geom_boxplot(width = 0.4, outlier.shape = NA, color = "gray40", fill = NA, na.rm = TRUE) +
  geom_jitter(width = 0.12, height = 0, size = 2.3, alpha = 0.85, na.rm = TRUE) +
  facet_wrap(~ genotype) +
  scale_color_manual(values = unname(pal[c("blue", "vermillion")])) +
  labs(x = "Treatment", y = "Distance traveled (m)") +
  theme(legend.position = "none")
print(fig_groups)

# 12d. Two numeric variables: a SCATTER plot (is there a relationship?)
fig_scatter <- ggplot(mice, aes(weight_g_wk3, distance_m, color = sex, shape = sex)) +
  geom_point(size = 2.5, alpha = 0.85, na.rm = TRUE) +
  scale_color_manual(values = unname(pal[c("pink", "green")])) +
  labs(x = "Body weight at week 3 (g)", y = "Distance traveled (m)",
       color = "Sex", shape = "Sex")
print(fig_scatter)

# 12e. Change over time: individual trajectories (thin lines) + group means
week_means <- weights_long |>
  group_by(genotype, treatment, week) |>
  summarize(mean = mean(weight_g, na.rm = TRUE),
            se = sd(weight_g, na.rm = TRUE) / sqrt(sum(!is.na(weight_g))),
            .groups = "drop")

fig_time <- ggplot(weights_long, aes(week, weight_g, color = treatment)) +
  geom_line(aes(group = mouse_id), alpha = 0.25, na.rm = TRUE) +
  geom_line(data = week_means, aes(y = mean, linetype = treatment), linewidth = 1.1) +
  geom_pointrange(data = week_means, aes(y = mean, ymin = mean - se, ymax = mean + se,
                                         shape = treatment), size = 0.45) +
  facet_wrap(~ genotype) +
  scale_color_manual(values = unname(pal[c("blue", "vermillion")])) +
  scale_x_continuous(breaks = 0:3) +
  labs(x = "Week", y = "Body weight (g)", color = "Treatment", linetype = "Treatment",
       shape = "Treatment")
print(fig_time)

# 12f. Misleading versus honest displays: the same data, three ways
bars_se <- mice |>
  group_by(treatment) |>
  summarize(mean = mean(distance_m, na.rm = TRUE),
            se = sd(distance_m, na.rm = TRUE) / sqrt(sum(!is.na(distance_m))))

bad1 <- ggplot(bars_se, aes(treatment, mean)) +
  geom_col(fill = "gray60", width = 0.5) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.15) +
  coord_cartesian(ylim = c(45, 60)) +            # TRUNCATED axis exaggerates the difference
  labs(title = "Truncated axis", x = NULL, y = "Distance traveled (m)")

bad2 <- ggplot(bars_se, aes(treatment, mean)) +
  geom_col(fill = "gray60", width = 0.5) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.15) +
  labs(title = "Bars hide the data", x = NULL, y = "Distance traveled (m)")

ci_tbl <- mice |>
  filter(!is.na(distance_m)) |>
  group_by(treatment) |>
  summarize(mean = mean(distance_m), n = n(),
            half = qt(0.975, n - 1) * sd(distance_m) / sqrt(n), .groups = "drop")

good <- ggplot(mice, aes(treatment, distance_m, color = treatment, shape = treatment)) +
  geom_jitter(width = 0.12, height = 0, size = 2, alpha = 0.7, na.rm = TRUE) +
  geom_pointrange(data = ci_tbl, aes(y = mean, ymin = mean - half, ymax = mean + half),
                  color = "black", shape = 16, size = 0.5, linewidth = 0.9,
                  position = position_nudge(x = 0.28)) +
  scale_color_manual(values = unname(pal[c("blue", "vermillion")])) +
  labs(title = "Show every point", x = NULL, y = "Distance traveled (m)") +
  theme(legend.position = "none")

fig_bad_good <- bad1 + bad2 + good + plot_annotation(tag_levels = "A")
print(fig_bad_good)

# 12g. A multi-panel figure for a report, saved at publication resolution
fig_report <- (fig_groups | fig_time) / (fig_dist) +
  plot_annotation(tag_levels = "A")
print(fig_report)

ggsave(here("output", "figures", "fig_distance_by_group.png"), fig_groups,
       width = 6.5, height = 4, dpi = 300)
ggsave(here("output", "figures", "fig_weight_over_time.png"), fig_time,
       width = 7.5, height = 4, dpi = 300)
ggsave(here("output", "figures", "fig_report_panels.png"), fig_report,
       width = 10, height = 8, dpi = 300)


# ---- 13. Record the computing environment -----------------------------------
sessionInfo()
