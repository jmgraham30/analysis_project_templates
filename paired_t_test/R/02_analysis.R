# =============================================================================
# 02_analysis.R
# Project:  Paired-samples t-test (open-field crossover: drug vs. vehicle)
# Purpose:  Complete, step-by-step analysis script. Run it from top to bottom
#           (or line by line with Ctrl/Cmd + Enter) to reproduce every number,
#           table and figure in the lab notebook.
#
# RESEARCH QUESTION
#   Do mice travel a different distance in an open field after a low dose of a
#   stimulant (Drug) than after a saline injection (Vehicle)? Each mouse was
#   tested under BOTH conditions (a counterbalanced crossover), so the data are
#   paired. A second part of the script shows why several neurons recorded from
#   the same mouse must NOT be treated as independent observations.
#
# BEFORE YOU RUN: open paired_t_test.Rproj so the working directory is the
#                 project folder.
# INPUT : Data/open_field_crossover.csv, Data/slice_recordings.csv
# OUTPUT: figures saved in output/figures/   (re-created each time you run this)
# =============================================================================


# ---- 1. Set up ---------------------------------------------------------------
source("R/00_setup.R")   # loads packages, plot theme, and helper functions
dir.create(here("output", "figures"), recursive = TRUE, showWarnings = FALSE)


# ---- 2. Decisions made BEFORE looking at the results -------------------------
alpha        <- 0.05   # significance level = accepted Type I error rate
target_power <- 0.80   # conventional minimum power
# Test: two-tailed PAIRED t-test (a one-sample t-test on the differences).
# Direction: Drug minus Vehicle. A POSITIVE difference means the mouse traveled
# FARTHER after the drug.


# ---- 3. Import the raw data --------------------------------------------------
of_raw <- read_csv(here("Data", "open_field_crossover.csv"), show_col_types = FALSE)

glimpse(of_raw)                  # one row per mouse per session ("long" format)
colSums(is.na(of_raw))           # missing values
count(of_raw, mouse_id) |> count(n)   # every mouse should appear exactly twice


# ---- 4. Prepare the data -----------------------------------------------------
of <- of_raw |>
  mutate(treatment = factor(treatment, levels = c("Vehicle", "Drug")),
         sex       = factor(sex, levels = c("F", "M"), labels = c("Female", "Male")))

# One row per mouse, the two measurements side by side, and the difference
of_wide <- of |>
  select(mouse_id, sex, sequence, treatment, distance_m) |>
  pivot_wider(names_from = treatment, values_from = distance_m) |>
  mutate(diff = Drug - Vehicle)

# A paired test needs BOTH values, so a mouse with one missing session is dropped
# from the test entirely (the other session is not used as an "unpaired" value).
lost_mice <- of_wide |> filter(is.na(diff))
pairs     <- of_wide |> drop_na(diff)
n_pairs   <- nrow(pairs)
lost_mice


# ---- 5. Describe the data ----------------------------------------------------
desc <- pairs |>
  select(mouse_id, Vehicle, Drug, Difference = diff) |>
  pivot_longer(-mouse_id, names_to = "Measure", values_to = "value") |>
  mutate(Measure = factor(Measure, levels = c("Vehicle", "Drug", "Difference"))) |>
  group_by(Measure) |>
  summarize(n = n(), M = mean(value), SD = sd(value), SE = SD / sqrt(n),
            Median = median(value), Min = min(value), Max = max(value),
            .groups = "drop")
desc

r_obs <- cor(pairs$Vehicle, pairs$Drug)   # correlation between a mouse's two values
r_obs

lims <- range(c(pairs$Vehicle, pairs$Drug)) + c(-2, 2)
p_scatter <- ggplot(pairs, aes(Vehicle, Drug)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray40") +
  geom_point(size = 2.5, color = pal["blue"]) +
  coord_equal(xlim = lims, ylim = lims) +
  labs(x = "Distance after vehicle (m)", y = "Distance after drug (m)")

p_bars <- pairs |>
  mutate(direction = if_else(diff >= 0, "Increase", "Decrease"),
         mouse_id  = fct_reorder(mouse_id, diff, .desc = TRUE)) |>
  ggplot(aes(mouse_id, diff, fill = direction)) +
  geom_col() + geom_hline(yintercept = 0) +
  scale_fill_manual(values = c(Increase = unname(pal["blue"]),
                               Decrease = unname(pal["vermillion"]))) +
  labs(x = "Mouse (sorted)", y = "Change in distance (m)", fill = NULL) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

fig_eda <- p_scatter + p_bars + plot_annotation(tag_levels = "A")
ggsave(here("output", "figures", "fig_eda.png"), fig_eda, width = 7, height = 3.8, dpi = 300)


# ---- 6. Check the assumptions (on the DIFFERENCES) ---------------------------
shapiro.test(pairs$diff)          # normality of the differences
pairs |> identify_outliers(diff)  # extreme differences (none = zero rows)

p_h <- ggplot(pairs, aes(diff)) +
  geom_histogram(binwidth = 4, fill = pal["sky"], color = "white") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x = "Drug - Vehicle (m)", y = "Number of mice")
p_q <- ggplot(pairs, aes(sample = diff)) +
  stat_qq(color = pal["blue"]) + stat_qq_line(color = "black") +
  labs(x = "Theoretical quantiles", y = "Sample quantiles")
ggsave(here("output", "figures", "fig_assumptions.png"),
       p_h + p_q + plot_annotation(tag_levels = "A"), width = 7, height = 3.4, dpi = 300)


# ---- 7. Power and sample size ------------------------------------------------
# For a paired test the effect size is dz = mean(D) / SD(D).
pwr.t.test(d = 0.8, sig.level = alpha, power = target_power, type = "paired")  # pairs needed
pwr.t.test(n = n_pairs, sig.level = alpha, power = target_power, type = "paired")  # smallest dz
pwr.t.test(n = n_pairs, d = 0.8, sig.level = alpha, type = "paired")$power
pwr.t.test(n = n_pairs, d = 0.5, sig.level = alpha, type = "paired")$power

# How the correlation between conditions changes power for the same raw effect
# (d = 0.5 in single-condition SD units): dz = d / sqrt(2 * (1 - r)).
tibble(r = c(0, 0.3, 0.6, 0.9)) |>
  mutate(dz = 0.5 / sqrt(2 * (1 - r)),
         power = map_dbl(dz, ~ pwr.t.test(n = n_pairs, d = .x, sig.level = alpha,
                                          type = "paired")$power),
         pairs_needed = map_dbl(dz, ~ ceiling(pwr.t.test(d = .x, power = target_power,
                                                          sig.level = alpha,
                                                          type = "paired")$n)))


# ---- 8. The paired t-test ----------------------------------------------------
# paired = TRUE matches the i-th Drug value with the i-th Vehicle value, so both
# columns must be in the same mouse order (they are: one row per mouse).
tt <- t.test(pairs$Drug, pairs$Vehicle, paired = TRUE,
             alternative = "two.sided", conf.level = 1 - alpha)
tt
tt_tidy   <- tidy(tt)
mean_diff <- tt_tidy$estimate
ci_low    <- tt_tidy$conf.low
ci_high   <- tt_tidy$conf.high

# Effect sizes
dz <- effectsize::cohens_d(pairs$diff)       # dz = mean(D) / SD(D), with 95% CI
gz <- effectsize::hedges_g(pairs$diff)       # small-sample-corrected dz
veh <- desc |> filter(Measure == "Vehicle"); drg <- desc |> filter(Measure == "Drug")
d_av <- mean_diff / sqrt((veh$SD^2 + drg$SD^2) / 2)   # comparable to a two-sample d
dz; gz; d_av

# What pairing bought us: the same data analyzed (incorrectly) as independent groups
t.test(pairs$Drug, pairs$Vehicle, paired = FALSE)       # Welch: wider CI, larger p
sqrt(veh$SD^2 + drg$SD^2 - 2 * r_obs * veh$SD * drg$SD) # equals SD of the differences
sd(pairs$diff)

# Robustness check: Wilcoxon signed-rank test (does not assume normal differences)
wilcox.test(pairs$Drug, pairs$Vehicle, paired = TRUE)


# ---- 9. Order effects in the crossover ---------------------------------------
# Counterbalancing spreads a session effect evenly across both treatments.
pairs |>
  group_by(sequence) |>
  summarize(pairs = n(), mean_diff = mean(diff), sd_diff = sd(diff))
# Half the gap between the two sequence means estimates the session effect.


# ---- 10. Main figure ----------------------------------------------------------
cond_long <- pairs |>
  select(mouse_id, Vehicle, Drug) |>
  pivot_longer(c(Vehicle, Drug), names_to = "treatment", values_to = "distance") |>
  mutate(treatment = factor(treatment, levels = c("Vehicle", "Drug")))

cond_sum <- cond_long |>
  group_by(treatment) |>
  summarize(m = mean(distance), se = sd(distance) / sqrt(n()), n = n(), .groups = "drop") |>
  mutate(lo = m - qt(0.975, n - 1) * se, hi = m + qt(0.975, n - 1) * se)

fig_a <- ggplot(cond_long, aes(treatment, distance)) +
  geom_line(aes(group = mouse_id), color = "gray60") +
  geom_point(aes(color = treatment), size = 2.4) +
  geom_pointrange(data = cond_sum, aes(y = m, ymin = lo, ymax = hi),
                  position = position_nudge(x = 0.2), color = "black",
                  shape = 16, size = 0.5, linewidth = 0.8) +
  scale_color_manual(values = unname(pal[c("blue", "vermillion")])) +
  labs(x = NULL, y = "Distance traveled (m)") +
  theme(legend.position = "none")

fig_b <- ggplot(tibble(m = mean_diff, lo = ci_low, hi = ci_high), aes(x = "Drug - Vehicle", y = m)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray30") +
  geom_pointrange(aes(ymin = lo, ymax = hi), color = pal["green"], size = 0.7, linewidth = 1) +
  labs(x = NULL, y = "Mean difference (m)") +
  theme(axis.ticks.x = element_blank())

fig_main <- fig_a + fig_b + plot_layout(widths = c(2, 1)) + plot_annotation(tag_levels = "A")
ggsave(here("output", "figures", "fig_main.png"), fig_main, width = 7, height = 4.5, dpi = 300)


# ---- 11. Type I error and power by simulation --------------------------------
simulate_paired_p <- function(n, dz, n_sims = 10000) {
  map_dbl(seq_len(n_sims), function(i) t.test(rnorm(n, mean = dz, sd = 1))$p.value)
}
set.seed(123)
p_h0 <- simulate_paired_p(n_pairs, dz = 0)      # H0 true
p_h1 <- simulate_paired_p(n_pairs, dz = 0.8)    # H1 true (large effect)
mean(p_h0 < alpha)    # Type I error rate (should be near .05)
mean(p_h1 < alpha)    # power at dz = 0.8 with our number of pairs


# =============================================================================
# PART B. MANY MEASUREMENTS FROM ONE ANIMAL (nested data)
# =============================================================================
slice_raw  <- read_csv(here("Data", "slice_recordings.csv"), show_col_types = FALSE)
slice_wide <- slice_raw |>
  pivot_wider(names_from = condition, values_from = firing_hz) |>
  mutate(change = Drug - Baseline)
stopifnot(!anyNA(slice_wide))

n_neurons    <- nrow(slice_wide)
n_slice_mice <- n_distinct(slice_wide$mouse_id)

# The unit of analysis is the MOUSE: average the neuron-level changes within mouse
mouse_means <- slice_wide |>
  group_by(mouse_id) |>
  summarize(neurons = n(), change = mean(change), .groups = "drop")
mouse_means

# WRONG: neurons treated as independent pairs (pseudoreplication)
t.test(slice_wide$Drug, slice_wide$Baseline, paired = TRUE)
# RIGHT: one value per mouse
correct <- t.test(mouse_means$change)
correct
effectsize::cohens_d(mouse_means$change)

# Intraclass correlation (ICC), design effect, and effective sample size
aov_tab <- anova(lm(change ~ mouse_id, data = slice_wide))
MSB <- aov_tab["mouse_id", "Mean Sq"]
MSW <- aov_tab["Residuals", "Mean Sq"]
sizes <- as.integer(table(slice_wide$mouse_id))
n0 <- (n_neurons - sum(sizes^2) / n_neurons) / (n_slice_mice - 1)
sd_between <- sqrt(max((MSB - MSW) / n0, 0))
sd_within  <- sqrt(MSW)
icc   <- sd_between^2 / (sd_between^2 + sd_within^2)
m_bar <- n_neurons / n_slice_mice
deff  <- 1 + (m_bar - 1) * icc
c(icc = icc, design_effect = deff, effective_n = n_neurons / deff)

# Simulate a world where the drug does nothing, with the same clustered structure
simulate_nested <- function(sizes, sd_b, sd_w, n_sims = 5000) {
  id <- rep(seq_along(sizes), times = sizes)
  map_dfr(seq_len(n_sims), function(i) {
    u   <- rnorm(length(sizes), 0, sd_b)
    chg <- rnorm(sum(sizes), mean = u[id], sd = sd_w)
    tibble(neurons = t.test(chg)$p.value,
           mice    = t.test(tapply(chg, id, mean))$p.value)
  })
}
set.seed(2468)
nested_sim <- simulate_nested(sizes, sd_between, sd_within)
colMeans(nested_sim < alpha)   # Type I error rate: neurons as units vs. mice as units

ggplot(slice_wide, aes(mouse_id, change)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray30") +
  geom_jitter(width = 0.1, height = 0, size = 2, alpha = 0.8, color = pal["blue"]) +
  geom_point(data = mouse_means, shape = 23, size = 4, fill = "white", color = "black") +
  labs(x = "Mouse", y = "Change in firing rate (Hz)")
ggsave(here("output", "figures", "fig_nested.png"), width = 7, height = 3.8, dpi = 300)
