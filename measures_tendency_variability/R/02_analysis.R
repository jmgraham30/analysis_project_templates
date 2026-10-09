# =============================================================================
# 02_analysis.R
# Project:  Measures of central tendency and variability (neuron firing rates)
# Purpose:  Complete, step-by-step analysis script. Run it from top to bottom
#           (or line by line with Ctrl/Cmd + Enter) to reproduce every number,
#           table and figure in the lab notebook.
#
# RESEARCH QUESTIONS (all descriptive; no hypothesis tests)
#   For 120 neurons recorded in three brain regions: What is a "typical" firing
#   rate, and how much do neurons differ? Which summary statistics describe the
#   data honestly, given that firing rates are skewed? How precisely do we know
#   the average?
#
# BEFORE YOU RUN: open measures_tendency_variability.Rproj so the working
#                 directory is the project folder.
# INPUT : Data/neuron_firing.csv   (raw data, never edited by hand)
# OUTPUT: figures saved in output/figures/   (re-created each time you run this)
# =============================================================================


# ---- 1. Set up ---------------------------------------------------------------
source("R/00_setup.R")   # loads packages, plot theme, and helper functions

dir.create(here("output", "figures"), recursive = TRUE, showWarnings = FALSE)


# ---- 2. Decisions made BEFORE looking at the results -------------------------
conf_level <- 0.95     # confidence level for intervals
trim       <- 0.10     # fraction trimmed from EACH end for the trimmed mean
# Rules: (1) choose each summary to match the variable's scale and shape;
#        (2) never delete extreme values only because they are extreme: we flag
#            them and report the results with and without them;
#        (3) always report n, and say whether error bars are SD, SE, or a CI.


# ---- 3. Import the raw data --------------------------------------------------
neurons_raw <- read_csv(here("Data", "neuron_firing.csv"), show_col_types = FALSE)

glimpse(neurons_raw)
summary(neurons_raw)
sum(is.na(neurons_raw))               # missing values? (0)
anyDuplicated(neurons_raw$neuron_id)  # duplicated IDs? (0)

neurons <- neurons_raw |>
  mutate(
    region    = factor(region, levels = c("Cortex", "Hippocampus", "Striatum")),
    cell_type = factor(cell_type, levels = c("Regular spiking", "Fast spiking")),
    # an ORDINAL variable: the categories have an order, but the gaps between them
    # are not necessarily equal
    tuning_rating = factor(tuning_rating, levels = 1:5, ordered = TRUE)
  )
N <- nrow(neurons)
count(neurons, region)
x <- neurons$firing_rate_hz     # shorthand for the main variable
w <- neurons$spike_width_ms


# ---- 4. NOMINAL and ORDINAL variables: counts, mode, median ------------------
# Nominal (categories with no order): the only sensible "average" is the MODE
# (the most common category). Never compute a mean of category codes.
neurons |> count(cell_type) |> mutate(percent = 100 * n / sum(n))
mode_cat(neurons$cell_type)

# Ordinal (categories with an order): the median (and quartiles) make sense; the
# mean is debatable because the steps between ratings may not be equal.
neurons |> count(tuning_rating) |> mutate(percent = 100 * n / sum(n))
median(as.integer(neurons$tuning_rating))        # converting to numbers to find the middle
quantile(as.integer(neurons$tuning_rating), c(0.25, 0.75))
mode_cat(neurons$tuning_rating)


# ---- 5. Look at the shape first: histogram, density, Q-Q ---------------------
# The right summary depends on the SHAPE, so always look before you summarize.
p_hist_rate <- ggplot(neurons, aes(firing_rate_hz)) +
  geom_histogram(binwidth = 2, fill = pal["blue"], color = "white") +
  labs(x = "Firing rate (spikes/s)", y = "Number of neurons")

p_hist_width <- ggplot(neurons, aes(spike_width_ms)) +
  geom_histogram(binwidth = 0.04, fill = pal["orange"], color = "white") +
  labs(x = "Spike width (ms)", y = "Number of neurons")

p_hist_rate + p_hist_width + plot_annotation(tag_levels = "A")
# A: firing rate has a long RIGHT tail (a few neurons fire very fast).
# B: spike width is much more symmetric (with a left shoulder from fast-spiking cells).


# ---- 6. Measures of CENTRAL TENDENCY -----------------------------------------
# MEAN   : the balance point; uses every value, so it is pulled toward extreme values.
# MEDIAN : the middle value; unaffected by how extreme the extremes are.
# TRIMMED MEAN: the mean after cutting a fraction from each end (a compromise).
# GEOMETRIC MEAN: exp(mean(log(x))); the natural average for positive, skewed data
#                 (it estimates the median of a log-normal distribution).
# MODE   : the most common value; for continuous data, the peak of the density.

center_stats <- function(v) {
  dens <- density(v, from = min(v), to = max(v))
  tibble(
    mean         = mean(v),
    median       = median(v),
    trimmed_mean = mean(v, trim = trim),
    geometric    = exp(mean(log(v))),
    mode_density = dens$x[which.max(dens$y)]      # peak of the smoothed histogram
  )
}
center_rate  <- center_stats(x)
center_width <- center_stats(w)
bind_rows(`Firing rate (spikes/s)` = center_rate, `Spike width (ms)` = center_width,
          .id = "variable")

# For skewed data, mean > median. For symmetric data, they are close together.
fig_center <- ggplot(neurons, aes(firing_rate_hz)) +
  geom_histogram(binwidth = 2, fill = "gray85", color = "white") +
  geom_vline(data = pivot_longer(center_rate, c(mean, median, geometric),
                                 names_to = "statistic", values_to = "value"),
             aes(xintercept = value, color = statistic, linetype = statistic), linewidth = 1.1) +
  scale_color_manual(values = c(mean = unname(pal["vermillion"]), median = unname(pal["blue"]),
                                geometric = unname(pal["green"]))) +
  labs(x = "Firing rate (spikes/s)", y = "Number of neurons", color = NULL, linetype = NULL)
print(fig_center)


# ---- 7. Measures of VARIABILITY ----------------------------------------------
# RANGE   : max - min (depends on only 2 values, so it grows with n and outliers)
# IQR     : the spread of the middle 50% (Q3 - Q1); robust
# VARIANCE: average squared distance from the mean, with n - 1 in the denominator
# SD      : square root of the variance; in the same units as the data
# MAD     : median absolute deviation; a robust version of the SD
# CV      : SD / mean; unit-free, for comparing variability of different variables
# SE      : SD / sqrt(n); precision of the MEAN (not variability of the DATA)

spread_stats <- function(v) {
  tibble(
    n     = length(v),
    min   = min(v), max = max(v), range = max(v) - min(v),
    q1    = quantile(v, 0.25), q3 = quantile(v, 0.75), iqr = IQR(v),
    var   = var(v), sd = sd(v),
    mad   = mad(v),            # scaled to equal the SD for normal data
    cv    = cv(v),
    se    = se(v)
  )
}
spread_rate  <- spread_stats(x)
spread_width <- spread_stats(w)
bind_rows(`Firing rate (spikes/s)` = spread_rate, `Spike width (ms)` = spread_width,
          .id = "variable") |>
  print(width = Inf)

# The SD by hand, to see what the function does:
dev     <- x - mean(x)                    # each neuron's deviation from the mean
ss      <- sum(dev^2)                     # sum of squared deviations
var_by_hand <- ss / (length(x) - 1)       # divide by n - 1 (see section 12a)
c(by_hand = sqrt(var_by_hand), sd_function = sd(x))

# The coefficient of variation lets us compare variability on different scales:
c(rate_cv = cv(x), width_cv = cv(w))


# ---- 8. Shape: skewness and kurtosis -----------------------------------------
shape_tbl <- tibble(
  variable = c("Firing rate (spikes/s)", "Spike width (ms)", "log(Firing rate)"),
  skewness = c(skew(x), skew(w), skew(log(x))),
  excess_kurtosis = c(kurt(x), kurt(w), kurt(log(x)))
)
shape_tbl
# Rough guide: |skewness| < 0.5 is about symmetric; > 1 is strongly skewed.
# After a log transformation, the firing rates are much more symmetric.

p_qq1 <- ggplot(neurons, aes(sample = firing_rate_hz)) +
  stat_qq(color = pal["blue"]) + stat_qq_line() +
  labs(title = "Firing rate", x = "Theoretical quantiles", y = "Sample quantiles")
p_qq2 <- ggplot(neurons, aes(sample = log(firing_rate_hz))) +
  stat_qq(color = pal["blue"]) + stat_qq_line() +
  labs(title = "log(Firing rate)", x = "Theoretical quantiles", y = "Sample quantiles")
p_qq1 + p_qq2 + plot_annotation(tag_levels = "A")

# Geometric mean and "multiplicative" spread for log-normal data
exp(mean(log(x)))          # geometric mean (back-transformed mean of the logs)
exp(sd(log(x)))            # geometric SD: a typical neuron is within a factor of this


# ---- 9. How outliers change each statistic -----------------------------------
# Flag potential outliers with the 1.5 x IQR rule (the same rule box plots use)
q <- quantile(x, c(0.25, 0.75))
fence_hi <- q[2] + 1.5 * IQR(x)
fence_lo <- q[1] - 1.5 * IQR(x)
neurons |> filter(firing_rate_hz > fence_hi | firing_rate_hz < fence_lo) |>
  select(neuron_id, region, cell_type, firing_rate_hz) |>
  arrange(desc(firing_rate_hz))
# Do NOT delete these. They are plausible values (fast-spiking neurons in a
# fast region). Delete a value only if you have a documented reason to believe it is
# an ERROR (see the data organization project).

# Sensitivity: what if a recording artifact added one extreme value? We add one
# fake neuron to the data with firing rates from 0 to 200 spikes/s and watch
# how each statistic responds.
artifact_values <- c(0, 10, 25, 50, 100, 200)
sens <- map(artifact_values, function(a) {
  v <- c(x, a)
  tibble(artifact = a, mean = mean(v), sd = sd(v), median = median(v), iqr = IQR(v))
}) |> list_rbind()
sens
# The mean and SD move a lot; the median and IQR barely change.

fig_sens <- sens |>
  pivot_longer(-artifact, names_to = "statistic", values_to = "value") |>
  mutate(group = if_else(statistic %in% c("mean", "median"), "Center", "Spread"),
         statistic = factor(statistic, levels = c("mean", "median", "sd", "iqr"),
                            labels = c("Mean", "Median", "SD", "IQR"))) |>
  ggplot(aes(artifact, value, color = statistic, linetype = statistic, shape = statistic)) +
  geom_line(linewidth = 1) + geom_point(size = 2.5) +
  facet_wrap(~ group, scales = "free_y") +
  scale_color_manual(values = unname(pal[c("vermillion", "blue", "vermillion", "blue")])) +
  scale_linetype_manual(values = c("solid", "solid", "dashed", "dashed")) +
  labs(x = "Value of one added extreme observation (spikes/s)", y = "Statistic (spikes/s)",
       color = NULL, linetype = NULL, shape = NULL)
print(fig_sens)


# ---- 10. Summaries BY GROUP ---------------------------------------------------
group_tbl <- neurons |>
  group_by(region) |>
  summarize(
    n = n(),
    mean = mean(firing_rate_hz), sd = sd(firing_rate_hz),
    median = median(firing_rate_hz), q1 = quantile(firing_rate_hz, .25),
    q3 = quantile(firing_rate_hz, .75),
    geo_mean = exp(mean(log(firing_rate_hz))), cv = cv(firing_rate_hz),
    sd_log = sd(log(firing_rate_hz))
  )
group_tbl
# The SD is bigger where the mean is bigger (typical of skewed, ratio data), so
# the CV is similar across regions and the SD of the LOGS is about the same.

fig_groups <- ggplot(neurons, aes(region, firing_rate_hz, color = region, shape = region)) +
  geom_jitter(width = 0.12, height = 0, size = 2, alpha = 0.7) +
  geom_boxplot(width = 0.35, outlier.shape = NA, color = "gray30", fill = NA) +
  # The diamonds are the ARITHMETIC means (computed from the raw values). Using
  # stat_summary(fun = mean) here would be wrong on the log panel, where ggplot
  # would average the logged values and show the geometric mean instead.
  geom_point(data = group_tbl, aes(y = mean), shape = 23, size = 3.5, fill = "white", color = "black") +
  scale_color_manual(values = unname(pal[c("blue", "orange", "green")])) +
  labs(x = "Brain region", y = "Firing rate (spikes/s)") +
  theme(legend.position = "none")

fig_groups_log <- fig_groups +
  scale_y_log10(breaks = c(0.5, 1, 2, 5, 10, 20, 40)) +
  labs(y = "Firing rate (spikes/s, log scale)")
fig_groups + fig_groups_log + plot_annotation(tag_levels = "A")

# Medians and IQRs by cell type, too
neurons |>
  group_by(cell_type) |>
  summarize(n = n(), median = median(firing_rate_hz), iqr = IQR(firing_rate_hz),
            mean = mean(firing_rate_hz), sd = sd(firing_rate_hz))


# ---- 11. SD versus SE, and confidence intervals --------------------------------
# SD describes how much NEURONS differ. SE describes how much the SAMPLE MEAN would
# differ from sample to sample (it shrinks as n grows; the SD does not).
# A 95% CI for the mean is mean +- t* x SE.
mean_ci(x)
group_ci <- neurons |>
  group_by(region) |>
  reframe(mean_ci(firing_rate_hz, conf_level), sd = sd(firing_rate_hz), se = se(firing_rate_hz))
group_ci

# The mean's CI is symmetric, but the data are skewed. A better approach for
# skewed data: the CI of the GEOMETRIC mean (t interval on the log scale, then
# back-transformed) or a bootstrap CI for the MEDIAN.
geo_ci <- mean_ci(log(x)) |> mutate(across(c(mean, lower, upper), exp))
geo_ci

# BOOTSTRAP: re-sample the data with replacement many times, compute the statistic
# each time, and take the middle 95% of the results.
set.seed(123)
B <- 5000
boot_medians <- replicate(B, median(sample(x, replace = TRUE)))
quantile(boot_medians, c(0.025, 0.975))
median(x)


# ---- 12. Simulations: why the formulas are what they are ---------------------

# 12a. Why do we divide by n - 1 when computing a sample variance?
# (Each simulation below sets its own seed, exactly as in the notebook, so the
# numbers match the notebook.)
set.seed(123)
# Draw many small samples (n = 5) from a population whose TRUE variance is 1.
# Dividing the sum of squares by n gives an average that is too small ("biased");
# dividing by n - 1 gives an average of 1.
sim_var <- replicate(10000, {
  s <- rnorm(5)
  c(divide_by_n = sum((s - mean(s))^2) / 5, divide_by_n_minus_1 = var(s))
}) |> t() |> as_tibble()
sim_var |> summarize(across(everything(), mean))

# 12b. The SAMPLING DISTRIBUTION of the mean. We treat a log-normal distribution
set.seed(123)
# fitted to our data as the "population" and take many samples of size n.
mu_hat    <- mean(log(x))
sigma_hat <- sd(log(x))
pop_mean  <- exp(mu_hat + sigma_hat^2 / 2)                                   # true mean
pop_sd    <- sqrt((exp(sigma_hat^2) - 1) * exp(2 * mu_hat + sigma_hat^2))     # true SD

sample_means <- function(n, n_sims = 10000) {
  colMeans(matrix(rlnorm(n * n_sims, mu_hat, sigma_hat), nrow = n))
}
sims_means <- map(c(5, 20, 80), ~ tibble(n = .x, mean = sample_means(.x))) |> list_rbind()

sims_means |>
  group_by(n) |>
  summarize(sd_of_means = sd(mean), theory_se = pop_sd / sqrt(first(n)), skewness = skew(mean))
# The SD of the sample means matches pop_sd / sqrt(n), and the distribution of the
# means becomes more symmetric (more normal) as n grows: the CENTRAL LIMIT THEOREM.

fig_clt <- sims_means |>
  mutate(n_lab = factor(n, labels = paste("n =", c(5, 20, 80)))) |>
  ggplot(aes(mean)) +
  geom_histogram(bins = 40, fill = pal["sky"], color = "white") +
  geom_vline(xintercept = pop_mean, linetype = "dashed") +
  facet_wrap(~ n_lab, ncol = 1, scales = "free_y") +
  coord_cartesian(xlim = c(0, 25)) +
  labs(x = "Sample mean firing rate (spikes/s)", y = "Number of simulated samples")
print(fig_clt)

# 12c. How often does the usual 95% t interval for the mean contain the TRUE mean?
set.seed(123)
# (If the method worked perfectly, 95% of the time.)
coverage <- function(n, n_sims = 10000) {
  m <- matrix(rlnorm(n * n_sims, mu_hat, sigma_hat), nrow = n)
  mean_x <- colMeans(m)
  se_x   <- sqrt(apply(m, 2, var) / n)
  half   <- qt(0.975, n - 1) * se_x
  mean(mean_x - half <= pop_mean & pop_mean <= mean_x + half)
}
coverage_tbl <- tibble(n = c(5, 10, 30, 100, 400)) |>
  mutate(coverage = map_dbl(n, coverage))
coverage_tbl
# For skewed data and small n, the interval "misses" more than 5% of the time.
# Larger samples fix this, which is the Central Limit Theorem again.


# ---- 13. Why variability matters (a preview of effect size and power) --------
# The same difference between two group means is easy to see when the spread is
# small and hard to see when the spread is large. Later projects turn this idea
# into an effect size (the difference divided by the SD) and into power.
overlap_demo <- expand_grid(sd = c(1, 3), group = c("A", "B"), x = seq(-8, 14, length.out = 400)) |>
  mutate(mean = if_else(group == "A", 3, 6),
         density = dnorm(x, mean, sd),
         panel = factor(sd, labels = c("Small spread (SD = 1)", "Large spread (SD = 3)")))

fig_overlap <- ggplot(overlap_demo, aes(x, density, color = group, linetype = group)) +
  geom_line(linewidth = 1.1) +
  facet_wrap(~ panel) +
  scale_color_manual(values = unname(pal[c("blue", "vermillion")])) +
  labs(x = "Measurement", y = "Density", color = "Group", linetype = "Group")
print(fig_overlap)


# ---- 14. Save figures ---------------------------------------------------------
ggsave(here("output", "figures", "fig_distributions.png"),
       p_hist_rate + p_hist_width + plot_annotation(tag_levels = "A"), width = 8, height = 3.6, dpi = 300)
ggsave(here("output", "figures", "fig_center.png"), fig_center, width = 6.5, height = 4, dpi = 300)
ggsave(here("output", "figures", "fig_outlier_sensitivity.png"), fig_sens, width = 8, height = 4, dpi = 300)
ggsave(here("output", "figures", "fig_by_region.png"),
       fig_groups + fig_groups_log + plot_annotation(tag_levels = "A"), width = 8, height = 4, dpi = 300)
ggsave(here("output", "figures", "fig_sampling_distribution.png"), fig_clt, width = 6, height = 6, dpi = 300)


# ---- 15. Write the result sentences ------------------------------------------
cat("Firing rate: median = ", fmt_num(median(x), 1), " spikes/s (IQR ",
    fmt_num(quantile(x, .25), 1), " to ", fmt_num(quantile(x, .75), 1), "); mean = ",
    fmt_num(mean(x), 1), ", SD = ", fmt_num(sd(x), 1), "\n", sep = "")


# ---- 16. Record the computing environment -----------------------------------
sessionInfo()
