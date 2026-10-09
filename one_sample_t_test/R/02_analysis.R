# =============================================================================
# 02_analysis.R
# Project:  One-sample t-test (novel object recognition)
# Purpose:  Complete, step-by-step analysis script. Run it from top to bottom
#           (or line by line with Ctrl/Cmd + Enter) to reproduce every number,
#           table and figure in the lab notebook.
#
# RESEARCH QUESTION
#   Do mice show memory for a familiar object? If they do, the discrimination
#   index (DI) for the novel object should be greater than 0 (= no preference).
#
# BEFORE YOU RUN: open one_sample_t_test.Rproj so the working directory is the
#                 project folder.
# INPUT : Data/novel_object_recognition.csv   (raw data, never edited by hand)
# OUTPUT: figures saved in output/figures/    (re-created each time you run this)
# =============================================================================


# ---- 1. Set up ---------------------------------------------------------------
source("R/00_setup.R")   # loads packages, plot theme, and helper functions


# ---- 2. Decisions made BEFORE looking at the results -------------------------
# Good practice: write down your decisions up front so they cannot be (even
# unintentionally) adjusted after seeing the data.

mu0         <- 0      # the null value: DI = 0 means "no preference"
alpha       <- 0.05   # significance level = the Type I error rate we will accept
min_explore <- 20     # exclusion rule: mice that explored the objects for less
                      # than 20 s in total are dropped (too little information)
target_power <- 0.80  # conventional minimum power


# ---- 3. Import the raw data --------------------------------------------------
nor_raw <- read_csv(here("Data", "novel_object_recognition.csv"),
                    show_col_types = FALSE)

# ALWAYS look at your data right after importing it.
glimpse(nor_raw)                  # variable names, types, first values
summary(nor_raw)                  # ranges: do any values look impossible?
sum(is.na(nor_raw))               # how many missing values? (here: 0)
anyDuplicated(nor_raw$mouse_id)   # duplicated IDs? (0 means none)


# ---- 4. Prepare the data (derive new variables, apply exclusions) ------------
# The raw file is left untouched; all changes happen here, in code, where they
# are documented and repeatable.

nor <- nor_raw |>
  mutate(
    total_s = time_novel_s + time_familiar_s,                 # total exploration
    DI      = (time_novel_s - time_familiar_s) / total_s,     # discrimination index
    sex     = factor(sex, levels = c("F", "M"), labels = c("Female", "Male"))
  )

# Report how many mice are excluded and why (journals expect this).
n_recorded <- nrow(nor)
excluded   <- filter(nor, total_s < min_explore)
nor_clean  <- filter(nor, total_s >= min_explore)
n          <- nrow(nor_clean)

cat("Mice recorded:", n_recorded, "\n")
cat("Mice excluded (total exploration <", min_explore, "s):", nrow(excluded), "\n")
print(select(excluded, mouse_id, total_s, DI))
cat("Mice analyzed:", n, "\n")


# ---- 5. Descriptive statistics ----------------------------------------------
desc <- nor_clean |>
  summarise(
    n      = n(),
    mean   = mean(DI),
    sd     = sd(DI),
    se     = sd / sqrt(n),
    median = median(DI),
    min    = min(DI),
    max    = max(DI)
  )
print(desc)


# ---- 6. Exploratory plots: look at the distribution BEFORE testing -----------
# (a) Histogram: is the distribution roughly symmetric and bell-shaped?
p_hist <- ggplot(nor_clean, aes(DI)) +
  geom_histogram(binwidth = 0.1, fill = pal["sky"], color = "white") +
  geom_vline(xintercept = mu0, linetype = "dashed") +
  labs(x = "Discrimination index", y = "Number of mice")

# (b) Normal Q-Q plot: points near the line suggest approximate normality.
p_qq <- ggplot(nor_clean, aes(sample = DI)) +
  stat_qq(color = pal["blue"]) +
  stat_qq_line(color = "black") +
  labs(x = "Theoretical quantiles", y = "Sample quantiles")

p_hist + p_qq + plot_annotation(tag_levels = "A")


# ---- 7. Check the assumptions of the one-sample t-test ----------------------
# 1. Independence  : each mouse contributes one value (design - cannot test).
# 2. Interval data : DI is a continuous measure (design).
# 3. Normality     : the DI values should come from an approximately normal
#                    distribution (matters most when n is small).
# 4. No extreme outliers.

shapiro <- tidy(shapiro.test(nor_clean$DI))   # H0: data are normally distributed
print(shapiro)                                # p.value > .05 -> no evidence against normality

# Outliers: flag values more than 3 SDs from the mean.
nor_clean <- nor_clean |> mutate(z = (DI - mean(DI)) / sd(DI))
nor_clean |> filter(abs(z) > 3)           # an empty result means none flagged


# ---- 8. The hypothesis test --------------------------------------------------
# H0: the population mean DI equals 0   (mu = 0; no memory / no preference)
# H1: the population mean DI differs from 0 (two-tailed)
#
# t = (sample mean - mu0) / (sd / sqrt(n))     with df = n - 1

tt <- t.test(nor_clean$DI, mu = mu0, alternative = "two.sided",
             conf.level = 1 - alpha)
print(tt)

# tidy() from the broom package turns the result into a one-row data frame
tt_tidy <- tidy(tt)
tt_tidy

# Pull out the pieces we need for reporting
t_value <- tt_tidy$statistic
t_df    <- tt_tidy$parameter      # degrees of freedom (n - 1)
p_value <- tt_tidy$p.value
ci_mean <- c(tt_tidy$conf.low, tt_tidy$conf.high)   # 95% CI for the mean DI

# The decision rule: compare p to alpha
if (p_value < alpha) {
  cat("p < alpha: REJECT H0. The mean DI differs from", mu0, "\n")
} else {
  cat("p >= alpha: FAIL TO REJECT H0. No evidence the mean DI differs from",
      mu0, "\n")
}


# ---- 9. Effect size ----------------------------------------------------------
# A p-value says whether an effect is detectable, NOT how big it is. Cohen's d
# expresses the difference in standard-deviation units:
#     d = (sample mean - mu0) / sd
# Rough benchmarks (use with care, context matters): 0.2 small, 0.5 medium,
# 0.8 large. Hedges' g is d with a small-sample bias correction.

d <- cohens_d(nor_clean$DI, mu = mu0)    # includes a 95% CI
g <- hedges_g(nor_clean$DI, mu = mu0)
print(d)
print(g)
interpret_cohens_d(d$Cohens_d, rules = "cohen1988")


# ---- 10. Publication-quality figure -----------------------------------------
# Individual mice (points), the distribution (box), and the mean with its 95% CI.
# The dashed line is the null value, so the reader can judge the result at a
# glance.

mean_df <- tibble(group = "All mice", m = mean(nor_clean$DI),
                  lo = ci_mean[1], hi = ci_mean[2])

fig_main <- ggplot(nor_clean, aes(x = "All mice", y = DI)) +
  geom_hline(yintercept = mu0, linetype = "dashed", color = "gray30") +
  geom_boxplot(width = 0.35, outlier.shape = NA, color = "gray40",
               fill = NA) +
  geom_jitter(width = 0.08, height = 0, size = 2.2, alpha = 0.8,
              color = pal["blue"]) +
  geom_pointrange(data = mean_df,
                  aes(x = 1.3, y = m, ymin = lo, ymax = hi),
                  color = pal["vermillion"], size = 0.7, linewidth = 1) +
  annotate("text", x = 1.3, y = 0.45, label = "Mean and 95% CI",
           color = pal["vermillion"], size = 3.5) +
  labs(x = NULL, y = "Discrimination index") +
  theme(axis.ticks.x = element_blank())
print(fig_main)

ggsave(here("output", "figures", "fig_discrimination_index.png"), fig_main,
       width = 4, height = 4.5, dpi = 300)


# ---- 11. Power analysis ------------------------------------------------------
# POWER = probability of rejecting H0 when it is really false (= 1 - beta).
# It depends on: effect size (d), sample size (n), alpha, and tails.
# Power analysis is best done BEFORE collecting data, to plan the sample size.

# (a) A priori: how many mice would we need to detect a medium effect (d = 0.5)
#     with 80% power at alpha = .05?
plan <- pwr.t.test(d = 0.5, sig.level = alpha, power = target_power,
                   type = "one.sample", alternative = "two.sided")
print(plan)   # look at n (round UP to the next whole mouse)

# (b) Sensitivity: given the n we actually have, what is the SMALLEST effect we
#     could detect with 80% power?
sens <- pwr.t.test(n = n, sig.level = alpha, power = target_power,
                   type = "one.sample", alternative = "two.sided")
print(sens)   # look at d

# (c) Power curve: how does power change with sample size?
power_curve <- expand_grid(n = 5:80, d = c(0.2, 0.5, 0.8)) |>
  mutate(power = map2_dbl(n, d, ~ pwr.t.test(n = .x, d = .y, sig.level = alpha,
                                             type = "one.sample")$power),
         d_label = factor(d, labels = c("Small (d = 0.2)",
                                        "Medium (d = 0.5)",
                                        "Large (d = 0.8)")))

fig_power <- ggplot(power_curve, aes(n, power, color = d_label)) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = target_power, linetype = "dashed") +
  geom_vline(xintercept = n, linetype = "dotted") +
  scale_color_manual(values = unname(pal[c("sky", "blue", "vermillion")])) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(x = "Sample size (n)", y = "Power", color = "True effect size")
print(fig_power)

ggsave(here("output", "figures", "fig_power_curve.png"), fig_power,
       width = 6, height = 4, dpi = 300)

# A caution about "observed power": computing power from the d you just
# observed adds no information beyond the p-value (they are mathematically
# linked). Prefer the sensitivity analysis and confidence intervals above.


# ---- 12. Type I and Type II errors: simulating many experiments -------------
# Imagine repeating this experiment 10,000 times in a world where we KNOW the
# truth, and counting how often the t-test reaches the wrong conclusion.
#
#                       | H0 is actually TRUE  | H0 is actually FALSE
#   --------------------+----------------------+----------------------------
#   Reject H0           | TYPE I error (alpha) | Correct (power = 1 - beta)
#   Fail to reject H0   | Correct (1 - alpha)  | TYPE II error (beta)

simulate_p_values <- function(true_mean, sd, n, mu0 = 0, n_sims = 10000) {
  map_dbl(seq_len(n_sims),
          ~ t.test(rnorm(n, mean = true_mean, sd = sd), mu = mu0)$p.value)
}

set.seed(123)
sd_est <- sd(nor_clean$DI)            # use our sample SD as the "true" SD
d_true <- 0.5                         # a medium true effect for scenario 2

p_h0 <- simulate_p_values(true_mean = mu0,              sd = sd_est, n = n)
p_h1 <- simulate_p_values(true_mean = mu0 + d_true * sd_est, sd = sd_est, n = n)

type1_rate <- mean(p_h0 < alpha)   # should be close to alpha (.05)
power_sim  <- mean(p_h1 < alpha)   # should be close to the pwr.t.test() answer
type2_rate <- 1 - power_sim

cat("Scenario 1 (H0 true):  Type I error rate =", type1_rate, "\n")
cat("Scenario 2 (H1 true, d = 0.5): power =", power_sim,
    "| Type II error rate =", type2_rate, "\n")
cat("Theoretical power from pwr.t.test():",
    pwr.t.test(n = n, d = d_true, sig.level = alpha, type = "one.sample")$power, "\n")

# Plot the distribution of p-values in each scenario
pvals <- bind_rows(
  tibble(p = p_h0, scenario = "H0 true (no real effect)"),
  tibble(p = p_h1, scenario = "H1 true (d = 0.5)")
)

fig_pvals <- ggplot(pvals, aes(p)) +
  geom_histogram(breaks = seq(0, 1, 0.05), fill = pal["sky"], color = "white") +
  geom_vline(xintercept = alpha, linetype = "dashed", color = pal["vermillion"]) +
  facet_wrap(~ scenario) +
  labs(x = "p-value", y = "Number of simulated experiments")
print(fig_pvals)

ggsave(here("output", "figures", "fig_pvalue_simulation.png"), fig_pvals,
       width = 7, height = 3.5, dpi = 300)


# ---- 13. Robustness check: non-parametric alternative ------------------------
# If normality looks doubtful, the Wilcoxon signed-rank test is a rank-based
# alternative. If it agrees with the t-test, the conclusion is more convincing.
wilcox <- wilcox.test(nor_clean$DI, mu = mu0, conf.int = TRUE, exact = FALSE) |>
  tidy()
print(wilcox)


# ---- 14. Write the result sentence -------------------------------------------
# Build the APA-style sentence from the saved objects so the numbers can never
# be mistyped.
result_sentence <- paste0(
  "Mice (N = ", n, ") showed a mean discrimination index of ",
  fmt_num(mean(nor_clean$DI)), " (SD = ", fmt_num(sd(nor_clean$DI)), "), ",
  "which ", ifelse(p_value < alpha, "differed", "did not differ"),
  " from zero, t(", t_df, ") = ", fmt_num(t_value), ", ", fmt_p(p_value),
  ", d = ", fmt_num(d$Cohens_d), ", 95% CI [", fmt_num(d$CI_low), ", ",
  fmt_num(d$CI_high), "]."
)
cat(result_sentence, "\n")


# ---- 15. Record the computing environment -----------------------------------
# Lets you (or a collaborator) see exactly which R and package versions made
# these results.
sessionInfo()
