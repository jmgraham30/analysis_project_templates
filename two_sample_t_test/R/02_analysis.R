# =============================================================================
# 02_analysis.R
# Project:  Two-sample t-test (elevated plus maze: stress vs. control)
# Purpose:  Complete, step-by-step analysis script. Run it from top to bottom
#           (or line by line with Ctrl/Cmd + Enter) to reproduce every number,
#           table and figure in the lab notebook.
#
# RESEARCH QUESTION
#   Do mice exposed to chronic mild stress spend a different percentage of time
#   in the open arms of an elevated plus maze than control mice? (Less open-arm
#   time is interpreted as higher anxiety-like behavior.)
#
# BEFORE YOU RUN: open two_sample_t_test.Rproj so the working directory is the
#                 project folder.
# INPUT : Data/elevated_plus_maze.csv   (raw data, never edited by hand)
# OUTPUT: figures saved in output/figures/   (re-created each time you run this)
# =============================================================================


# ---- 1. Set up ---------------------------------------------------------------
source("R/00_setup.R")   # loads packages, plot theme, and helper functions

# Create a folder for the figures this script saves (does nothing if it exists)
dir.create(here("output", "figures"), recursive = TRUE, showWarnings = FALSE)


# ---- 2. Decisions made BEFORE looking at the results -------------------------
# Good practice: write down your decisions up front so they cannot be (even
# unintentionally) adjusted after seeing the data.

alpha        <- 0.05   # significance level = accepted Type I error rate
target_power <- 0.80   # conventional minimum power
# Test: Welch's two-sample t-test, two-tailed (explained in section 8).
# Direction of the difference: Control minus Stress. A POSITIVE difference
# therefore means the control mice spent MORE time in the open arms.


# ---- 3. Import the raw data --------------------------------------------------
epm_raw <- read_csv(here("Data", "elevated_plus_maze.csv"),
                    show_col_types = FALSE)

# ALWAYS look at your data right after importing it.
glimpse(epm_raw)                      # variable names, types, first values
summary(epm_raw)                      # ranges: do any values look impossible?
colSums(is.na(epm_raw))               # missing values in each variable
anyDuplicated(epm_raw$mouse_id)       # duplicated IDs? (0 means none)
count(epm_raw, group)                 # how many mice were tested per group?


# ---- 4. Prepare the data -----------------------------------------------------
# The raw file is left untouched; all changes happen here, in code.

epm <- epm_raw |>
  mutate(
    # Put the groups in a meaningful order (Control first)
    group = factor(group, levels = c("Control", "Stress")),
    sex   = factor(sex, levels = c("F", "M"), labels = c("Female", "Male"))
  )

# Missing data: one mouse (S09) jumped off the maze, so its outcome is NA.
# We use a "complete-case" analysis: that mouse is left out of the test.
# ALWAYS report how many animals were lost and why.
missing_mice <- epm |> filter(is.na(open_arm_pct))
print(missing_mice)

epm_clean <- epm |> drop_na(open_arm_pct)

# Sample sizes recorded vs. analyzed, by group
n_table <- count(epm, group, name = "recorded") |>
  left_join(count(epm_clean, group, name = "analyzed"), by = "group")
print(n_table)


# ---- 5. Descriptive statistics -----------------------------------------------
desc <- epm_clean |>
  group_by(group) |>
  summarize(
    n      = n(),
    mean   = mean(open_arm_pct),
    sd     = sd(open_arm_pct),
    se     = sd / sqrt(n),
    median = median(open_arm_pct),
    min    = min(open_arm_pct),
    max    = max(open_arm_pct)
  )
print(desc)

# Check the other variables too: were the groups comparable on sex and general
# activity (total arm entries)? A big activity difference would make the
# open-arm result harder to interpret.
epm_clean |> count(group, sex)
epm_clean |>
  group_by(group) |>
  summarize(mean_entries = mean(total_entries), sd_entries = sd(total_entries))


# ---- 6. Exploratory plots: look at the distributions BEFORE testing ----------
# (a) Histograms by group: shape, spread, and overlap of the two groups.
p_hist <- ggplot(epm_clean, aes(open_arm_pct, fill = group)) +
  geom_histogram(binwidth = 5, color = "white") +
  facet_wrap(~ group, ncol = 1) +
  scale_fill_manual(values = unname(pal[c("blue", "vermillion")])) +
  labs(x = "Time in open arms (%)", y = "Number of mice") +
  theme(legend.position = "none")

# (b) Normal Q-Q plots by group: points near the line suggest normality.
p_qq <- ggplot(epm_clean, aes(sample = open_arm_pct, color = group)) +
  stat_qq() +
  stat_qq_line(color = "black") +
  facet_wrap(~ group, ncol = 1) +
  scale_color_manual(values = unname(pal[c("blue", "vermillion")])) +
  labs(x = "Theoretical quantiles", y = "Sample quantiles") +
  theme(legend.position = "none")

p_hist + p_qq + plot_annotation(tag_levels = "A")


# ---- 7. Check the assumptions of the two-sample t-test ----------------------
# 1. Independence : each mouse is in only one group and is measured once, and
#                   mice do not influence each other (design - cannot test;
#                   cage-mates can violate this, so consider housing).
# 2. Continuous   : open-arm percentage is a continuous measure (design).
# 3. Normality    : the outcome is approximately normal WITHIN EACH GROUP.
# 4. Equal spread : classic Student's t-test assumes equal variances; Welch's
#                   test (used here) does NOT need this assumption.
# 5. No extreme outliers.

# Normality within each group (H0: the data are normally distributed)
epm_clean |>
  group_by(group) |>
  reframe(tidy(shapiro.test(open_arm_pct)))

# Spread: compare the SDs, and (optionally) Levene's test of equal variances
desc |> summarize(sd_ratio = max(sd) / min(sd))   # ratios under ~2 are mild
epm_clean |> levene_test(open_arm_pct ~ group)    # rstatix; p > .05: no evidence

# Outliers within each group (rstatix flags values beyond 1.5 x IQR, and marks
# "extreme" ones beyond 3 x IQR)
epm_clean |>
  group_by(group) |>
  identify_outliers(open_arm_pct)


# ---- 8. The hypothesis test --------------------------------------------------
# H0: the two population means are equal      (mu_Control - mu_Stress = 0)
# H1: the two population means are different  (mu_Control - mu_Stress != 0)
#
# WELCH's t-test does not assume equal variances and is nearly as powerful as
# Student's test when variances ARE equal, so it is a safe default. In R,
# t.test() uses Welch's test unless you say var.equal = TRUE.
#
# Formula interface:  outcome ~ grouping variable.  R subtracts the SECOND level
# from the FIRST, so the difference is Control - Stress.

tt <- t.test(open_arm_pct ~ group, data = epm_clean,
             alternative = "two.sided", var.equal = FALSE,
             conf.level = 1 - alpha)
print(tt)

# tidy() from broom gives a one-row data frame; estimate1/estimate2 are the
# group means.
tt_tidy <- tidy(tt)
tt_tidy

t_value   <- tt_tidy$statistic
t_df      <- tt_tidy$parameter     # Welch df is usually NOT a whole number
p_value   <- tt_tidy$p.value
mean_diff <- tt_tidy$estimate      # Control mean minus Stress mean
ci_diff   <- c(tt_tidy$conf.low, tt_tidy$conf.high)

if (p_value < alpha) {
  cat("p < alpha: REJECT H0. The group means differ.\n")
} else {
  cat("p >= alpha: FAIL TO REJECT H0. No evidence the group means differ.\n")
}

# For comparison only: the classic Student's t-test (assumes equal variances).
# With similar spreads the answers are nearly identical.
student <- t.test(open_arm_pct ~ group, data = epm_clean, var.equal = TRUE) |>
  tidy()
student


# ---- 9. Effect size ----------------------------------------------------------
# The raw effect size is the mean difference in the ORIGINAL UNITS (percentage
# points of time in open arms): easy to interpret, so always report it.
# The standardized effect size, Cohen's d, expresses the difference in
# standard-deviation units so effects can be compared across studies:
#     d = (mean1 - mean2) / pooled SD
# Rough benchmarks (use with care): 0.2 small, 0.5 medium, 0.8 large.
# Hedges' g applies a small-sample correction (useful when groups are small).
# (We call effectsize:: explicitly because rstatix has a function of the same
# name.)

d <- effectsize::cohens_d(open_arm_pct ~ group, data = epm_clean)   # with 95% CI
g <- effectsize::hedges_g(open_arm_pct ~ group, data = epm_clean)
print(d)
print(g)
effectsize::interpret_cohens_d(d$Cohens_d, rules = "cohen1988")


# ---- 10. Publication-quality figure -----------------------------------------
# (A) Every mouse, the group means and 95% CIs (which show how precisely each
#     mean is estimated).
# (B) The estimated difference between groups with its 95% CI. If the interval
#     excludes 0 (dashed line), that agrees with a significant test.

group_ci <- epm_clean |>
  group_by(group) |>
  summarize(m = mean(open_arm_pct), se = sd(open_arm_pct) / sqrt(n()),
            n = n(), .groups = "drop") |>
  mutate(lo = m - qt(0.975, n - 1) * se,
         hi = m + qt(0.975, n - 1) * se)

fig_a <- ggplot(epm_clean, aes(group, open_arm_pct, color = group, shape = group)) +
  geom_boxplot(width = 0.35, outlier.shape = NA, color = "gray40", fill = NA) +
  geom_jitter(width = 0.1, height = 0, size = 2.2, alpha = 0.8) +
  geom_pointrange(data = group_ci,
                  aes(x = group, y = m, ymin = lo, ymax = hi),
                  position = position_nudge(x = 0.32), color = "black",
                  shape = 16, size = 0.5, linewidth = 0.8, inherit.aes = FALSE) +
  scale_color_manual(values = unname(pal[c("blue", "vermillion")])) +
  labs(x = NULL, y = "Time in open arms (%)") +
  theme(legend.position = "none")

diff_df <- tibble(m = mean_diff, lo = ci_diff[1], hi = ci_diff[2])

fig_b <- ggplot(diff_df, aes(x = "Control - Stress", y = m)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray30") +
  geom_pointrange(aes(ymin = lo, ymax = hi), color = pal["green"],
                  size = 0.7, linewidth = 1) +
  labs(x = NULL, y = "Mean difference (percentage points)") +
  theme(axis.ticks.x = element_blank())

fig_main <- fig_a + fig_b + plot_layout(widths = c(2, 1)) +
  plot_annotation(tag_levels = "A")
print(fig_main)

ggsave(here("output", "figures", "fig_open_arm_by_group.png"), fig_main,
       width = 7, height = 4.5, dpi = 300)


# ---- 11. Power analysis ------------------------------------------------------
# POWER = probability of rejecting H0 when it is really false (= 1 - beta).
# It depends on: effect size (d), sample sizes, alpha, and tails.

n1 <- sum(epm_clean$group == "Control")
n2 <- sum(epm_clean$group == "Stress")

# (a) A priori: mice needed PER GROUP to detect a large effect (d = 0.8) with
#     80% power (equal group sizes).
plan <- pwr.t.test(d = 0.8, sig.level = alpha, power = target_power,
                   type = "two.sample", alternative = "two.sided")
print(plan)   # n is per group; round UP

# (b) Sensitivity: with the (unequal) group sizes we actually have, what is the
#     SMALLEST effect we could detect with 80% power? (pwr.t2n.test handles
#     unequal n.)
sens <- pwr.t2n.test(n1 = n1, n2 = n2, sig.level = alpha, power = target_power,
                     alternative = "two.sided")
print(sens)   # look at d

# (c) Power to detect a large (d = 0.8) and a medium (d = 0.5) effect with our n
pwr.t2n.test(n1 = n1, n2 = n2, d = 0.8, sig.level = alpha)$power
pwr.t2n.test(n1 = n1, n2 = n2, d = 0.5, sig.level = alpha)$power

# (d) Power curve: how does power change with the number of mice per group?
power_curve <- expand_grid(n = 5:80, d = c(0.2, 0.5, 0.8)) |>
  mutate(power = map2_dbl(n, d, ~ pwr.t.test(n = .x, d = .y, sig.level = alpha,
                                             type = "two.sample")$power),
         d_label = factor(d, labels = c("Small (d = 0.2)", "Medium (d = 0.5)",
                                        "Large (d = 0.8)")))

# Lines differ in both color AND line type so they can be told apart in grayscale.
fig_power <- ggplot(power_curve, aes(n, power, color = d_label, linetype = d_label)) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = target_power, linetype = "dashed") +
  scale_color_manual(values = unname(pal[c("sky", "blue", "vermillion")])) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(x = "Mice per group", y = "Power",
       color = "True effect size", linetype = "True effect size")
print(fig_power)

ggsave(here("output", "figures", "fig_power_curve.png"), fig_power,
       width = 6, height = 4, dpi = 300)

# A caution about "observed power": computing power from the d you just
# observed adds no information beyond the p-value. Prefer the sensitivity
# analysis and confidence intervals above.


# ---- 12. Type I and Type II errors: simulating many experiments -------------
#                       | H0 is actually TRUE  | H0 is actually FALSE
#   --------------------+----------------------+----------------------------
#   Reject H0           | TYPE I error (alpha) | Correct (power = 1 - beta)
#   Fail to reject H0   | Correct (1 - alpha)  | TYPE II error (beta)
#
# We repeat the experiment 10,000 times in worlds where we KNOW the truth and
# count how often the t-test reaches the wrong conclusion.

simulate_p_values <- function(n1, n2, diff, sd1, sd2, var_equal = FALSE,
                              n_sims = 10000) {
  map_dbl(seq_len(n_sims), function(i) {
    x <- rnorm(n1, mean = diff, sd = sd1)   # "Control" group
    y <- rnorm(n2, mean = 0,    sd = sd2)   # "Stress" group
    t.test(x, y, var.equal = var_equal)$p.value
  })
}

set.seed(123)
sd_pooled <- sqrt(((n1 - 1) * sd(epm_clean$open_arm_pct[epm_clean$group == "Control"])^2 +
                   (n2 - 1) * sd(epm_clean$open_arm_pct[epm_clean$group == "Stress"])^2) /
                  (n1 + n2 - 2))
d_true <- 0.8   # a large true effect for scenario 2

# Scenario 1: H0 true (no real difference)   Scenario 2: H1 true (d = 0.8)
p_h0 <- simulate_p_values(n1, n2, diff = 0,                sd1 = sd_pooled, sd2 = sd_pooled)
p_h1 <- simulate_p_values(n1, n2, diff = d_true * sd_pooled, sd1 = sd_pooled, sd2 = sd_pooled)

type1_rate <- mean(p_h0 < alpha)    # should be close to alpha (.05)
power_sim  <- mean(p_h1 < alpha)    # should be close to pwr.t2n.test()
type2_rate <- 1 - power_sim

cat("Scenario 1 (H0 true): Type I error rate =", type1_rate, "\n")
cat("Scenario 2 (H1 true, d = 0.8): power =", power_sim,
    "| Type II error rate =", type2_rate, "\n")
cat("Theoretical power from pwr.t2n.test():",
    pwr.t2n.test(n1 = n1, n2 = n2, d = d_true, sig.level = alpha)$power, "\n")

pvals <- bind_rows(
  tibble(p = p_h0, scenario = "H0 true (no real difference)"),
  tibble(p = p_h1, scenario = "H1 true (d = 0.8)")
)

fig_pvals <- ggplot(pvals, aes(p)) +
  geom_histogram(breaks = seq(0, 1, 0.05), fill = pal["sky"], color = "white") +
  geom_vline(xintercept = alpha, linetype = "dashed", color = pal["vermillion"]) +
  facet_wrap(~ scenario) +
  labs(x = "p-value", y = "Number of simulated experiments")
print(fig_pvals)

ggsave(here("output", "figures", "fig_pvalue_simulation.png"), fig_pvals,
       width = 7, height = 3.5, dpi = 300)


# ---- 13. Why Welch's test is the safer default -------------------------------
# Student's test assumes equal variances. What if that is wrong AND the group
# sizes are unequal, with the SMALLER group more variable? Simulate a world
# with NO true difference (H0 is true) and see how often each test falsely
# rejects H0. Both tests should reject about 5% of the time.

set.seed(456)
p_welch   <- simulate_p_values(n1 = 10, n2 = 30, diff = 0, sd1 = 15, sd2 = 5,
                               var_equal = FALSE)
p_student <- simulate_p_values(n1 = 10, n2 = 30, diff = 0, sd1 = 15, sd2 = 5,
                               var_equal = TRUE)

tibble(test = c("Welch", "Student"),
       type1_rate = c(mean(p_welch < alpha), mean(p_student < alpha)))
# Student's Type I error rate is far above 5% here; Welch's stays near 5%.


# ---- 14. Robustness check: non-parametric alternative ------------------------
# The Wilcoxon rank-sum (Mann-Whitney) test compares the two groups using ranks
# and does not assume normality. If it agrees with the t-test, the conclusion is
# more convincing.
wilcox <- wilcox.test(open_arm_pct ~ group, data = epm_clean, exact = FALSE) |>
  tidy()
print(wilcox)


# ---- 15. Write the result sentence -------------------------------------------
# Build the APA-style sentence from the saved objects so the numbers can never
# be mistyped.
result_sentence <- paste0(
  "Control mice (n = ", n1, ", M = ", fmt_num(desc$mean[1]), ", SD = ",
  fmt_num(desc$sd[1]), ") spent ",
  ifelse(p_value < alpha, "significantly more", "not significantly more"),
  " time in the open arms than stressed mice (n = ", n2, ", M = ",
  fmt_num(desc$mean[2]), ", SD = ", fmt_num(desc$sd[2]), "), Welch's t(",
  fmt_num(t_df), ") = ", fmt_num(t_value), ", ", fmt_p(p_value),
  ", mean difference = ", fmt_num(mean_diff), " percentage points, 95% CI [",
  fmt_num(ci_diff[1]), ", ", fmt_num(ci_diff[2]), "], d = ", fmt_num(d$Cohens_d),
  ", 95% CI [", fmt_num(d$CI_low), ", ", fmt_num(d$CI_high), "]."
)
cat(result_sentence, "\n")


# ---- 16. Record the computing environment -----------------------------------
# Lets you (or a collaborator) see exactly which R and package versions made
# these results.
sessionInfo()
