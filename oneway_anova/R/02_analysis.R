# =============================================================================
# 02_analysis.R
# Project:  One-way ANOVA (drug dose and locomotor activity)
# Purpose:  Complete, step-by-step analysis script. Run it from top to bottom
#           (or line by line with Ctrl/Cmd + Enter) to reproduce every number,
#           table and figure in the lab notebook.
#
# RESEARCH QUESTION
#   Does the distance mice travel in an open field (locomotor activity) depend
#   on the dose of a stimulant drug (Vehicle, Low, Medium, High)?
#
# BEFORE YOU RUN: open oneway_anova.Rproj so the working directory is the
#                 project folder.
# INPUT : Data/locomotor_dose.csv   (raw data, never edited by hand)
# OUTPUT: figures saved in output/figures/   (re-created each time you run this)
# =============================================================================


# ---- 1. Set up ---------------------------------------------------------------
source("R/00_setup.R")   # loads packages, plot theme, and helper functions

# Create a folder for the figures this script saves (does nothing if it exists)
dir.create(here("output", "figures"), recursive = TRUE, showWarnings = FALSE)


# ---- 2. Decisions made BEFORE looking at the results -------------------------
alpha        <- 0.05   # significance level = accepted Type I error rate
target_power <- 0.80   # conventional minimum power
# Omnibus test: one-way ANOVA (4 independent groups).
# Follow-up: Tukey HSD comparisons of all 6 pairs, ONLY if the ANOVA is
# significant. (Dunnett's test, comparing each dose with the vehicle, is shown
# as an alternative when that is the only question of interest.)


# ---- 3. Import the raw data --------------------------------------------------
loco_raw <- read_csv(here("Data", "locomotor_dose.csv"), show_col_types = FALSE)

# ALWAYS look at your data right after importing it.
glimpse(loco_raw)                  # variable names, types, first values
summary(loco_raw)                  # ranges: do any values look impossible?
sum(is.na(loco_raw))               # how many missing values? (here: 0)
anyDuplicated(loco_raw$mouse_id)   # duplicated IDs? (0 means none)
count(loco_raw, dose_group, dose_mgkg)   # how many mice per group?


# ---- 4. Prepare the data -----------------------------------------------------
# The raw file is left untouched; all changes happen here, in code.
# The grouping variable must be a FACTOR with the levels in a sensible order
# (otherwise R sorts them alphabetically: High, Low, Medium, Vehicle).

loco <- loco_raw |>
  mutate(
    dose_group = factor(dose_group, levels = c("Vehicle", "Low", "Medium", "High")),
    sex        = factor(sex, levels = c("F", "M"), labels = c("Female", "Male"))
  )

k <- nlevels(loco$dose_group)   # number of groups
N <- nrow(loco)                 # total number of mice


# ---- 5. Descriptive statistics -----------------------------------------------
desc <- loco |>
  group_by(dose_group) |>
  summarize(
    n      = n(),
    mean   = mean(distance_m),
    sd     = sd(distance_m),
    se     = sd / sqrt(n),
    median = median(distance_m),
    min    = min(distance_m),
    max    = max(distance_m)
  )
print(desc)

# Check the other variable: is sex balanced across groups?
count(loco, dose_group, sex)


# ---- 6. Exploratory plot: look at the data BEFORE testing -------------------
# Every mouse (points), plus a box for the median and quartiles. viridis colors
# are colorblind-friendly and ordered from low to high dose.
fig_eda <- ggplot(loco, aes(dose_group, distance_m, color = dose_group,
                            shape = dose_group)) +
  geom_boxplot(width = 0.4, outlier.shape = NA, color = "gray40", fill = NA) +
  geom_jitter(width = 0.12, height = 0, size = 2.3, alpha = 0.85) +
  scale_color_viridis_d(end = 0.85) +
  labs(x = "Dose group", y = "Distance traveled (m)") +
  theme(legend.position = "none")
print(fig_eda)


# ---- 7. Fit the model and check the assumptions ------------------------------
# aov() fits the one-way ANOVA:   outcome ~ grouping variable
fit <- aov(distance_m ~ dose_group, data = loco)

# ANOVA assumptions:
# 1. Independence : each mouse is in one group and measured once (design).
# 2. Normality    : the RESIDUALS (each mouse's distance from its group mean)
#                   are approximately normal.
# 3. Equal spread : the groups have similar variances ("homogeneity of variance").
# 4. No extreme outliers.

res_df <- tibble(.fitted = fitted(fit), .resid = resid(fit))   # per-mouse fitted value and residual

# Normality of residuals: Shapiro-Wilk test (H0: normal) and Q-Q plot
shapiro.test(res_df$.resid) |> tidy()

p_qq <- ggplot(res_df, aes(sample = .resid)) +
  stat_qq(color = pal["blue"]) +
  stat_qq_line(color = "black") +
  labs(x = "Theoretical quantiles", y = "Residuals")

# Equal spread: residuals vs. fitted values should show a similar vertical
# spread in every group (no "funnel" shape).
p_res <- ggplot(res_df, aes(.fitted, .resid)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(color = pal["blue"], size = 2, alpha = 0.8) +
  labs(x = "Fitted values (group means)", y = "Residuals")

p_qq + p_res + plot_annotation(tag_levels = "A")

desc |> summarize(sd_ratio = max(sd) / min(sd))   # ratios under ~2 are mild
loco |> levene_test(distance_m ~ dose_group)       # rstatix; p > .05: no evidence

# Outliers within each group (rstatix flags values beyond 1.5 x IQR; "extreme"
# values are beyond 3 x IQR)
loco |>
  group_by(dose_group) |>
  identify_outliers(distance_m)


# ---- 8. The omnibus F test ---------------------------------------------------
# H0: all group population means are equal   (mu_V = mu_L = mu_M = mu_H)
# H1: at least one group mean is different
#
# F = (variability BETWEEN group means) / (variability WITHIN groups)
#   = MS_between / MS_within,   with df = (k - 1, N - k)
# A large F means the group means differ more than chance variation within
# groups would produce.

summary(fit)                  # the classic ANOVA table
anova_tbl <- tidy(fit)        # the same table as a data frame
anova_tbl

F_value <- anova_tbl$statistic[1]
df1     <- anova_tbl$df[1]
df2     <- anova_tbl$df[2]
p_anova <- anova_tbl$p.value[1]

if (p_anova < alpha) {
  cat("p < alpha: REJECT H0. At least one group mean differs.\n")
} else {
  cat("p >= alpha: FAIL TO REJECT H0. Do not run post hoc tests.\n")
}


# ---- 9. Effect size ----------------------------------------------------------
# Eta-squared (eta^2): the proportion of ALL variability in the outcome that is
# explained by group. Omega-squared (omega^2) is a less biased version.
# Rough benchmarks (use with care): eta^2 of .01 small, .06 medium, .14 large.
# Cohen's f is another version used in power analysis: .10 small, .25 medium,
# .40 large.
# For eta^2, a 90% CI (not 95%) is conventional because F tests are one-sided.

eta2   <- eta_squared(fit, ci = 0.90, alternative = "two.sided")
omega2 <- omega_squared(fit, ci = 0.90, alternative = "two.sided")
f_eff  <- cohens_f(fit, ci = 0.90, alternative = "two.sided")
print(eta2)
print(omega2)
print(f_eff)
interpret_eta_squared(eta2$Eta2, rules = "field2013")


# ---- 10. Post hoc comparisons: which groups differ? -------------------------
# A significant F only says "not all means are equal". To find WHERE the
# differences are we compare groups in pairs. But 4 groups make 6 pairs, and
# doing 6 ordinary t-tests would inflate the chance of at least one false
# positive (see section 13). TUKEY'S HSD adjusts for this, keeping the
# "family-wise" Type I error rate at alpha across all 6 comparisons.
#
# emmeans() estimates each group mean; contrast(..., "revpairwise") compares
# every pair as (later group - earlier group), so a positive difference means
# the higher-dose group traveled farther.

emm <- emmeans(fit, ~ dose_group)
emm

tukey <- contrast(emm, method = "revpairwise", adjust = "tukey")
tukey_tbl <- summary(tukey, infer = TRUE) |> as_tibble()   # estimate, CI, adj. p
tukey_tbl

# Standardized pairwise effect sizes (Cohen's d, using the pooled within-group SD)
d_tbl <- eff_size(emm, sigma = sigma(fit), edf = df.residual(fit),
                  method = "revpairwise") |>
  as_tibble()
d_tbl

# ALTERNATIVE: Dunnett's test compares each dose with ONE control (the vehicle).
# It has more power than Tukey when those 3 comparisons are all you need.
dunnett_tbl <- contrast(emm, method = "trt.vs.ctrl", ref = 1) |>   # adjust = Dunnett
  summary(infer = TRUE) |>
  as_tibble()
dunnett_tbl


# ---- 11. Publication-quality figure -----------------------------------------
# (A) Every mouse, with each group's mean and 95% CI.
# (B) The 6 pairwise differences with TUKEY-adjusted 95% CIs. A difference is
#     significant at the family-wise alpha if its interval excludes 0.

group_ci <- as_tibble(summary(emm))   # emmean, SE, df, lower.CL, upper.CL

fig_a <- ggplot(loco, aes(dose_group, distance_m, color = dose_group,
                          shape = dose_group)) +
  geom_boxplot(width = 0.4, outlier.shape = NA, color = "gray40", fill = NA) +
  geom_jitter(width = 0.1, height = 0, size = 2.2, alpha = 0.8) +
  geom_pointrange(data = group_ci,
                  aes(x = dose_group, y = emmean, ymin = lower.CL, ymax = upper.CL),
                  position = position_nudge(x = 0.33), color = "black",
                  shape = 16, size = 0.5, linewidth = 0.8, inherit.aes = FALSE) +
  scale_color_viridis_d(end = 0.85) +
  labs(x = "Dose group", y = "Distance traveled (m)") +
  theme(legend.position = "none")

fig_b <- tukey_tbl |>
  mutate(contrast = fct_rev(fct_inorder(contrast))) |>
  ggplot(aes(estimate, contrast)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray30") +
  geom_pointrange(aes(xmin = lower.CL, xmax = upper.CL), color = pal["green"],
                  size = 0.5, linewidth = 0.9) +
  labs(x = "Mean difference (m)\nwith Tukey-adjusted 95% CI", y = NULL)

fig_main <- fig_a + fig_b + plot_layout(widths = c(1.2, 1)) +
  plot_annotation(tag_levels = "A")
print(fig_main)

ggsave(here("output", "figures", "fig_distance_by_dose.png"), fig_main,
       width = 9, height = 4.5, dpi = 300)


# ---- 12. Power analysis ------------------------------------------------------
# POWER = probability of rejecting H0 when it is really false (= 1 - beta).
# For one-way ANOVA it depends on: effect size (Cohen's f), number of groups,
# n per group, and alpha.

n_per_group <- unique(desc$n)   # 12 (the design is balanced)

# (a) A priori: mice needed PER GROUP for 80% power at a medium effect (f = 0.25)
plan <- pwr.anova.test(k = k, f = 0.25, sig.level = alpha, power = target_power)
print(plan)   # round n UP

# (b) Sensitivity: with our n per group, what is the SMALLEST effect we could
#     detect with 80% power?
sens <- pwr.anova.test(k = k, n = n_per_group, sig.level = alpha,
                       power = target_power)
print(sens)   # look at f

# (c) Power to detect a medium (f = 0.25) and a large (f = 0.40) effect
pwr.anova.test(k = k, n = n_per_group, f = 0.25, sig.level = alpha)$power
pwr.anova.test(k = k, n = n_per_group, f = 0.40, sig.level = alpha)$power

# (d) Power curve
power_curve <- expand_grid(n = 5:80, f = c(0.10, 0.25, 0.40)) |>
  mutate(power = map2_dbl(n, f, ~ pwr.anova.test(k = k, n = .x, f = .y,
                                                 sig.level = alpha)$power),
         f_label = factor(f, labels = c("Small (f = 0.10)", "Medium (f = 0.25)",
                                        "Large (f = 0.40)")))

fig_power <- ggplot(power_curve, aes(n, power, color = f_label, linetype = f_label)) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = target_power, linetype = "dashed") +
  geom_vline(xintercept = n_per_group, linetype = "dotted") +
  scale_color_manual(values = unname(pal[c("sky", "blue", "vermillion")])) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(x = "Mice per group", y = "Power",
       color = "True effect size", linetype = "True effect size")
print(fig_power)

ggsave(here("output", "figures", "fig_power_curve.png"), fig_power,
       width = 6, height = 4, dpi = 300)

# (Avoid "observed power": computing power from the effect you just observed
# adds no information beyond the p-value.)


# ---- 13. Type I and Type II errors, and the multiple-comparisons problem ----
#                       | H0 is actually TRUE  | H0 is actually FALSE
#   --------------------+----------------------+----------------------------
#   Reject H0           | TYPE I error (alpha) | Correct (power = 1 - beta)
#   Fail to reject H0   | Correct (1 - alpha)  | TYPE II error (beta)
#
# We repeat the experiment 10,000 times in worlds where we KNOW the truth. In
# each simulated experiment (k groups of n mice) we record three p-values:
#   - the omnibus ANOVA p-value
#   - the SMALLEST of the 6 ordinary (unadjusted) pairwise t-test p-values
#   - the SMALLEST of the 6 Tukey-adjusted p-values
# For speed, the statistics are computed directly from group means and the pooled
# within-group variance (MSE), which is exactly what ANOVA and Tukey use.

simulate_anova <- function(means, sd, n, n_sims = 10000) {
  k   <- length(means)
  N   <- k * n
  idx <- combn(k, 2)                      # all pairs of groups (2 x 6 matrix)

  results <- map(seq_len(n_sims), function(i) {
    x   <- matrix(rnorm(N, mean = rep(means, each = n), sd = sd), nrow = n)
    m   <- colMeans(x)                    # group means
    mse <- mean(apply(x, 2, var))         # pooled within-group variance
    F   <- n * var(m) / mse               # F statistic (balanced design)
    big  <- max(abs(m[idx[1, ]] - m[idx[2, ]]))   # the LARGEST pairwise difference
    # The smallest p-value among all pairs always belongs to the largest difference.
    t    <- big / sqrt(2 * mse / n)       # its t statistic
    q    <- big / sqrt(mse / n)           # its studentized range statistic (Tukey)
    c(p_anova     = pf(F, k - 1, N - k, lower.tail = FALSE),
      p_min_unadj = 2 * pt(t, N - k, lower.tail = FALSE),
      p_min_tukey = ptukey(q, k, N - k, lower.tail = FALSE))
  })

  as_tibble(do.call(rbind, results))      # one row per simulated experiment
}

# Means that give a chosen Cohen's f when the SD is sd (equally spaced means)
means_for_f <- function(f, sd, k, base = 40) {
  step <- f * sd / sqrt((k^2 - 1) / 12)
  base + (0:(k - 1)) * step
}

set.seed(123)
sd_within <- sqrt(anova_tbl$meansq[2])   # our pooled within-group SD

sim_h0     <- simulate_anova(rep(40, k),                      sd_within, n_per_group)
sim_medium <- simulate_anova(means_for_f(0.25, sd_within, k), sd_within, n_per_group)
sim_large  <- simulate_anova(means_for_f(0.56, sd_within, k), sd_within, n_per_group)

sim_summary <- bind_rows(
  `H0 true (no real effect)`   = summarize(sim_h0,     across(everything(), ~ mean(.x < alpha))),
  `H1 true, medium (f = 0.25)` = summarize(sim_medium, across(everything(), ~ mean(.x < alpha))),
  `H1 true, large (f = 0.56)`  = summarize(sim_large,  across(everything(), ~ mean(.x < alpha))),
  .id = "world"
)
print(sim_summary)
# Row 1 is the Type I error rate for each approach: the omnibus ANOVA and
# the Tukey-adjusted comparisons stay near .05, but running all pairwise
# t-tests without adjustment finds a "significant" pair far more often.
# Rows 2-3 are power: the chance of detecting a real effect.

cat("Theoretical ANOVA power, f = 0.25:",
    pwr.anova.test(k = k, n = n_per_group, f = 0.25, sig.level = alpha)$power, "\n")

# Figure: false-positive rates when H0 is true
fig_fwer <- sim_summary |>
  filter(world == "H0 true (no real effect)") |>
  pivot_longer(-world, names_to = "approach", values_to = "rate") |>
  mutate(approach = factor(approach,
           levels = c("p_anova", "p_min_tukey", "p_min_unadj"),
           labels = c("ANOVA\n(omnibus F test)", "Tukey-adjusted\npairwise tests",
                      "Unadjusted pairwise\nt-tests"))) |>
  ggplot(aes(approach, rate, fill = approach)) +
  geom_col(width = 0.6) +
  geom_hline(yintercept = alpha, linetype = "dashed") +
  scale_fill_manual(values = unname(pal[c("blue", "green", "vermillion")])) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 0.3)) +
  labs(x = NULL, y = "Chance of at least one false positive") +
  theme(legend.position = "none")
print(fig_fwer)

ggsave(here("output", "figures", "fig_familywise_error.png"), fig_fwer,
       width = 6, height = 4, dpi = 300)


# ---- 14. Robustness checks ---------------------------------------------------
# (a) Welch's ANOVA does not assume equal variances.
welch <- oneway.test(distance_m ~ dose_group, data = loco, var.equal = FALSE) |>
  tidy()
print(welch)

# (b) Games-Howell post hoc test: the pairwise follow-up that goes with Welch's
#     ANOVA (does not assume equal variances).
loco |> games_howell_test(distance_m ~ dose_group)

# (c) Kruskal-Wallis test: a rank-based alternative that does not assume
#     normality.
kruskal <- kruskal.test(distance_m ~ dose_group, data = loco) |> tidy()
print(kruskal)


# ---- 15. Write the result sentence -------------------------------------------
# Build the APA-style sentence from the saved objects so the numbers can never
# be mistyped.
result_sentence <- paste0(
  "A one-way ANOVA ", ifelse(p_anova < alpha, "showed", "did not show"),
  " an effect of dose on distance traveled, F(", df1, ", ", df2, ") = ",
  fmt_num(F_value), ", ", fmt_p(p_anova), ", eta-squared = ",
  fmt_apa(eta2$Eta2), ", 90% CI [", fmt_apa(eta2$CI_low), ", ",
  fmt_apa(eta2$CI_high), "]."
)
cat(result_sentence, "\n")


# ---- 16. Record the computing environment -----------------------------------
sessionInfo()
