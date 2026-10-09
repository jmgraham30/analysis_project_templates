# =============================================================================
# 02_analysis.R
# Project:  Two-way ANOVA (housing, stress, and hippocampal BDNF)
# Purpose:  Complete, step-by-step analysis script. Run it from top to bottom
#           (or line by line with Ctrl/Cmd + Enter) to reproduce every number,
#           table and figure in the lab notebook.
#
# RESEARCH QUESTION
#   Chronic stress is expected to lower hippocampal BDNF. Does living in an
#   enriched environment (instead of a standard cage) change BDNF, and does it
#   change how much stress lowers BDNF (that is, is there an interaction)?
#
# BEFORE YOU RUN: open twoway_anova.Rproj so the working directory is the
#                 project folder.
# INPUT : Data/bdnf_housing_stress.csv   (raw data, never edited by hand)
# OUTPUT: figures saved in output/figures/   (re-created each time you run this)
# =============================================================================


# ---- 1. Set up ---------------------------------------------------------------
source("R/00_setup.R")   # loads packages, plot theme, and helper functions

# Create a folder for the figures this script saves (does nothing if it exists)
dir.create(here("output", "figures"), recursive = TRUE, showWarnings = FALSE)


# ---- 2. Decisions made BEFORE looking at the results -------------------------
alpha        <- 0.05   # significance level = accepted Type I error rate
target_power <- 0.80   # conventional minimum power
# Test: 2 x 2 between-subjects ANOVA with three questions, each with its own F test:
#   (1) main effect of housing, (2) main effect of stress, (3) housing x stress
#   interaction.
# Plan: look at the INTERACTION first. If it is significant, the main effects can
# be misleading, so we follow up with SIMPLE EFFECTS (the effect of stress within
# each housing condition), adjusted for the 2 comparisons with the Holm method.


# ---- 3. Import the raw data --------------------------------------------------
bdnf_raw <- read_csv(here("Data", "bdnf_housing_stress.csv"), show_col_types = FALSE)

# ALWAYS look at your data right after importing it.
glimpse(bdnf_raw)                  # variable names, types, first values
summary(bdnf_raw)                  # ranges: do any values look impossible?
sum(is.na(bdnf_raw))               # how many missing values? (here: 0)
anyDuplicated(bdnf_raw$mouse_id)   # duplicated IDs? (0 means none)
count(bdnf_raw, housing, stress)   # how many mice per cell?


# ---- 4. Prepare the data -----------------------------------------------------
# The raw file is left untouched; all changes happen here, in code. Both factors
# must be FACTORS with the levels in a sensible order (the first level is the
# "reference": Standard housing and Control).
bdnf <- bdnf_raw |>
  mutate(
    housing = factor(housing, levels = c("Standard", "Enriched")),
    stress  = factor(stress,  levels = c("Control", "Stress")),
    sex     = factor(sex, levels = c("F", "M"), labels = c("Female", "Male")),
    cell    = interaction(housing, stress, sep = " / ")   # one label per cell
  )

# Sum-to-zero ("effect") coding, the usual setup for factorial designs. With equal
# n per cell, as here, it changes nothing. If the cells ever have different n, you
# must ALSO request Type III tests (for example car::Anova(fit, type = 3)); aov()
# alone gives sequential (Type I) tests that depend on the order of the terms.
options(contrasts = c("contr.sum", "contr.poly"))

n_cell <- unique(count(bdnf, housing, stress)$n)   # 10 (balanced)
N      <- nrow(bdnf)


# ---- 5. Descriptive statistics -----------------------------------------------
# Cell means (the 4 groups), then the marginal means (each factor's levels
# averaged over the other factor).
desc_cells <- bdnf |>
  group_by(housing, stress) |>
  summarize(n = n(), mean = mean(bdnf_pgmg), sd = sd(bdnf_pgmg),
            se = sd / sqrt(n), median = median(bdnf_pgmg),
            min = min(bdnf_pgmg), max = max(bdnf_pgmg), .groups = "drop")
print(desc_cells)

bdnf |> group_by(housing) |> summarize(n = n(), mean = mean(bdnf_pgmg), sd = sd(bdnf_pgmg))
bdnf |> group_by(stress)  |> summarize(n = n(), mean = mean(bdnf_pgmg), sd = sd(bdnf_pgmg))

count(bdnf, housing, stress, sex)   # is sex balanced across cells?


# ---- 6. Exploratory plot: look at the data BEFORE testing -------------------
# Color AND shape both show housing, so the figure works in grayscale and for
# readers with color blindness.
dodge <- position_dodge(width = 0.6)

fig_eda <- ggplot(bdnf, aes(stress, bdnf_pgmg, color = housing, shape = housing)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.6),
             size = 2.3, alpha = 0.85) +
  stat_summary(fun = median, geom = "crossbar", width = 0.45, linewidth = 0.4,
               position = dodge, show.legend = FALSE) +
  scale_color_manual(values = unname(pal[c("blue", "vermillion")])) +
  labs(x = "Stress condition", y = "Hippocampal BDNF (pg/mg protein)",
       color = "Housing", shape = "Housing")
print(fig_eda)


# ---- 7. Fit the model and check the assumptions ------------------------------
# aov(outcome ~ housing * stress) fits BOTH main effects AND their interaction.
# (housing * stress is shorthand for housing + stress + housing:stress.)
fit <- aov(bdnf_pgmg ~ housing * stress, data = bdnf)

# ANOVA assumptions (the same as one-way ANOVA, applied to the 4 cells):
# 1. Independence : each mouse is in one cell and measured once (design).
# 2. Normality    : the RESIDUALS are approximately normal.
# 3. Equal spread : the 4 cells have similar variances.
# 4. No extreme outliers.

res_df <- tibble(.fitted = fitted(fit), .resid = resid(fit))

shapiro.test(res_df$.resid) |> tidy()

p_qq <- ggplot(res_df, aes(sample = .resid)) +
  stat_qq(color = pal["blue"]) +
  stat_qq_line(color = "black") +
  labs(x = "Theoretical quantiles", y = "Residuals")

p_res <- ggplot(res_df, aes(.fitted, .resid)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(color = pal["blue"], size = 2, alpha = 0.8) +
  labs(x = "Fitted values (cell means)", y = "Residuals")

p_qq + p_res + plot_annotation(tag_levels = "A")

desc_cells |> summarize(sd_ratio = max(sd) / min(sd))   # ratios under ~2 are mild
bdnf |> levene_test(bdnf_pgmg ~ cell)                   # p > .05: no evidence of unequal spread

bdnf |>
  group_by(housing, stress) |>
  identify_outliers(bdnf_pgmg)   # beyond 1.5 x IQR; "extreme" beyond 3 x IQR


# ---- 8. The ANOVA table: three F tests ---------------------------------------
# H0 (housing)     : mean BDNF is the same in standard and enriched housing
#                    (averaged over stress).
# H0 (stress)      : mean BDNF is the same in control and stressed mice
#                    (averaged over housing).
# H0 (interaction) : the effect of stress is the same in both housing conditions.
#
# Each effect has its own F = MS_effect / MS_error with df = (1, N - 4).

summary(fit)
anova_tbl <- tidy(fit)     # term, df, sumsq, meansq, statistic, p.value
anova_tbl

df_error <- anova_tbl$df[anova_tbl$term == "Residuals"]
get_row  <- function(term) filter(anova_tbl, .data$term == !!term)
a_housing <- get_row("housing")
a_stress  <- get_row("stress")
a_inter   <- get_row("housing:stress")


# ---- 9. Effect sizes ---------------------------------------------------------
# PARTIAL eta-squared: the proportion of the variability left over after removing
# the other effects that is explained by this effect. (In a one-way ANOVA it equals
# plain eta-squared.) Omega-squared is a less biased estimate.
# Rough benchmarks (use with care): partial eta^2 of .01 small, .06 medium,
# .14 large. A 90% CI is conventional: eta^2 cannot be
# negative and the F test is one-sided, so it matches a two-sided test at alpha = .05.

eta2   <- eta_squared(fit,   partial = TRUE, ci = 0.90, alternative = "two.sided")
omega2 <- omega_squared(fit, partial = TRUE, ci = 0.90, alternative = "two.sided")
f_eff  <- cohens_f(fit,      partial = TRUE, ci = 0.90, alternative = "two.sided")
print(eta2)
print(omega2)
print(f_eff)


# ---- 10. Interaction and simple effects -------------------------------------
# If the interaction is real, "the effect of stress" has no single answer: it
# depends on housing. So we look at the SIMPLE EFFECTS: the stress effect (Stress
# minus Control) separately within each housing condition.

emm_cells <- emmeans(fit, ~ housing * stress)      # the 4 cell means with 95% CIs
emm_cells

emm_by_housing <- emmeans(fit, ~ stress | housing)

simple_stress <- contrast(emm_by_housing, method = "revpairwise", adjust = "none") |>
  rbind(adjust = "holm")                       # Holm adjustment across the 2 tests
simple_tbl <- summary(simple_stress, infer = TRUE) |> as_tibble()
simple_tbl

# Standardized simple effects (Cohen's d, pooled within-cell SD)
simple_d <- eff_size(emm_by_housing, sigma = sigma(fit), edf = df.residual(fit),
                     method = "revpairwise") |> as_tibble()
simple_d

# The interaction contrast itself: the DIFFERENCE between the two stress effects.
# (Stress - Control in Enriched) - (Stress - Control in Standard)
interaction_contrast <- contrast(emmeans(fit, ~ housing * stress),
                                 interaction = c("revpairwise", "revpairwise")) |>
  summary(infer = TRUE) |> as_tibble()
interaction_contrast

# Marginal means (for reporting main effects, if they are of interest)
emmeans(fit, ~ housing)
emmeans(fit, ~ stress)


# ---- 11. Publication-quality figure -----------------------------------------
# (A) Every mouse, with each cell's mean and 95% CI joined by lines (an
#     "interaction plot"). Non-parallel lines suggest an interaction.
# (B) The stress effect (Stress - Control) within each housing condition, with
#     Bonferroni-adjusted 95% CIs (emmeans gives Bonferroni intervals whenever
#     p-values are Holm-adjusted; Holm adjusts the p-values only).

cell_ci <- as_tibble(summary(emm_cells))   # emmean, SE, df, lower.CL, upper.CL

fig_a <- ggplot(bdnf, aes(stress, bdnf_pgmg, color = housing, shape = housing,
                          linetype = housing, group = housing)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.5),
             size = 1.8, alpha = 0.45) +
  geom_line(data = cell_ci, aes(y = emmean), position = position_dodge(width = 0.5),
            linewidth = 0.9) +
  geom_pointrange(data = cell_ci, aes(y = emmean, ymin = lower.CL, ymax = upper.CL),
                  position = position_dodge(width = 0.5), size = 0.6, linewidth = 0.9) +
  scale_color_manual(values = unname(pal[c("blue", "vermillion")])) +
  labs(x = "Stress condition", y = "Hippocampal BDNF (pg/mg protein)",
       color = "Housing", shape = "Housing", linetype = "Housing")

fig_b <- simple_tbl |>
  mutate(housing = factor(housing, levels = c("Standard", "Enriched"))) |>
  ggplot(aes(estimate, housing, color = housing, shape = housing)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray30") +
  geom_pointrange(aes(xmin = lower.CL, xmax = upper.CL), size = 0.6, linewidth = 0.9) +
  scale_color_manual(values = unname(pal[c("blue", "vermillion")])) +
  scale_y_discrete(limits = rev) +
  labs(x = "Effect of stress on BDNF (pg/mg)\nwith Bonferroni-adjusted 95% CI", y = "Housing") +
  theme(legend.position = "none")

fig_main <- fig_a + fig_b + plot_layout(widths = c(1.2, 1)) +
  plot_annotation(tag_levels = "A")
print(fig_main)

ggsave(here("output", "figures", "fig_bdnf_housing_stress.png"), fig_main,
       width = 9, height = 4.5, dpi = 300)


# ---- 12. Power analysis ------------------------------------------------------
# POWER = probability of rejecting H0 when it is really false (= 1 - beta).
# For each effect in a 2 x 2 ANOVA (numerator df = 1), power depends on Cohen's
# f-squared (= partial eta^2 / (1 - partial eta^2)), the total N, and alpha.
# We compute it directly from the noncentral F distribution:
#   noncentrality = f^2 * N,   df = (1, N - 4)

power_2x2 <- function(f2, n, alpha = 0.05, cells = 4) {
  N   <- n * cells
  df2 <- N - cells
  pf(qf(1 - alpha, 1, df2), 1, df2, ncp = f2 * N, lower.tail = FALSE)
}
f2_medium <- 0.25^2    # Cohen's "medium" f = .25
f2_large  <- 0.40^2    # Cohen's "large"  f = .40

# (a) A priori: mice needed PER CELL for 80% power to detect a medium effect
n_needed <- function(f2) {
  n <- 3
  while (power_2x2(f2, n, alpha) < target_power) n <- n + 1
  n
}
n_medium <- n_needed(f2_medium)
n_medium
n_medium * 4   # total mice

# (b) Sensitivity: the smallest f^2 detectable with 80% power at our n per cell
f2_min <- uniroot(function(f2) power_2x2(f2, n_cell, alpha) - target_power,
                  interval = c(0.001, 5))$root
sqrt(f2_min)   # as Cohen's f

# (c) Power to detect a medium and a large effect with our n
power_2x2(f2_medium, n_cell, alpha)
power_2x2(f2_large,  n_cell, alpha)

# (d) Power curve
power_curve <- expand_grid(n = 3:80, f = c(0.10, 0.25, 0.40)) |>
  mutate(power   = map2_dbl(n, f, ~ power_2x2(.y^2, .x, alpha)),
         f_label = factor(f, labels = c("Small (f = 0.10)", "Medium (f = 0.25)",
                                        "Large (f = 0.40)")))

fig_power <- ggplot(power_curve, aes(n, power, color = f_label, linetype = f_label)) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = target_power, linetype = "dashed") +
  geom_vline(xintercept = n_cell, linetype = "dotted") +
  scale_color_manual(values = unname(pal[c("sky", "blue", "vermillion")])) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(x = "Mice per cell", y = "Power",
       color = "True effect size", linetype = "True effect size")
print(fig_power)

ggsave(here("output", "figures", "fig_power_curve.png"), fig_power,
       width = 6, height = 4, dpi = 300)

# Interactions are the hardest effects to detect: an interaction of the same
# size as a main effect needs a larger sample to detect it, and interactions
# are often smaller than main effects. A study sized for the main effects is
# usually underpowered for the interaction.
# (Avoid "observed power": computing power from the effect you just observed
# adds no information beyond the p-value.)


# ---- 13. Type I and Type II errors in a factorial design ---------------------
#                       | H0 is actually TRUE  | H0 is actually FALSE
#   --------------------+----------------------+----------------------------
#   Reject H0           | TYPE I error (alpha) | Correct (power = 1 - beta)
#   Fail to reject H0   | Correct (1 - alpha)  | TYPE II error (beta)
#
# A 2 x 2 ANOVA runs THREE F tests at once, each with its own Type I error rate of
# alpha. So if NOTHING is going on, the chance that at least one of the three is
# "significant" is about 1 - (1 - .05)^3 = 14%. We check this by simulation, and
# also estimate the power of each test if the TRUE pattern is the one we used to
# generate the data. In each of 10,000 simulated experiments we compute the three F
# statistics directly from the 4 cell means and the pooled variance (for a
# balanced 2 x 2 design this gives exactly the ANOVA F tests).

simulate_2x2 <- function(means, sd, n, n_sims = 10000) {
  # means: cell means in this order: Standard/Control, Standard/Stress,
  #        Enriched/Control, Enriched/Stress
  df2 <- 4 * (n - 1)
  results <- map(seq_len(n_sims), function(i) {
    x   <- matrix(rnorm(4 * n, mean = rep(means, each = n), sd = sd), nrow = n)
    m   <- colMeans(x)
    mse <- mean(apply(x, 2, var))
    # Each effect is a contrast of the cell means; SS = n * contrast^2 / 4
    ss_housing <- n * (m[3] + m[4] - m[1] - m[2])^2 / 4
    ss_stress  <- n * (m[2] + m[4] - m[1] - m[3])^2 / 4
    ss_inter   <- n * (m[1] - m[2] - m[3] + m[4])^2 / 4
    c(housing     = pf(ss_housing / mse, 1, df2, lower.tail = FALSE),
      stress      = pf(ss_stress  / mse, 1, df2, lower.tail = FALSE),
      interaction = pf(ss_inter   / mse, 1, df2, lower.tail = FALSE))
  })
  as_tibble(do.call(rbind, results))   # one row per simulated experiment
}

set.seed(123)
sd_within <- sqrt(anova_tbl$meansq[anova_tbl$term == "Residuals"])  # our pooled SD

sim_h0   <- simulate_2x2(rep(100, 4),           sd_within, n_cell)   # nothing is going on
sim_true <- simulate_2x2(c(100, 70, 108, 102), sd_within, n_cell)    # the true pattern

rates <- function(sim) {
  summarize(sim,
            any_effect = mean(housing < alpha | stress < alpha | interaction < alpha),
            across(c(housing, stress, interaction), ~ mean(.x < alpha))) |>
    relocate(any_effect, .after = last_col())
}
sim_summary <- bind_rows(`H0 true (no effects at all)` = rates(sim_h0),
                         `True pattern (interaction present)` = rates(sim_true),
                         .id = "world")
print(sim_summary)
# Row 1: each test rejects about 5% of the time (Type I error), but at least one
# of the three does so about 14% of the time.
# Row 2: power of each test (the chance of detecting the real effect).

# Theoretical power for the true pattern (compare with row 2). Using the
# population values: partial eta^2 = effect variance / (effect variance + SD^2),
# with the same SD that the simulation used.
effect_var <- c(housing = 100, stress = 81, interaction = 36)   # from the true cell means
true_pe    <- effect_var / (effect_var + sd_within^2)
power_2x2(true_pe / (1 - true_pe), n_cell, alpha)

fig_error <- sim_summary |>
  filter(world == "H0 true (no effects at all)") |>
  pivot_longer(-world, names_to = "test", values_to = "rate") |>
  mutate(test = factor(test,
           levels = c("housing", "stress", "interaction", "any_effect"),
           labels = c("Housing", "Stress", "Interaction", "At least one\nof the three"))) |>
  ggplot(aes(test, rate, fill = test)) +
  geom_col(width = 0.6) +
  geom_hline(yintercept = alpha, linetype = "dashed") +
  scale_fill_manual(values = unname(pal[c("blue", "orange", "green", "vermillion")])) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 0.2)) +
  labs(x = NULL, y = "Chance of a false positive") +
  theme(legend.position = "none")
print(fig_error)

ggsave(here("output", "figures", "fig_type1_error.png"), fig_error,
       width = 6, height = 4, dpi = 300)


# ---- 14. Robustness checks ---------------------------------------------------
# (a) Welch t-tests for the stress effect within each housing condition. These do
#     not assume equal variances. (t.test compares Control - Stress, so the sign
#     is flipped to match "Stress - Control".)
welch_simple <- bdnf |>
  group_by(housing) |>
  reframe(tidy(t.test(bdnf_pgmg ~ stress, var.equal = FALSE))) |>
  mutate(effect = -estimate, p_holm = p.adjust(p.value, method = "holm"))
welch_simple

# (b) Wilcoxon rank-sum tests within each housing condition (no normality needed)
bdnf |>
  group_by(housing) |>
  reframe(tidy(wilcox.test(bdnf_pgmg ~ stress, exact = FALSE))) |>
  mutate(p_holm = p.adjust(p.value, method = "holm"))

# (c) A rank-based check on the interaction: run the same ANOVA on the RANKS of
#     the data. (This is a rough check, not a replacement for the main analysis.)
summary(aov(rank(bdnf_pgmg) ~ housing * stress, data = bdnf))


# ---- 15. Write the result sentences ------------------------------------------
# Build the APA-style sentences from the saved objects so the numbers cannot be
# mistyped.
fmt_F <- function(row, pe) {
  paste0("F(", row$df, ", ", df_error, ") = ", fmt_num(row$statistic), ", ",
         fmt_p(row$p.value), ", partial eta-squared = ", fmt_apa(pe))
}
pe <- set_names(eta2$Eta2_partial, eta2$Parameter)
cat("Housing:     ", fmt_F(a_housing, pe["housing"]), "\n")
cat("Stress:      ", fmt_F(a_stress,  pe["stress"]),  "\n")
cat("Interaction: ", fmt_F(a_inter,   pe["housing:stress"]), "\n")


# ---- 16. Record the computing environment -----------------------------------
sessionInfo()
