# =============================================================================
# 02_analysis.R
# Project:  Chi-square tests (genotype ratios and seizure outcomes)
# Purpose:  Complete, step-by-step analysis script. Run it from top to bottom
#           (or line by line with Ctrl/Cmd + Enter) to reproduce every number,
#           table and figure in the lab notebook.
#
# TWO RESEARCH QUESTIONS (one chi-square test each)
#   Part A (goodness of fit): Do the genotypes of pups from heterozygous
#     crosses follow the Mendelian 1 : 2 : 1 ratio, or is the knockout genotype
#     under-represented?
#   Part B (test of independence): Is having a seizure after a convulsant drug
#     related to pretreatment (anticonvulsant drug vs. vehicle)?
#
# BEFORE YOU RUN: open chisquare_test.Rproj so the working directory is the
#                 project folder.
# INPUT : Data/genotype_counts.csv, Data/seizure_treatment.csv  (raw, never edited)
# OUTPUT: figures saved in output/figures/   (re-created each time you run this)
# =============================================================================


# ---- 1. Set up ---------------------------------------------------------------
source("R/00_setup.R")   # loads packages, plot theme, and helper functions

dir.create(here("output", "figures"), recursive = TRUE, showWarnings = FALSE)


# ---- 2. Decisions made BEFORE looking at the results -------------------------
alpha        <- 0.05   # significance level = accepted Type I error rate
target_power <- 0.80   # conventional minimum power
# Part A: chi-square GOODNESS-OF-FIT test against the 1 : 2 : 1 Mendelian ratio.
# Part B: chi-square TEST OF INDEPENDENCE on the 2 x 2 table (treatment x
#         seizure), WITHOUT the continuity correction (see section 12 for why),
#         with Fisher's exact test as a check. Both tests are chosen to be valid
#         only if every EXPECTED count is at least 5 (checked below).
mendel <- c(WT = 0.25, Het = 0.50, KO = 0.25)


# =============================================================================
# PART A: GOODNESS OF FIT  (one categorical variable, compared with a theory)
# =============================================================================

# ---- 3. Import and prepare the genotype data ---------------------------------
geno_raw <- read_csv(here("Data", "genotype_counts.csv"), show_col_types = FALSE)

glimpse(geno_raw)
sum(is.na(geno_raw))                # missing values?
anyDuplicated(geno_raw$pup_id)      # duplicated IDs?

# The raw file has one row per pup. For a chi-square test we need the COUNT in
# each category. The factor levels fix the order (WT, Het, KO).
geno <- geno_raw |>
  mutate(genotype = factor(genotype, levels = c("WT", "Het", "KO")))

obs_geno <- count(geno, genotype, .drop = FALSE)   # .drop = FALSE keeps zero counts
obs_geno
N_geno <- nrow(geno)


# ---- 4. The goodness-of-fit test --------------------------------------------
# H0: the genotype probabilities are 25% / 50% / 25%.
# H1: at least one probability is different.
#
# chi-square = sum over categories of (Observed - Expected)^2 / Expected
# with df = (number of categories - 1) = 2.

gof <- chisq.test(x = obs_geno$n, p = mendel)   # p must sum to 1
gof
gof_tbl <- tidy(gof)

# Expected counts and residuals. The test is only trustworthy if EVERY expected
# count is at least 5 (a common rule of thumb).
gof_detail <- obs_geno |>
  mutate(
    expected  = N_geno * mendel,
    residual  = n - expected,
    std_resid = as.numeric(gof$stdres)   # (O - E) / sqrt(E (1 - p)): about +-2 is notable
  )
gof_detail
min(gof_detail$expected)

# Effect size: Cohen's w = sqrt(chi-square / N). Benchmarks: .10 small, .30
# medium, .50 large. (For a 2 x 2 table, w equals phi.)
w_geno <- cohens_w(obs_geno$n, p = mendel, ci = 0.95, alternative = "two.sided")
w_geno

# Plot: observed counts with the expected counts from the 1 : 2 : 1 ratio
fig_geno <- gof_detail |>
  ggplot(aes(genotype)) +
  geom_col(aes(y = n, fill = genotype), width = 0.6) +
  geom_point(aes(y = expected, shape = "Expected (1 : 2 : 1)"),
             size = 4, color = "black") +
  scale_fill_manual(values = unname(pal[c("sky", "blue", "vermillion")]), guide = "none") +
  scale_shape_manual(values = 18, name = NULL) +
  labs(x = "Genotype", y = "Number of pups")
print(fig_geno)


# =============================================================================
# PART B: TEST OF INDEPENDENCE  (two categorical variables)
# =============================================================================

# ---- 5. Import and prepare the seizure data ---------------------------------
seiz_raw <- read_csv(here("Data", "seizure_treatment.csv"), show_col_types = FALSE)

glimpse(seiz_raw)
sum(is.na(seiz_raw))
anyDuplicated(seiz_raw$mouse_id)

# Factor levels: put the control first and "Yes" (the event of interest) first
# in the outcome, so the table reads naturally.
seiz <- seiz_raw |>
  mutate(
    treatment = factor(treatment, levels = c("Vehicle", "Drug")),
    seizure   = factor(seizure,   levels = c("Yes", "No")),
    sex       = factor(sex, levels = c("F", "M"), labels = c("Female", "Male"))
  )
N <- nrow(seiz)

# The CONTINGENCY TABLE: counts for every combination of the two variables
tab <- table(Treatment = seiz$treatment, Seizure = seiz$seizure)
tab
addmargins(tab)

# Percent with a seizure in each group (the "row proportions")
prop_tbl <- seiz |>
  count(treatment, seizure) |>
  group_by(treatment) |>
  mutate(percent = 100 * n / sum(n)) |>
  ungroup()
prop_tbl

count(seiz, treatment, sex)   # sex balanced across groups?


# ---- 6. The test of independence --------------------------------------------
# H0: seizure and treatment are INDEPENDENT (the probability of a seizure is the
#     same after vehicle and after the drug).
# H1: they are not independent (the probability differs).
#
# For each cell: Expected = (row total x column total) / grand total, that is,
# the count we would see if treatment and outcome were unrelated.
# chi-square = sum of (Observed - Expected)^2 / Expected,
# df = (rows - 1) x (columns - 1) = 1 for a 2 x 2 table.

pearson <- chisq.test(tab, correct = FALSE)   # the standard Pearson test
pearson
pearson$expected                              # CHECK: all expected counts >= 5?
round(pearson$stdres, 2)                      # adjusted standardized residuals

# Cell-by-cell contributions: which cells depart most from independence?
(tab - pearson$expected)^2 / pearson$expected

# R's default for 2 x 2 tables applies Yates' continuity correction; compare:
yates <- chisq.test(tab, correct = TRUE)
yates
# Fisher's exact test does not use the chi-square approximation at all:
fisher <- fisher.test(tab)
fisher


# ---- 7. Effect sizes ---------------------------------------------------------
# (a) phi (= Cramer's V for a 2 x 2 table): the strength of association, 0 to 1.
#     Benchmarks: .10 small, .30 medium, .50 large.
phi_es <- phi(tab, ci = 0.95, alternative = "two.sided", adjust = FALSE)
phi_es

# (b) Risk difference, risk ratio, and odds ratio.
# Rows are Vehicle, Drug; columns are Yes, No. We compare the DRUG group with
# VEHICLE (the reference), so values below 1 mean the drug lowers the risk.
risk_vehicle <- tab["Vehicle", "Yes"] / sum(tab["Vehicle", ])
risk_drug    <- tab["Drug", "Yes"]    / sum(tab["Drug", ])
risk_vehicle; risk_drug

rd <- prop.test(c(tab["Drug", "Yes"], tab["Vehicle", "Yes"]),
                c(sum(tab["Drug", ]), sum(tab["Vehicle", ])), correct = FALSE) |> tidy()
rd   # estimate1 - estimate2 = risk difference (Drug - Vehicle), with 95% CI

# Risk ratio and odds ratio, with 95% confidence intervals computed on the log
# scale (the standard "Wald" method).
drug_yes <- tab["Drug", "Yes"];    drug_no <- tab["Drug", "No"]
veh_yes  <- tab["Vehicle", "Yes"]; veh_no  <- tab["Vehicle", "No"]
z <- qnorm(0.975)

rr_hat <- (drug_yes / (drug_yes + drug_no)) / (veh_yes / (veh_yes + veh_no))  # risk: Drug / Vehicle
rr_se  <- sqrt(1/drug_yes - 1/(drug_yes + drug_no) + 1/veh_yes - 1/(veh_yes + veh_no))  # SE of log(RR)
rr_ci  <- exp(log(rr_hat) + c(-1, 1) * z * rr_se)
c(RR = rr_hat, lower = rr_ci[1], upper = rr_ci[2])

or_hat <- (drug_yes * veh_no) / (drug_no * veh_yes)       # odds: Drug / Vehicle
or_se  <- sqrt(1/drug_yes + 1/drug_no + 1/veh_yes + 1/veh_no)   # SE of log(OR)
or_ci  <- exp(log(or_hat) + c(-1, 1) * z * or_se)
c(OR = or_hat, lower = or_ci[1], upper = or_ci[2])
# Both are below 1: a seizure is LESS likely (and has lower odds) after the drug.
# (Reminder: an odds ratio is not a risk ratio. They are similar only when the
# event is rare. Here seizures are common, so the OR is further from 1 than the RR.)

# (c) Cohen's h compares two proportions on a transformed scale (.2 small, .5 medium,
#     .8 large). It is used for power analysis of two-proportion tests.
h_es <- cohens_h(tab, ci = 0.95, alternative = "two.sided")
h_es


# ---- 8. Publication-quality figure ------------------------------------------
# (A) Percent of mice with a seizure in each group, with 95% confidence
#     intervals (computed from the binomial distribution).
# (B) The risk difference (Drug minus Vehicle) with its 95% CI.

group_ci <- seiz |>
  group_by(treatment) |>
  summarize(n = n(), n_yes = sum(seizure == "Yes"), .groups = "drop") |>
  mutate(ci = map2(n_yes, n, ~ tidy(prop.test(.x, .y, correct = FALSE)))) |>
  unnest(ci) |>
  select(treatment, n, n_yes, estimate, conf.low, conf.high)
group_ci

fig_a <- ggplot(group_ci, aes(treatment, 100 * estimate, fill = treatment)) +
  geom_col(width = 0.55) +
  geom_errorbar(aes(ymin = 100 * conf.low, ymax = 100 * conf.high), width = 0.15,
                linewidth = 0.8) +
  geom_text(aes(label = paste0(n_yes, "/", n), y = 3), color = "white",
            fontface = "bold", size = 4) +
  scale_fill_manual(values = unname(pal[c("gray", "orange")]), guide = "none") +
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
  labs(x = "Pretreatment", y = "Mice with a seizure (%)")

fig_b <- tibble(estimate = rd$estimate1 - rd$estimate2,
                lower = rd$conf.low, upper = rd$conf.high) |>
  ggplot(aes(100 * estimate, "Drug minus Vehicle")) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray30") +
  geom_pointrange(aes(xmin = 100 * lower, xmax = 100 * upper), color = pal["green"],
                  size = 0.6, linewidth = 0.9) +
  labs(x = "Difference in seizure rate\n(percentage points, with 95% CI)", y = NULL)

fig_main <- fig_a + fig_b + plot_layout(widths = c(1, 1.2)) +
  plot_annotation(tag_levels = "A")
print(fig_main)

ggsave(here("output", "figures", "fig_seizure_by_treatment.png"), fig_main,
       width = 8, height = 4, dpi = 300)
ggsave(here("output", "figures", "fig_genotype_counts.png"), fig_geno,
       width = 5.5, height = 4, dpi = 300)


# ---- 9. Power analysis (Part B) ----------------------------------------------
# POWER = probability of rejecting H0 when it is really false (= 1 - beta).
# For chi-square tests, power depends on the effect size w, the total N, the
# degrees of freedom, and alpha. pwr.chisq.test() does the calculation.

# (a) A priori: total N needed for 80% power to detect a medium effect (w = 0.30)
plan <- pwr.chisq.test(w = 0.30, df = 1, sig.level = alpha, power = target_power)
plan          # round N UP

# (b) Sensitivity: with our N, what is the smallest w we could detect with 80% power?
sens <- pwr.chisq.test(N = N, df = 1, sig.level = alpha, power = target_power)
sens

# (c) Power to detect a medium (w = 0.30) and small (w = 0.10) effect with our N
pwr.chisq.test(w = 0.30, N = N, df = 1, sig.level = alpha)$power
pwr.chisq.test(w = 0.10, N = N, df = 1, sig.level = alpha)$power

# (d) Power curve
power_curve <- expand_grid(N = seq(10, 400, by = 5), w = c(0.10, 0.30, 0.50)) |>
  mutate(power   = map2_dbl(N, w, ~ pwr.chisq.test(w = .y, N = .x, df = 1,
                                                   sig.level = alpha)$power),
         w_label = factor(w, labels = c("Small (w = 0.10)", "Medium (w = 0.30)",
                                        "Large (w = 0.50)")))

fig_power <- ggplot(power_curve, aes(N, power, color = w_label, linetype = w_label)) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = target_power, linetype = "dashed") +
  geom_vline(xintercept = N, linetype = "dotted") +
  scale_color_manual(values = unname(pal[c("sky", "blue", "vermillion")])) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(x = "Total number of mice (N)", y = "Power",
       color = "True effect size", linetype = "True effect size")
print(fig_power)

ggsave(here("output", "figures", "fig_power_curve.png"), fig_power,
       width = 6, height = 4, dpi = 300)

# Power for Part A (goodness of fit, df = 2): our N of pups and w = 0.30
pwr.chisq.test(w = 0.30, N = N_geno, df = 2, sig.level = alpha)$power
pwr.chisq.test(N = N_geno, df = 2, sig.level = alpha, power = target_power)$w

# (Avoid "observed power": computing power from the effect you just observed adds
# no information beyond the p-value.)


# ---- 10. Type I and Type II errors by simulation -----------------------------
#                       | H0 is actually TRUE  | H0 is actually FALSE
#   --------------------+----------------------+----------------------------
#   Reject H0           | TYPE I error (alpha) | Correct (power = 1 - beta)
#   Fail to reject H0   | Correct (1 - alpha)  | TYPE II error (beta)
#
# We simulate 10,000 experiments for two groups (n mice each) in worlds where we
# KNOW the true seizure probabilities, and run three tests on each 2 x 2 table:
# Pearson chi-square, chi-square with Yates' correction, and Fisher's exact test.
# A table is determined by the two counts of seizures, so there are only (n+1)^2
# different tables; we compute each test once and look up the p-values.

simulate_2x2 <- function(p1, p2, n, n_sims = 10000) {
  x1 <- rbinom(n_sims, n, p1)   # seizures in group 1
  x2 <- rbinom(n_sims, n, p2)   # seizures in group 2

  # Compute all three p-values for each distinct table
  tables <- distinct(tibble(x1, x2)) |>
    mutate(
      p_pearson = map2_dbl(x1, x2, \(a, b) {
        m <- matrix(c(a, n - a, b, n - b), 2, byrow = TRUE)
        if (any(rowSums(m) == 0) || any(colSums(m) == 0)) return(1)   # nothing to test
        suppressWarnings(chisq.test(m, correct = FALSE)$p.value)
      }),
      p_yates = map2_dbl(x1, x2, \(a, b) {
        m <- matrix(c(a, n - a, b, n - b), 2, byrow = TRUE)
        if (any(rowSums(m) == 0) || any(colSums(m) == 0)) return(1)
        suppressWarnings(chisq.test(m, correct = TRUE)$p.value)
      }),
      p_fisher = map2_dbl(x1, x2, \(a, b) {
        fisher.test(matrix(c(a, n - a, b, n - b), 2, byrow = TRUE))$p.value
      })
    )
  tibble(x1, x2) |> left_join(tables, by = c("x1", "x2")) |> select(starts_with("p_"))
}

set.seed(123)
p_null <- mean(seiz$seizure == "Yes")     # overall seizure rate, used when H0 is true

sim_pilot_h0 <- simulate_2x2(p_null, p_null, n = 8)       # small study, H0 true
sim_main_h0  <- simulate_2x2(p_null, p_null, n = 50)      # our study size, H0 true
sim_main_h1  <- simulate_2x2(0.60, 0.36, n = 50)          # our study size, true effect
sim_pilot_h1 <- simulate_2x2(0.60, 0.36, n = 8)           # small study, true effect

rate <- function(sim) summarize(sim, across(everything(), ~ mean(.x < alpha)))
sim_summary <- bind_rows(
  `n = 8 per group, H0 true`    = rate(sim_pilot_h0),
  `n = 50 per group, H0 true`   = rate(sim_main_h0),
  `n = 8 per group, H1 true`    = rate(sim_pilot_h1),
  `n = 50 per group, H1 true`   = rate(sim_main_h1),
  .id = "scenario"
)
print(sim_summary)
# Rows 1-2 are Type I error rates (should be near .05); rows 3-4 are power.
# With tiny samples, Fisher's test (and Yates) are conservative: their Type I error
# is well BELOW .05, and they have less power.

# Theoretical power for the true effect at n = 50 per group (compare with row 4)
true_phi <- abs(0.60 - 0.36) / (2 * sqrt(0.48 * 0.52))
pwr.chisq.test(w = true_phi, N = 100, df = 1, sig.level = alpha)$power

fig_error <- sim_summary |>
  filter(str_detect(scenario, "H0 true")) |>
  pivot_longer(-scenario, names_to = "test", values_to = "rate") |>
  mutate(test = factor(test, levels = c("p_pearson", "p_yates", "p_fisher"),
                       labels = c("Pearson", "Yates-corrected", "Fisher's exact"))) |>
  ggplot(aes(test, rate, fill = test)) +
  geom_col(width = 0.6) +
  geom_hline(yintercept = alpha, linetype = "dashed") +
  facet_wrap(~ scenario) +
  scale_fill_manual(values = unname(pal[c("blue", "orange", "green")])) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 0.10)) +
  labs(x = NULL, y = "Type I error rate") +
  theme(legend.position = "none")
print(fig_error)

ggsave(here("output", "figures", "fig_type1_error.png"), fig_error,
       width = 7, height = 4, dpi = 300)


# ---- 11. Small samples: when expected counts are below 5 --------------------
# To see the problem, pretend we only had a small PILOT study: the first 8 mice in
# each group. (This is a demonstration; in a real study you would not throw
# data away.)
pilot <- seiz |>
  group_by(treatment) |>
  slice_head(n = 8) |>
  ungroup()

tab_pilot <- table(Treatment = pilot$treatment, Seizure = pilot$seizure)
tab_pilot

chisq_pilot <- suppressWarnings(chisq.test(tab_pilot, correct = FALSE))
chisq_pilot$expected             # some expected counts are below 5!
chisq_pilot$p.value              # NOT trustworthy here
fisher_pilot <- fisher.test(tab_pilot)
fisher_pilot$p.value             # exact; use this one
fisher_pilot$estimate            # conditional MLE of the odds ratio (can differ from the sample odds ratio)


# ---- 12. Robustness checks for the main 2 x 2 table --------------------------
# Compare the four ways of testing the same table. If they all agree, the
# conclusion does not depend on which one we choose. (Yates is generally considered
# too conservative today; we show it because R applies it by default.)

# Likelihood-ratio (G) test: another large-sample test, G = 2 * sum(O * log(O / E))
obs <- as.vector(tab)
exp_counts <- as.vector(pearson$expected)
G <- 2 * sum(obs * log(obs / exp_counts))
p_G <- pchisq(G, df = 1, lower.tail = FALSE)

robust_tbl <- tibble(
  Test = c("Pearson chi-square", "Yates-corrected chi-square",
           "Likelihood-ratio (G) test", "Fisher's exact test"),
  p    = c(pearson$p.value, yates$p.value, p_G, fisher$p.value)
)
robust_tbl

# Part A: the goodness-of-fit p-value is close to .05, so check it two other ways.
# (a) A Monte Carlo p-value: simulate 10,000 samples of the same size from the
#     1 : 2 : 1 ratio and see how often chi-square is at least as large as ours.
set.seed(123)
gof_mc <- chisq.test(obs_geno$n, p = mendel, simulate.p.value = TRUE, B = 10000)
gof_mc$p.value
# (b) The likelihood-ratio (G) test for the same counts.
G_geno <- 2 * sum(obs_geno$n * log(obs_geno$n / (N_geno * mendel)))
pchisq(G_geno, df = 2, lower.tail = FALSE)


# ---- 13. Write the result sentences ------------------------------------------
fmt_chisq <- function(test, N, es_label, es) {
  paste0("chi-square(", test$parameter, ", N = ", N, ") = ", fmt_num(test$statistic),
         ", ", fmt_p(test$p.value), ", ", es_label, " = ", fmt_apa(es))
}
cat("Part A:", fmt_chisq(gof, N_geno, "w", w_geno$Cohens_w), "\n")
cat("Part B:", fmt_chisq(pearson, N, "phi", phi_es$phi), "\n")


# ---- 14. Record the computing environment -----------------------------------
sessionInfo()
