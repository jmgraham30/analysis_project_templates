# =============================================================================
# 02_analysis.R
# Project:  Correlation and regression (running and adult neurogenesis)
# Purpose:  Complete, step-by-step analysis script. Run it from top to bottom
#           (or line by line with Ctrl/Cmd + Enter) to reproduce every number,
#           table and figure in the lab notebook.
#
# RESEARCH QUESTION
#   In 60 mice, do animals that run more on a wheel have more new neurons
#   (DCX+ cells) in the dentate gyrus? How strong is the relationship, how well
#   can we predict cell counts from running distance, and does the relationship
#   hold after accounting for age?
#
# HYPOTHESES (two-sided)
#   H0: the population correlation between running distance and DCX+ cell count
#       is zero (rho = 0); equivalently the regression slope is zero.
#   H1: rho is not zero (the slope is not zero).
#
# BEFORE YOU RUN: open correlation_regression.Rproj so the working directory is
#                 the project folder.
# INPUT : Data/neurogenesis_running.csv   (raw data, never edited by hand)
# OUTPUT: figures saved in output/figures/   (re-created each time you run this)
# =============================================================================


# ---- 1. Set up ---------------------------------------------------------------
source("R/00_setup.R")   # loads packages, plot theme, and helper functions

dir.create(here("output", "figures"), recursive = TRUE, showWarnings = FALSE)


# ---- 2. Decisions made BEFORE looking at the results -------------------------
alpha      <- 0.05    # Type I error rate we are willing to accept
conf_level <- 0.95    # confidence level for intervals
power_goal <- 0.80    # conventional target for statistical power
# Rules: (1) X (predictor) = wheel_km, Y (outcome) = dcx_cells;
#        (2) look at the scatterplot before computing anything;
#        (3) Pearson r is our primary statistic; Spearman is a check;
#        (4) check the residuals before trusting the regression.


# ---- 3. Import and inspect the data -------------------------------------------
neuro <- read_csv(here("Data", "neurogenesis_running.csv"), show_col_types = FALSE) |>
  mutate(sex = factor(sex, levels = c("F", "M")))

glimpse(neuro)
stopifnot(!anyNA(neuro), !anyDuplicated(neuro$mouse_id))   # no missing values or repeated IDs

neuro |>
  summarize(across(c(age_weeks, wheel_km, dcx_cells),
                   list(mean = mean, sd = sd, min = min, max = max)))


# ---- 4. Explore: LOOK at the relationship first ------------------------------
# A correlation is only a good summary if the relationship is roughly a straight
# line, so always plot first.
ggplot(neuro, aes(wheel_km, dcx_cells)) +
  geom_point(size = 2.2, alpha = 0.8, color = pal[["blue"]]) +
  labs(x = "Running distance (km per night)", y = "DCX+ cells in dentate gyrus")

# Distributions of the two variables (are they roughly symmetric? any outliers?)
p_x <- ggplot(neuro, aes(wheel_km)) + geom_histogram(binwidth = 0.75, fill = "gray80", color = "white") +
  labs(x = "Running distance (km)", y = "Mice")
p_y <- ggplot(neuro, aes(dcx_cells)) + geom_histogram(binwidth = 20, fill = "gray80", color = "white") +
  labs(x = "DCX+ cells", y = "Mice")
p_x + p_y

# All pairwise correlations among the numeric variables
neuro |> select(age_weeks, wheel_km, dcx_cells) |> cor() |> round(2)


# ---- 5. Correlation ----------------------------------------------------------
n_mice <- nrow(neuro)

pearson  <- cor_row(neuro$wheel_km, neuro$dcx_cells, "pearson")
spearman <- cor_row(neuro$wheel_km, neuro$dcx_cells, "spearman")
kendall  <- cor_row(neuro$wheel_km, neuro$dcx_cells, "kendall")
bind_rows(pearson, spearman, kendall)

# r-squared: the proportion of variability in Y "shared" with X.
pearson$r^2

# Cohen's rough benchmarks for r: .10 small, .30 medium, .50 large.
interpret_r(pearson$r, rules = "cohen1988")


# ---- 6. Power for a correlation -----------------------------------------------
# (a) Power for a "medium" true correlation, r = .30, with our 60 mice.
# (We do NOT compute "power" from the r we observed: that number is just a
# re-statement of the p-value and tells us nothing new.)
pwr.r.test(n = n_mice, r = 0.30, sig.level = alpha)$power

# (b) Sample size needed for 80% power at different true correlations
tibble(r = c(0.1, 0.3, 0.5, 0.7)) |>
  mutate(n_needed = map_dbl(r, ~ ceiling(pwr.r.test(r = .x, sig.level = alpha,
                                                    power = power_goal)$n)))

# (c) Sensitivity: the smallest r we could detect with 80% power at n = 60
pwr.r.test(n = n_mice, sig.level = alpha, power = power_goal)$r


# ---- 7. Simple linear regression ---------------------------------------------
# Model: dcx_cells = intercept + slope * wheel_km + error
fit1 <- lm(dcx_cells ~ wheel_km, data = neuro)

tidy(fit1, conf.int = TRUE, conf.level = conf_level)   # intercept and slope with CIs
glance(fit1)                                           # R-squared, residual SE, F test

# The slope test is identical to the correlation test: t^2 = F, same p-value.
tidy(fit1) |> filter(term == "wheel_km") |> pull(statistic)
pearson$r * sqrt((n_mice - 2) / (1 - pearson$r^2))     # same t from r

# The intercept is the predicted count at 0 km, which is outside the range of
# the data (the least active mouse ran 1 km). To make the intercept meaningful,
# CENTER the predictor so that 0 means "an average runner":
fit1c <- lm(dcx_cells ~ I(wheel_km - mean(wheel_km)), data = neuro)
tidy(fit1c)

# Predictions with two kinds of uncertainty
new_x <- tibble(wheel_km = c(2, 4, 6))
bind_cols(new_x,
          predict(fit1, new_x, interval = "confidence") |> as_tibble() |> rename_with(~ paste0("ci_", .x)),
          predict(fit1, new_x, interval = "prediction") |> as_tibble() |> select(-fit) |> rename_with(~ paste0("pi_", .x)))
# confidence interval: where the AVERAGE count for mice running that far lies
# prediction interval: where a SINGLE new mouse's count is likely to fall (wider)

# Extrapolation: the model will happily predict for 15 km, but no mouse ran that far.
predict(fit1, tibble(wheel_km = 15))
range(neuro$wheel_km)


# ---- 8. Check the regression assumptions (residual diagnostics) --------------
aug <- augment(fit1, data = neuro)    # adds .fitted, .resid, .hat, .cooksd, .std.resid

# (a) Linearity and equal spread: residuals vs fitted should be a flat, even band
ggplot(aug, aes(.fitted, .resid)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(color = pal[["blue"]], alpha = 0.8) +
  geom_smooth(method = "loess", formula = y ~ x, se = FALSE, color = pal[["vermillion"]], linewidth = 0.8)

# (b) Normality of residuals: points should follow the line
ggplot(aug, aes(sample = .std.resid)) + stat_qq() + stat_qq_line()
shapiro.test(aug$.resid)   # informational only; look at the plot first

# (c) Influential points: Cook's distance > 4/n is a common flag
n_flag <- 4 / n_mice
aug |> filter(.cooksd > n_flag) |> select(mouse_id, wheel_km, dcx_cells, .std.resid, .hat, .cooksd)
max(aug$.cooksd)

# (d) Does the conclusion survive the non-normal residuals? Check with a
# bootstrap: resample mice with replacement 2000 times and refit the line.
set.seed(2024)
boot_slopes <- replicate(2000, {
  b <- neuro[sample(n_mice, replace = TRUE), ]
  coef(lm(dcx_cells ~ wheel_km, data = b))[["wheel_km"]]
})
quantile(boot_slopes, c(0.025, 0.975))   # compare with the t-based CI above

# (e) Independence: each row must be a different, independent mouse. This is a
# design question; no test can fix it. (Here it holds: one row = one mouse.)


# ---- 9. Multiple regression: account for age ---------------------------------
# Neurogenesis declines with age. If older mice also ran less, age could confound
# the running effect. Add age as a second predictor.
cor(neuro$wheel_km, neuro$age_weeks)

fit2 <- lm(dcx_cells ~ wheel_km + age_weeks, data = neuro)
tidy(fit2, conf.int = TRUE, conf.level = conf_level)
glance(fit2)

# Did age add anything? Compare the nested models.
anova(fit1, fit2)

# Standardized coefficients put the predictors on a common scale (SD units)
standardize_parameters(fit2)

# Unique contribution of each predictor = how much R-squared drops if removed
r2_full <- summary(fit2)$r.squared
c(wheel = r2_full - summary(lm(dcx_cells ~ age_weeks, neuro))$r.squared,
  age   = r2_full - summary(fit1)$r.squared)

# Residuals of the two-predictor model (the omitted-age pattern should be gone)
shapiro.test(resid(fit2))
augment(fit2, data = neuro) |> summarize(max_cooks = max(.cooksd), max_abs_std_resid = max(abs(.std.resid)))

# Collinearity check (variance inflation factor, for two predictors = 1/(1-r^2))
1 / (1 - cor(neuro$wheel_km, neuro$age_weeks)^2)


# ---- 10. Pitfalls ------------------------------------------------------------
# (a) Anscombe's quartet: four data sets, the same r and regression line.
anscombe_long <- anscombe |>
  pivot_longer(everything(), names_to = c(".value", "set"), names_pattern = "(.)(.)")
anscombe_long |>
  group_by(set) |>
  summarize(mean_x = mean(x), mean_y = mean(y), r = cor(x, y),
            slope = coef(lm(y ~ x))[2])

# (b) Restricted range: look only at mice that ran between 3 and 5 km
restricted <- neuro |> filter(wheel_km >= 3, wheel_km <= 5)
cor_row(restricted$wheel_km, restricted$dcx_cells)   # r shrinks even though the biology is unchanged

# (c) One outlier can create or destroy a correlation (add one extreme mouse)
with_outlier <- bind_rows(neuro, tibble(mouse_id = "M61", sex = "M", age_weeks = 14,
                                        wheel_km = 12, dcx_cells = 60))
cor_row(with_outlier$wheel_km, with_outlier$dcx_cells)


# ---- 11. Simulations: what happens over many repeated experiments? ------------
set.seed(123)
n_sims <- 10000

# Draw n pairs (X, Y) with a chosen TRUE correlation rho and return the sample r's
sim_r <- function(n, rho, n_sims) {
  replicate(n_sims, {
    x <- rnorm(n); y <- rho * x + sqrt(1 - rho^2) * rnorm(n)
    cor(x, y)
  })
}
r_to_p <- function(r, n) 2 * pt(-abs(r * sqrt((n - 2) / (1 - r^2))), n - 2)

# (a) The sampling distribution of r when rho = 0.6 (the SD shrinks as n grows)
sizes <- c(10, 30, 60, 200)
sim_dist <- map_dfr(sizes, ~ tibble(n = .x, r = sim_r(.x, 0.6, n_sims)))
sim_dist |> group_by(n) |> summarize(mean_r = mean(r), sd_r = sd(r),
                                     lo = quantile(r, 0.025), hi = quantile(r, 0.975))

# (b) Type I error: when rho = 0, how often is p < .05? (should be about 5%)
r0 <- sim_r(n_mice, 0, n_sims)
mean(r_to_p(r0, n_mice) < alpha)

# (c) Power and Type II error: how often do we detect a REAL correlation?
power_grid <- expand_grid(n = c(10, 20, 30, 40, 60, 80, 100), rho = c(0.3, 0.6)) |>
  mutate(power_sim = map2_dbl(n, rho, ~ mean(r_to_p(sim_r(.x, .y, 4000), .x) < alpha)),
         power_pwr = map2_dbl(n, rho, ~ pwr.r.test(n = .x, r = .y, sig.level = alpha)$power),
         type2     = 1 - power_sim)
power_grid

# (d) Selection and the "winner's curse": among only the SIGNIFICANT results from
# small studies (n = 20, true rho = 0.3), the average observed r is inflated.
r_small <- sim_r(20, 0.3, n_sims)
sig <- r_to_p(r_small, 20) < alpha
c(mean_all = mean(r_small), mean_significant_only = mean(r_small[sig]), power = mean(sig))
