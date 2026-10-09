# =============================================================================
# 00_setup.R
# Project:  Measures of central tendency and variability (neuron firing rates)
# Purpose:  Load the packages, plotting theme, and small helper functions that
#           every other script and the Quarto lab notebook rely on.
#
# HOW TO USE: You do not usually run this file by itself. The analysis script
#             (02_analysis.R) and the lab notebook (measures_tendency_variability.qmd) both
#             call source("R/00_setup.R") at the top. Keeping this "boilerplate"
#             in one place means you fix or change it once, not in every file.
#
# NOTE: Paths are written relative to the PROJECT FOLDER (the folder that holds
#       the .Rproj file). Always open the .Rproj file first so RStudio sets that
#       folder as the working directory.
# =============================================================================


# ---- 1. Packages ------------------------------------------------------------
# Packages are add-on toolboxes. Each one is used for a specific job:
#   tidyverse  - data import (readr), tidying (tidyr), wrangling (dplyr),
#                text (stringr), dates (lubridate), plotting (ggplot2)
#   here       - builds file paths that work on any computer
#   knitr, kableExtra - nicely formatted tables
#   patchwork  - combine several plots into one figure

required_packages <- c("tidyverse", "here", "knitr", "kableExtra", "patchwork")

# Install any package that is not already on this computer. This needs an
# internet connection and only happens the first time. (If you are offline, or
# the install fails, install the packages above from the Packages pane first.)
# These projects were written for R 4.1 or newer and a current tidyverse
# (dplyr 1.1 or newer). If a function is "not found", update your packages.
missing_packages <- setdiff(required_packages, rownames(installed.packages()))
if (length(missing_packages) > 0) {
  message("Installing missing packages: ", paste(missing_packages, collapse = ", "))
  install.packages(missing_packages, repos = "https://cloud.r-project.org")
}

# Load the packages (quietly, so the notebook is not cluttered with startup text)
suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(knitr)
  library(kableExtra)
  library(patchwork)
})


# ---- 2. Color palette and plot theme --------------------------------------
# Okabe-Ito colors are distinguishable by people with common forms of color
# blindness and still work when printed in grayscale.
pal <- c(blue = "#0072B2", orange = "#E69F00", green = "#009E73",
         vermillion = "#D55E00", sky = "#56B4E9", pink = "#CC79A7",
         yellow = "#F0E442", gray = "#666666")

# A clean, journal-style theme: white background, no grid lines, readable text.
theme_pub <- function(base_size = 12) {
  theme_classic(base_size = base_size) +
    theme(
      axis.text        = element_text(color = "black"),
      axis.title       = element_text(face = "plain"),
      plot.title       = element_text(face = "bold", size = base_size + 1),
      legend.position  = "bottom",
      strip.background = element_blank(),
      strip.text       = element_text(face = "bold")
    )
}
theme_set(theme_pub())   # make it the default for every ggplot in the project


# ---- 3. Helper functions for reporting results (APA style) ------------------

# Format a p-value the APA way: no leading zero, and "< .001" for tiny values.
#   fmt_p(0.0234)           ->  "p = .023"
#   fmt_p(0.00002)          ->  "p < .001"
#   fmt_p(0.0234, md = TRUE) -> "*p* = .023"   (italic p for the Quarto notebook)
fmt_p <- function(p, digits = 3, md = FALSE) {
  lab <- if (md) "*p*" else "p"
  if (is.na(p)) return(paste(lab, "= NA"))
  if (p < 0.001) return(paste(lab, "< .001"))
  paste0(lab, " = ", sub("^0", "", formatC(p, format = "f", digits = digits)))
}

# Format a number to a fixed number of decimals (keeps trailing zeros, e.g. 0.50).
fmt_num <- function(x, digits = 2) formatC(x, format = "f", digits = digits)

# Same, but drops the leading zero. APA style does this for statistics that
# cannot exceed 1 (correlations r, proportions, p-values, eta-squared).
fmt_apa <- function(x, digits = 2) sub("^(-?)0\\.", "\\1.", fmt_num(x, digits))

# Vectorized p-value formatter for table columns: "< .001" or ".023" (no leading
# zero). Unlike fmt_p(), it works on a whole column and omits the "p =" label.
fmt_p_col <- function(p, digits = 3) {
  if_else(p < 0.001, "< .001", sub("^0", "", formatC(p, format = "f", digits = digits)))
}

# A consistent, clean table style for every table in the notebook (journal-like:
# horizontal rules only, no shading). Extra arguments go to knitr::kbl(), e.g.
# col.names = c(...), digits = 2, align = "lrr".
pub_table <- function(x, ...) {
  x |>
    kbl(...) |>
    kable_classic(full_width = FALSE, html_font = "inherit", position = "center")
}


# ---- 4. Helper functions for describing data -----------------------------------
# Base R has mean(), median(), sd(), var(), IQR(), mad(), and range(), but no
# built-in functions for the standard error, coefficient of variation, mode,
# skewness, or kurtosis. We define small, readable versions here so you can see
# exactly what each statistic is. All of them ignore missing values (NA).

# Standard error of the mean: how precisely the sample mean estimates the
# population mean. SE = SD / sqrt(n).
se <- function(x) {
  x <- x[!is.na(x)]
  sd(x) / sqrt(length(x))
}

# Coefficient of variation: the SD as a fraction of the mean (unit-free). Only
# meaningful for ratio-scale data with a true zero and a positive mean.
cv <- function(x) sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE)

# Mode of a CATEGORICAL variable: the most common category (returns all of them
# if there is a tie).
mode_cat <- function(x) {
  tab <- table(x)
  names(tab)[tab == max(tab)]
}

# Skewness: 0 for symmetric data, positive when the right tail is longer, negative
# when the left tail is longer. (This is the moment coefficient, g1.)
skew <- function(x) {
  x <- x[!is.na(x)]
  m <- mean(x)
  mean((x - m)^3) / mean((x - m)^2)^1.5
}

# Excess kurtosis: 0 for a normal distribution; positive means heavier tails
# (more extreme values) than a normal distribution.
kurt <- function(x) {
  x <- x[!is.na(x)]
  m <- mean(x)
  mean((x - m)^4) / mean((x - m)^2)^2 - 3
}

# The mean with a t-based confidence interval (default 95%), as a one-row table.
mean_ci <- function(x, level = 0.95) {
  x <- x[!is.na(x)]
  n <- length(x)
  half <- qt(1 - (1 - level) / 2, df = n - 1) * sd(x) / sqrt(n)
  tibble(n = n, mean = mean(x), lower = mean(x) - half, upper = mean(x) + half)
}
