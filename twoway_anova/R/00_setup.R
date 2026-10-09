# =============================================================================
# 00_setup.R
# Project:  Two-way ANOVA (housing, stress, and hippocampal BDNF)
# Purpose:  Load the packages, plotting theme, and small helper functions that
#           every other script and the Quarto lab notebook rely on.
#
# HOW TO USE: You do not usually run this file by itself. The analysis script
#             (02_analysis.R) and the lab notebook (twoway_anova.qmd) both
#             call source("R/00_setup.R") at the top. Keeping this "boilerplate"
#             in one place means you fix or change it once, not in every file.
#
# NOTE: Paths are written relative to the PROJECT FOLDER (the folder that holds
#       the .Rproj file). Always open the .Rproj file first so RStudio sets that
#       folder as the working directory.
# =============================================================================


# ---- 1. Packages ------------------------------------------------------------
# Packages are add-on toolboxes. Each one is used for a specific job:
#   tidyverse  - data import, wrangling (dplyr) and plotting (ggplot2)
#   here       - builds file paths that work on any computer
#   broom      - turns test results into tidy data frames (a tidymodels package)
#   rstatix    - tidy-friendly helpers for outlier and variance checks
#   emmeans    - estimated group means and post hoc (multiple) comparisons
#   effectsize - effect sizes (Cohen's d, Hedges' g) with confidence intervals
#   pwr        - power analysis for common tests
#   knitr, kableExtra - nicely formatted tables
#   patchwork  - combine several plots into one figure

required_packages <- c("tidyverse", "here", "broom", "rstatix", "emmeans",
                       "effectsize", "pwr", "knitr", "kableExtra", "patchwork")

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
  library(broom)
  library(rstatix)
  library(emmeans)
  library(effectsize)
  library(pwr)
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
