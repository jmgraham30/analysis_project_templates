# =============================================================================
# 01_simulate_data.R
# Project:  Chi-square tests (genotype ratios and seizure outcomes)
# Purpose:  Document HOW the two teaching data sets were created.
#
# WHY SIMULATED DATA?
#   These data are computer-generated (simulated); no animals were tested. When
#   we simulate, we know the "truth" (the real probabilities), so we can check
#   how well our statistics recover it and demonstrate power and error rates
#   honestly. Real raw data are never created by a script like this; they come
#   from your lab notebook or instrument and are saved untouched.
#
# ABOUT THE RANDOM SEED
#   A random sample can, by chance, look very different from the population it
#   came from. For a clear teaching example we wanted REPRESENTATIVE samples, so
#   we used this rule, decided before looking at any p-values: use the first seed
#   (counting 1, 2, 3, ...) for which BOTH samples have an effect size within 0.05
#   of the TRUE effect size (Cohen's w of the genotype counts, and phi of the
#   seizure-by-treatment table). That was seed 8. (Our first attempt used 96
#   pups; the genotype test then landed almost exactly on p = .05, so we increased
#   the sample to 120 pups and applied the same rule again. The result was still
#   close to the cutoff, p = .049, and we kept it on purpose: the notebook uses
#   it to teach why p-values near .05 deserve caution.) Choosing a seed and a
#   sample size after seeing results is a form of selection, so we disclose it here; a real experiment gets one sample and
#   no chance to choose.
#
# YOU DO NOT NEED TO RUN THIS FILE to do the analysis. The data already exist in
# Data/. Running it again re-creates identical files because the seed is set
# (same seed = same "random" numbers).
# =============================================================================

library(tidyverse)
library(here)

# ---- 1. The scenarios ---------------------------------------------------------
# DATA SET 1: genotype counts (a "goodness-of-fit" question).
#   Heterozygous (+/-) mice are bred together. If the gene is not needed for
#   survival, Mendel predicts offspring genotypes in a 1 : 2 : 1 ratio
#   (25% wild type, 50% heterozygous, 25% knockout). Here the knockout is
#   partly lethal, so knockouts are UNDER-represented. We genotype 120 pups
#   from 15 litters of 8.
#
# DATA SET 2: seizure outcome by treatment (a "test of independence").
#   100 mice receive a convulsant drug. Half were pretreated with an
#   anticonvulsant drug and half with vehicle (a control injection). Outcome:
#   did the mouse have a seizure (Yes/No)?

# ---- 2. The "true" population values used to generate the data ---------------
true_geno <- c(WT = 0.33, Het = 0.52, KO = 0.15)   # true genotype probabilities
mendel    <- c(WT = 0.25, Het = 0.50, KO = 0.25)   # the H0 ratio (1 : 2 : 1)
true_seiz <- c(Vehicle = 0.60, Drug = 0.36)        # true probability of a seizure
n_pups    <- 120    # 15 litters x 8 pups
n_group   <- 50     # mice per treatment group

# Effect sizes implied by the truth (used by the seed rule)
true_w   <- sqrt(sum((true_geno - mendel)^2 / mendel))               # Cohen's w
pbar     <- mean(true_seiz)
true_phi <- abs(diff(true_seiz)) / (2 * sqrt(pbar * (1 - pbar)))     # phi for a 2 x 2 table

# ---- 3. A function that generates both data sets for a given seed ------------
make_data <- function(seed) {
  set.seed(seed)

  geno <- tibble(
    pup_id   = sprintf("P%03d", 1:n_pups),
    litter   = rep(sprintf("L%02d", 1:15), each = 8),
    genotype = sample(names(true_geno), n_pups, replace = TRUE, prob = true_geno)
  )

  seiz <- tibble(
    mouse_id  = sprintf("M%03d", 1:(2 * n_group)),
    treatment = rep(c("Vehicle", "Drug"), each = n_group),
    sex       = unlist(map(1:2, ~ sample(rep(c("F", "M"), length.out = n_group)))),
    seizure   = if_else(rbinom(2 * n_group, 1, rep(true_seiz, each = n_group)) == 1,
                        "Yes", "No")
  )

  list(geno = geno, seiz = seiz)
}

# ---- 4. Generate and save (the seed rule is explained at the top) -----------
seed <- 8
dat  <- make_data(seed)

write_csv(dat$geno, here("Data", "genotype_counts.csv"))
write_csv(dat$seiz, here("Data", "seizure_treatment.csv"))
message("Saved Data/genotype_counts.csv (", nrow(dat$geno), " rows) and ",
        "Data/seizure_treatment.csv (", nrow(dat$seiz), " rows)")
