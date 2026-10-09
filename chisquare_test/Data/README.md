# Data README: genotype and seizure data

This project uses **two** small data sets, one for each chi-square test.

## At a glance

| | `genotype_counts.csv` | `seizure_treatment.csv` |
|:--|:--|:--|
| **Used for** | Part A: chi-square goodness-of-fit test | Part B: chi-square test of independence |
| **Topic** | Offspring genotypes from heterozygous crosses | Seizure after a convulsant drug, with or without an anticonvulsant pretreatment |
| **Data type** | **Simulated** (computer-generated; no animals were tested) | **Simulated** (computer-generated; no animals were tested) |
| **Rows / columns** | 120 rows (one per pup) × 3 columns | 100 rows (one per mouse) × 4 columns |
| **Created by** | `R/01_simulate_data.R` (random seed 8) | `R/01_simulate_data.R` (random seed 8) |
| **Last updated** | 2026-10-09 | 2026-10-09 |

## Background

**Genotype counts.** Heterozygous (+/−) mice were bred together. If the gene is not needed for survival, Mendel predicts offspring genotypes in a 1 : 2 : 1 ratio (25% wild type, 50% heterozygous, 25% knockout). The data are the genotypes of 120 pups from 15 litters of 8. Whether the knockout genotype is under-represented is the question.

**Seizure data.** One hundred mice received a convulsant drug. Fifty were pretreated with an anticonvulsant drug and 50 with vehicle (a control injection), by random assignment, with 25 females and 25 males in each group. Each mouse was scored for whether it had a seizure (Yes or No).

## Variables (data dictionary)

### `genotype_counts.csv`

| Variable | Type | Description | Allowed values |
|:--|:--|:--|:--|
| `pup_id` | Text | Unique pup identifier | `P001`–`P120` |
| `litter` | Text | Litter the pup was born in | `L01`–`L15` (8 pups each) |
| `genotype` | Text (category) | Genotype of the pup | `WT` (wild type), `Het` (heterozygous), `KO` (knockout) |

### `seizure_treatment.csv`

| Variable | Type | Description | Allowed values |
|:--|:--|:--|:--|
| `mouse_id` | Text | Unique animal identifier | `M001`–`M100` |
| `treatment` | Text (category) | Pretreatment | `Vehicle`, `Drug` |
| `sex` | Text (category) | Sex of the mouse | `F` = female, `M` = male |
| `seizure` | Text (category) | Did the mouse have a seizure? | `Yes`, `No` |

Both files are in **"tidy" form: one row per animal.** A chi-square test needs *counts*, so the analysis code counts the rows in each category (see the notebook). Many real data sets are instead stored as a table of counts; either form is fine as long as it is documented.

## Missing values

There are no missing values and no duplicated IDs in either file. (When data are missing, code them as `NA`, never as a category such as "0" or "unknown," unless "unknown" is a real category you want to analyze, and document them.)

## How the data were generated

**Genotypes** were drawn at random from a population with these true probabilities (the knockout is partly lethal, so it is under-represented):

| Genotype | True probability | Probability under the Mendelian 1 : 2 : 1 ratio |
|:--|:--:|:--:|
| WT | .33 | .25 |
| Het | .52 | .50 |
| KO | .15 | .25 |

The true effect is about Cohen's *w* = .26 (between small and medium).

**Seizures** were drawn from independent yes/no (Bernoulli) trials with a true seizure probability of **.60 after vehicle** and **.36 after the drug**, a true phi of about .24. Sex was assigned at random within each group and has **no** effect in the simulation.

**A note on the seed and the sample size.** A random sample can by chance look quite different from the population it came from. To get representative teaching examples, we used a rule decided before looking at any *p*-values: use the first seed (1, 2, 3, ...) for which *both* samples had an effect size within 0.05 of its true value. That was seed 8. Our first attempt used 96 pups and seed 7, and the genotype test then landed almost exactly on *p* = .05, so we **increased the sample to 120 pups** and applied the same rule again. The result was still close to the cutoff (*p* = .049), and we **kept it on purpose**: the notebook uses it to teach why p-values near .05 deserve caution. Choosing a seed and a sample size after seeing results is a form of selection, so we disclose it here. A real experiment gets one sample and no chance to choose. See `R/01_simulate_data.R` for the exact code. Re-running it reproduces these same files.

## Good data habits shown here

- **One row per observation** (one animal), **one column per variable**, and a column that identifies the group.
- **Meaningful labels** (`Yes`/`No`, `Vehicle`/`Drug`) instead of unexplained codes such as 0/1.
- **A litter identifier**, so the non-independence of littermates can be examined.
- **Raw data are never edited by hand.** All counting and calculations happen in code.
- **Documentation travels with the data** (this file).

## Using your own data

Save your data as a CSV with a header row, one row per subject, and one column for each categorical variable. Put it in this folder and update the file names and variable names in `R/02_analysis.R` and the notebook. Write a README like this one.

## Citation

These are teaching data sets with no real-world scientific content. They should not be cited as evidence about any gene or drug.
