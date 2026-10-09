# Data README: `open_field_raw.csv` and `animals.csv`

**These two files are messy on purpose.** They imitate the data-entry sheet and the animal records of a real lab, and the project teaches you to find and fix the problems in code. **Do not edit them by hand.** The "clean" version is created by the analysis and saved to `output/data/`.

## At a glance

| | `open_field_raw.csv` | `animals.csv` |
|:--|:--|:--|
| **What it is** | The lab's data-entry sheet (behavior and weights) | The colony manager's animal records |
| **Used in** | `data_organization_and_display` project | `data_organization_and_display` project |
| **Data type** | **Simulated** (computer-generated; no animals were tested) | **Simulated** (computer-generated; no animals were tested) |
| **Rows / columns** | 41 rows (40 mice + 1 duplicated row) × 10 columns | 40 rows (one per mouse) × 5 columns |
| **Created by** | `R/01_simulate_data.R` (random seed 2026) | `R/01_simulate_data.R` (random seed 2026) |
| **Last updated** | 2026-10-09 | 2026-10-09 |

## Background

Forty adult mice from two genotypes (wild type, `WT`, and a knockout, `KO`) received either vehicle (a control injection) or a drug every day for three weeks. Body weight was measured at the start (week 0) and then weekly (weeks 1 to 3). After the last weight measurement, each mouse spent 10 minutes in an open-field arena, and a video system recorded how far it traveled and how long it spent in the center.

Mice lived in 8 cages of 5, each cage holding mice of one genotype and one sex. Treatments were assigned at random within each cage. The design has 10 mice in each of the four genotype × treatment groups (half female and half male in each genotype).

## Variables: `open_field_raw.csv` (as typed by the lab)

| Column | Meaning | Units | Notes |
|:--|:--|:--|:--|
| `Animal ID` | Unique mouse identifier | none | `M01`–`M40` |
| `Test Date` | Date of the open-field test | YYYY-MM-DD | |
| `Treatment` | Treatment group | none | Should be `Vehicle` or `Drug` |
| `Weight Wk0 (g)` … `Weight Wk3 (g)` | Body weight at weeks 0, 1, 2, and 3 | grams | Wide format: one column per week |
| `Distance (m)` | Distance traveled in the open field | meters | |
| `Center Time (s)` | Time spent in the center of the open field | seconds | |
| `Notes` | Free-text comments | none | |

## Variables: `animals.csv`

| Column | Meaning | Allowed values |
|:--|:--|:--|
| `animal_id` | Unique mouse identifier | `M01`–`M40` |
| `genotype` | Genotype | `WT`, `KO` |
| `cage` | Cage the mouse lives in | `C01`–`C08` (5 mice each) |
| `sex` | Sex | Should be `Female` or `Male` |
| `birth_date` | Date of birth | YYYY-MM-DD |

## The problems added on purpose

A real data set seldom arrives clean. These problems were added to the simulated data, and the notebook finds and fixes each one:

| Problem | Where | How many |
|:--|:--|:--:|
| Treatment spelled several ways (`Vehicle`, `vehicle`, `VEH`, `Vehicle ` with a trailing space, `Drug`, `drug`, `DRUG`, `Drug `) | `Treatment` | 8 spellings as typed (6 after `read_csv()` trims spaces) of 2 categories |
| Sex spelled several ways (`F`, `f`, `Female`, `M`, `m`, `Male`) | `animals.csv`, `sex` | 6 spellings of 2 categories |
| An animal ID typed in lowercase (`m14`) | `Animal ID` | 1 |
| A unit typed inside a number cell (`20.7 g`) | `Weight Wk0 (g)` | 1 |
| A missing-value code (`-999`) | `Weight Wk2 (g)` | 1 |
| The text `NA` in a number column | `Weight Wk1 (g)` | 1 |
| A blank cell with no explanation (a lost measurement) | `Weight Wk3 (g)` | 1 |
| A decimal-point typo (26.0 recorded as 260) | `Weight Wk2 (g)` | 1 |
| An exact duplicate row (mouse M27 entered twice) | whole row | 1 |
| A mouse that died (no later weights, no open-field test) | M33 | 1 |
| Tracking software failure (no open-field values) | M12 | 1 |
| Column names with spaces, capitals, and symbols | `open_field_raw.csv` | all |

## Missing values

In the raw file, missing values appear as blank cells, `NA`, or `-999`. After cleaning, all become `NA`, and **no value is ever replaced by 0.** Five mice have at least one missing measurement. Of these, M33 died (day 12) and M12 had a tracking failure (both recorded in `Notes`). The other three have a single missing weight with no recorded reason.

## How the data were generated

First a clean data set was created from a simple model, and then the problems above were added to a copy. The clean model was:

- **Body weight** at week 0: about 20.5 g for females and 25.5 g for males (1 g lower for KO), with a standard deviation of 1.3 g. Mice then gain about 0.45 g per week (0.55 g for males), and the drug lowers weight gain by 0.7 g per week.
- **Distance traveled**: 40 m, plus 6 m for KO, 12 m for drug, and 5 m for males, with a standard deviation of 10 m.
- **Time in the center**: right-skewed (gamma distribution) with a mean of 45 s (36 s for drug mice).

The random seed (2026) was the first one tried; it was **not** chosen to produce any particular pattern. This project runs no hypothesis tests, so there is no result to select for. See `R/01_simulate_data.R` for the exact code. Re-running it reproduces the same files.

## Good data habits shown (and violated) here

The raw files are a catalog of **what not to do** when entering data. The cleaned data show what to do instead:

- **One spelling per category**, and a data dictionary that lists the allowed values.
- **One row per observation, one column per variable** ("tidy" data). Weight at four weeks is *one* variable, so it belongs in one column, with the week in another.
- **Dates as YYYY-MM-DD.**
- **Blank cells or `NA` for missing values**, never 0 or -999, and comments in a separate `Notes` column, never in a number cell.
- **Raw data are never edited by hand.** All cleaning is done in code and logged, so it can be checked and repeated.
- **Documentation travels with the data** (this file).

## Using your own data

Save your data as a CSV with a header row and one row per subject (or per measurement). Put it in this folder and update the file names and variable names in `R/02_analysis.R` and the notebook. Messy data are fine, but write the cleaning rules *before* you start, and log every change. Write a README like this one.

## Citation

These are teaching data sets with no real-world scientific content. They should not be cited as evidence about any gene or drug.
