# Data README: `locomotor_dose.csv`

## At a glance

| | |
|:--|:--|
| **File** | `locomotor_dose.csv` |
| **Used in** | `oneway_anova` project |
| **Topic** | Locomotor activity in mice after a stimulant drug at four doses |
| **Data type** | **Simulated** (computer-generated for teaching; no animals were tested) |
| **Rows / columns** | 48 rows (one per mouse) × 5 columns |
| **Created by** | `R/01_simulate_data.R` (random seed 1; see the note on seed selection below) |
| **Last updated** | 2026-10-09 |

## Background

Stimulant drugs often increase locomotor activity (movement) in an open field. In this example, 48 mice each received one injection and were then placed in an open-field chamber. The distance each mouse traveled was recorded for 30 minutes.

The study compares **four independent groups** that differ in dose:

| Group | Dose (mg/kg) | n |
|:--|:--:|:--:|
| Vehicle (control injection, no drug) | 0 | 12 |
| Low | 2 | 12 |
| Medium | 5 | 12 |
| High | 10 | 12 |

Mice were randomly assigned to groups. Each mouse is in only one group and was tested once.

## Variables (data dictionary)

| Variable | Type | Units | Description | Allowed values |
|:--|:--|:--|:--|:--|
| `mouse_id` | Text | none | Unique animal identifier | `M01`–`M48` |
| `dose_group` | Text (category) | none | Dose group | `Vehicle`, `Low`, `Medium`, `High` |
| `dose_mgkg` | Number | mg/kg | Dose given | 0, 2, 5, 10 |
| `sex` | Text (category) | none | Sex of the mouse | `F` = female, `M` = male |
| `distance_m` | Number (1 decimal) | meters in 30 minutes | Total distance traveled | ≥ 0 |

## Missing values

There are no missing values and no duplicated IDs. (When data are missing, code them as `NA`, never as 0, and document them.)

## How the data were generated

Distances were simulated from normal distributions with the same standard deviation (10 m) in every group and these true means:

| Group | True mean (m) |
|:--|:--:|
| Vehicle | 40 |
| Low | 45 |
| Medium | 50 |
| High | 55 |

The true effect is a large one (Cohen's *f* of about 0.56; η² of about 0.24). Sex was assigned at random within each group and has **no** effect in the simulation. Values were rounded to 0.1 and kept at or above 0.

**A note on the seed.** A random sample can by chance look quite different from the population it came from. To get a representative teaching example, we used a rule decided before looking at any *p*-values: use the first seed (1, 2, 3, ...) whose sample η² was within 0.04 of the true value. That was seed 1. Choosing a seed is a form of selection, so we disclose it. A real experiment gets only one sample and no chance to choose. The notebook uses simulation to show how often samples like this one, and unlucky ones, occur. See `R/01_simulate_data.R` for the exact code. Re-running it reproduces this same file.

## Good data habits shown here

- **One row per observation** (one mouse), **one column per variable** ("tidy" data).
- **Both a label and a number for dose** (`dose_group` and `dose_mgkg`), so the group can be treated as a category (ANOVA) or as a number (a possible follow-up regression).
- **Raw data are never edited by hand.** All cleaning and calculations happen in code.
- **Documentation travels with the data** (this file).

## Using your own data

Save your data as a CSV with a header row, one row per subject, and one column giving each subject's group. Put it in this folder and update the file name and variable names in `R/02_analysis.R` and the notebook. Write a README like this one.

## Citation

This is a teaching data set with no real-world scientific content. It should not be cited as evidence about any drug.
