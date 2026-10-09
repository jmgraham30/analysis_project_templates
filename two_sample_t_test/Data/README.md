# Data README: `elevated_plus_maze.csv`

## At a glance

| | |
|:--|:--|
| **File** | `elevated_plus_maze.csv` |
| **Used in** | `two_sample_t_test` project |
| **Topic** | Anxiety-like behavior in mice on the elevated plus maze: chronic stress vs. control |
| **Data type** | **Simulated** (computer-generated for teaching; no animals were tested) |
| **Rows / columns** | 29 rows (one per mouse) × 5 columns |
| **Created by** | `R/01_simulate_data.R` (random seed 2026) |
| **Last updated** | 2026-10-09 |

## Background

The elevated plus maze is a plus-shaped maze raised off the floor, with two open arms and two closed (walled) arms. Mice prefer the enclosed arms, and more anxious mice tend to spend a smaller percentage of the test in the open arms. Each mouse explored the maze once for 5 minutes.

The study compares two **independent groups**:

- **Control** (15 mice): handled normally.
- **Stress** (14 mice): exposed to two weeks of chronic mild stress.

Mice were randomly assigned to groups, and each mouse is in only one group and was tested once.

## Variables (data dictionary)

| Variable | Type | Units | Description | Allowed values |
|:--|:--|:--|:--|:--|
| `mouse_id` | Text | none | Unique animal identifier | `C01`–`C15` (control), `S01`–`S14` (stress) |
| `group` | Text (category) | none | Experimental group | `Control`, `Stress` |
| `sex` | Text (category) | none | Sex of the mouse | `F` = female, `M` = male |
| `open_arm_pct` | Number (1 decimal) | % of the 5-minute test | Percent of the test spent in the open arms | 0–100; `NA` if the test was not completed |
| `total_entries` | Whole number | count | Total number of arm entries (a rough measure of general activity) | ≥ 0; `NA` if the test was not completed |

## Missing values

**One mouse (`S09`, stress group) has `NA` for both `open_arm_pct` and `total_entries`.** This mouse jumped off the maze, so its test could not be completed. Missing values are recorded as `NA` (not as 0, which would falsely claim the mouse spent 0% of its time in the open arms). The analysis uses a complete-case approach: this mouse is left out of the test and the loss is reported. There are no other missing values and no duplicated IDs.

## How the data were generated

The values were simulated in R from two normal distributions (means and standard deviations in percent of time in the open arms; values were kept between 0 and 100 and rounded to 0.1):

| Group | n recorded | True mean | True SD |
|:--|:--:|:--:|:--:|
| Control | 15 | 34 | 9 |
| Stress | 14 | 25 | 11 |

The true standardized difference is about *d* = 0.9 (a large effect), and the spreads are slightly unequal on purpose. Total arm entries were drawn from a Poisson distribution with a mean of 24 in both groups. The seed (2026) was the first one tried; it was **not** chosen to produce a particular result. Because these are random samples, the sample means differ from the true values, which is exactly the sampling variability the analysis is designed to handle. See `R/01_simulate_data.R` for the exact code. You do not need to run it to do the analysis; re-running it reproduces this same file.

## Good data habits shown here

- **One row per observation** (one mouse), **one column per variable** ("tidy" data), with a column that identifies the group.
- **Short, consistent variable names** with units where helpful (`_pct`).
- **Missing values coded as `NA`**, never as 0 or a blank, and documented.
- **Raw data are never edited by hand.** All cleaning, exclusions, and calculations happen in code.
- **Documentation travels with the data** (this file).

## Using your own data

To adapt this project for your own study, save your data as a CSV with the same layout (a header row of variable names, then one row per subject, with one column giving each subject's group), put it in this folder, and update the file name and variable names in `R/02_analysis.R` and the lab notebook. Write a README like this one for your file.

## Citation

This is a teaching data set with no real-world scientific content. It should not be cited as evidence about stress or anxiety-like behavior in mice.
