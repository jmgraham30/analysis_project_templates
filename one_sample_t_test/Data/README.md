# Data README: `novel_object_recognition.csv`

## At a glance

| | |
|:--|:--|
| **File** | `novel_object_recognition.csv` |
| **Used in** | `one_sample_t_test` project |
| **Topic** | Novel object recognition (NOR) memory test in mice |
| **Data type** | **Simulated** (computer-generated for teaching; no animals were tested) |
| **Rows / columns** | 26 rows (one per mouse) × 5 columns |
| **Created by** | `R/01_simulate_data.R` (random seed 5) |
| **Last updated** | 2026-10-09 |

## Background

In the novel object recognition test, a mouse first becomes familiar with an object. Later, the familiar object is paired with a **novel** object, and an observer records how long the mouse actively explores each one. Mice that remember the familiar object tend to spend more time with the novel one.

Each row of this file holds the exploration times for one mouse during the **test phase**.

## Variables (data dictionary)

| Variable | Type | Units | Description | Allowed values |
|:--|:--|:--|:--|:--|
| `mouse_id` | Text | none | Unique animal identifier | `M01`–`M26` |
| `sex` | Text (category) | none | Sex of the mouse | `F` = female, `M` = male |
| `age_weeks` | Whole number | weeks | Age at testing | 10–14 |
| `time_novel_s` | Number (1 decimal) | seconds | Time spent actively exploring the **novel** object | ≥ 0 |
| `time_familiar_s` | Number (1 decimal) | seconds | Time spent actively exploring the **familiar** object | ≥ 0 |

**Missing values:** none. **Duplicated IDs:** none.

## Derived variables (calculated in the analysis, *not* stored in the file)

The raw file contains only what an experimenter would record. These values are calculated in `R/02_analysis.R` and the lab notebook:

- **Total exploration time** = `time_novel_s + time_familiar_s`
- **Discrimination index (DI)** = (`time_novel_s` − `time_familiar_s`) / (`time_novel_s` + `time_familiar_s`). DI ranges from −1 to +1; 0 means no preference.

## Exclusion rule

Mice with **total exploration time under 20 seconds** are excluded from the analysis because their DI is based on too little information. Two mice (`M07` and `M19`) meet this rule. They remain in the file, so the raw data are complete, and they are removed **in code**.

## How the data were generated

The values were simulated in R from a population with a true mean DI of 0.15 and a true standard deviation of 0.28 (a true Cohen's *d* of about 0.54, a medium effect). Total exploration time was drawn from a right-skewed (lognormal) distribution centered near 40 s. Times were rounded to 0.1 s. See `R/01_simulate_data.R` for the exact code. You do not need to run it to do the analysis; re-running it reproduces this same file.

## Good data habits shown here

- **One row per observation** (one mouse), **one column per variable** ("tidy" data).
- **Short, consistent variable names** with units included (`_s` for seconds, `_weeks`).
- **Raw data are never edited by hand.** All cleaning, exclusions, and calculations happen in code.
- **Documentation travels with the data** (this file).

## Using your own data

To adapt this project for your own study, save your data as a CSV with the same layout (a header row of variable names, then one row per subject), put it in this folder, and update the file name and variable names in `R/02_analysis.R` and the lab notebook. Write a README like this one for your file.

## Citation

This is a teaching data set with no real-world scientific content. It should not be cited as evidence about mouse behavior.
