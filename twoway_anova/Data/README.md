# Data README: `bdnf_housing_stress.csv`

## At a glance

| | |
|:--|:--|
| **File** | `bdnf_housing_stress.csv` |
| **Used in** | `twoway_anova` project |
| **Topic** | Hippocampal BDNF in mice: housing (standard vs. enriched) crossed with chronic stress (control vs. stress) |
| **Data type** | **Simulated** (computer-generated for teaching; no animals were tested) |
| **Rows / columns** | 40 rows (one per mouse) × 5 columns |
| **Created by** | `R/01_simulate_data.R` (random seed 14; see the note on seed selection below) |
| **Last updated** | 2026-10-09 |

## Background

Chronic stress lowers brain-derived neurotrophic factor (BDNF), a protein that supports neuron survival and plasticity, in the hippocampus. Environmental enrichment (larger cages with a running wheel, tunnels, and nesting material) is thought to protect against this. The study asks whether enrichment *buffers* the effect of stress, which is a question about an **interaction**.

The design is a **2 × 2 factorial**: two factors, each with two levels, with every combination present.

| | Control | Stress |
|:--|:--:|:--:|
| **Standard housing** | 10 mice | 10 mice |
| **Enriched housing** | 10 mice | 10 mice |

Each mouse is in exactly one cell and was measured once. Hippocampal tissue was assayed for BDNF protein by ELISA.

## Variables (data dictionary)

| Variable | Type | Units | Description | Allowed values |
|:--|:--|:--|:--|:--|
| `mouse_id` | Text | none | Unique animal identifier | `M01`–`M40` |
| `housing` | Text (category) | none | Housing condition | `Standard`, `Enriched` |
| `stress` | Text (category) | none | Stress condition | `Control`, `Stress` |
| `sex` | Text (category) | none | Sex of the mouse | `F` = female, `M` = male |
| `bdnf_pgmg` | Number (1 decimal) | pg BDNF per mg total protein | Hippocampal BDNF concentration | > 0 |

## Missing values

There are no missing values and no duplicated IDs. (When data are missing, code them as `NA`, never as 0, and document them.)

## How the data were generated

BDNF values were simulated from normal distributions with the same standard deviation (15 pg/mg) in every cell and these true means:

| Cell | True mean (pg/mg) |
|:--|:--:|
| Standard / Control | 100 |
| Standard / Stress | 70 |
| Enriched / Control | 108 |
| Enriched / Stress | 102 |

Stress lowers BDNF by 30 units in standard cages but only by 6 units in enriched cages. That difference is the interaction (true partial η² of about .14, a large effect). The true main effects are also large (partial η² of about .31 for housing and .26 for stress). Sex was assigned at random within each cell (5 females and 5 males) and has **no** effect in the simulation.

**A note on the seed.** A random sample can by chance look quite different from the population it came from. To get a representative teaching example, we used a rule decided before looking at any *p*-values: use the first seed (1, 2, 3, ...) for which the sample partial η² of each of the three effects was within 0.05 of its true value. That was seed 14. Choosing a seed is a form of selection, so we disclose it. A real experiment gets only one sample and no chance to choose. The notebook uses simulation to show how often experiments like this one detect the interaction (about two in three). See `R/01_simulate_data.R` for the exact code. Re-running it reproduces this same file.

## Good data habits shown here

- **One row per observation** (one mouse), **one column per variable** ("tidy" data), with one column for each factor. Do not spread the four cells across four columns.
- **Short, consistent variable names** with units where helpful (`_pgmg`).
- **Raw data are never edited by hand.** All cleaning and calculations happen in code.
- **Documentation travels with the data** (this file).

## Using your own data

Save your data as a CSV with a header row, one row per subject, one column for each of the two grouping variables, and one column for the outcome. Put it in this folder and update the file name and variable names in `R/02_analysis.R` and the notebook. Write a README like this one.

## Citation

This is a teaching data set with no real-world scientific content. It should not be cited as evidence about stress, enrichment, or BDNF.
