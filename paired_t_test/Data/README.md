# Data README: `open_field_crossover.csv` and `slice_recordings.csv`

## At a glance

| | `open_field_crossover.csv` | `slice_recordings.csv` |
|:--|:--|:--|
| **Used in** | Main paired t-test analysis | "Many measurements from one animal" lesson |
| **Topic** | Open-field locomotion in mice: drug vs. vehicle (crossover) | Neuron firing rates in brain slices: baseline vs. drug |
| **Data type** | **Simulated** (computer-generated for teaching; no animals were tested) | **Simulated** |
| **Rows / columns** | 28 rows (14 mice × 2 sessions) × 6 columns | 68 rows (34 neurons × 2 conditions) × 4 columns |
| **Created by** | `R/01_simulate_data.R` (random seed 5) | `R/01_simulate_data.R` (random seed 5) |
| **Last updated** | 2026-10-09 | 2026-10-09 |

Both files are in tidy **long** format (one row per measurement). The analysis reshapes them to **wide** format (one row per mouse, or per neuron) with `pivot_wider()` so each subject's two values sit side by side.

---

## `open_field_crossover.csv`

### Background

Fourteen mice were each tested in an open-field arena twice, one week apart: once after a low dose of a stimulant (**Drug**) and once after a saline injection (**Vehicle**). Half of the mice received the drug first and half the vehicle first (a **counterbalanced crossover design**). The outcome is the total distance traveled in a 10-minute session.

### Variables (data dictionary)

| Variable | Type | Units | Description | Allowed values |
|:--|:--|:--|:--|:--|
| `mouse_id` | Text | none | Unique animal identifier (each appears in two rows) | `M01`–`M14` |
| `sex` | Text (category) | none | Sex of the mouse | `F` = female, `M` = male |
| `sequence` | Text (category) | none | Order of treatments | `Drug first`, `Vehicle first` |
| `session` | Whole number | none | Session number (one week apart) | `1`, `2` |
| `treatment` | Text (category) | none | Treatment given before that session | `Drug`, `Vehicle` |
| `distance_m` | Number (1 decimal) | meters | Distance traveled in the 10-minute session | > 0; `NA` if tracking failed |

### Missing values

**One value is missing: mouse `M07`, session 2, `distance_m` is `NA`.** The video tracking failed in that session. Missing values are recorded as `NA` (not 0, which would claim the mouse did not move). A paired analysis needs both measurements, so the paired test uses a **complete-pairs** approach: this mouse is left out of the test (leaving 13 pairs) and the loss is reported. There are no other missing values and no duplicated mouse-by-session rows.

### How the data were generated

Each mouse has its own baseline activity level (mean 30 m, SD 6 m) and its own response to the drug (mean +5 m, SD 4 m). Each session adds a little random noise (SD 3 m), and the second session is 2 m lower on average than the first (habituation to the arena). Values were rounded to 0.1 m. The true standardized effect is therefore moderate to large, and the two measurements from a mouse are positively correlated because of the shared baseline.

---

## `slice_recordings.csv`

### Background

Brain slices from **8 mice** were used to record neurons (3 to 6 neurons per mouse, 34 in all). Each neuron's firing rate was measured at baseline and again after bath application of the drug. The data are paired at the level of the **neuron**, but neurons from one mouse are not independent of each other, so the **mouse** is the unit of analysis. This file exists to demonstrate that point.

### Variables (data dictionary)

| Variable | Type | Units | Description | Allowed values |
|:--|:--|:--|:--|:--|
| `mouse_id` | Text | none | Unique animal identifier | `S01`–`S08` |
| `neuron_id` | Text | none | Unique neuron identifier (mouse ID plus neuron number) | e.g. `S01-N1` |
| `condition` | Text (category) | none | Recording condition | `Baseline`, `Drug` |
| `firing_hz` | Number (1 decimal) | Hz (spikes per second) | Firing rate | ≥ 0 |

There are **no missing values** in this file.

### How the data were generated

Each mouse has its own baseline firing level (SD 1.0 Hz) and its own shift in response to the drug (SD 1.2 Hz); the neurons within a mouse then vary around those values (SD 1.5 Hz at baseline, SD 1.5 Hz in the change). The true average drug effect is +1.3 Hz. Neurons per mouse were set by hand to 4, 3, 5, 4, 6, 3, 4, 5. The estimated intraclass correlation in the resulting sample is about .28.

---

## About the random seed

Both files come from one script run with a single seed (**5**). The seed was chosen by a rule fixed before looking at any results: the first seed (trying 1, 2, 3, ...) for which (a) the mean of the analyzed paired differences (open field) and the mean neuron-level change (slices) were each within 10% of the true effect built into the simulation, and (b) a Shapiro–Wilk test did not flag the analyzed differences as non-normal (*p* > .20). The rule keeps the examples typical of the true situation. The sample sizes were chosen for a classroom example, not by a power analysis. A different seed would give somewhat different results. See `R/01_simulate_data.R` for the exact code; you do not need to run it to do the analysis.

## Good data habits shown here

- **One row per observation** and **one column per variable** ("tidy" data), with ID columns that link the rows belonging to the same animal.
- **Short, consistent variable names** with units where helpful (`_m`, `_hz`).
- **Missing values coded as `NA`**, never as 0 or a blank, and documented.
- **Raw data are never edited by hand.** All reshaping, exclusions, and calculations happen in code.
- **Documentation travels with the data** (this file).

## Using your own data

Save your data as a CSV in the same long layout (one row per subject per condition, with an ID column that identifies the subject), put it in this folder, and update the file names and variable names in `R/02_analysis.R` and the lab notebook. If you have several measurements per animal, decide on your unit of analysis first. Write a README like this one for your file.

## Citation

These are teaching data sets with no real-world scientific content. They should not be cited as evidence about drug effects on locomotion or neuronal firing.
