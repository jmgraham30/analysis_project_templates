# Data: `neurogenesis_running.csv`

**These data are simulated.** No animals were tested. The script `R/01_simulate_data.R` documents exactly how they were generated.

## Scenario

Sixty adult mice (8 to 20 weeks old) lived with a running wheel for two weeks. For each mouse we recorded the average distance run per night, and afterward counted doublecortin-positive (DCX+) cells, a marker of newly born neurons, in the dentate gyrus of the hippocampus. Each row is one mouse.

## Variables

| Column | Type | Description |
|---|---|---|
| `mouse_id` | text | Unique mouse label (M01–M60) |
| `sex` | categorical | `F` (30 mice) or `M` (30 mice) |
| `age_weeks` | numeric | Age in weeks (8–20) |
| `wheel_km` | numeric | Average distance run per night, in kilometers (predictor) |
| `dcx_cells` | numeric (count) | DCX+ cells counted in the dentate gyrus (outcome) |

There are no missing values.

## How the data were generated

In the population, the expected cell count is `110 + 12 × wheel_km − 3 × (age_weeks − 14)`, with random noise (SD 25). The true correlation between running and cell count is about r = .62. Sex has no effect in the simulation.

## Seed

The script uses `set.seed(4)`. Seed 4 was the first seed (counting from 1) for which the sample correlation was within 0.03 of the true value *and* the random noise was not noticeably non-normal. Seed 1 met only the first rule and had strongly skewed noise, so a second rule was added before any analysis was run. This is disclosed in the script header.

Re-running `R/01_simulate_data.R` re-creates an identical file.
