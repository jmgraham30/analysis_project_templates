# Data: `neuron_firing.csv`

**These data are simulated.** No animals were recorded. The script `R/01_simulate_data.R` documents exactly how they were generated.

## Scenario

Extracellular recordings from 120 neurons in three brain regions (Cortex, Hippocampus, Striatum; 40 neurons each), collected from 12 mice (10 neurons per mouse). Each row is one neuron.

## Variables

| Column | Type | Description |
|---|---|---|
| `neuron_id` | text | Unique neuron label (N001–N120) |
| `mouse_id` | text | Mouse the neuron was recorded from (M01–M12) |
| `region` | nominal | Cortex, Hippocampus, or Striatum |
| `cell_type` | nominal | Regular spiking (RS) or fast spiking (FS); about 20% FS |
| `firing_rate_hz` | ratio (continuous) | Average firing rate, spikes per second. Right-skewed (lognormal) |
| `spike_width_ms` | ratio (continuous) | Trough-to-peak spike width in milliseconds. Roughly symmetric; narrower in FS cells |
| `tuning_rating` | ordinal | Rated stimulus selectivity, 1 (none) to 5 (very sharp) |

There are no missing values.

## Important notes

- **Pseudoreplication.** Neurons from the same mouse are not independent. The notebook treats the 120 neurons as a simple sample for teaching descriptive statistics and flags where this matters for inference.
- **Seed.** The script uses `set.seed(2026)`, the first seed tried; it was not selected for how the sample looks.
- Re-running `R/01_simulate_data.R` re-creates an identical file.
