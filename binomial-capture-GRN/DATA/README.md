# DATA/

Simulation output read by the `ggplot_*` scripts. The `.gitkeep` files only
preserve the folder structure — **replace them with the original CSVs** (too
large to include here). The exact file names from the original commit are below;
keep these names, since the scripts read them literally.

## Files at `DATA/` root

- `gg_df_1S_S1_Ru12_Rb52_Su1000_B4_n3.csv` — single-gene auto-regulatory loop, tidy data frame read by `scripts/ggplot_1S_2S_result.R` (variable `gg_df_1S`).
- `gg_df_2S_L40_H200_B3_n2.csv`
- `gg_df_2S_L40_H200_B3_n2_case2.csv`
- `gg_df_2S_L40_H200_B3_n2_case3.csv` — two-gene toggle switch, read by `scripts/ggplot_1S_2S_result.R` (variable `gg_df_2S`; the `write.csv` line in `simulation_result_2S.R` points at `_case3`).
- `gg_df_3S_B2_dt0_2_n3.csv` — three-gene repressilator, tidy data frame.
- `3S_TS_B1_Tmax50_5000sim.csv` — repressilator time-series (5000 sims).

## `3S_periodogram/` — three-gene repressilator

Actual-process files at the folder root:

- `act_for_density.csv`      — actual-process protein distributions
- `act_periodogram_df.csv`   — actual-process power spectra
- `act_traj_df.csv`          — actual-process trajectories (large)

Per capture-probability subfolders `p_0.5/`, `p_0.25/`, `p_0.1/`, each containing:

- `obs_for_density.csv`, `obs_periodogram_df.csv`   — observed (true-observed)
- `map_for_density.csv`, `map_periodogram_df.csv`   — mapped (predicted-observed)
- (`obs_traj_df.csv`, `map_traj_df.csv` — trajectories, large; optional)

`unnormalized/` — un-normalized periodogram variants (not read by the active
plotting scripts; kept for completeness).

## `discrete_works/` and `discrete_breaks/`

CSVs read via `list.files()` by `scripts/ggplot_discrete_result.R` and
`scripts/ggplot_discrete_breaks.R` (all `*.csv` in each folder are loaded):

- `discrete_works/`  — discrete CME regime where the mapping holds
- `discrete_breaks/` — discrete CME regime where the mapping breaks down

## Regenerating

Run the matching `scripts/simulation_result_*` script with its `write.csv(...)`
lines uncommented, then the corresponding `ggplot_*` script.
