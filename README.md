# Imperfect molecular detection can renormalize apparent kinetic rates in stochastic gene regulatory networks

Simulation code accompanying the paper by I. Zabaikina and R. Grima.

The code implements a hybrid Gillespie (PDMP) simulator for stochastic gene
regulatory networks under a binomial/Gaussian molecular-capture observation
model, and reproduces the simulation figures in the manuscript (auto-regulatory
loop, two-gene toggle switch, and three-gene repressilator, together with the
discrete chemical-master-equation results).

## Requirements

- R (>= 4.1 recommended)
- The packages listed in `install_packages.R`. Install them with:

```r
source("install_packages.R")
```

Key packages: `data.table`, `ggplot2`/`cowplot`/`ggpubr`/`wesanderson` (plots),
`future`/`future.apply`/`parallel` (parallel simulation), `deSolve` (deterministic
repressilator), `pracma`, `RConics`, `logspline`, `here`.

## How to run

Open `binomial-capture-GRN.Rproj` in RStudio (this sets the project root, which
all paths resolve against via the `here` package), then run scripts from
`scripts/`. Each simulation script sources `R/0preambule.R`, which loads the
libraries and all function files in `R/`.

Parallelism uses `parallel::detectCores() - 1` cores by default.

## Repository structure

```
R/         function definitions and shared plotting helpers (sourced, not run directly)
scripts/   top-level scripts that run simulations and produce the figures
DATA/      simulation output (CSV) read by the plotting scripts (see DATA/README.md)
```

- `R/simulation_functions.R` — core hybrid-Gillespie / PDMP routines
- `R/simulation_parallel.R` — parallel wrappers over many runs
- `R/simulation_trajectories_correct.R`, `R/simulation_trajectories_on_grid.R` — trajectory samplers
- `R/simulation_periodogram_functions.R` — periodogram / power-spectrum estimation (repressilator)
- `R/deterministic_3S_solution.R` — deterministic repressilator ODE (via `deSolve`)
- `R/ggplot_ggbase0.R`, `R/ggplot_line_logspline.R`, `R/ggplot_3S_pdf_csv_helper.R` — plotting helpers
- `R/0preambule.R` — libraries + `source()` of the files above


| Script | Produces |
|---|---|
| `scripts/simulation_result_1S.R` | single-gene auto-regulatory loop — simulation data |
| `scripts/simulation_result_2S.R` | two-gene toggle switch — simulation data |
| `scripts/simulation_result_3S.R` | three-gene repressilator — trajectories + periodograms |
| `scripts/ggplot_1S_2S_result.R` | figure: auto-regulatory loop and toggle switch distributions (Fig. 3-4) |
| `scripts/ggplot_3S_pdf_result.R` | figure: repressilator protein distributions (Fig. 5, panels A–C) |
| `scripts/ggplot_3S_periodogram_result.R` | figure: repressilator power spectra (Fig. 5, panels D–F) |
| `scripts/ggplot_discrete_result.R` | figure: discrete CME results (Fig. 2A) |
| `scripts/ggplot_discrete_breaks.R` | figure: discrete CME breakdown at small p (Fig. 2B) |

Simulation scripts write CSVs into `DATA/` (the write lines are commented out by
default to avoid overwriting); the `ggplot_*` scripts read those CSVs back to
render the figures. To reproduce a figure from scratch: run the corresponding
`simulation_result_*` script (uncommenting the `write.csv` lines), then the
matching `ggplot_*` script.

## Data availability

All figures can be generated using this project.
The aggregated data needed to reproduce every figure (the `*_for_density.csv`
and `*_periodogram_df.csv` files, and the `gg_df_*` frames) are included under
`DATA/`. Some raw simulation outputs — the large trajectory files
(`*_traj_df.csv`, `3S_TS_B1_Tmax50_5000sim.csv`) — are omitted here due to their
size and are **available from the authors on request**; they can also
be regenerated with the `scripts/simulation_result_*` scripts. See
`DATA/README.md` for the full file listing.

## Citation

If you use this code, please cite the paper and this repository (Zenodo DOI).
