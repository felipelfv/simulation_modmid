# README

This repository contains the code, materials, and manuscript for the simulation
study **"Latent Variable Moderated Mediation Using Correct and Misspecified
Models"** by Felipe Fontana Vieira and Yves Rosseel.

The study compares three estimators of a latent moderated-mediation model under
varying sample size, reliability, exogenous-predictor distribution, and structural
misspecification:

- **LSAM** — local structural-after-measurement, `lavaan::sam(sam.method = "local", se = "local")`
- **LMS** — latent moderated structural equations, `modsem::modsem_da(method = "lms", robust.se = TRUE)`
- **UPI** — extended unconstrained product indicators (double mean-centering), `modsem::modsem(method = "dblcent", estimator = "MLR")`

The design is a full factorial of 4 sample sizes × 3 interaction effects × 2
reliabilities × 5 distributions × 9 misspecifications = 1080 conditions, each run
for 1,000 replications.

<sub>The code in this repository is licensed under the [MIT License](LICENSE). The
manuscript, figures, and other non-code content are licensed under a [Creative
Commons Attribution 4.0 International License (CC BY 4.0)](https://creativecommons.org/licenses/by/4.0/).
This project and the preregistration for it are archived at
<https://doi.org/10.5281/zenodo.20703257>.</sub>

## Repository Structure

The project consists of a set of R scripts forming a three-stage pipeline (simulation,
analysis and figures, manuscript), plus a Quarto manuscript and a Nix environment
definition. Calibration is a one-off prerequisite whose output is committed.

### 1. Simulation

```
.
├── gen_env.R          # generates default.nix via {rix} (pinned R + packages + Quarto + TeX)
├── default.nix        # generated environment definition (built with nix-shell)
├── calibration/
│   ├── calibrate.R                 # calibrates misspecification coefficients & residual variances
│   ├── calibration_explanation.qmd # full derivation of the calibration scheme
│   └── calibration_results.rds     # output; copy to results/calibration.rds before Step 1
├── sim_mc.R           # main driver: builds the design grid and runs all conditions in parallel
└── helpers_mc.R       # data generation, the three estimators, MC/Wald inference, admissibility
                       #   checks, and the per-replication / per-condition workers (sourced by sim_mc.R)
```

The simulation writes its outputs to `results/`:

```
results/
├── calibration.rds          # calibration lookup consumed by sim_mc.R
├── design.rds               # the 1080-row design grid
├── estimates_mc.rds         # combined per-fit estimates (~1.8 GB, ~58M rows)
├── metrics_mc.rds           # combined per-replication fit metrics / timings
├── progress.log             # per-condition progress written during the run
├── summary.rds              # performance metrics      (← analyze_mc.R)
├── rejection_imm_a3.rds     # Type I error / power     (← analyze_mc.R)
└── convergence_report.rds   # convergence & outlier counts (← analyze_mc.R)
```

The per-condition files `condition_mc_0001.rds … condition_mc_1080.rds` (one
per condition, written by `sim_mc.R` for resumability) are intermediate
artifacts combined into `estimates_mc.rds` at the end of Step 1, and are *not
included* in this repository.

`results/` is not committed. The large combined outputs are archived on Zenodo
(DOI above); download them into `results/` to re-run `analyze_mc.R`.

### 2. Analysis and figures

```
.
├── analyze_mc.R       # results/estimates_mc.rds → summary / rejection / convergence tables
└── figures.R          # results/{summary,rejection_imm_a3}.rds → plots/*.png (16 figures, 14 used in the paper)
```

### 3. Manuscript

```
manuscript/
├── article.qmd        # paper source (APA 7 via the apaquarto extension)
├── references.bib     # bibliography
├── article.pdf        # rendered output
├── article.tex        # kept LaTeX source of the rendered output
└── _extensions/       # vendored apaquarto extension
docs/
├── article.pdf        # copy of the rendered PDF served via GitHub Pages
└── index.html         # redirect to article.pdf
```

`article.pdf`, `article.tex`, and `docs/` are rebuilt and committed on every
push to `main` by the GitHub Actions workflow (`.github/workflows/render.yaml`).

The manuscript embeds the figures from `plots/` directly (via
`knitr::include_graphics`).

## Pipeline

The pipeline runs in three stages, each consuming the previous stage's output.
Calibration (`calibration/calibrate.R`) is a prerequisite whose result is
already provided as `calibration/calibration_results.rds`; `sim_mc.R` reads it
from `results/calibration.rds`, so copy it there first.

### Step 1: Run the simulation

```sh
nix-shell --pure --run 'OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 Rscript sim_mc.R'
```

`sim_mc.R` builds the design grid (`results/design.rds`), reads
`results/calibration.rds`, and runs every condition in parallel via
`parallel::mclapply` with reproducible `L'Ecuyer-CMRG` streams (master seed
`1234`). Each condition is saved to `results/condition_mc_<id>.rds` as it
completes (so the run is **resumable**), then all are combined into
`results/estimates_mc.rds` and `results/metrics_mc.rds`.

### Step 2: Process results and build figures

```r
Rscript analyze_mc.R    # → results/{summary,rejection_imm_a3,convergence_report}.rds
Rscript figures.R       # → plots/*.png
```

`analyze_mc.R` applies the two convergence criteria (soft / strict), removes
outliers, and computes performance metrics with Monte Carlo standard errors.
`figures.R` turns those tables into 16 figures, 14 of which appear in the paper.

### Step 3: Render the manuscript

```sh
quarto render manuscript/article.qmd    # → manuscript/article.pdf
```

### Data flow

```
calibrate.R ──> calibration.rds
                     │
helpers_mc.R ──┐     ▼
               └─> sim_mc.R ──> condition_mc_*.rds ──> estimates_mc.rds
                                                            │
                                            analyze_mc.R ──┤──> summary.rds
                                                            │    rejection_imm_a3.rds
                                                            │    convergence_report.rds
                                                            ▼
                                               figures.R ──> plots/*.png
                                                            │
                                            article.qmd ───> article.pdf
```

## Cloning this repository

```sh
git clone https://github.com/felipelfv/simulation_modmid
cd simulation_modmid
nix-shell        # builds the pinned environment from default.nix
```

The environment is fully pinned with [Nix](https://nixos.org/). `default.nix` is
generated by `gen_env.R` (using the `{rix}` package) and fixes the R version, all
R packages, the Quarto CLI, and the LaTeX dependencies to a single snapshot
(nixpkgs date `2026-05-11`), with `lavaan` and `modsem` built from specific GitHub
commits. To regenerate the environment definition:

```r
Rscript gen_env.R    # requires the {rix} package; rewrites default.nix
```

The same environment is used by the GitHub Actions workflow
(`.github/workflows/render.yaml`), which rebuilds it and re-renders the manuscript
on every push to `main`.

## Reproducing the simulation

Running `sim_mc.R` (Step 1) requires the following packages, all pinned in
`default.nix`. Parallelism uses base R's `parallel`.

| Package | Version | Citation |
|---|---|---|
| modsem | 1.0.22 (git `5c36547`) | Slupphaug, Mehmetoglu & Mittner (2025) |
| lavaan | 0.7.2.3207 (git `34f69c7`) | Rosseel, Jorgensen & De Wilde (2026) |
| covsim | 1.1.0 | Grønneberg, Foldnes & Marcoulides (2022) |
| rvinecopulib | 0.7.3.1.0 | Nagler & Vatter (2025) |

## Reproducing the results

Running `analyze_mc.R` and `figures.R` (Step 2) requires:

| Package | Version | Citation |
|---|---|---|
| dplyr | 1.2.1 | Wickham, François, Henry, Müller & Vaughan (2026) |
| simhelpers | 0.3.1 | Joshi & Pustejovsky (2025) |
| ggplot2 | 4.0.3 | Wickham (2016) |
| tidyr | 1.3.2 | Wickham, Vaughan & Girlich (2025) |

The environment itself is generated with `{rix}` 0.18.2 (Rodrigues & Baumann,
2026). Versions above are those resolved by the pinned environment under R 4.6.0.
