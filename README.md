# README

This repository contains the code, materials, and manuscript for the simulation
study **"Latent Variable Moderated Mediation Using Correct and Misspecified
Models"** by Felipe Fontana Vieira and Yves Rosseel.

The study compares four estimators of a latent moderated-mediation model under
varying sample size, reliability, exogenous-predictor distribution, and structural
misspecification:

- **LSAM** — local structural-after-measurement, `lavaan::sam(sam.method = "local", se = "local")`
- **LMS** — latent moderated structural equations, `modsem::modsem_da(method = "lms", robust.se = TRUE)`
- **UPI** — extended unconstrained product indicators (double mean-centering), `modsem::modsem(method = "dblcent", estimator = "MLR")`
- **FSR** — uncorrected factor score regression with Bartlett scores (Devlieger, Mayer & Rosseel, 2016; Ng & Chan, 2020): one CFA of the measurement model, Bartlett scores, product of the X and W scores, path model on the scores

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
definition. Two one-off prerequisites, the exogenous populations and the
misspecification calibration, have their outputs committed.

### 1. Simulation

```
.
├── gen_env.R          # generates default.nix via {rix} (pinned R + packages + Quarto + TeX)
├── default.nix        # generated environment definition (built with nix-shell)
├── calibration/
│   ├── copulas.R                   # builds the fixed (X, W) populations: one vine per distribution x copula family
│   ├── vines.rds                   # output of copulas.R, read by calibrate.R and sim_mc.R
│   ├── calibrate.R                 # calibrates misspecification coefficients & residual variances
│   ├── calibration_explanation.qmd # full derivation of the calibration scheme
│   └── calibration_results.rds     # output; copy to results_main/calibration.rds before Step 1
├── sim_mc.R           # main driver: builds the design grid and runs all conditions in parallel
└── helpers_mc.R       # data generation, the four estimators, MC/Wald inference, admissibility
                       #   checks, and the per-replication / per-condition workers (sourced by sim_mc.R)
```

The simulation writes its outputs to `results_main/`:

```
results_main/
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

`results_main/` is not committed. The large combined outputs are archived on Zenodo
(DOI above); download them into `results_main/` to re-run `analyze_mc.R`.

### 2. Analysis and figures

```
.
├── analyze_mc.R       # results_main/estimates_mc.rds → summary / rejection / convergence tables
└── figures.R          # results_main/{summary,rejection_imm_a3}.rds → plots/*.png (16 figures, 14 used in the paper)
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
Two prerequisites come first and their results are committed. `calibration/copulas.R`
builds the exogenous populations once: for every distribution and copula family it
solves the copula parameter so that Cor(X, W) = 0.3 exactly (the correlation is
computed by numerical integration over the copula, not by simulation), verifies the
saved object on an independent draw, and stores the vines in `calibration/vines.rds`.
Every replication of the simulation then draws from these fixed populations, so a
cell has one population rather than a fresh calibration per replication.
`calibration/calibrate.R` reads the same vines, computes the population moments of
X and W exactly by the same integration, derives the misspecification coefficients
and the Y-variance ratios from them in closed form, and writes
`calibration/calibration_results.rds`, which `sim_mc.R` reads from
`results_main/calibration.rds`, so copy it there first. Neither prerequisite
involves random draws, so both are deterministic. If `rho`, the degrees of
freedom or the margins change, rerun `copulas.R`, then `calibrate.R`, then the
simulation.

### Step 1: Run the simulation

```sh
nix-shell --pure --run 'OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 Rscript sim_mc.R'
```

`sim_mc.R` builds the design grid (`results_main/design.rds`), reads
`results_main/calibration.rds`, and runs every condition in parallel via
`parallel::mclapply` with reproducible `L'Ecuyer-CMRG` streams (master seed
`1234`). Each condition is saved to `results_main/condition_mc_<id>.rds` as it
completes (so the run is **resumable**), then all are combined into
`results_main/estimates_mc.rds` and `results_main/metrics_mc.rds`.

### Step 2: Process results and build figures

```r
Rscript analyze_mc.R    # → results_main/{summary,rejection_imm_a3,convergence_report}.rds
Rscript figures.R       # → plots/*.png
```

`analyze_mc.R` applies the two convergence criteria (soft / strict), removes
outliers, and computes performance metrics with Monte Carlo standard errors.
`figures.R` turns those tables into 16 figures, 14 of which appear in the paper.

### Robustness run: Gaussian copula for the nonnormal margins

In the main run the four nonnormal exogenous distributions use the Clayton
copula while the normal condition uses a Gaussian copula, so copula and margins
change together across the distribution factor.
The robustness run repeats the 864 nonnormal conditions with a Gaussian copula.
Design, seeds, and estimators are unchanged. A single flag at the top of
`calibration/calibrate.R`, `sim_mc.R`, and `analyze_mc.R` controls it:

```r
copula_robustness <- TRUE    # FALSE (the committed value) is the main run
```

With the flag set to `TRUE` the scripts use the Gaussian copula for the four
nonnormal margins, drop the normal condition from the design grid, and read and
write `results_robustness/` instead of `results_main/`. The commands are then
the same as for the main run:

```sh
# 1. recalibrate the misspecification coefficients under the gaussian copula.
#    needed: the c3 coefficient depends on Var(XW) and Cov(XW, Y), which differ
#    by ~20% between copulas for the t5 and same-skew chi-square margins.
#    writes calibration/calibration_results_robustness.rds
nix-shell --pure --run 'cd calibration && Rscript calibrate.R'
mkdir -p results_robustness
cp calibration/calibration_results_robustness.rds results_robustness/calibration.rds

# 2. simulate the four nonnormal distributions (864 conditions)
nix-shell --pure --run 'OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 Rscript sim_mc.R'

# 3. summarise into results_robustness/{summary,rejection_imm_a3,convergence_report}.rds
nix-shell --pure --run 'Rscript analyze_mc.R'
```

Set the flag back to `FALSE` afterwards so the committed scripts describe the
main run.

### Step 3: Render the manuscript

```sh
quarto render manuscript/article.qmd    # → manuscript/article.pdf
```

### Data flow

```
copulas.R ──> vines.rds ──> calibrate.R ──> calibration.rds
                  │                               │
                  └───┬───────────────────────────┘
helpers_mc.R ──┐      ▼
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
(nixpkgs date `2026-10-13`), which ships the CRAN releases of `lavaan` (0.7-3) and
of the other packages, with `modsem` built from the GitHub commit that is
source-identical to its CRAN 1.0.24 release (published after the snapshot date).
To regenerate the environment definition:

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
| modsem | 1.0.24 (git `c723119`, = CRAN 1.0.24) | Slupphaug, Mehmetoglu & Mittner (2025) |
| lavaan | 0.7-3 (CRAN) | Rosseel, Jorgensen & De Wilde (2026) |
| rvinecopulib | 1.0.0.1.0 | Nagler & Vatter (2025) |

## Reproducing the results

Running `analyze_mc.R` and `figures.R` (Step 2) requires:

| Package | Version | Citation |
|---|---|---|
| dplyr | 1.2.1 | Wickham, François, Henry, Müller & Vaughan (2026) |
| simhelpers | 0.3.1 | Joshi & Pustejovsky (2025) |
| ggplot2 | 4.0.3 | Wickham (2016) |
| tidyr | 1.3.2 | Wickham, Vaughan & Girlich (2025) |

The environment itself is generated with `{rix}` 0.18.2 (Rodrigues & Baumann,
2026). Versions above are those resolved by the pinned environment under R 4.6.1.
