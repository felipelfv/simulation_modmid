# code for the simulation: 
# moderated-mediation simulation comparing three estimators.

# lsam    : lavaan::sam(sam.method = "local", se = "local")
# lms     : modsem_da(method = "lms", robust.se = TRUE)
# dblcent : modsem(method = "dblcent", estimator = "MLR")

# 4 n x 3 a3 x 2 rel x 5 distr_exo x 9 misspec = 1080 conditions

# note to myself: probably keep model fit values also. that implies retaining 
# other values in the final output
# note to myself: if checking cfa indices for the first step as in brandt et al. 2020, 
# then we should have 1 block for all latent variables

library(modsem); library(lavaan)
library(covsim); library(rvinecopulib)

master_seed <- 1234L; RNGkind("L'Ecuyer-CMRG"); set.seed(master_seed)

r_per_condition  <- 1000L # (!) 
m_mc <- 20000L # justification for the indirect effect mc (?)
n_cores <- max(1L, parallel::detectCores(logical = FALSE) - 1L)
direc <- "."
results_dir <- file.path(direc, "results")
dir.create(results_dir, showWarnings = FALSE)

# population and measurement model
a1 <- 0.4; a2 <- 0.3; b <- 0.5; cp <- 0.15
rho <- 0.3
loadings <- c(1, 0.8, 0.7, 0.6)
rel_levels <- list(low = 0.5, high = 0.7)

res_m_var <- 0.8531
res_y_var <- 0.9350
var_m_at <- function(a3v) a1^2 + a2^2 + 2*a1*a2*rho + a3v^2*(1 + rho^2) + res_m_var
var_y_at <- function(a3v) b^2*var_m_at(a3v) + cp^2 + 2*b*cp*(a1 + a2*rho) + res_y_var

err_sd_for <- function(target_rel, var_eta = 1) {
  loadings * sqrt(var_eta * (1 - target_rel) / target_rel)
}

# 12 free loadings (first indicator of each factor is fixed to 1 for identification);
# true free-loading values are 0.8 / 0.7 / 0.6 per factor (loadings[2:4]).
load_names <- c("lx2","lx3","lx4","lw2","lw3","lw4","lm2","lm3","lm4","ly2","ly3","ly4")
load_true  <- setNames(rep(loadings[2:4], 4), load_names)
pars <- c("a1", "a2", "a3", "b", "cp", "imm", load_names)

model_lv <- "
  X =~ x1 + lx2*x2 + lx3*x3 + lx4*x4
  W =~ w1 + lw2*w2 + lw3*w3 + lw4*w4
  M =~ m1 + lm2*m2 + lm3*m3 + lm4*m4
  Y =~ y1 + ly2*y2 + ly3*y3 + ly4*y4
  M ~ a1*X + a2*W + a3*X:W
  Y ~ b*M  + cp*X
  imm := a3 * b
"

# structural misspecifications
# note: we always estimate the correctly-specified structural model, so
# (est - true_value) captures the misspec-induced bias.
# c2: true Y eq adds c2 * W (omitted W -> Y direct path)
# c3: true Y eq adds c3 * X*W (omitted X*W -> Y direct path)
# cl: indicator m3 also loads on X (omitted cross-loading)
# rcov: cor(d_m3, e_y3) = rcov (omitted residual covariance)
# note the c2 and c3 nomenclature comes from X
misspec_specs <- list(
  none             = list(c2 = 0,    c3 = 0,    cl = 0,    rcov = 0),
  c2_small         = list(c2 = 0.25, c3 = 0,    cl = 0,    rcov = 0),
  c2_medium        = list(c2 = 0.50, c3 = 0,    cl = 0,    rcov = 0),
  c3_small         = list(c2 = 0,    c3 = 0.25, cl = 0,    rcov = 0),
  c3_medium        = list(c2 = 0,    c3 = 0.50, cl = 0,    rcov = 0),
  crossload_small  = list(c2 = 0,    c3 = 0,    cl = 0.25, rcov = 0),
  crossload_medium = list(c2 = 0,    c3 = 0,    cl = 0.50, rcov = 0),
  rescov_small     = list(c2 = 0,    c3 = 0,    cl = 0,    rcov = 0.25),
  rescov_medium    = list(c2 = 0,    c3 = 0,    cl = 0,    rcov = 0.50)
)

# design grid (full factorial)
design <- expand.grid(
  n = c(100, 200, 300, 500),
  a3 = c(0, 0.20, 0.40),
  rel = c("low", "high"),
  distr_exo = c("normal", "uniform", "t5", "chisq_same", "chisq_diff"),
  misspec = names(misspec_specs),
  stringsAsFactors = FALSE
)
saveRDS(design, file.path(results_dir, "design.rds"))

chi_df <- 2L; t_df <- 5L

# calibration (from calibrate_final.R): misspec_coefs = beta coefs per cell (pure lookup);
# vy_ratios = var(Y) inflation ratios to hold Y reliability at its `none` value under c2/c3.
cal           <- readRDS(file.path(results_dir, "calibration.rds"))
misspec_coefs <- cal$misspec_coefs
vy_ratios     <- cal$vy_ratios

# machinery: data generation, inference, admissibility checks, estimators,
# and the per-rep / per-condition workers. sourced after all config above so
# its functions can read these globals at call time.
source(file.path(direc, "helpers_mc.R"))

t_start <- Sys.time()
# TRUE (default, what we use) for mc.set.seed: 
# each child gets a different rng state. The exact mechanism depends on RNGkind():
# (i) with the default Mersenne-Twister: each child seeds itself from Sys.time() + Sys.getpid().
# should be "different enough" but not formally independent
# (ii) with RNGkind("L'Ecuyer-CMRG") 
# (what we set at the top of sim_mc.R): each child calls nextRNGStream() to advance to a
# guaranteed-independent substream (substreams are 2^127 draws apart and provably non-overlapping)
parallel::mclapply(seq_len(nrow(design)), run_condition,
                   mc.cores = n_cores, mc.preschedule = FALSE,
                   mc.set.seed = TRUE)

message(sprintf("done in %.1f min",
                as.numeric(difftime(Sys.time(), t_start, units = "mins"))))

# aggregate per-condition files into estimates_mc.rds (per-fit estimates) and
# metrics_mc.rds (per-rep dataset metrics). uncomment to run.
all_files <- list.files(results_dir, "^condition_mc_\\d+\\.rds$", full.names = TRUE)
all_data  <- lapply(all_files, readRDS)
estimates <- do.call(rbind, lapply(all_data, `[[`, "rows"))
metrics   <- do.call(rbind, lapply(all_data, `[[`, "metrics"))
saveRDS(estimates, file.path(results_dir, "estimates_mc.rds"))
saveRDS(metrics,   file.path(results_dir, "metrics_mc.rds"))