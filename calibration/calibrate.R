# calibration of the misspecification sizes: for every (distribution x a3) cell,
# the raw coefficient that gives each omitted path a standardized size of 0.25
# ("small") or 0.50 ("medium"), plus the var(Y) inflation ratios that gen_data()
# uses to hold Y reliability at its baseline under c2 / c3. the derivation of the
# standardized-size definition is in calibration_explanation.qmd.

# the population moments are computed by numerical integration of the saved
# populations (vines.rds, built once by copulas.R with vita). the three mixed
# moments of the standardized exogenous variables, E[X^2 W], E[X W^2] and
# E[X^2 W^2], are integrals over the unit square of the quantile functions
# against the pair-copula density, evaluated by nested adaptive quadrature.
# every variance and covariance of M and Y then follows from the structural
# equations. no simulation is involved here, so for a given vines.rds the
# result is deterministic. the normal cell has closed forms, which the script
# checks against the gaussian parameter vita stored in that vine.

library(rvinecopulib)

# copula robustness run (see README). FALSE calibrates the main run (clayton
# populations), TRUE the gaussian-copula populations, each to its own output file.
copula_robustness <- FALSE
vines      <- readRDS("vines.rds")$vines
exo_family <- if (copula_robustness) "gaussian" else "clayton"
out_file   <- if (copula_robustness) "calibration_results_robustness.rds" else "calibration_results.rds"

a1 <- 0.4; a2 <- 0.3; b <- 0.5; cp <- 0.15
rho <- 0.3
lambda3 <- 0.7                     # loading of m3, the indicator that carries the cross-loading
rel_levels <- c(low = 0.5, high = 0.7)
res_m_var <- 0.8531                # disturbance psi_M
res_y_var <- 0.9350                # disturbance psi_Y
chi_df <- 2L; t_df <- 5L
sizes <- c(0.25, 0.50)             # standardized misspecification sizes, small and medium
quad_tol <- 1e-9                   # relative tolerance of each one-dimensional integration

# error variance of m3 as gen_data() sets it: from the normal-theory var(M),
# calibrated under normality and held fixed across distributions. the
# cross-loading calibration must use the same value.
var_m_normal <- function(a3) a1^2 + a2^2 + 2*a1*a2*rho + a3^2*(1 + rho^2) + res_m_var
theta3_for   <- function(rel, a3) lambda3^2 * var_m_normal(a3) * (1 - rel) / rel

# nutshell explanation (details in calibration_explanation.qmd):
# we add an omitted predictor p to a baseline outcome y0 with a raw coefficient
# beta, producing y0 + beta p. we do not specify beta directly, we specify how big
# the omission is in standardized terms,
# t = beta·sd(p) / sd(y0 + beta p),
# because t is what is interpretable across cells ("small" 0.25 or "medium" 0.50)
# and beta is not. so we fix t and solve for the beta that produces it.
# mom = list(Vp = Var(p), C = Cov(p, y0), V0 = Var(y0)) describes the clean
# baseline, before the misspecification term is added.
solve_std <- function(t, mom) {
  A <- mom$Vp*(1-t^2); B <- -2*t^2*mom$C; D <- -t^2*mom$V0
  (-B + sqrt(B^2 - 4*A*D)) / (2*A)
}
solve_sizes <- function(mom) setNames(sapply(sizes, solve_std, mom = mom), as.character(sizes))

# exact moments of the standardized exogenous variables
# X = qx(U), W = qw(V) with (U, V) from the cell's copula. the quantile functions
# are standardized, so E[X] = E[W] = 0 and E[X^2] = E[W^2] = 1. the chisq_diff
# cell negates W, exactly as gen_exo() does.
qstd <- list(
  normal  = function(u) qnorm(u),
  uniform = function(u) qunif(u, -sqrt(3), sqrt(3)),
  t5      = function(u) qt(u, t_df) / sqrt(t_df / (t_df - 2)),
  chisq   = function(u) (qchisq(u, chi_df) - chi_df) / sqrt(2 * chi_df))

moment <- function(bc, qx, qw, px, pw) {
  inner <- function(u) vapply(u, function(ui) integrate(function(v)
      qx(ui)^px * qw(v)^pw * dbicop(cbind(rep(ui, length(v)), v), bc),
      0, 1, rel.tol = quad_tol, subdivisions = 2000L, stop.on.error = FALSE)$value, numeric(1))
  integrate(inner, 0, 1, rel.tol = quad_tol, subdivisions = 2000L, stop.on.error = FALSE)$value
}

exo_moments <- function(cell) {
  vd <- if (cell == "normal") vines$normal else vines[[exo_family]][[cell]]
  bc <- vd$copula$pair_copulas[[1]][[1]]
  qx <- qstd[[sub("chisq_.*", "chisq", cell)]]
  qw <- if (cell == "chisq_diff") function(u) -qx(u) else qx
  c(xw = moment(bc, qx, qw, 1, 1), x2w = moment(bc, qx, qw, 2, 1),
    xw2 = moment(bc, qx, qw, 1, 2), x2w2 = moment(bc, qx, qw, 2, 2))
}

# calibration per cell 
# with P = XW:  M = a1 X + a2 W + a3 P + zeta_M,  Y0 = b M + cp X + eps_Y
misspec_coefs <- list(); vy_ratios <- list(); moments <- list()
for (de in c("normal", "uniform", "t5", "chisq_same", "chisq_diff")) {
  mo <- exo_moments(de)
  if (de == "normal") { # closed forms under bivariate normality check the integration
    r_vine <- vines$normal$copula$pair_copulas[[1]][[1]]$parameters[1]
    stopifnot(abs(mo[["xw"]] - r_vine) < 1e-6, abs(mo[["x2w"]]) < 1e-6,
              abs(mo[["xw2"]]) < 1e-6, abs(mo[["x2w2"]] - (1 + 2*r_vine^2)) < 1e-6)
  }
  r <- mo[["xw"]]; var_p <- mo[["x2w2"]] - r^2; cov_xp <- mo[["x2w"]]; cov_wp <- mo[["xw2"]]
  moments[[de]] <- c(mo, var_p = var_p)
  cat(sprintf("%-10s E[XW] %.6f  E[X2W] %+.5f  E[XW2] %+.5f  E[X2W2] %.5f\n",
              de, r, cov_xp, cov_wp, mo[["x2w2"]]))
  for (a3 in c(0, 0.20, 0.40)) {
    var_m  <- a1^2 + a2^2 + a3^2*var_p + 2*a1*a2*r + 2*a1*a3*cov_xp + 2*a2*a3*cov_wp + res_m_var
    cov_xm <- a1 + a2*r + a3*cov_xp
    cov_wm <- a1*r + a2 + a3*cov_wp
    cov_pm <- a1*cov_xp + a2*cov_wp + a3*var_p
    var_y0 <- b^2*var_m + cp^2 + 2*b*cp*cov_xm + res_y_var
    cov_wy <- b*cov_wm + cp*r
    cov_py <- b*cov_pm + cp*cov_xp
    c2 <- solve_sizes(list(Vp = 1,     C = cov_wy, V0 = var_y0))
    c3 <- solve_sizes(list(Vp = var_p, C = cov_py, V0 = var_y0))
    # var(Y0 + beta p) / var(Y0): the inflation that gen_data() cancels in the
    # Y indicator error variance so Y reliability stays at its `none` value
    vy_ratios[[de]][[as.character(a3)]] <- list(
      c2 = sapply(c2, function(beta) (var_y0 + 2*beta*cov_wy + beta^2) / var_y0),
      c3 = sapply(c3, function(beta) (var_y0 + 2*beta*cov_py + beta^2*var_p) / var_y0))
    # cross-loading of m3 on X: m3 = lambda3 M + delta3. the only coefficient
    # keyed by reliability level, through theta3.
    cl <- lapply(rel_levels, function(rel)
      solve_sizes(list(Vp = 1, C = lambda3*cov_xm, V0 = lambda3^2*var_m + theta3_for(rel, a3))))
    misspec_coefs[[de]][[as.character(a3)]] <- list(c2 = c2, c3 = c3, cl = cl)
  }
}

saveRDS(list(misspec_coefs = misspec_coefs, vy_ratios = vy_ratios,
             moments = moments, exo_family = exo_family), out_file)