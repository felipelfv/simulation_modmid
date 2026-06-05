library(covsim); library(rvinecopulib)

set.seed(20240501L) # clarification to myself: seed does L internally
n_cal <- 1e6L
results_dir <- "results"; dir.create(results_dir, showWarnings = FALSE)

a1 <- 0.4; a2 <- 0.3; b <- 0.5; cp <- 0.15 
rho <- 0.3
loadings <- c(1, 0.8, 0.7, 0.6)
rel_levels <- list(low = 0.5, high = 0.7)

res_m_var <- 0.8531 # disturbance psi_M
res_y_var <- 0.9350 # disturbance psi_Y

# var_m_at(a3v) is literally var(S_M) + psi_M
# the signal variance var(S_M) = a1^2 + a2^2 + 2*a1*a2*rho + a3v^2*(1+rho^2) 
# R^2_M is 0.30 only at that reference (0.2) and is allowed to drift as a_3
# and the distribution change
var_m_at <- function(a3v) { # a3v to avoid clashing with a3 (0,0.2,0.4)
  a1^2 + a2^2 + 2*a1*a2*rho + a3v^2*(1 + rho^2) + res_m_var
}

# reliability of indicator j:
# w_j = signal / (signal + noise) = lambda_j^2·V / (lambda_j^2·V + theta_j)
# var_eta = 1 as default but not always (!): 
# mediator M and outcome Y are not standardized
err_sd_for <- function(target_rel, var_eta = 1) { 
  loadings * sqrt(var_eta * (1 - target_rel) / target_rel)
}

# detailed explained behind this function/formula is given in the file:
# calibration_explanation.qmd
targets <- c(0.25, 0.50)

# nutshell explanation:
# we want to add an omitted predictor p to a baseline outcome y_o with some 
# raw coefficient beta, producing y_o + beta_p. But we dont want to specify beta 
# directly, we want to specify how big the omission is in standardized terms:
# t = beta·sd(p) / sd(y_o + beta_p)  
# thus,  t is what's interpretable across cells ("small" 0.25 or "medium" 0.50), 
# beta is not. the same beta means wildly different things depending on the variances.
# So we fix t and solve for the beta that produces it. 
# That inversion is solve_std
solve_std <- function(t, m) if (t == 0) 0 else {
  A <- m$Vp*(1-t^2); B <- -2*t^2*m$C; D <- -t^2*m$V0
  (-B + sqrt(B^2 - 4*A*D)) / (2*A)
}

# mom is a named list of three moments describing one misspecification in one cell:
# mom = list(Vp = <Var(p)>, C = <Cov(p, y0)>, V0 = <Var(y0)>)
# hence, the mom values are computed on the clean baseline, before the 
# misspecification term is added
solve_targets <- function(mom) { 
  setNames(sapply(targets, solve_std, m = mom), as.character(targets)) 
} # note to myself: as.character drops the 0 in 0.50 (i.e., 0.5)

chi_df <- 2L; t_df <- 5L

# data generation choices (at the code level) are explained in the main document 
# see sim_mc.R for more information
gen_exo <- function(n, distr_exo) {
  if (distr_exo == "normal") {
    sigma <- matrix(c(1, rho, rho, 1), 2, 2)
    margin <- list(distr = "norm", mean = 0, sd = 1)
    vd <- covsim::vita(list(margin, margin), sigma, verbose = FALSE, 
                       Nmax = 10^6, family_set = "gauss", cores = 1)
    e <- rvinecopulib::rvine(n, vine = vd, cores = 1)
    return(list(x = e[, 1], w = e[, 2]))
  }
  if (distr_exo == "uniform") {
    sigma <- matrix(c(1, rho, rho, 1), 2, 2)
    margin <- list(distr = "unif", min = -sqrt(3), max = sqrt(3))
    vd <- covsim::vita(list(margin, margin), sigma, verbose = FALSE, 
                       Nmax = 10^6, cores = 1)
    e <- rvinecopulib::rvine(n, vine = vd, cores = 1)
    return(list(x = e[, 1], w = e[, 2]))
  }
  if (distr_exo == "t5") {
    nat_var <- t_df / (t_df - 2)
    sigma <- nat_var * matrix(c(1, rho, rho, 1), 2, 2)
    margin <- list(distr = "t", df = t_df)
    vd <- covsim::vita(list(margin, margin), sigma, verbose = FALSE, 
                       Nmax = 10^6, cores = 1)
    e <- rvinecopulib::rvine(n, vine = vd, cores = 1)
    return(list(x = e[, 1] / sqrt(nat_var), w = e[, 2] / sqrt(nat_var)))
  }
  if (distr_exo %in% c("chisq_same", "chisq_diff")) {
    rho_in <- if (distr_exo == "chisq_same") rho else -rho
    nat_var <- 2 * chi_df
    sigma <- nat_var * matrix(c(1, rho_in, rho_in, 1), 2, 2)
    margin <- list(distr = "chisq", df = chi_df)
    vd <- covsim::vita(list(margin, margin), sigma, verbose = FALSE, 
                       Nmax = 10^6, cores = 1)
    e <- rvinecopulib::rvine(n, vine = vd, cores = 1)
    x <- (e[, 1] - chi_df) / sqrt(nat_var)
    w <- (e[, 2] - chi_df) / sqrt(nat_var)
    if (distr_exo == "chisq_diff") w <- -w
    return(list(x = x, w = w))
  }
  stop("not specified: ", distr_exo)
}

# might change to a data.frame or something similar as currently things are 
# arguably disorganized (05/06/26)
misspec_coefs <- list(); vy_ratios <- list()
for (de in c("normal", "uniform", "t5", "chisq_same", "chisq_diff")) {
  for (a3 in c(0, 0.20, 0.40)) {
    e <- gen_exo(n_cal, de); x <- e$x; w <- e$w; xw <- x * w
    m  <- a1*x + a2*w + a3*xw + rnorm(n_cal, sd = sqrt(res_m_var))
    y0 <- b*m + cp*x+ rnorm(n_cal, sd = sqrt(res_y_var))
    rc2 <- solve_targets(list(Vp = var(w),  C = cov(w,  y0), V0 = var(y0)))
    rc3 <- solve_targets(list(Vp = var(xw), C = cov(xw, y0), V0 = var(y0)))
    # For each beta in rc2, compute k = Var(y0 + βW)/Var(y0). this is the factor by which
    # adding the omitted path inflates Y's variance. Likewise for rc3 with XW.
    vy_ratios[[de]][[as.character(a3)]] <- list(
      c2 = sapply(rc2, function(rc) var(y0 + rc * w)  / var(y0)),
      c3 = sapply(rc3, function(rc) var(y0 + rc * xw) / var(y0)))
    # adding c2/c3 raises Var(Y), which would raise Y's indicator reliability 
    # as a side effect. The simulation multiplies the Y indicator noise 
    # by k to cancel that, holding Y reliability at its no-misspecification value.
    # (k is reliability-independent, so it's stored once per cell; like the comment
    # says)
    
    # cl is the only rel-dependent coefficient, so it is the only one keyed by
    # reliability level; c2/c3 (rel-independent) are stored once, like vy_ratios.
    rcl <- lapply(rel_levels, function(rel) {
      m3 <- loadings[3]*m + rnorm(n_cal, sd = err_sd_for(rel, var_m_at(a3))[3])
      solve_targets(list(Vp = var(x), C = cov(x, m3), V0 = var(m3)))
    })
    misspec_coefs[[de]][[as.character(a3)]] <- list(c2 = rc2, c3 = rc3, cl = rcl)
  }
}

saveRDS(list(misspec_coefs = misspec_coefs, vy_ratios = vy_ratios), 
        file.path(results_dir, "calibration.rds"))
