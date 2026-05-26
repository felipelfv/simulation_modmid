# self-contained moderated-mediation simulation comparing three estimators.
# imm: MC quantile CI via lav_mvrnorm draws on vcov(fit)[c("a3","b"), .].
# other params: Wald CI from model vcov.
#
# lsam    : lavaan::sam(sam.method = "local", se = "local")
# lms     : modsem_da(method = "lms", robust.se = TRUE)
# dblcent : modsem(method = "dblcent", estimator = "MLR")
#
# 4 n x 3 a3 x 2 rel x 5 distr_exo x 9 misspec = 1080 cells.

library(modsem); library(lavaan)

master_seed <- 7000L; RNGkind("L'Ecuyer-CMRG"); set.seed(master_seed)

r_per_cell  <- 100L # (!)
m_mc        <- 20000L
n_cores     <- max(1L, parallel::detectCores(logical = FALSE) - 1L)
direc       <- "."
results_dir <- file.path(direc, "results")
dir.create(results_dir, showWarnings = FALSE)

# population and measurement model
a1 <- 0.4; a2 <- 0.3; b <- 0.5; cp <- 0.15
rho <- 0.3
loadings   <- c(1, 0.8, 0.7, 0.6)
rel_levels <- list(low = 0.5, high = 0.7)

res_m_var <- 0.7341
res_y_var <- 0.8214
var_m_emp <- c("0" = 1.0576, "0.2" = 1.0718, "0.4" = 1.1402)
var_y_emp <- c("0" = 1.1824, "0.2" = 1.1804, "0.4" = 1.1916)
var_m_at  <- function(a3v) unname(var_m_emp[as.character(a3v)])
var_y_at  <- function(a3v) unname(var_y_emp[as.character(a3v)])

err_sd_for <- function(target_rel, var_eta = 1) {
  loadings * sqrt(var_eta * (1 - target_rel) / target_rel)
}

pars <- c("a1", "a2", "a3", "b", "cp", "imm")

model_lv <- "
  X =~ x1 + x2 + x3 + x4
  W =~ w1 + w2 + w3 + w4
  M =~ m1 + m2 + m3 + m4
  Y =~ y1 + y2 + y3 + y4
  M ~ a1*X + a2*W + a3*X:W
  Y ~ b*M  + cp*X
  imm := a3 * b
"

# structural misspecifications
# note we always estimate the correctly-specified structural model, so
# (est - true_value) captures the misspec-induced bias.
# c2: true Y eq adds c2 * W (omitted W -> Y direct path)
# c3: true Y eq adds c3 * X*W (omitted X*W -> Y direct path)
# cl: indicator m3 also loads on X (omitted cross-loading)
# rcov: cor(d_m3, e_y3) = rcov (omitted residual covariance)
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

# design grid
grid <- expand.grid(
  n         = c(100, 200, 300, 500),
  a3        = c(0, 0.20, 0.40),
  rel       = c("low", "high"),
  distr_exo = c("normal", "uniform", "t5", "chisq_same", "chisq_diff"),
  misspec   = names(misspec_specs),
  stringsAsFactors = FALSE
)
saveRDS(grid, file.path(results_dir, "grid.rds"))

chi_df <- 2L; t_df <- 5L

# data generation 
gen_exo <- function(n, distr_exo) {
  if (distr_exo == "normal") {
    sigma  <- matrix(c(1, rho, rho, 1), 2, 2)
    margin <- list(distr = "norm", mean = 0, sd = 1)
    vd <- covsim::vita(list(margin, margin), sigma, verbose = FALSE, Nmax = 10^6)
    e  <- rvinecopulib::rvine(n, vine = vd)
    return(list(x = e[, 1], w = e[, 2]))
  }
  if (distr_exo == "uniform") {
    sigma  <- matrix(c(1, rho, rho, 1), 2, 2)
    margin <- list(distr = "unif", min = -sqrt(3), max = sqrt(3))
    vd <- covsim::vita(list(margin, margin), sigma, verbose = FALSE, Nmax = 10^6)
    e  <- rvinecopulib::rvine(n, vine = vd)
    return(list(x = e[, 1], w = e[, 2]))
  }
  if (distr_exo == "t5") {
    nat_var <- t_df / (t_df - 2)
    sigma   <- nat_var * matrix(c(1, rho, rho, 1), 2, 2)
    margin  <- list(distr = "t", df = t_df)
    vd <- covsim::vita(list(margin, margin), sigma, verbose = FALSE, Nmax = 10^6)
    e  <- rvinecopulib::rvine(n, vine = vd)
    return(list(x = e[, 1] / sqrt(nat_var), w = e[, 2] / sqrt(nat_var)))
  }
  if (distr_exo %in% c("chisq_same", "chisq_diff")) {
    rho_in  <- if (distr_exo == "chisq_same") rho else -rho
    nat_var <- 2 * chi_df
    sigma   <- nat_var * matrix(c(1, rho_in, rho_in, 1), 2, 2)
    margin  <- list(distr = "chisq", df = chi_df)
    vd <- covsim::vita(list(margin, margin), sigma, verbose = FALSE, Nmax = 10^6)
    e  <- rvinecopulib::rvine(n, vine = vd)
    x <- (e[, 1] - chi_df) / sqrt(nat_var)
    w <- (e[, 2] - chi_df) / sqrt(nat_var)
    if (distr_exo == "chisq_diff") w <- -w
    return(list(x = x, w = w))
  }
  stop("unknown distr_exo: ", distr_exo)
}

gen_data <- function(n, a3_true, rel_key, distr_exo, misspec_key = "none") {
  msp <- misspec_specs[[misspec_key]]
  exo <- gen_exo(n, distr_exo); x <- exo$x; w <- exo$w
  m   <- a1*x + a2*w + a3_true*x*w + rnorm(n, sd = sqrt(res_m_var))
  # y eq adds omitted c2*W and c3*X*W direct paths when misspec is non-zero
  y   <- b*m + cp*x + msp$c2*w + msp$c3*(x*w) + rnorm(n, sd = sqrt(res_y_var))
  
  target_rel <- rel_levels[[rel_key]]
  e_sd_xw    <- err_sd_for(target_rel, var_eta = 1)
  e_sd_m     <- err_sd_for(target_rel, var_eta = var_m_at(a3_true))
  e_sd_y     <- err_sd_for(target_rel, var_eta = var_y_at(a3_true))
  mk <- function(eta, prefix, e_sd) {
    z <- sapply(seq_along(loadings),
                function(j) loadings[j]*eta + rnorm(n, sd = e_sd[j]))
    colnames(z) <- paste0(prefix, seq_along(loadings))
    as.data.frame(z)
  }
  X_ind <- mk(x, "x", e_sd_xw)
  W_ind <- mk(w, "w", e_sd_xw)
  M_ind <- mk(m, "m", e_sd_m)
  Y_ind <- mk(y, "y", e_sd_y)
  
  # omitted cross-loading: m3 also loads on X with lambda = msp$cl
  if (msp$cl != 0) M_ind$m3 <- M_ind$m3 + msp$cl * x
  
  # omitted residual covariance: replace independent residuals for m3 and y3
  # with a correlated pair (corr = msp$rcov, marginal SDs preserved)
  if (msp$rcov != 0) {
    sd_m3 <- e_sd_m[3]; sd_y3 <- e_sd_y[3]
    Sigma <- matrix(c(sd_m3^2,             msp$rcov*sd_m3*sd_y3,
                      msp$rcov*sd_m3*sd_y3, sd_y3^2), 2, 2) 
    # we want correlation, not covariance (!)
    # otherwise the meaning would change with reliability level
    # at low rel the residuals are bigger, so the same covariance number 
    # would imply a different correlation than at high rel.
    pair  <- lavaan:::lav_mvrnorm(n, mu = c(0, 0), sigma_1 = Sigma)
    M_ind$m3 <- loadings[3]*m + pair[, 1]
    Y_ind$y3 <- loadings[3]*y + pair[, 2]
  }
  
  cbind(X_ind, W_ind, M_ind, Y_ind)
}

# mc + wald inference from coef + vcov
mc_inference <- function(co, V, M = m_mc, alpha = 0.05) {
  na6 <- setNames(rep(NA_real_, length(pars)), pars)
  if (!is.matrix(V) ||
      !all(c("a1","a2","a3","b","cp") %in% rownames(V)) ||
      !all(c("a1","a2","a3","b","cp") %in% names(co)))
    return(list(est = na6, se = na6, lo = na6, hi = na6))
  
  pn        <- c("a1","a2","a3","b","cp")
  z         <- qnorm(1 - alpha/2)
  se_struct <- sqrt(diag(V)[pn])
  lo_struct <- co[pn] - z * se_struct
  hi_struct <- co[pn] + z * se_struct
  
  V_ab     <- V[c("a3","b"), c("a3","b")]
  sims     <- lavaan:::lav_mvrnorm(n = M, mu = c(co["a3"], co["b"]),
                                   sigma_1 = V_ab)
  imm_sims <- sims[, 1] * sims[, 2]
  imm_se   <- sd(imm_sims)
  imm_lo   <- as.numeric(quantile(imm_sims, alpha/2))
  imm_hi   <- as.numeric(quantile(imm_sims, 1 - alpha/2))
  
  est <- c(co[pn], imm = unname(co["a3"] * co["b"]))
  se  <- c(se_struct, imm = imm_se)
  lo  <- c(lo_struct, imm = imm_lo)
  hi  <- c(hi_struct, imm = imm_hi)
  
  list(est = est[pars], se = se[pars], lo = lo[pars], hi = hi[pars])
}

# admissibility helpers (heywood cases)
# neg_theta: any residual variance (diag of theta) <= 0
# npd_psi: factor covariance matrix psi has any non-positive eigenvalue
neg_theta_diag <- function(theta_diag_vec) {
  if (!length(theta_diag_vec)) return(NA)
  any(theta_diag_vec <= 0, na.rm = TRUE)
}
npd_psi_mat <- function(psi_mat) {
  if (!length(psi_mat) || !nrow(psi_mat)) return(NA)
  if (any(!is.finite(psi_mat))) return(TRUE)
  min(eigen(psi_mat, symmetric = TRUE, only.values = TRUE)$values) <= 0
}
admis_from_pe <- function(pe, latents) {
  pe <- as.data.frame(pe)
  vv <- pe[pe$op == "~~", , drop = FALSE]
  th <- vv[vv$lhs == vv$rhs & !(vv$lhs %in% latents), , drop = FALSE]
  ps <- vv[vv$lhs %in% latents & vv$rhs %in% latents, , drop = FALSE]
  L_names <- intersect(latents, unique(c(ps$lhs, ps$rhs)))
  psi <- matrix(0, length(L_names), length(L_names),
                dimnames = list(L_names, L_names))
  for (k in seq_len(nrow(ps))) {
    psi[ps$lhs[k], ps$rhs[k]] <- ps$est[k]
    psi[ps$rhs[k], ps$lhs[k]] <- ps$est[k]
  }
  list(neg_theta = neg_theta_diag(th$est),
       npd_psi   = npd_psi_mat(psi))
}
safe_admis <- function(expr) {
  res <- try(expr, silent = TRUE)
  if (inherits(res, "try-error") || !is.list(res) ||
      !all(c("neg_theta", "npd_psi") %in% names(res)))
    return(list(neg_theta = NA, npd_psi = NA))
  res
}

# per-method fitters: return list(co, V, neg_theta, npd_psi) or NULL on failure
fit_lsam <- function(data) {
  fit <- try(lavaan::sam(model_lv, data = data,
                         sam.method = "local", se = "local"),
             silent = TRUE)
  if (inherits(fit, "try-error")) return(NULL)
  V  <- try(vcov(fit), silent = TRUE)
  co <- try(coef(fit), silent = TRUE)
  pe <- try(lavaan::parameterEstimates(fit), silent = TRUE)
  if (inherits(V, "try-error") || inherits(co, "try-error")) return(NULL)
  adm <- if (inherits(pe, "try-error")) list(neg_theta = NA, npd_psi = NA)
  else safe_admis(admis_from_pe(pe, latents = c("X","W","M","Y")))
  list(co = co, V = V, neg_theta = adm$neg_theta, npd_psi = adm$npd_psi)
}

fit_lms <- function(data) {
  fit <- try(modsem_da(model_lv, data = data, method = "lms",
                       robust.se = TRUE),
             silent = TRUE)
  if (inherits(fit, "try-error")) return(NULL)
  V  <- try(vcov(fit), silent = TRUE)
  co <- try(coef(fit), silent = TRUE)
  pe <- try(parameter_estimates(fit), silent = TRUE)
  if (inherits(V, "try-error") || inherits(co, "try-error")) return(NULL)
  adm <- if (inherits(pe, "try-error")) list(neg_theta = NA, npd_psi = NA)
  else safe_admis(admis_from_pe(pe, latents = c("X","W","M","Y")))
  list(co = co, V = V, neg_theta = adm$neg_theta, npd_psi = adm$npd_psi)
}

fit_dblcent <- function(data) {
  fit <- try(modsem(model_lv, data = data, method = "dblcent",
                    estimator = "MLR"),
             silent = TRUE)
  if (inherits(fit, "try-error")) return(NULL)
  V  <- try(vcov(fit), silent = TRUE)
  co <- try(coef(fit), silent = TRUE)
  pe <- try(parameter_estimates(fit), silent = TRUE)
  if (inherits(V, "try-error") || inherits(co, "try-error")) return(NULL)
  # dblcent introduces an extra XW latent from product indicators
  adm <- if (inherits(pe, "try-error")) list(neg_theta = NA, npd_psi = NA)
  else safe_admis(admis_from_pe(pe, latents = c("X","W","XW","M","Y")))
  list(co = co, V = V, neg_theta = adm$neg_theta, npd_psi = adm$npd_psi)
}

fitters <- list(lsam    = fit_lsam,
                lms     = fit_lms,
                dblcent = fit_dblcent)
methods <- names(fitters)

# per-rep wrapper
true_value_for <- function(a3v) setNames(c(a1, a2, a3v, b, cp, b * a3v), pars)

work_for_cell <- function(i, r) {
  # capture l'ecuyer substream before any draw, for per-rep replay
  rng_stream <- paste(get(".Random.seed", envir = .GlobalEnv), collapse = ":")
  
  cell  <- grid[i, ]
  dat   <- gen_data(cell$n, cell$a3, cell$rel, cell$distr_exo, cell$misspec)
  truth <- true_value_for(cell$a3)
  na6   <- setNames(rep(NA_real_, length(pars)), pars)
  
  # helper: run expr while collecting any warnings into a single string
  with_warnings <- function(expr) {
    warns <- character(0)
    out <- withCallingHandlers(expr,
      warning = function(w) {
        warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning")
      })
    list(value = out,
         warnings = if (length(warns)) paste(unique(warns), collapse = " | ")
                    else NA_character_)
  }

  rows_all <- lapply(methods, function(m) {
    ww <- with_warnings({
      fitted <- fitters[[m]](dat)
      list(fitted = fitted,
           res = if (is.null(fitted)) list(est = na6, se = na6, lo = na6, hi = na6)
                 else                 mc_inference(fitted$co, fitted$V))
    })
    fitted <- ww$value$fitted
    res    <- ww$value$res
    adm    <- if (is.null(fitted)) list(neg_theta = NA, npd_psi = NA)
              else list(neg_theta = fitted$neg_theta, npd_psi = fitted$npd_psi)
    data.frame(parameter = pars,
               est       = as.numeric(res$est),
               se        = unname(res$se),
               ci_lo     = unname(res$lo),
               ci_hi     = unname(res$hi),
               neg_theta = adm$neg_theta,
               npd_psi   = adm$npd_psi,
               warnings  = ww$warnings,
               method    = m,
               stringsAsFactors = FALSE)
  })
  rows <- do.call(rbind, rows_all)

  rows$cell       <- sprintf("c%04d", i)
  rows$n          <- cell$n
  rows$a3         <- cell$a3
  rows$rel        <- cell$rel
  rows$distr_exo  <- cell$distr_exo
  rows$misspec    <- cell$misspec
  rows$rep        <- r
  rows$true_value <- truth[rows$parameter]
  rows$rng_stream <- rng_stream
  rows[, c("cell","n","a3","rel","distr_exo","misspec","rep","method","parameter",
           "true_value","est","se","ci_lo","ci_hi","neg_theta","npd_psi",
           "warnings","rng_stream")]
}

# parallel run over cells
run_cell <- function(i) {
  cell_file <- file.path(results_dir, sprintf("cell_mc_%04d.rds", i))
  if (file.exists(cell_file)) return(invisible(NULL))
  
  t0        <- Sys.time()
  rows_list <- lapply(seq_len(r_per_cell),
                      function(r) try(work_for_cell(i, r), silent = TRUE))
  ok        <- vapply(rows_list, is.data.frame, logical(1))
  cell_raw  <- do.call(rbind, rows_list[ok])
  # guard!: don't write a NULL/empty file (would lock the cell out of future
  # resume-skip retries even though it has no usable data)
  if (!is.null(cell_raw) && nrow(cell_raw) > 0)
    saveRDS(cell_raw, cell_file)
  
  dt <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  cat(sprintf("c%04d  %.1fs  %d/%d\n", i, dt, sum(ok), r_per_cell),
      file = file.path(results_dir, "progress.log"), append = TRUE)
  invisible(NULL)
}

t_start <- Sys.time()
parallel::mclapply(seq_len(nrow(grid)), run_cell,
                   mc.cores = n_cores, mc.preschedule = FALSE,
                   mc.set.seed = TRUE)
message(sprintf("done in %.1f min",
                as.numeric(difftime(Sys.time(), t_start, units = "mins"))))

# aggregate raw cell files into raw_mc.rds
raw <- do.call(rbind, lapply(
  list.files(results_dir, "^cell_mc_\\d+\\.rds$", full.names = TRUE),
  readRDS))
saveRDS(raw, file.path(results_dir, "raw_mc.rds"))