# code for the simulation: 
# moderated-mediation simulation comparing three estimators.
# 
# lsam    : lavaan::sam(sam.method = "local", se = "local")
# lms     : modsem_da(method = "lms", robust.se = TRUE)
# dblcent : modsem(method = "dblcent", estimator = "MLR")
#
# 4 n x 3 a3 x 2 rel x 5 distr_exo x 9 misspec = 1080 conditions

# note to myself: write in the first article taht all have non normal copulas
# althought the calibration was done empirically under the non normal copula with marginal gaussians
# in this paper, we should then calibrate with everything normal (maginal and copula)

# note to myself: probably keep model fit values also. that implies retaining other values in the final output
# note to myself: also keep track of the r2 and reliability under each dataset generated

library(modsem); library(lavaan)
library(covsim); library(rvinecopulib)

master_seed <- 1234L; RNGkind("L'Ecuyer-CMRG"); set.seed(master_seed)

r_per_condition  <- 100L # (!)
m_mc <- 20000L
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
var_m_emp <- c("0" = 1.1746, "0.2" = 1.2188, "0.4" = 1.3498)
var_y_emp <- c("0" = 1.3249, "0.2" = 1.3350, "0.4" = 1.3682)
var_m_at <- function(a3v) unname(var_m_emp[as.character(a3v)])
var_y_at <- function(a3v) unname(var_y_emp[as.character(a3v)])

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

# data generation
# every distr_exo (normal, uniform, t5, chisq) is standardized to
# mean = 0 and variance = 1, so that differences across branches come from
# distribution shape, not from scale or location.
# normal condition uses a forced Gaussian copula (clean reference); the four
# non-normal conditions use vita's default (clayton-first) — so they vary the
# marginals while keeping a consistent non-normal joint dependence.
gen_exo <- function(n, distr_exo) {
  if (distr_exo == "normal") {
    # truly bivariate normal baseline: force gaussian copula
    sigma <- matrix(c(1, rho, rho, 1), 2, 2)
    margin <- list(distr = "norm", mean = 0, sd = 1)
    vd <- covsim::vita(list(margin, margin), sigma, verbose = FALSE, Nmax = 10^6,
                       family_set = "gauss") 
    e <- rvinecopulib::rvine(n, vine = vd)
    return(list(x = e[, 1], w = e[, 2]))
  }
  if (distr_exo == "uniform") {
    sigma <- matrix(c(1, rho, rho, 1), 2, 2)
    # note: var(U(a,b)) = (b-a)^2 / 12. we have a symmetry constrain: a = -b
    # this for the mean 0, var 1
    # thus U{-sqrt(3),sqrt(3)} has var = 1
    margin <- list(distr = "unif", min = -sqrt(3), max = sqrt(3)) # bounded -+ sqrt(3)
    vd <- covsim::vita(list(margin, margin), sigma, verbose = FALSE, Nmax = 10^6)
    e <- rvinecopulib::rvine(n, vine = vd)
    return(list(x = e[, 1], w = e[, 2]))
  }
  if (distr_exo == "t5") {
    # t is symmetric so mean = 0 is free, but var(t) = t_df / (t_df - 2)
    # for t5 that is 5/3, not 1
    # strategy: 
    # target the covariance at the native scale 
    # (diag = nat_var, off-diag = nat_var * rho,
    # so the implied correlation is still rho), then divide by
    # sqrt(nat_var) at the end. correlation is scale-invariant so cor stays at
    # rho, and the heavy tails are preserved
    nat_var <- t_df / (t_df - 2)
    sigma <- nat_var * matrix(c(1, rho, rho, 1), 2, 2)
    margin <- list(distr = "t", df = t_df)
    vd <- covsim::vita(list(margin, margin), sigma, verbose = FALSE, Nmax = 10^6)
    e <- rvinecopulib::rvine(n, vine = vd)
    return(list(x = e[, 1] / sqrt(nat_var), w = e[, 2] / sqrt(nat_var)))
  }
  if (distr_exo %in% c("chisq_same", "chisq_diff")) {
    # chisq has mean = chi_df and var = 2*chi_df (= nat_var), so unlike unif and t
    # we need both centering and scaling
    # strategy: sample at native scale, then subtract chi_df (mean) and divide 
    # by sqrt(nat_var) (sd) to standardize
    rho_in <- if (distr_exo == "chisq_same") rho else -rho
    nat_var <- 2 * chi_df
    sigma <- nat_var * matrix(c(1, rho_in, rho_in, 1), 2, 2)
    margin <- list(distr = "chisq", df = chi_df)
    vd <- covsim::vita(list(margin, margin), sigma, verbose = FALSE, Nmax = 10^6)
    e <- rvinecopulib::rvine(n, vine = vd)
    x <- (e[, 1] - chi_df) / sqrt(nat_var)
    w <- (e[, 2] - chi_df) / sqrt(nat_var)
    if (distr_exo == "chisq_diff") w <- -w
    return(list(x = x, w = w))
  }
  stop("unknown distr_exo: ", distr_exo)
}

gen_data <- function(n, a3_true, rel_key, distr_exo, misspec_key = "none") {
  # easiest way thus far to get misspecification from the list (to be decided); 
  # ask claude:
  msp <- misspec_specs[[misspec_key]] 
  exo <- gen_exo(n, distr_exo); x <- exo$x; w <- exo$w
  m <- a1*x + a2*w + a3_true*x*w + rnorm(n, sd = sqrt(res_m_var))
  # y eq adds omitted c2*W and c3*X*W direct paths when misspec is non-zero (!)
  y <- b*m + cp*x + msp$c2*w + msp$c3*(x*w) + rnorm(n, sd = sqrt(res_y_var))
  
  target_rel <- rel_levels[[rel_key]]
  # not be confused: both are the same and we use the population var = 1 (!):
  e_sd_xw <- err_sd_for(target_rel, var_eta = 1) 
  e_sd_m <- err_sd_for(target_rel, var_eta = var_m_at(a3_true))
  e_sd_y <- err_sd_for(target_rel, var_eta = var_y_at(a3_true))
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
    Sigma <- matrix(c(sd_m3^2, msp$rcov*sd_m3*sd_y3,
                      msp$rcov*sd_m3*sd_y3, sd_y3^2), 2, 2) 
    # we want correlation, not covariance (!)
    # otherwise the meaning would change with reliability level
    # at low rel the residuals are bigger, so the same covariance number 
    # would imply a different correlation than at high rel.
    pair <- lavaan:::lav_mvrnorm(n, mu = c(0, 0), sigma_1 = Sigma)
    M_ind$m3 <- loadings[3]*m + pair[, 1]
    Y_ind$y3 <- loadings[3]*y + pair[, 2]
  }
  
  cbind(X_ind, W_ind, M_ind, Y_ind)
}

# mc + wald inference from coef + vcov
# note:
# imm: mc quantile CI via lav_mvrnorm draws on vcov(fit)[c("a3","b"), .]
# other parameters just use the: wald CI from model vcov
# this was validated with the lavaan implementation (commit d3f423c)
# i use my own implementation for consistency throughout all methods though
mc_inference <- function(co, V, M = m_mc, alpha = 0.05) {
  na6 <- setNames(rep(NA_real_, length(pars)), pars)
  if (!is.matrix(V) ||
      !all(c("a1","a2","a3","b","cp") %in% rownames(V)) ||
      !all(c("a1","a2","a3","b","cp") %in% names(co)))
    return(list(est = na6, se = na6, lo = na6, hi = na6))
  
  pn <- c("a1","a2","a3","b","cp")
  z <- qnorm(1 - alpha/2)
  se_struct <- sqrt(diag(V)[pn])
  lo_struct <- co[pn] - z * se_struct
  hi_struct <- co[pn] + z * se_struct
  
  V_ab <- V[c("a3","b"), c("a3","b")]
  sims <- lavaan:::lav_mvrnorm(n = M, mu = c(co["a3"], co["b"]),
                                   sigma_1 = V_ab)
  imm_sims <- sims[, 1] * sims[, 2]
  imm_se <- sd(imm_sims)
  imm_lo <- as.numeric(quantile(imm_sims, alpha/2))
  imm_hi <- as.numeric(quantile(imm_sims, 1 - alpha/2))
  
  est <- c(co[pn], imm = unname(co["a3"] * co["b"]))
  se <- c(se_struct, imm = imm_se)
  lo <- c(lo_struct, imm = imm_lo)
  hi <- c(hi_struct, imm = imm_hi)
  
  list(est = est[pars], se = se[pars], lo = lo[pars], hi = hi[pars])
}

# admissibility/implausible helpers (e.g., heywood cases)
# note that lavaan has all of this already, but not modsem (from my knowledge)
# therefore, also for consistency, the same exact check throughout all appraoches
# neg_theta: any residual variance (diag of theta) <= 0
# npd_psi: factor covariance matrix psi has any non-positive eigenvalue
# maybe move to a helpers file. not sure. i generally dislike an extra helpers.R 
# note also this means we include term x:w in the PD check for lsam
neg_theta_diag <- function(v) any(v <= 0, na.rm = TRUE)

npd_psi_mat <- function(psi_mat) {
  if (!length(psi_mat) || !nrow(psi_mat)) return(NA)
  if (any(!is.finite(psi_mat))) return(TRUE)
  min(eigen(psi_mat, symmetric = TRUE, only.values = TRUE)$values) <= 0
}

admis_from_pe <- function(pe, latents) {
  vv <- pe[pe$op == "~~", ]
  th <- vv[vv$lhs == vv$rhs & !(vv$lhs %in% latents), ]
  ps <- vv[vv$lhs %in% latents & vv$rhs %in% latents, ]
  L_names <- intersect(latents, unique(c(ps$lhs, ps$rhs)))
  psi <- matrix(0, length(L_names), length(L_names),
                dimnames = list(L_names, L_names))
  psi[cbind(ps$lhs, ps$rhs)] <- ps$est
  psi[cbind(ps$rhs, ps$lhs)] <- ps$est
  list(neg_theta = neg_theta_diag(th$est),
       npd_psi   = npd_psi_mat(psi))
}

safe_admis <- function(expr) {
  res <- try(expr, silent = TRUE)
  if (inherits(res, "try-error")) return(list(neg_theta = NA, npd_psi = NA))
  res
}

# estimation per method
# any error bubbles up to try_with_warnings in work_for_condition, where the
# message is captured into the `warnings` column and an NA row is produced.
fit_lsam <- function(data) {
  # mm.list: X and W (interacting variables) share one block;
  # M and Y get their own blocks. three measurement blocks total.
  fit <- lavaan::sam(model_lv, data = data,
                     sam.method = "local", se = "local",
                     mm.list = list(c("X", "W"), "M", "Y"))
  V <- vcov(fit)
  co <- coef(fit)
  # remove.step1 as FALSE is important for the implausible checks
  pe <- lavaan::parameterEstimates(fit, remove.step1 = FALSE) 
  adm <- safe_admis(admis_from_pe(pe, latents = c("X","W","M","Y","X:W")))
  list(co = co, V = V, neg_theta = adm$neg_theta, npd_psi = adm$npd_psi)
}

fit_lms <- function(data) {
  fit <- modsem_da(model_lv, data = data, method = "lms",
                   robust.se = TRUE) # robust SEs
  V <- vcov(fit)
  co <- coef(fit)
  pe <- parameter_estimates(fit)
  adm <- safe_admis(admis_from_pe(pe, latents = c("X","W","M","Y")))
  list(co = co, V = V, neg_theta = adm$neg_theta, npd_psi = adm$npd_psi)
}

fit_dblcent <- function(data) {
  fit <- modsem(model_lv, data = data, method = "dblcent",
                estimator = "MLR") # robust SEs but also statistics diff (?) 
  V <- vcov(fit)
  co <- coef(fit)
  # dblcent introduces an extra XW latent from product indicators
  pe <- parameter_estimates(fit)
  adm <- safe_admis(admis_from_pe(pe, latents = c("X","W","XW","M","Y")))
  list(co = co, V = V, neg_theta = adm$neg_theta, npd_psi = adm$npd_psi)
}

methods <- list(lsam = fit_lsam,
                lms = fit_lms,
                dblcent = fit_dblcent)

# wrapper per replication
true_value_for <- function(a3v) setNames(c(a1, a2, a3v, b, cp, b * a3v), pars)

work_for_condition <- function(i, r) {
  condition  <- design[i, ]
  truth <- true_value_for(condition$a3)

  # helper: run expr capturing BOTH errors and warnings
  # returns list(value, error, warnings) — error = NA on success
  # tidyverse has functions for this, but i personally want to avoid extra large packages
  try_with_warnings <- function(expr) {
    warns <- character(0); err <- NA_character_
    val <- tryCatch(
      withCallingHandlers(expr,
        warning = function(w) {
          warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning")
        }),
      error = function(e) { err <<- conditionMessage(e); NULL })
    notes <- c(if (!is.na(err)) paste0("ERROR: ", err), unique(warns))
    list(value = val, error = err,
         warnings = if (length(notes)) paste(notes, collapse = " | ")
                    else NA_character_)
  }

  # gen_data with up to 5 retries. vita can fail on rare rng states
  dat <- NULL
  for (attempt in seq_len(5)) {
    attempt_stream <- paste(get(".Random.seed", envir = .GlobalEnv), collapse = ":")
    dat <- tryCatch(
      gen_data(condition$n, condition$a3, condition$rel,
               condition$distr_exo, condition$misspec),
      error = function(e) NULL)
    if (!is.null(dat)) { rng_stream <- attempt_stream; break }
  }
  if (is.null(dat)) {
    cat(sprintf("c%04d  r%04d  data failure\n", i, r),
        file = file.path(results_dir, "progress.log"), append = TRUE)
    return(NULL)
  }

  # per-method fit + mc_inference, each wrapped independently so one
  # method's error does NOT block the others. Method failures produce
  # NA rows for that method only, with the error captured in 'warnings'
  rows_all <- lapply(names(methods), function(m) {
    mw <- try_with_warnings({
      fitted <- methods[[m]](dat)
      list(fitted = fitted, res = mc_inference(fitted$co, fitted$V))
    })

    if (!is.na(mw$error)) {
      # method level failure — NA row + the error message stored in warnings
      return(data.frame(parameter = pars,
                        est = NA_real_, se = NA_real_,
                        ci_lo = NA_real_, ci_hi = NA_real_,
                        neg_theta = NA, npd_psi = NA,
                        warnings = mw$warnings,
                        method = m, stringsAsFactors = FALSE))
    }

    fitted <- mw$value$fitted
    res <- mw$value$res
    data.frame(parameter = pars,
               est = as.numeric(res$est),
               se = unname(res$se),
               ci_lo = unname(res$lo),
               ci_hi = unname(res$hi),
               neg_theta = fitted$neg_theta,
               npd_psi = fitted$npd_psi,
               warnings = mw$warnings,
               method = m,
               stringsAsFactors = FALSE)
  })
  rows <- do.call(rbind, rows_all)

  rows$condition <- sprintf("c%04d", i)
  rows$n <- condition$n
  rows$a3 <- condition$a3
  rows$rel <- condition$rel
  rows$distr_exo <- condition$distr_exo
  rows$misspec <- condition$misspec
  rows$rep <- r
  rows$true_value <- truth[rows$parameter]
  rows$rng_stream <- rng_stream
  rows[, c("condition","n","a3","rel","distr_exo","misspec","rep","method","parameter",
           "true_value","est","se","ci_lo","ci_hi","neg_theta","npd_psi",
           "warnings","rng_stream")]
}

# parallel run over conditions
run_condition <- function(i) {
  condition_file <- file.path(results_dir, sprintf("condition_mc_%04d.rds", i))
  if (file.exists(condition_file)) return(invisible(NULL))
  
  t0 <- Sys.time()
  rows_list <- lapply(seq_len(r_per_condition),
                      function(r) try(work_for_condition(i, r), silent = TRUE))
  ok <- vapply(rows_list, is.data.frame, logical(1))
  condition_raw  <- do.call(rbind, rows_list[ok])
  # guard!: don't write a NULL/empty file (would lock the condition out of future
  # resume skip retries even though it has no usable data)
  if (!is.null(condition_raw) && nrow(condition_raw) > 0)
    saveRDS(condition_raw, condition_file)
  
  dt <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  cat(sprintf("c%04d  %.1fs  %d/%d\n", i, dt, sum(ok), r_per_condition),
      file = file.path(results_dir, "progress.log"), append = TRUE)
  invisible(NULL)
}

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

# aggregate raw condition files into raw_mc.rds
raw <- do.call(rbind, lapply(
  list.files(results_dir, "^condition_mc_\\d+\\.rds$", full.names = TRUE),
  readRDS))
saveRDS(raw, file.path(results_dir, "raw_mc.rds"))