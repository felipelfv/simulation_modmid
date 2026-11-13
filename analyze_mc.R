# post-processing for the moderated-mediation MC sim: builds the result tables.
# figures are made separately by figures.R.
# in : results/estimates_mc.rds
# out: results/{summary,rejection_imm_a3,convergence_report,se_extremes}.rds

# convergence is judged on all estimated parameters: the structural paths of
# interest (poi = a1, a2, a3, b, cp, imm) plus the 12 free loadings. three
# passes, each reported and each producing its own tables:
# soft   : the optimizer reported nonconvergence, OR any est/se/ci (any
# parameter) is NA/NaN/Inf.
# strict : soft, OR an inadmissible solution flagged at fit time in sim_mc.R
# (neg_theta = residual variance <= 0, npd_psi = non-PD factor cov).
# unscreened : soft convergence with NO outlier removal. a sensitivity check for
# the supplement, which reports the largest change in each measure when the
# screened replications are retained.
# warnings are counted (n_warn_method) but do NOT exclude. so are converged fits
# whose hessian at the solution is not positive definite (n_npd_info).

# exclusion is per-method, then global:
# - per-method counts (n_conv_method, n_outlier_method) are kept BEFORE the
# global step, so we can see e.g. lms converging more than lsam in a cell.
# - global step: a rep is dropped from ALL methods unless every method passed.
# so after global convergence + global outlier removal every method in a
# condition runs on the identical rep set (n_final is shared).
# outliers: IQR +- 3 fences on poi est per (condition, method, parameter), plus
# a degenerate SE (> se_max) on any poi parameter. a rep flagged on any
# method/parameter is dropped condition-wide.
# metrics carry jackknife MCSEs (*_mcse) from simhelpers, per (condition, method,
# parameter): rel_bias = mean(est)/true - 1, rel_rmse = rmse/|true|,
# coverage = P(lo <= true <= hi). rel_bias_var = E[se^2]/var(est) and
# rel_rmse_var = relative RMSE of se^2 about var(est). relative measures fall
# back to absolute at true == 0.
# se_extremes: sensitivity of the SE/SD ratio of a3 and the imm to extreme
# standard errors, those more than 3 IQR above the upper quartile of their cell.

library(dplyr); library(simhelpers)

# TRUE summarises the gaussian-copula robustness run instead of the main run.
copula_robustness <- FALSE
results_dir <- if (copula_robustness) "results_robustness" else "results_main"

estimates <- readRDS(file.path(results_dir, "estimates_mc.rds"))

cond_keys   <- c("condition", "n", "a3", "rel", "distr_exo", "misspec")
methods_all <- c("lsam", "lms", "dblcent", "fsr")
poi         <- c("a1", "a2", "a3", "b", "cp", "imm")  # structural paths of interest
se_max      <- 10  # a larger SE on any poi parameter is a degenerate fit

is_bad <- function(x) is.na(x) | is.infinite(x)

# per (condition, method, rep) convergence flags
conv <- estimates |>
  group_by(across(all_of(c(cond_keys, "method", "rep")))) |>
  summarise(
    soft_nonconv = any(converged %in% FALSE | is_bad(est) | is_bad(se) | is_bad(ci_lo) | is_bad(ci_hi)),
    npd_info     = !soft_nonconv & any(min_eig_info <= 0, na.rm = TRUE),   # reported, not excluded
    inadmissible = any(neg_theta %in% TRUE) | any(npd_psi %in% TRUE),
    warned       = any(!is.na(warnings)),
    .groups = "drop") |>
  mutate(strict_nonconv = soft_nonconv | inadmissible)

# conv is the only consumer of warnings/neg_theta/npd_psi; everything downstream
# needs just the metric columns. project to those and free the full frame before
# processing (warnings holds unique error strings -- the biggest real-heap cost).
est_slim <- estimates |>
  select(all_of(c(cond_keys, "method", "parameter", "rep")),
         true_value, est, se, ci_lo, ci_hi)
rm(estimates); gc()

# per-rep global convergence: a rep survives only if ALL methods converged
rep_conv <- conv |>
  group_by(across(all_of(c(cond_keys, "rep")))) |>
  summarise(
    all_present   = n_distinct(method) == length(methods_all),
    soft_global   = all_present & !any(soft_nonconv),
    strict_global = all_present & !any(strict_nonconv),
    .groups = "drop")

# summary helpers
# one row per (condition x method x parameter). rel_bias subtracts 1 from
# simhelpers' ratio (0-centered, MCSE unchanged); relative measures fall back to
# absolute at true == 0.
calc_one <- function(df) {
  true   <- df$true_value[[1]]
  df$se2 <- df$se^2
  ab  <- calc_absolute(df, estimates = est, true_param = true_value)
  cov <- calc_coverage(df, lower_bound = ci_lo, upper_bound = ci_hi, true_param = true_value)
  rv  <- calc_relative_var(df, estimates = est, var_estimates = se2)
  if (true != 0) {
    rel <- calc_relative(df, estimates = est, true_param = true_value)
    rb <- rel$rel_bias - 1; rb_se <- rel$rel_bias_mcse
    rr <- rel$rel_rmse;     rr_se <- rel$rel_rmse_mcse
  } else {
    rb <- ab$bias; rb_se <- ab$bias_mcse
    rr <- ab$rmse; rr_se <- ab$rmse_mcse
  }
  tibble(
    n_used   = ab$K_absolute,
    true     = true,
    mean_est = mean(df$est),
    sd_est   = sd(df$est),
    mean_se  = mean(df$se),
    bias = ab$bias, bias_mcse = ab$bias_mcse,
    rmse = ab$rmse, rmse_mcse = ab$rmse_mcse,
    rel_bias = rb, rel_bias_mcse = rb_se,
    rel_rmse = rr, rel_rmse_mcse = rr_se,
    coverage = cov$coverage, coverage_mcse = cov$coverage_mcse,
    width = cov$width, width_mcse = cov$width_mcse,
    rel_bias_var = rv$rel_bias_var, rel_bias_var_mcse = rv$rel_bias_var_mcse,
    rel_rmse_var = rv$rel_rmse_var, rel_rmse_var_mcse = rv$rel_rmse_var_mcse,
    se_sd = mean(df$se) / sd(df$est)
  )
}

# reject = 1 - coverage-of-zero (CI excludes 0). Type I error at a3 = 0, power at
# a3 = 0.2, for imm and a3.
calc_reject <- function(df) {
  df$zero <- 0
  cov <- calc_coverage(df, lower_bound = ci_lo, upper_bound = ci_hi, true_param = zero)
  tibble(reject = 1 - cov$coverage, reject_mcse = cov$coverage_mcse, n_used = nrow(df))
}

# per-criterion pipeline
# builds, for one criterion: the convergence/outlier report (counts), the metric
# summary, and the rejection table. exclusion is global; per-method counts are
# captured before the global step so they stay informative.
process <- function(criterion, screen_outliers = TRUE, label = criterion) {
  nonconv_col <- paste0(criterion, "_nonconv")
  global_col  <- paste0(criterion, "_global")

  # per-method convergence within a condition (pre-global)
  conv_counts <- conv |>
    group_by(across(all_of(c(cond_keys, "method")))) |>
    summarise(n_total       = n(),
              n_conv_method  = sum(!.data[[nonconv_col]]),
              n_npd_info     = sum(npd_info),
              n_warn_method  = sum(warned),
              .groups = "drop") |>
    mutate(conv_rate_method = n_conv_method / n_total)

  # global convergence (shared across methods in a condition)
  global_counts <- rep_conv |>
    group_by(across(all_of(cond_keys))) |>
    summarise(n_reps        = n(),
              n_global_conv  = sum(.data[[global_col]]),
              .groups = "drop") |>
    mutate(global_conv_rate = n_global_conv / n_reps)

  # rows from globally-converged reps (every method passed, so no per-method
  # filtering needed)
  kept_reps  <- rep_conv |> filter(.data[[global_col]]) |>
    select(all_of(c(cond_keys, "rep")))

  # outliers on the globally-converged set (slim frame: keys + est + se).
  # computed straight off est_slim so we never materialise a full-width copy.
  rep_out <- est_slim |>
    semi_join(kept_reps, by = c(cond_keys, "rep")) |>
    select(all_of(c(cond_keys, "method", "parameter", "rep")), est, se) |>
    filter(parameter %in% poi) |>
    group_by(across(all_of(c(cond_keys, "method", "parameter")))) |>
    mutate(q1 = quantile(est, 0.25), q3 = quantile(est, 0.75), iqr = q3 - q1,
           is_out = screen_outliers &
                    (est < q1 - 3 * iqr | est > q3 + 3 * iqr | se > se_max)) |>
    group_by(across(all_of(c(cond_keys, "method", "rep")))) |>
    summarise(method_out = any(is_out), .groups = "drop")

  # per-method outlier counts (pre-global)
  out_counts_method <- rep_out |>
    group_by(across(all_of(c(cond_keys, "method")))) |>
    summarise(n_outlier_method    = sum(method_out),
              outlier_rate_method = mean(method_out),
              .groups = "drop")

  # a rep flagged on any method is dropped condition-wide
  rep_out_global <- rep_out |>
    group_by(across(all_of(c(cond_keys, "rep")))) |>
    summarise(rep_out = any(method_out), .groups = "drop")

  out_counts_global <- rep_out_global |>
    group_by(across(all_of(cond_keys))) |>
    summarise(n_final             = sum(!rep_out),
              global_outlier_rate = mean(rep_out),
              .groups = "drop")

  final_reps <- rep_out_global |> filter(!rep_out) |>
    select(all_of(c(cond_keys, "rep")))
  clean <- est_slim |> semi_join(final_reps, by = c(cond_keys, "rep"))

  # extreme standard errors and the SE/SD ratio, for a3 and the imm. pct_ext is
  # the share of reps whose SE lies above the fence, se_sd the ratio as in the
  # summary, and se_sd_trim the ratio with those reps dropped from the mean SE
  # only. the SD of the estimates keeps every rep, since it is the quantity the
  # SEs are meant to track.
  se_ext_tbl <- clean |>
    filter(parameter %in% c("a3", "imm")) |>
    group_by(across(all_of(c(cond_keys, "method", "parameter")))) |>
    mutate(ext = se > quantile(se, 0.75) + 3 * IQR(se)) |>
    summarise(pct_ext    = 100 * mean(ext),
              se_sd      = mean(se) / sd(est),
              se_sd_trim = mean(se[!ext]) / sd(est),
              .groups = "drop") |>
    mutate(criterion = label)

  # report: one row per (condition x method)
  report <- conv_counts |>
    left_join(select(global_counts, all_of(cond_keys), n_reps, n_global_conv, global_conv_rate),
              by = cond_keys) |>
    left_join(out_counts_method, by = c(cond_keys, "method")) |>
    left_join(out_counts_global, by = cond_keys) |>
    mutate(n_outlier_method = coalesce(n_outlier_method, 0L),
           n_final          = coalesce(n_final, 0L),
           criterion        = label)

  # metrics on the globally-clean set (all parameters, incl. loadings)
  summary_tbl <- clean |>
    group_by(across(all_of(c(cond_keys, "method", "parameter")))) |>
    group_modify(~ calc_one(.x)) |>
    ungroup() |>
    mutate(criterion = label)

  rejection_tbl <- clean |>
    filter(parameter %in% c("imm", "a3"), a3 %in% c(0, 0.2)) |>
    group_by(across(all_of(c(cond_keys, "method", "parameter")))) |>
    group_modify(~ calc_reject(.x)) |>
    ungroup() |>
    mutate(quantity  = if_else(a3 == 0, "Type I error", "Power"),
           criterion = label)

  list(summary = summary_tbl, rejection = rejection_tbl, report = report,
       se_ext = se_ext_tbl)
}

soft   <- process("soft")
gc()
strict <- process("strict")
gc()
unscreened <- process("soft", screen_outliers = FALSE, label = "unscreened")

summary_tbl   <- bind_rows(soft$summary,   strict$summary,   unscreened$summary)
rejection_tbl <- bind_rows(soft$rejection, strict$rejection, unscreened$rejection)
conv_report   <- bind_rows(soft$report,    strict$report,    unscreened$report)
se_extremes   <- bind_rows(soft$se_ext,    strict$se_ext,    unscreened$se_ext)

saveRDS(summary_tbl, file.path(results_dir, "summary.rds"))
saveRDS(rejection_tbl, file.path(results_dir, "rejection_imm_a3.rds"))
saveRDS(conv_report, file.path(results_dir, "convergence_report.rds"))
saveRDS(se_extremes, file.path(results_dir, "se_extremes.rds"))
