# post-processing for the moderated-mediation MC sim: builds the result tables.
# figures are made separately by figures.R.
# in : results/estimates_mc.rds
# out: results/{summary,rejection_imm_a3,convergence_report}.rds

# convergence is judged on all estimated parameters: the structural paths of
# interest (poi = a1, a2, a3, b, cp, imm) plus the 12 free loadings. two
# criteria, each reported and each producing its own metrics:
#   soft   : any est/se/ci (any parameter) is NA/NaN/Inf.
#   strict : soft, OR an inadmissible solution flagged at fit time in sim_mc.R
#            (neg_theta = residual variance <= 0; npd_psi = non-PD factor cov).
# warnings are counted (n_warn_method) but do NOT exclude.

# exclusion follows the per-method-then-global pattern of Results(1).R:
# - per-method counts (n_conv_method, n_outlier_method) are kept BEFORE the
# global step, so we can see e.g. lms converging more than lsam in a cell.
# - global step: a rep is dropped from ALL methods unless every method passed.
# so after global convergence + global outlier removal every method in a
# condition runs on the identical rep set (n_final is shared).
# outliers: IQR +- 3 fences on poi est per (condition, method, parameter); a rep
# flagged on any method/parameter is dropped condition-wide.
# metrics carry jackknife MCSEs (*_mcse) from simhelpers, per (condition, method,
# parameter): rel_bias = mean(est)/true - 1, rel_rmse = rmse/|true|,
# coverage = P(lo <= true <= hi). rel_bias_var = E[se^2]/var(est) and
# rel_rmse_var = relative RMSE of se^2 about var(est). relative measures fall
# back to absolute at true == 0.

library(dplyr); library(simhelpers)

results_dir <- "results"

estimates <- readRDS(file.path(results_dir, "estimates_mc.rds"))

cond_keys     <- c("condition", "n", "a3", "rel", "distr_exo", "misspec")
methods_all   <- c("lsam", "lms", "dblcent")
n_methods_all <- length(methods_all)
poi           <- c("a1", "a2", "a3", "b", "cp", "imm")  # structural paths of interest

is_bad <- function(x) is.na(x) | is.nan(x) | is.infinite(x)

# --- per (condition, method, rep) convergence flags ---------------------------
conv <- estimates |>
  group_by(across(all_of(c(cond_keys, "method", "rep")))) |>
  summarise(
    soft_nonconv = any(is_bad(est) | is_bad(se) | is_bad(ci_lo) | is_bad(ci_hi)),
    inadmissible = any(neg_theta %in% TRUE) | any(npd_psi %in% TRUE),
    warned       = any(!is.na(warnings)),
    .groups = "drop") |>
  mutate(strict_nonconv = soft_nonconv | inadmissible)

# per-rep global convergence: a rep survives only if ALL methods converged
rep_conv <- conv |>
  group_by(across(all_of(c(cond_keys, "rep")))) |>
  summarise(
    all_present   = n_distinct(method) == n_methods_all,
    soft_global   = all_present & !any(soft_nonconv),
    strict_global = all_present & !any(strict_nonconv),
    .groups = "drop")

# --- summary helpers (unchanged) ---------------------------------------------
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

# --- per-criterion pipeline ---------------------------------------------------
# builds, for one criterion: the convergence/outlier report (counts), the metric
# summary, and the rejection table. exclusion is global; per-method counts are
# captured before the global step so they stay informative.
process <- function(criterion) {
  nonconv_col <- paste0(criterion, "_nonconv")
  global_col  <- paste0(criterion, "_global")

  # per-method convergence within a condition (pre-global)
  conv_counts <- conv |>
    group_by(across(all_of(c(cond_keys, "method")))) |>
    summarise(n_total       = n(),
              n_conv_method  = sum(!.data[[nonconv_col]]),
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
  clean_conv <- estimates |> semi_join(kept_reps, by = c(cond_keys, "rep"))

  # outliers on the globally-converged set (slim frame: keys + est only)
  rep_out <- clean_conv |>
    select(all_of(c(cond_keys, "method", "parameter", "rep")), est) |>
    filter(parameter %in% poi) |>
    group_by(across(all_of(c(cond_keys, "method", "parameter")))) |>
    mutate(q1 = quantile(est, 0.25, na.rm = TRUE),
           q3 = quantile(est, 0.75, na.rm = TRUE),
           iqr = q3 - q1,
           is_out = est < q1 - 3 * iqr | est > q3 + 3 * iqr) |>
    group_by(across(all_of(c(cond_keys, "method", "rep")))) |>
    summarise(method_out = any(is_out, na.rm = TRUE), .groups = "drop")

  # per-method outlier counts (pre-global)
  out_counts_method <- rep_out |>
    group_by(across(all_of(c(cond_keys, "method")))) |>
    summarise(n_eval_outlier   = n(),
              n_outlier_method = sum(method_out),
              .groups = "drop") |>
    mutate(outlier_rate_method = n_outlier_method / n_eval_outlier)

  # a rep flagged on any method is dropped condition-wide
  rep_out_global <- rep_out |>
    group_by(across(all_of(c(cond_keys, "rep")))) |>
    summarise(rep_out = any(method_out), .groups = "drop")

  out_counts_global <- rep_out_global |>
    group_by(across(all_of(cond_keys))) |>
    summarise(n_eval_global = n(),
              n_final        = sum(!rep_out),
              .groups = "drop") |>
    mutate(global_outlier_rate = (n_eval_global - n_final) / n_eval_global) |>
    select(all_of(cond_keys), n_final, global_outlier_rate)

  final_reps <- rep_out_global |> filter(!rep_out) |>
    select(all_of(c(cond_keys, "rep")))
  clean <- clean_conv |> semi_join(final_reps, by = c(cond_keys, "rep"))

  # report: one row per (condition x method)
  report <- conv_counts |>
    left_join(select(global_counts, all_of(cond_keys), n_reps, n_global_conv, global_conv_rate),
              by = cond_keys) |>
    left_join(out_counts_method, by = c(cond_keys, "method")) |>
    left_join(out_counts_global, by = cond_keys) |>
    mutate(n_outlier_method = coalesce(n_outlier_method, 0L),
           n_final          = coalesce(n_final, 0L),
           criterion        = criterion)

  # metrics on the globally-clean set (all parameters, incl. loadings)
  summary_tbl <- clean |>
    group_by(across(all_of(c(cond_keys, "method", "parameter")))) |>
    group_modify(~ calc_one(.x)) |>
    ungroup() |>
    mutate(criterion = criterion)

  rejection_tbl <- clean |>
    filter(parameter %in% c("imm", "a3"), a3 %in% c(0, 0.2)) |>
    group_by(across(all_of(c(cond_keys, "method", "parameter")))) |>
    group_modify(~ calc_reject(.x)) |>
    ungroup() |>
    mutate(quantity  = if_else(a3 == 0, "Type I error", "Power"),
           criterion = criterion)

  list(summary = summary_tbl, rejection = rejection_tbl, report = report)
}

soft   <- process("soft")
strict <- process("strict")

summary_tbl   <- bind_rows(soft$summary,   strict$summary)
rejection_tbl <- bind_rows(soft$rejection, strict$rejection)
conv_report   <- bind_rows(soft$report,    strict$report)

saveRDS(summary_tbl, file.path(results_dir, "summary.rds"))
saveRDS(rejection_tbl, file.path(results_dir, "rejection_imm_a3.rds"))
saveRDS(conv_report, file.path(results_dir, "convergence_report.rds"))
