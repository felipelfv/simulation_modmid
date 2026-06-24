# paper figures from the analyze_mc.R tables (results/summary.rds,
# results/rejection_imm_a3.rds) to plots/. methods are distinguished by
# point shape + line type, not colour. titles are omitted; captions come from the
# quarto chunk.

library(dplyr); library(ggplot2)

# summary.rds / rejection carry both soft and strict criteria; the paper uses soft.
summary_tbl <- readRDS(file.path("results", "summary.rds")) |> filter(criterion == "soft")
rej_tbl     <- readRDS(file.path("results", "rejection_imm_a3.rds")) |> filter(criterion == "soft")
plot_dir    <- "plots"
dir.create(plot_dir, showWarnings = FALSE)
unlink(list.files(plot_dir, pattern = "\\.png$", full.names = TRUE))

mis_lv  <- c("none","c2_small","c2_medium","c3_small","c3_medium",
             "crossload_small","crossload_medium","rescov_small","rescov_medium")
# plotmath labels: c gets a subscript, cl/rc/none stay upright. rendered on the
# x-axis via x_misspec (scale_x_discrete + parse) below.
mis_lab <- c("n",
             "c[2]*'.s'", "c[2]*'.m'", "c[3]*'.s'", "c[3]*'.m'",
             "'cl.s'", "'cl.m'", "'rc.s'", "'rc.m'")
dis_lv  <- c("normal","uniform","t5","chisq_same","chisq_diff")

# distribution strip labels (plain text so strip.text bold applies; plotmath
# symbols do not honour fontface)
dist_labeller <- as_labeller(
  c(normal = "Normal", uniform = "Uniform", t5 = "t(5)",
    chisq_same = "Chi-square (s)", chisq_diff = "Chi-square (d)"))

n_labeller <- as_labeller(function(x) paste0("N = ", x))

# parameter strip labels as plotmath (subscripts for a_k, prime for c'); a3 column
# (the design interaction value) rendered as a[3] == value. plotmath strips are not
# bold (fontface is not honoured), unlike the plain-text distribution strips.
param_labeller <- as_labeller(
  c(a1 = "a[1]", a2 = "a[2]", a3 = "a[3]", b = "b", cp = "c*\"'\"", imm = "imm"),
  label_parsed)
a3_labeller <- as_labeller(function(x) paste0("a[3] == ", x), label_parsed)

# x-axis misspecification labels: parse the plotmath strings in mis_lab
x_misspec <- scale_x_discrete(labels = function(l) parse(text = l))

# factor coding shared by the metric and rejection plots
prep_factors <- function(d) {
  d |> mutate(
    misspec   = factor(misspec, levels = mis_lv, labels = mis_lab),
    distr_exo = factor(distr_exo, levels = dis_lv),
    rel       = factor(rel, levels = c("low", "high"),
                       labels = c("Reliability: low", "Reliability: high")),
    method    = factor(method, levels = c("lsam", "lms", "dblcent"), labels = c("LSAM", "LMS", "UPI")))
}

method_shapes <- c(LSAM = 16, LMS = 17, UPI = 15)
method_lines  <- c(LSAM = "solid", LMS = "dashed", UPI = "dotdash")
method_fills  <- c(LSAM = "grey25", LMS = "grey55", UPI = "grey85")

# shared facet + scale + theme layers; caller may override the facet
apa_grid <- function(facets = facet_grid(distr_exo ~ rel + n,
                       labeller = labeller(distr_exo = dist_labeller, n = n_labeller))) list(
  facets,
  scale_shape_manual(values = method_shapes),
  scale_linetype_manual(values = method_lines),
  theme_bw(base_size = 10),
  theme(legend.position = "bottom", panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold")))

# per metric: axis label, MCSE column (NULL = none), acceptable band, reference line
metric_spec <- list(
  rel_bias = list(lab = "Relative bias", mcse = "rel_bias_mcse", band = c(-0.10, 0.10),   ref = 0),
  rel_rmse = list(lab = "Relative RMSE", mcse = "rel_rmse_mcse", band = NULL,             ref = 0),
  coverage = list(lab = "Coverage",      mcse = "coverage_mcse", band = c(0.925, 0.975),  ref = 0.95),
  se_sd    = list(lab = "SE/SD",         mcse = NULL,            band = c(0.90, 1.10),    ref = 1),
  rel_bias_var = list(lab = "Rel. bias of variance", mcse = "rel_bias_var_mcse", band = c(0.90, 1.10), ref = 1),
  typeI    = list(lab = "Type I error",  mcse = "reject_mcse",   band = c(0.025, 0.075),  ref = 0.05),
  power    = list(lab = "Power",          mcse = "reject_mcse",   band = NULL,             ref = 0.80)
)

# band + ref line + +/-1 MCSE bars + line + points for a metric (caller adds the
# facet via apa_grid); shared by the per-condition slice figures.
slice_geoms <- function(d, metric) {
  m <- metric_spec[[metric]]
  p <- ggplot(d, aes(misspec, .data[[metric]], shape = method, linetype = method, group = method))
  if (!is.null(m$band))
    p <- p + annotate("rect", xmin = -Inf, xmax = Inf, ymin = m$band[1], ymax = m$band[2],
                      fill = "grey60", alpha = 0.25)
  p <- p + geom_hline(yintercept = m$ref, linetype = "dotted")
  if (!is.null(m$mcse))
    p <- p + geom_errorbar(aes(ymin = .data[[metric]] - .data[[m$mcse]],
                               ymax = .data[[metric]] + .data[[m$mcse]]),
                           width = .25, position = position_dodge(.5), linewidth = .3,
                           linetype = "solid", color = "grey30")
  p +
    geom_line(position = position_dodge(.5), linewidth = .5) +
    geom_point(position = position_dodge(.5), size = 1.6) +
    x_misspec +
    labs(x = "Misspecification", y = m$lab, shape = NULL, linetype = NULL)
}

# one metric for one parameter at fixed (a3, n): distribution x reliability grid,
# x = misspecification, line + point per method, +/- 1 MCSE bars, acceptable band.
fig_slice <- function(param = "imm", a3v = 0.2, metric = "rel_bias", nvals = c(200, 500)) {
  d <- summary_tbl |>
    filter(parameter == param, a3 == a3v, n %in% nvals) |>
    prep_factors()
  # drop rel. variance bias > 10: a few cells blow up and swamp the y-axis (flagged for review)
  if (metric == "rel_bias_var") d <- filter(d, rel_bias_var <= 10)
  slice_geoms(d, metric) + apa_grid()
}

# baseline relative bias for the imm and a3 under the normal, correctly specified
# condition, across ALL sample sizes, both reliabilities, and both nonzero
# interaction values: parameter (rows) x reliability + a3 (cols), x = sample size,
# line + point per method, +/- 1 MCSE bars.
fig_baseline_bias <- function(a3v = c(0.2, 0.4)) {
  m <- metric_spec[["rel_bias"]]
  d <- summary_tbl |>
    filter(parameter %in% c("imm", "a3"), distr_exo == "normal",
           misspec == "none", a3 %in% a3v) |>
    prep_factors() |>
    mutate(parameter = factor(parameter, levels = c("imm", "a3")),
           n = factor(n, levels = c(100, 200, 300, 500)))
  ggplot(d, aes(n, rel_bias, shape = method, linetype = method, group = method)) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = m$band[1], ymax = m$band[2],
             fill = "grey60", alpha = 0.25) +
    geom_hline(yintercept = m$ref, linetype = "dotted") +
    geom_errorbar(aes(ymin = rel_bias - rel_bias_mcse, ymax = rel_bias + rel_bias_mcse),
                  width = .25, position = position_dodge(.5), linewidth = .3,
                  linetype = "solid", color = "grey30") +
    geom_line(position = position_dodge(.5), linewidth = .5) +
    geom_point(position = position_dodge(.5), size = 1.8) +
    facet_grid(parameter ~ rel + a3,
               labeller = labeller(parameter = param_labeller, a3 = a3_labeller)) +
    scale_shape_manual(values = method_shapes) +
    scale_linetype_manual(values = method_lines) +
    theme_bw(base_size = 10) +
    theme(legend.position = "bottom", panel.grid.minor = element_blank(),
          strip.text = element_text(face = "bold")) +
    labs(x = "Sample size (N)", y = m$lab, shape = NULL, linetype = NULL)
}

# rel. bias of variance, imm and a3 together at one reliability level: distribution
# (rows) x parameter x N (cols). high reliability goes in the manuscript, low in the
# appendix (low-reliability cells blow up under heavy tails and swamp the y-axis).
fig_relvar_both <- function(a3v = 0.2, rel_lvl = "high") {
  d <- summary_tbl |>
    filter(parameter %in% c("imm", "a3"), a3 == a3v, n %in% c(200, 500),
           rel == rel_lvl, rel_bias_var <= 10) |>
    prep_factors() |>
    mutate(parameter = factor(parameter, levels = c("imm", "a3")))
  slice_geoms(d, "rel_bias_var") +
    apa_grid(facet_grid(distr_exo ~ parameter + n,
                        labeller = labeller(distr_exo = dist_labeller,
                                            parameter = param_labeller, n = n_labeller)))
}

# imm and a3 together as greyscale grouped bars: parameter x distribution (rows),
# reliability x N (cols), one bar per method with +/- 1 MCSE. d2 supplies an
# alternate table (e.g. the rejection table) with its own value/MCSE columns.
fig_bar_both <- function(metric = "rel_rmse", a3v = 0.2, nvals = c(200, 500),
                         d2 = NULL, value = metric, mcse = NULL, ylab = NULL,
                         band = NULL, ref = NULL) {
  m <- metric_spec[[metric]]
  if (is.null(band)) band <- m$band
  if (is.null(ref))  ref  <- m$ref
  if (is.null(mcse)) mcse <- m$mcse
  if (is.null(ylab)) ylab <- m$lab
  src <- if (is.null(d2)) summary_tbl else d2
  d <- src |>
    filter(parameter %in% c("imm", "a3"), a3 == a3v, n %in% nvals) |>
    prep_factors() |>
    mutate(parameter = factor(parameter, levels = c("imm", "a3")))
  if (metric == "rel_bias_var") d <- filter(d, rel_bias_var <= 10)
  p <- ggplot(d, aes(misspec, .data[[value]], fill = method, group = method))
  if (!is.null(band))
    p <- p + annotate("rect", xmin = -Inf, xmax = Inf, ymin = band[1], ymax = band[2],
                      fill = "grey60", alpha = 0.25)
  p <- p + geom_hline(yintercept = ref, linetype = "dotted") +
    geom_col(position = position_dodge(.8), width = .7, colour = "grey30", linewidth = .2)
  if (!is.null(mcse))
    p <- p + geom_errorbar(aes(ymin = .data[[value]] - .data[[mcse]],
                               ymax = .data[[value]] + .data[[mcse]]),
                           position = position_dodge(.8), width = .3, linewidth = .25)
  p +
    facet_grid(parameter + distr_exo ~ rel + n,
               labeller = labeller(distr_exo = dist_labeller,
                                   parameter = param_labeller, n = n_labeller)) +
    scale_fill_manual(values = method_fills) +
    x_misspec +
    theme_bw(base_size = 9) +
    theme(legend.position = "bottom", panel.grid.minor = element_blank(),
          strip.text = element_text(face = "bold"),
          axis.text.x = element_text(size = 6)) +
    labs(x = "Misspecification", y = ylab, fill = NULL)
}

save_bar_both <- function(metric = "rel_rmse", a3v = 0.2, nvals = c(200, 500),
                          w = 10, h = 12, dpi = 150, ...) {
  f <- file.path(plot_dir, sprintf("both_%s_a3-%s.png", metric, a3v))
  ggsave(f, fig_bar_both(metric, a3v, nvals, ...), width = w, height = h, dpi = dpi)
}

# imm and a3 together in one panel: imm as greyscale bars (method by shade), a3 as
# a slate line + points (method by shape/line type). suits bias, where imm is the
# larger, more structured effect and a3 the small one. distr (rows) x rel x N (cols).
fig_combo_both <- function(metric = "rel_bias", a3v = 0.2, nvals = c(200, 500)) {
  m <- metric_spec[[metric]]
  base   <- summary_tbl |> filter(a3 == a3v, n %in% nvals)
  d_bar  <- base |> filter(parameter == "imm") |> prep_factors()
  d_line <- base |> filter(parameter == "a3")  |> prep_factors()
  p <- ggplot(mapping = aes(misspec))
  if (!is.null(m$band))
    p <- p + annotate("rect", xmin = -Inf, xmax = Inf, ymin = m$band[1], ymax = m$band[2],
                      fill = "grey60", alpha = 0.22)
  p <- p + geom_hline(yintercept = m$ref, linetype = "dotted") +
    geom_col(data = d_bar, aes(y = .data[[metric]], fill = method, group = method),
             position = position_dodge(.8), width = .7, colour = "grey35", linewidth = .15) +
    geom_line(data = d_line, aes(y = .data[[metric]], linetype = method, group = method),
              position = position_dodge(.4), colour = "#4C72B0", linewidth = .5) +
    geom_point(data = d_line, aes(y = .data[[metric]], shape = method, group = method),
               position = position_dodge(.4), colour = "#4C72B0", size = 1.3)
  p +
    facet_grid(distr_exo ~ rel + n,
               labeller = labeller(distr_exo = dist_labeller, n = n_labeller)) +
    scale_fill_manual(values = method_fills, name = "imm") +
    scale_shape_manual(values = method_shapes, name = expression(a[3])) +
    scale_linetype_manual(values = method_lines, name = expression(a[3])) +
    x_misspec +
    theme_bw(base_size = 10) +
    theme(legend.position = "bottom", panel.grid.minor = element_blank(),
          strip.text = element_text(face = "bold"), axis.text.x = element_text(size = 6)) +
    labs(x = "Misspecification", y = m$lab)
}

# structural-coefficient overview for one metric at fixed reliability and n:
# parameter (rows) x distribution (cols), x = misspecification, per method.
fig_param <- function(metric = "coverage", params = c("a1", "a2", "b", "cp"),
                      a3v = 0.2, nval = 500, rel_lvl = "high", xangle = 0) {
  m <- metric_spec[[metric]]
  d <- summary_tbl |>
    filter(parameter %in% params, a3 == a3v, n == nval, rel == rel_lvl) |>
    mutate(misspec   = factor(misspec, levels = mis_lv, labels = mis_lab),
           distr_exo = factor(distr_exo, levels = dis_lv),
           parameter = factor(parameter, levels = params),
           method    = factor(method, levels = c("lsam", "lms", "dblcent"), labels = c("LSAM", "LMS", "UPI")))
  if (metric == "rel_bias_var") d <- filter(d, rel_bias_var <= 10)
  p <- ggplot(d, aes(misspec, .data[[metric]], shape = method, linetype = method, group = method))
  if (!is.null(m$band))
    p <- p + annotate("rect", xmin = -Inf, xmax = Inf, ymin = m$band[1], ymax = m$band[2],
                      fill = "grey60", alpha = 0.25)
  p <- p + geom_hline(yintercept = m$ref, linetype = "dotted")
  if (!is.null(m$mcse))
    p <- p + geom_errorbar(aes(ymin = .data[[metric]] - .data[[m$mcse]],
                               ymax = .data[[metric]] + .data[[m$mcse]]),
                           width = .25, position = position_dodge(.5), linewidth = .3,
                           linetype = "solid", color = "grey30")
  p +
    geom_line(position = position_dodge(.5), linewidth = .5) +
    geom_point(position = position_dodge(.5), size = 1.4) +
    facet_grid(parameter ~ distr_exo,
               labeller = labeller(distr_exo = dist_labeller, parameter = param_labeller)) +
    scale_shape_manual(values = method_shapes) +
    scale_linetype_manual(values = method_lines) +
    x_misspec +
    labs(x = "Misspecification", y = m$lab, shape = NULL, linetype = NULL) +
    theme_bw(base_size = 10) +
    theme(legend.position = "bottom", panel.grid.minor = element_blank(),
          strip.text = element_text(face = "bold"),
          axis.text.x = element_text(angle = xangle,
                                     hjust = if (xangle != 0) 1 else 0.5,
                                     vjust = if (xangle == 90) 0.5 else 1))
}

# all four performance measures for a SINGLE loading (default lm3, the only one
# the cross-loading moves) in one figure: measure (rows, free y) x distribution
# (cols), x = misspecification, line + point per method, per-measure band + ref.
fig_loading_measures <- function(loading = "lm3", a3v = 0.2, nval = 500, rel_lvl = "high") {
  meas_lv  <- c("rel_bias", "coverage", "se_sd", "rel_bias_var")
  meas_lab <- c("Relative bias", "Coverage", "SE/SD", "Rel. bias of variance")
  d <- summary_tbl |>
    filter(parameter == loading, a3 == a3v, n == nval, rel == rel_lvl) |>
    mutate(misspec   = factor(misspec, levels = mis_lv, labels = mis_lab),
           distr_exo = factor(distr_exo, levels = dis_lv),
           method    = factor(method, levels = c("lsam", "lms", "dblcent"), labels = c("LSAM", "LMS", "UPI"))) |>
    tidyr::pivot_longer(all_of(meas_lv), names_to = "measure", values_to = "value") |>
    filter(!(measure == "rel_bias_var" & value > 10)) |>
    mutate(measure = factor(measure, levels = meas_lv, labels = meas_lab))
  refdf <- data.frame(
    measure = factor(meas_lab, levels = meas_lab),
    ref = c(0, 0.95, 1, 1),
    blo = c(-0.10, 0.925, 0.90, 0.90),
    bhi = c( 0.10, 0.975, 1.10, 1.10))
  ggplot(d, aes(misspec, value, shape = method, linetype = method, group = method)) +
    geom_rect(data = refdf, inherit.aes = FALSE, xmin = -Inf, xmax = Inf,
              aes(ymin = blo, ymax = bhi), fill = "grey60", alpha = 0.25) +
    geom_hline(data = refdf, inherit.aes = FALSE, aes(yintercept = ref), linetype = "dotted") +
    geom_line(position = position_dodge(.5), linewidth = .5) +
    geom_point(position = position_dodge(.5), size = 1.4) +
    facet_grid(measure ~ distr_exo, scales = "free_y",
               labeller = labeller(distr_exo = dist_labeller)) +
    scale_shape_manual(values = method_shapes) +
    scale_linetype_manual(values = method_lines) +
    x_misspec +
    theme_bw(base_size = 10) +
    theme(legend.position = "bottom", panel.grid.minor = element_blank(),
          strip.text = element_text(face = "bold"),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6)) +
    labs(x = "Misspecification", y = NULL, shape = NULL, linetype = NULL)
}

save_loading_measures <- function(loading = "lm3", a3v = 0.2, nval = 500, rel_lvl = "high",
                                  w = 9, h = 8, dpi = 150) {
  f <- file.path(plot_dir, sprintf("load_%s_allmeasures_a3-%s_n%d_%s.png",
                                   loading, a3v, nval, rel_lvl))
  ggsave(f, fig_loading_measures(loading, a3v, nval, rel_lvl), width = w, height = h, dpi = dpi)
}

# draw a (landscape) plot rotated 90 degrees onto a portrait canvas, so the saved
# PNG is sideways and sits large on a portrait end-page. w, h are the plot's
# native landscape dimensions; the canvas is h (wide) x w (tall).
save_rotated <- function(plot, file, w, h, dpi = 150) {
  grDevices::png(file, width = h, height = w, units = "in", res = dpi)
  on.exit(grDevices::dev.off())
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(width = grid::unit(w, "in"),
                                    height = grid::unit(h, "in"), angle = 90))
  grid::grid.draw(ggplot2::ggplotGrob(plot))
}

# --- create all figures ---------------------------------------------------------
# wide multi-panel figures are saved rotated 90 degrees so they sit large on the
# portrait end-pages; the tall bar figures and lm3 stay upright.

# relative bias: imm bars + a3 slate line, one panel (wide -> rotated)
save_rotated(fig_combo_both("rel_bias", 0.2),
             file.path(plot_dir, "combo_rel_bias_a3-0.2.png"), w = 9.5, h = 8)

# baseline relative bias under the normal, correctly specified condition across all
# sample sizes, both reliabilities, and both nonzero interaction values
ggsave(file.path(plot_dir, "baseline_rel_bias.png"),
       fig_baseline_bias(c(0.2, 0.4)), width = 10, height = 5.5, dpi = 150)

# relative RMSE, Type I error (a3 = 0) and power (a3 = 0.2): tall bar panels, upright
save_bar_both("rel_rmse", 0.2)
save_bar_both("typeI", 0, d2 = rej_tbl, value = "reject")
save_bar_both("power", 0.2, d2 = rej_tbl, value = "reject")

# coverage and SE/SD as line slices (wide -> rotated); bars compress these
# (coverage near .95, SE/SD near 1)
for (p in c("imm", "a3"))
  for (m in c("coverage", "se_sd"))
    save_rotated(fig_slice(p, 0.2, m),
                 file.path(plot_dir, sprintf("%s_%s_a3-0.2_byrel.png", p, m)),
                 w = 9.5, h = 8)

# rel. bias of variance: imm + a3 in one panel (wide -> rotated); high reliability
# in the manuscript, low reliability in the appendix (the low-reliability cells blow
# up under heavy-tailed predictors and would swamp the shared y-axis)
for (rl in c("high", "low"))
  save_rotated(fig_relvar_both(0.2, rl),
               file.path(plot_dir, sprintf("rel_bias_var_a3-0.2_%s.png", rl)),
               w = 9.5, h = 8)

# structural-coefficient overviews (wide -> rotated)
for (m in c("rel_bias", "coverage", "se_sd"))
  save_rotated(fig_param(m, a3v = 0.2, nval = 500, rel_lvl = "high"),
               file.path(plot_dir, sprintf("param_%s_a3-0.2_n500_high.png", m)),
               w = 12, h = 10)

# loadings: all four measures for the cross-loaded indicator lm3, the only loading
# the misspecifications move (the other eleven are flat across every measure)
save_loading_measures("lm3", 0.2, 500, "high")
