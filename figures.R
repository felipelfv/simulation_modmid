# figures from the analyze_mc.R tables (results_main/summary.rds and
# results_main/rejection_imm_a3.rds), written to plots/. methods are
# distinguished by point shape and line type or by grey shade, not by colour.
# the figures carry no titles because the captions come from the quarto chunks.

library(dplyr); library(ggplot2)

# the tables carry the soft, strict and unscreened criteria. the paper uses soft.
# the rejection rate is the type I error at a3 = 0 and the power otherwise, so it
# is stored under both names, which lets metric_spec address it like any measure.
summary_tbl <- readRDS(file.path("results_main", "summary.rds")) |> filter(criterion == "soft")
rej_tbl     <- readRDS(file.path("results_main", "rejection_imm_a3.rds")) |>
  filter(criterion == "soft") |> mutate(typeI = reject, power = reject)
plot_dir    <- "plots"
dir.create(plot_dir, showWarnings = FALSE)
unlink(list.files(plot_dir, pattern = "\\.png$", full.names = TRUE))

# plotmath labels of the misspecifications, parsed on the x-axis by x_misspec
mis_lab <- c(none = "n", c2_small = "c[2]*'.s'", c2_medium = "c[2]*'.m'",
             c3_small = "c[3]*'.s'", c3_medium = "c[3]*'.m'",
             crossload_small = "'cl.s'", crossload_medium = "'cl.m'",
             rescov_small = "'rc.s'", rescov_medium = "'rc.m'")
dis_lab <- c(normal = "Normal", uniform = "Uniform", t5 = "t(5)",
             chisq_same = "Chi-square (s)", chisq_diff = "Chi-square (d)")
x_misspec <- scale_x_discrete(labels = function(l) parse(text = l))

# factor coding shared by all figures
prep_factors <- function(d) d |> mutate(
  misspec   = factor(misspec, levels = names(mis_lab), labels = mis_lab),
  distr_exo = factor(distr_exo, levels = names(dis_lab), labels = dis_lab),
  rel       = factor(rel, levels = c("low", "high"),
                     labels = c("Reliability: low", "Reliability: high")),
  method    = factor(method, levels = c("lsam", "lms", "dblcent", "fsr"),
                     labels = c("LSAM", "LMS", "UPI", "FSR")))

# rows of one figure: the given parameters at one interaction value and one
# reliability level
cell <- function(d, params, a3v, rel_lvl) d |>
  filter(parameter %in% params, a3 == a3v, rel == rel_lvl) |>
  prep_factors() |>
  mutate(parameter = factor(parameter, levels = params))

# strip labels. sample size reads N = value, the interaction value and the
# parameters are plotmath. plotmath strips are not bold, plain text ones are.
strip_labels <- labeller(
  n         = function(x) paste0("N = ", x),
  a3        = as_labeller(function(x) paste0("a[3] == ", x), label_parsed),
  parameter = as_labeller(c(a1 = "a[1]", a2 = "a[2]", a3 = "a[3]", b = "b",
                            cp = "c*\"'\"", imm = "imm"), label_parsed))

# methods as point shapes and line types, and as grey shades for bars
method_marks <- list(
  scale_shape_manual(values = c(LSAM = 16, LMS = 17, UPI = 15, FSR = 18)),
  scale_linetype_manual(values = c(LSAM = "solid", LMS = "dashed", UPI = "dotdash", FSR = "dotted")))
method_fills <- scale_fill_manual(values = c(LSAM = "grey25", LMS = "grey55", UPI = "grey85", FSR = "white"))

# theme shared by the line and bar figures
apa <- list(theme_bw(base_size = 10),
            theme(legend.position = "bottom", panel.grid.minor = element_blank(),
                  strip.text = element_text(face = "bold")))
# smaller axis text for figures with nine misspecification labels per panel
small_x <- theme(axis.text.x = element_text(size = 6.5))

# per measure: axis label, MCSE column (NULL if none), acceptable band (NULL if
# none) and reference line
metric_spec <- list(
  rel_bias     = list(lab = "Relative bias", mcse = "rel_bias_mcse",     band = c(-0.10, 0.10),  ref = 0),
  rel_rmse     = list(lab = "Relative RMSE", mcse = "rel_rmse_mcse",     band = NULL,            ref = 0),
  coverage     = list(lab = "Coverage",      mcse = "coverage_mcse",     band = c(0.925, 0.975), ref = 0.95),
  se_sd        = list(lab = "SE/SD",         mcse = NULL,                band = c(0.90, 1.10),   ref = 1),
  rel_bias_var = list(lab = "SE²/VAR ratio", mcse = "rel_bias_var_mcse", band = c(0.90, 1.10),   ref = 1),
  typeI        = list(lab = "Type I error",  mcse = "reject_mcse",       band = c(0.025, 0.075), ref = 0.05),
  power        = list(lab = "Power",         mcse = "reject_mcse",       band = NULL,            ref = 0.80))

# acceptable band and reference line of a measure
band_ref <- function(metric) {
  m <- metric_spec[[metric]]
  list(if (!is.null(m$band))
         annotate("rect", xmin = -Inf, xmax = Inf, ymin = m$band[1], ymax = m$band[2],
                  fill = "grey60", alpha = 0.25),
       geom_hline(yintercept = m$ref, linetype = "dotted"))
}

# plus and minus one MCSE around a measure, for marks dodged by `dodge`
mcse_bars <- function(metric, dodge) {
  mcse <- metric_spec[[metric]]$mcse
  if (!is.null(mcse))
    geom_errorbar(aes(ymin = .data[[metric]] - .data[[mcse]], ymax = .data[[metric]] + .data[[mcse]]),
                  position = position_dodge(dodge), width = .25, linewidth = .3,
                  linetype = "solid", colour = "grey30")
}

# one line with points per method for a measure against the misspecification or
# against the sample size. the caller adds the facets.
fig_lines <- function(d, metric, x = "misspec") {
  ggplot(d, aes(.data[[x]], .data[[metric]], shape = method, linetype = method, group = method)) +
    band_ref(metric) + mcse_bars(metric, .5) +
    geom_line(position = position_dodge(.5), linewidth = .5) +
    geom_point(position = position_dodge(.5), size = 1.6) +
    (if (x == "misspec") x_misspec) + method_marks + apa +
    labs(x = c(misspec = "Misspecification", n = "Sample size (N)")[[x]],
         y = metric_spec[[metric]]$lab, shape = NULL, linetype = NULL)
}

# one measure under the correctly specified model against the sample size:
# distribution (rows) by the two variables in `cols` (columns). ylim cuts the
# axis where one method would otherwise flatten the others.
fig_baseline <- function(params, a3v, metric, cols, ylim = NULL) {
  d <- summary_tbl |> filter(parameter %in% params, misspec == "none", a3 %in% a3v) |>
    prep_factors() |> mutate(n = factor(n), parameter = factor(parameter, levels = params))
  fig_lines(d, metric, x = "n") + coord_cartesian(ylim = ylim) +
    facet_grid(reformulate(cols, "distr_exo"), labeller = strip_labels)
}

# type I error or power of the imm under the correctly specified model against
# the sample size: distribution (rows) by reliability (columns)
fig_baseline_reject <- function(metric, a3v) {
  d <- rej_tbl |> filter(parameter == "imm", misspec == "none", a3 == a3v) |>
    prep_factors() |> mutate(n = factor(n))
  fig_lines(d, metric, x = "n") + facet_grid(distr_exo ~ rel)
}

# one measure against the misspecification at the sample sizes in `ns`:
# distribution (rows) by the two variables in `cols` (columns). d is the table
# that holds the measure.
fig_misspec <- function(d, params, a3v, metric, cols, ns) {
  d <- d |> filter(parameter %in% params, a3 == a3v, n %in% ns) |>
    prep_factors() |> mutate(parameter = factor(parameter, levels = params))
  fig_lines(d, metric) + small_x +
    facet_grid(reformulate(cols, "distr_exo"), labeller = strip_labels)
}

# one measure of one parameter at one interaction value and one reliability
# level: distribution (rows) by sample size (columns)
fig_slice <- function(param, a3v, metric, rel_lvl)
  fig_lines(cell(summary_tbl, param, a3v, rel_lvl), metric) +
    facet_grid(distr_exo ~ n, labeller = strip_labels)

# the imm and a3 together as grey bars, one per method: parameter and
# distribution (rows) by sample size (columns). used where the two parameters
# share a scale, that is, for relative RMSE, type I error and power. d is the
# table that holds the measure.
fig_bars <- function(d, metric, a3v, rel_lvl) {
  ggplot(cell(d, c("imm", "a3"), a3v, rel_lvl),
         aes(misspec, .data[[metric]], fill = method, group = method)) +
    band_ref(metric) +
    geom_col(position = position_dodge(.8), width = .7, colour = "grey30", linewidth = .2) +
    mcse_bars(metric, .8) +
    facet_grid(parameter + distr_exo ~ n, labeller = strip_labels) +
    x_misspec + method_fills + apa +
    theme(axis.text.x = element_text(size = 6)) +
    labs(x = "Misspecification", y = metric_spec[[metric]]$lab, fill = NULL)
}

# SE²/VAR ratio of the imm and a3 at a3 = 0.2 and one reliability level:
# distribution (rows) by parameter and sample size (columns). ratios above 10
# occur at low reliability under heavy tails and are left out because they would
# swamp the shared axis. the axis stops at 2.5 so that the ratios close to 1 can
# be read. the ratios above it are not shown: one for LMS at high reliability,
# and those of UPI at N = 200 with chi-square predictors at low reliability.
fig_relvar <- function(rel_lvl) {
  d <- cell(summary_tbl, c("imm", "a3"), 0.2, rel_lvl) |>
    filter(n %in% c(200, 500), rel_bias_var <= 10)
  fig_lines(d, "rel_bias_var") + facet_grid(distr_exo ~ parameter + n, labeller = strip_labels) +
    coord_cartesian(ylim = c(0.5, 2.5))
}

# one measure of the remaining structural coefficients at a3 = 0.2, N = 500 and
# one reliability level: parameter (rows) by distribution (columns). scales =
# "free_y" gives each parameter its own axis.
fig_param <- function(metric, rel_lvl, scales = "fixed") {
  d <- cell(summary_tbl, c("a1", "a2", "b", "cp"), 0.2, rel_lvl) |> filter(n == 500)
  fig_lines(d, metric) +
    facet_grid(parameter ~ distr_exo, labeller = strip_labels, scales = scales)
}

# the measures in `meas` of the loading of the cross-loaded indicator m3, the
# only loading that the misspecifications move, at a3 = 0.2, N = 500 and one
# reliability level: measure (rows, each with its own axis, band and reference
# line) by distribution (columns). SE²/VAR ratios above 10 are left out.
fig_load <- function(rel_lvl, meas) {
  lab <- sapply(metric_spec[meas], `[[`, "lab")
  d   <- cell(summary_tbl, "lm3", 0.2, rel_lvl) |> filter(n == 500) |>
    tidyr::pivot_longer(all_of(meas), names_to = "measure") |>
    filter(!(measure == "rel_bias_var" & value > 10)) |>
    mutate(measure = factor(measure, levels = meas, labels = lab))
  ref <- data.frame(measure = factor(lab, levels = lab),
                    ref = sapply(metric_spec[meas], `[[`, "ref"),
                    t(sapply(metric_spec[meas], `[[`, "band")))
  ggplot(d, aes(misspec, value, shape = method, linetype = method, group = method)) +
    geom_rect(data = ref, inherit.aes = FALSE, xmin = -Inf, xmax = Inf,
              aes(ymin = X1, ymax = X2), fill = "grey60", alpha = 0.25) +
    geom_hline(data = ref, aes(yintercept = ref), linetype = "dotted") +
    geom_line(position = position_dodge(.5), linewidth = .5) +
    geom_point(position = position_dodge(.5), size = 1.4) +
    facet_grid(measure ~ distr_exo, scales = "free_y") + x_misspec + method_marks + apa +
    labs(x = "Misspecification", y = NULL, shape = NULL, linetype = NULL)
}

# saves a figure to plots/ at 150 dpi, w by h inches. rotate = TRUE draws a
# landscape figure sideways on a portrait canvas, so that it fills a portrait
# page.
save_fig <- function(plot, name, w, h, rotate = FALSE) {
  file <- file.path(plot_dir, paste0(name, ".png"))
  if (!rotate) return(ggsave(file, plot, width = w, height = h, dpi = 150))
  png(file, width = h, height = w, units = "in", res = 150)
  on.exit(dev.off())
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(width = grid::unit(w, "in"),
                                    height = grid::unit(h, "in"), angle = 90))
  grid::grid.draw(ggplotGrob(plot))
}
save_bars <- function(d, metric, a3v, rel_lvl)
  save_fig(fig_bars(d, metric, a3v, rel_lvl),
           sprintf("both_%s_a3-%s_%s", metric, a3v, rel_lvl), 10, 12)

# the imm and a3 under the correctly specified model. relative bias has one
# figure per parameter with both nonzero interaction values. SE/SD and coverage
# show the two parameters side by side at a3 = 0.2, and the coverage axis stops
# at 0.6 because the coverage of FSR falls far below that of the other methods.
for (p in c("imm", "a3"))
  save_fig(fig_baseline(p, c(0.2, 0.4), "rel_bias", c("rel", "a3")),
           paste0("baseline_rel_bias_", p), 10, 9)
save_fig(fig_baseline(c("a3", "imm"), 0.2, "se_sd", c("parameter", "rel")),
         "baseline_se_sd", 10, 9)
save_fig(fig_baseline(c("a3", "imm"), 0.2, "coverage", c("parameter", "rel"), ylim = c(0.6, 1)),
         "baseline_coverage", 10, 9)

# type I error (a3 = 0) and power (a3 = 0.2) of the imm as two blocks, which are
# separate plots because the two rates need different axes. they are drawn side
# by side above one shared legend, with the device opened first because the
# grobs are sized for the open device.
blocks <- list(
  fig_baseline_reject("typeI", 0) +
    theme(strip.text.y = element_blank(), strip.background.y = element_blank()),
  fig_baseline_reject("power", 0.2) + coord_cartesian(ylim = c(0, 1)))
png(file.path(plot_dir, "baseline_reject_imm.png"), width = 10, height = 9, units = "in", res = 150)
legend <- ggplotGrob(blocks[[1]])
legend <- legend$grobs[[which(legend$layout$name == "guide-box-bottom")]]
grid::pushViewport(grid::viewport(layout = grid::grid.layout(2, 2, heights = grid::unit(c(8.6, 0.4), "in"))))
for (i in 1:2)
  print(blocks[[i]] + theme(legend.position = "none"),
        vp = grid::viewport(layout.pos.row = 1, layout.pos.col = i))
grid::pushViewport(grid::viewport(layout.pos.row = 2))
grid::grid.draw(legend)
dev.off()

# the imm and a3 under misspecification. the bias induced by an omission barely
# changes with the sample size, so relative bias shows both parameters at
# N = 500. SE/SD and the type I error depend on the sample size, so they show the
# imm at the smallest and the largest one.
save_fig(fig_misspec(summary_tbl, c("a3", "imm"), 0.2, "rel_bias", c("parameter", "rel"), 500),
         "misspec_rel_bias", 10, 9)
save_fig(fig_misspec(summary_tbl, "imm", 0.2, "se_sd", c("rel", "n"), c(100, 500)),
         "misspec_se_sd_imm", 10, 9)
save_fig(fig_misspec(rej_tbl, "imm", 0, "typeI", c("rel", "n"), c(100, 500)),
         "misspec_typeI_imm", 10, 9)

# the imm and a3 in every condition, one figure per reliability level and
# interaction value. relative bias, coverage and SE/SD are lines with one figure
# per parameter, because bars compress coverage and SE/SD and because the bias of
# the two parameters differs too much in scale to share an axis. power is
# available at a3 = 0.2 only.
for (rl in c("low", "high")) {
  for (a in c(0.2, 0.4)) {
    for (p in c("imm", "a3")) for (m in c("rel_bias", "coverage", "se_sd"))
      save_fig(fig_slice(p, a, m, rl), sprintf("%s_%s_a3-%s_%s", p, m, a, rl), 9.5, 8, rotate = TRUE)
    save_bars(summary_tbl, "rel_rmse", a, rl)
  }
  save_bars(rej_tbl, "typeI", 0, rl)
  save_bars(rej_tbl, "power", 0.2, rl)
  save_fig(fig_relvar(rl), paste0("rel_bias_var_a3-0.2_", rl), 9.5, 8, rotate = TRUE)
}

# the remaining structural coefficients and the m3 loading. the figures at high
# reliability go to the supplement. those at low reliability, where the
# approaches differ most, go to the article: the relative bias of the
# coefficients, each on its own axis because that of c' is far larger, and the
# loading without the SE²/VAR ratio. the coverage and the SE/SD ratio of the
# coefficients at low reliability go to the supplement as well.
for (m in c("rel_bias", "coverage", "se_sd"))
  save_fig(fig_param(m, "high"), sprintf("param_%s_a3-0.2_n500_high", m), 12, 10, rotate = TRUE)
for (m in c("coverage", "se_sd"))
  save_fig(fig_param(m, "low"), sprintf("param_%s_a3-0.2_n500_low", m), 12, 10, rotate = TRUE)
save_fig(fig_param("rel_bias", "low", "free_y") + small_x, "param_rel_bias_a3-0.2_n500_low", 10, 9)
save_fig(fig_load("high", c("rel_bias", "coverage", "se_sd", "rel_bias_var")) +
           theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6)),
         "load_lm3_allmeasures_a3-0.2_n500_high", 9, 8)
save_fig(fig_load("low", c("rel_bias", "coverage", "se_sd")) + small_x,
         "load_lm3_a3-0.2_n500_low", 10, 8)

# copula robustness: each measure of the imm or a3 under the gaussian copula
# (robustness run) against its value under the clayton copula (main run), one
# point per condition. points on the identity line are unaffected by the copula.
# estimators are columns, measures are rows, point shape is the distribution.
# the two runs number their conditions differently, so cells are matched on the
# design factors.
cop_pair <- function(f) {
  read <- function(dir) readRDS(file.path(dir, f)) |>
    filter(criterion == "soft", parameter %in% c("imm", "a3"))
  inner_join(read("results_robustness"), read("results_main"), suffix = c("_g", "_c"),
             by = c("n", "a3", "rel", "distr_exo", "misspec", "method", "parameter")) |>
    prep_factors()
}
cop_sum <- cop_pair("summary.rds")
cop_rej <- cop_pair("rejection_imm_a3.rds")
cop_fmt <- function(v) sub("^(-?)0\\.", "\\1.", format(v, drop0trailing = TRUE, trim = TRUE))

# one row of panels for one measure. the axes share one range so that the
# identity line is the diagonal, and the acceptable band is shaded on both.
cop_row <- function(d, value, lab, band = NULL, strips = FALSE, legend = FALSE) {
  d <- d |> mutate(x = .data[[paste0(value, "_c")]], y = .data[[paste0(value, "_g")]])
  lim <- range(d$x, d$y); lim <- lim + c(-1, 1) * 0.04 * diff(lim)
  ggplot(d, aes(x, y, shape = distr_exo)) +
    (if (!is.null(band)) list(
      annotate("rect", xmin = band[1], xmax = band[2], ymin = -Inf, ymax = Inf,
               fill = "grey60", alpha = 0.2),
      annotate("rect", xmin = -Inf, xmax = Inf, ymin = band[1], ymax = band[2],
               fill = "grey60", alpha = 0.2))) +
    geom_abline(linewidth = .3) +
    geom_point(size = .9, stroke = .3, alpha = .7) +
    facet_grid(. ~ method) +
    scale_shape_manual(values = c(Uniform = 1, `t(5)` = 2, `Chi-square (s)` = 0, `Chi-square (d)` = 5)) +
    scale_x_continuous(labels = cop_fmt, n.breaks = 4) +
    scale_y_continuous(labels = cop_fmt, n.breaks = 4) +
    coord_fixed(xlim = lim, ylim = lim, expand = FALSE) +
    labs(x = if (legend) "Clayton copula (main design)", y = paste0(lab, "\nGaussian copula"),
         shape = NULL) +
    guides(shape = guide_legend(override.aes = list(size = 2.2, alpha = 1))) +
    theme_bw(base_size = 9) +
    # the reduced top and bottom margins leave each row enough height for its
    # square panels to reach full width, which keeps the rows aligned
    theme(panel.grid.minor = element_blank(), strip.text = element_text(face = "bold"),
          legend.position = if (legend) "bottom" else "none",
          plot.margin = margin(2, 6, 2, 6)) +
    if (!strips) theme(strip.text = element_blank(), strip.background = element_blank())
}

# the five rows are separate plots because each measure needs its own axis
# range. they are stacked with grid and given a common left margin so that the
# panels line up. relative measures are undefined at a3 = 0.
for (p in c("imm", "a3")) {
  s  <- filter(cop_sum, parameter == p)
  s1 <- filter(s, a3 != 0)
  h  <- c(1.93, 1.68, 1.68, 1.68, 2.43)
  # the device is opened first because the grobs are sized for the open device
  png(file.path(plot_dir, sprintf("copula_agreement_%s.png", p)),
      width = 6.5, height = sum(h), units = "in", res = 300)
  gs <- lapply(list(
    cop_row(s1, "rel_bias", "Relative bias", c(-.10, .10), strips = TRUE),
    cop_row(s1, "rel_rmse", "Relative RMSE"),
    cop_row(s, "coverage", "Coverage", c(.925, .975)),
    cop_row(s, "se_sd", "SE/SD", c(.90, 1.10)),
    cop_row(filter(cop_rej, parameter == p, a3 == 0), "reject", "Type I error",
            c(.025, .075), legend = TRUE)), ggplotGrob)
  w  <- do.call(grid::unit.pmax, lapply(gs, function(g) g$widths))
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(5, 1, heights = grid::unit(h, "in"))))
  for (i in seq_along(gs)) {
    gs[[i]]$widths <- w
    grid::pushViewport(grid::viewport(layout.pos.row = i, layout.pos.col = 1))
    grid::grid.draw(gs[[i]])
    grid::popViewport()
  }
  dev.off()
}
