# fixed exogenous populations for the simulation: one vine distribution per
# (distribution x copula family), built once here and saved to vines.rds.
# gen_exo() (helpers_mc.R) and calibrate.R use these saved objects, so every
# replication in a cell samples the same population.

# for each nonnormal margin the copula parameter is solved so that the pearson
# correlation of the transformed variables equals rho (0.3 in magnitude). the
# correlation is computed for each candidate parameter as the integral of
# q(u) q(v) c(u, v) over the unit square by nested adaptive quadrature, and
# uniroot finds the parameter. no simulation is involved in the definition of
# the population. the saved object is then checked on an independent monte
# carlo draw, which is reported in the summary table and more commonly used.

# families: "clayton" is the main run, "gaussian" the copula robustness run.
# the normal condition is the bivariate normal (gaussian copula, parameter rho)
# in both runs and needs no solve.

# rotations: clayton can only express positive dependence at rotation 0. the
# opposite-skew chi-square condition needs Cor(X, W') = -rho before W = -W' is
# applied in gen_exo(), so it uses the 90-degree rotation. 
# the gaussian copula has no rotation.

library(rvinecopulib)

rho    <- 0.3
t_df   <- 5L
chi_df <- 2L

quad_tol <- 1e-9   # relative tolerance of each one-dimensional integration
root_tol <- 1e-8   # tolerance of uniroot on the copula parameter
n_check  <- 2e7    # monte carlo draws used only to check the saved object

# raw margins as stored in the vine (what rvine() draws from)
margins <- list(
  normal  = list(distr = "norm",  mean = 0, sd = 1),
  uniform = list(distr = "unif",  min = -sqrt(3), max = sqrt(3)),
  t5      = list(distr = "t",     df = t_df),
  chisq   = list(distr = "chisq", df = chi_df))

# standardized quantile functions (mean 0, variance 1), so that the integral in
# cor_quad() is the pearson correlation itself. correlation is invariant to the
# standardization, so the solved parameter applies to the raw margins above.
qstd <- list(
  uniform = function(u) qunif(u, -sqrt(3), sqrt(3)),
  t5      = function(u) qt(u, t_df) / sqrt(t_df / (t_df - 2)),
  chisq   = function(u) (qchisq(u, chi_df) - chi_df) / sqrt(2 * chi_df))

# the four nonnormal design cells and the correlation the copula has to produce
# between the two margins before any negation
targets <- c(uniform = rho, t5 = rho, chisq_same = rho, chisq_diff = -rho)
margin_of <- function(cell) sub("chisq_.*", "chisq", cell)

# pearson correlation of (q(U), q(V)) under the copula with parameter theta,
# by nested adaptive quadrature over the unit square
cor_quad <- function(theta, family, rotation, qf) {
  bc <- bicop_dist(family, rotation = rotation, parameters = theta)
  inner <- function(u) vapply(u, function(ui) integrate(function(v)
      qf(ui) * qf(v) * dbicop(cbind(rep(ui, length(v)), v), bc),
      0, 1, rel.tol = quad_tol, subdivisions = 2000L, stop.on.error = FALSE)$value, numeric(1))
  integrate(inner, 0, 1, rel.tol = quad_tol, subdivisions = 2000L, stop.on.error = FALSE)$value
}

solve_theta <- function(family, rotation, qf, target) {
  interval <- if (family == "clayton") c(0.02, 10) else if (target > 0) c(0.01, 0.99) else c(-0.99, -0.01)
  uniroot(function(th) cor_quad(th, family, rotation, qf) - target,
          interval = interval, tol = root_tol)$root
}

build_vine <- function(family, rotation, theta, m) {
  vine_dist(list(m, m),
            pair_copulas = list(list(bicop_dist(family, rotation = rotation, parameters = theta))),
            structure = dvine_structure(1:2))
}

# check of the saved object on an independent draw: correlation with a monte
# carlo standard error from 20 blocks, and the first two moments of each margin
check_vine <- function(vd) {
  set.seed(20261110L)
  e <- rvine(n_check, vd)
  blocks <- split(seq_len(n_check), rep(1:20, each = n_check / 20))
  cors <- vapply(blocks, function(i) cor(e[i, 1], e[i, 2]), numeric(1))
  c(cor = cor(e[, 1], e[, 2]), cor_mcse = sd(cors) / sqrt(20),
    mean1 = mean(e[, 1]), var1 = var(e[, 1]), mean2 = mean(e[, 2]), var2 = var(e[, 2]))
}

row_for <- function(cell, family, rotation, theta, target, chk) {
  data.frame(cell = cell, family = family, rotation = rotation, parameter = theta,
             ktau = par_to_ktau(family, rotation, theta), target = target,
             check_cor = unname(chk["cor"]), check_mcse = unname(chk["cor_mcse"]),
             mean1 = unname(chk["mean1"]), var1 = unname(chk["var1"]),
             mean2 = unname(chk["mean2"]), var2 = unname(chk["var2"]))
}

vines <- list(clayton = list(), gaussian = list())
tbl   <- list()

for (family in c("clayton", "gaussian")) {
  for (cell in names(targets)) {
    target   <- targets[[cell]]
    rotation <- if (family == "clayton" && target < 0) 90L else 0L
    theta <- solve_theta(family, rotation, qstd[[margin_of(cell)]], target)
    vd    <- build_vine(family, rotation, theta, margins[[margin_of(cell)]])
    chk   <- check_vine(vd)
    vines[[family]][[cell]] <- vd
    tbl[[length(tbl) + 1]] <- row_for(cell, family, rotation, theta, target, chk)
    cat(sprintf("%-10s %-8s rotation %3d  theta %.6f  check cor %.4f\n",
                cell, family, rotation, theta, chk["cor"]))
  }
}

# the normal condition: bivariate normal, identical in both runs
vines$normal <- build_vine("gaussian", 0L, rho, margins$normal)
tbl[[length(tbl) + 1]] <- row_for("normal", "gaussian", 0L, rho, rho, check_vine(vines$normal))

tbl <- do.call(rbind, tbl)
rownames(tbl) <- NULL
print(tbl, digits = 5)

saveRDS(list(vines = vines, summary = tbl, rho = rho, t_df = t_df, chi_df = chi_df), 
        "vines.rds")