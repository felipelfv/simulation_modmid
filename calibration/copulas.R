# fixed exogenous populations for the simulation: one vine distribution per
# (distribution x copula family), built once here with the VITA algorithm
# (Gronneberg & Foldnes, 2017; covsim::vita) and saved to vines.rds.
# gen_exo() (helpers_mc.R) and calibrate.R use these saved objects, so every
# replication in a cell samples the same population.
#
# vita takes the margins and the target covariance matrix and searches for the
# pair-copula parameter that reproduces the target. the search evaluates the
# correlation on Nmax simulated draws, so the result depends on the seed. the
# seed is fixed here and the resulting family, rotation and parameter are
# printed and checked, so that the saved population is reproducible and known.
#
# families: "clayton" is the main run, "gaussian" the copula robustness run.
# the normal condition is the bivariate normal (gaussian copula, normal
# margins) in both runs.
#
# rotations: clayton expresses positive dependence at rotation 0. the
# opposite-skew chi-square condition needs Cor(X, W') = -rho before W = -W' is
# applied in gen_exo(), which clayton reaches at the 90-degree rotation. vita
# can also reach a positive target with the 180-degree rotation (upper instead
# of lower tail dependence), so the rotation of every saved vine is checked
# against the expected one below. the gaussian copula has no rotation.

library(covsim); library(rvinecopulib)

rho    <- 0.3
t_df   <- 5L
chi_df <- 2L

Nmax    <- 1e7      # draws per evaluation inside vita's search
seed    <- 20261111L
n_check <- 2e7      # monte carlo draws used only to check the saved object

# raw margins as stored in the vine (what rvine() draws from). gen_exo()
# centers and scales t and chisq afterwards. correlation is invariant to that,
# so the target is set at the native scale of each margin.
margins <- list(
  normal  = list(distr = "norm",  mean = 0, sd = 1),
  uniform = list(distr = "unif",  min = -sqrt(3), max = sqrt(3)),
  t5      = list(distr = "t",     df = t_df),
  chisq   = list(distr = "chisq", df = chi_df))
nat_var <- c(normal = 1, uniform = 1, t5 = t_df / (t_df - 2), chisq = 2 * chi_df)

# the four nonnormal design cells, the correlation the copula has to produce
# between the two margins before any negation, and the expected rotation
cells <- data.frame(
  cell     = c("uniform", "t5", "chisq_same", "chisq_diff"),
  margin   = c("uniform", "t5", "chisq",      "chisq"),
  target   = c(rho,       rho,  rho,          -rho),
  stringsAsFactors = FALSE)
expected_rotation <- list(clayton = c(0L, 0L, 0L, 90L), gaussian = c(0L, 0L, 0L, 0L))

fit_vita <- function(margin, target, family) {
  m <- margins[[margin]]
  sigma <- nat_var[[margin]] * matrix(c(1, target, target, 1), 2, 2)
  set.seed(seed)
  vita(list(m, m), sigma, family_set = family, Nmax = Nmax, verbose = FALSE, cores = 1)
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

row_for <- function(cell, family, bc, target, chk) {
  data.frame(cell = cell, family = bc$family, rotation = bc$rotation, parameter = bc$parameters[1],
             ktau = par_to_ktau(bc$family, bc$rotation, bc$parameters), target = target,
             check_cor = unname(chk["cor"]), check_mcse = unname(chk["cor_mcse"]),
             mean1 = unname(chk["mean1"]), var1 = unname(chk["var1"]),
             mean2 = unname(chk["mean2"]), var2 = unname(chk["var2"]))
}

vines <- list(clayton = list(), gaussian = list())
tbl   <- list()

for (family in c("clayton", "gaussian")) {
  for (i in seq_len(nrow(cells))) {
    cell <- cells$cell[i]; target <- cells$target[i]
    vd <- fit_vita(cells$margin[i], target, if (family == "gaussian") "gauss" else family)
    bc <- vd$copula$pair_copulas[[1]][[1]]
    chk <- check_vine(vd)
    cat(sprintf("%-10s %-8s -> %s rotation %3d  parameter %.6f  check cor %.4f\n",
                cell, family, bc$family, bc$rotation, bc$parameters[1], chk["cor"]))
    if (bc$rotation != expected_rotation[[family]][i])
      stop(sprintf("%s / %s: vita selected rotation %d, expected %d. change the seed and rerun.",
                   cell, family, bc$rotation, expected_rotation[[family]][i]))
    vines[[family]][[cell]] <- vd
    tbl[[length(tbl) + 1]] <- row_for(cell, family, bc, target, chk)
  }
}

# the normal condition: bivariate normal, identical in both runs
vines$normal <- fit_vita("normal", rho, "gauss")
bc <- vines$normal$copula$pair_copulas[[1]][[1]]
chk <- check_vine(vines$normal)
cat(sprintf("%-10s %-8s -> %s rotation %3d  parameter %.6f  check cor %.4f\n",
            "normal", "gaussian", bc$family, bc$rotation, bc$parameters[1], chk["cor"]))
tbl[[length(tbl) + 1]] <- row_for("normal", "gaussian", bc, rho, chk)

tbl <- do.call(rbind, tbl)
rownames(tbl) <- NULL
print(tbl, digits = 5)

saveRDS(list(vines = vines, summary = tbl, rho = rho, t_df = t_df, chi_df = chi_df,
             vita = list(Nmax = Nmax, seed = seed)), "vines.rds")
