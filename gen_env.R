library(rix)

rix(date = "2026-05-25", # available_dates()
    r_pkgs = c("modsem", "covsim", "rvinecopulib", "ggplot2"),
    git_pkgs = list(
      package_name = "lavaan",
      repo_url = "https://github.com/yrosseel/lavaan",
      commit = "52d303c"
    ),
    project_path = ".",
    overwrite = TRUE
)

# rix bug workaround (per the official docs: "manually inspect and fix the file")
# https://docs.ropensci.org/rix/articles/remote-dependencies.html
# lavaan's DESCRIPTION has `Depends: R (>= 3.4)`, which rix mis-parses as an
# R-package dependency. R lives at pkgs.R, not pkgs.rPackages.R 
# current solution (28/05/26): strip the bad line so the nix build succeeds
nix <- readLines("default.nix")
nix <- nix[!grepl("^\\s+R\\s*$", nix)]
writeLines(nix, "default.nix")

