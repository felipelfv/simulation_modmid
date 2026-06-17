library(rix)

rix(date = "2026-05-11", # available_dates(); 1.8.x quarto (1.9.37 has a pandoc bug)
    r_pkgs = c("covsim", "rvinecopulib",
               "knitr"),                   # knitr: R engine quarto uses for code chunks
                                           # modsem comes from our fork (git_pkgs below)
    system_pkgs = "quarto",                # quarto CLI to render the paper
    # LaTeX deps for the Quarto PDF; extend as the build reports missing .sty
    tex_pkgs = c("amsmath", "framed", "fvextra", "fancyvrb", "booktabs",
                 "caption", "etoolbox", "xcolor", "geometry", "hyperref",
                 "float"),
    git_pkgs = list(
      list(package_name = "lavaan",
           repo_url = "https://github.com/yrosseel/lavaan",
           commit = "8d99872"),
      # modsem from our fork: MplusAutomation moved Imports -> Suggests, which
      # drops the gsubfn -> xvfb -> mesa build chain from the closure (it was
      # only pulled in for modsem's optional Mplus interface, never used here).
      list(package_name = "modsem",
           repo_url = "https://github.com/felipelfv/modsem",
           commit = "b189bfec2a743ccfcca4c197127861d17ba88fd3")
    ),
    project_path = ".",
    overwrite = TRUE
)

# rix bug workaround (per the official docs: "manually inspect and fix the file")
# https://docs.ropensci.org/rix/articles/remote-dependencies.html
# lavaan's and modsem's DESCRIPTION have `Depends: R (>= ...)`, which rix
# mis-parses as an R-package dependency. R lives at pkgs.R, not pkgs.rPackages.R
# current solution (28/05/26): strip the bad line so the nix build succeeds
nix <- readLines("default.nix")
nix <- nix[!grepl("^\\s+R\\s*$", nix)]
writeLines(nix, "default.nix")

