library(rix)

rix(date = "2026-10-13", # available_dates(); snapshot with CRAN lavaan 0.7-3
    r_pkgs = c("lavaan",                   # CRAN 0.7-3 (as shipped by the snapshot)
               "covsim", "rvinecopulib",
               "tidyr", "ggplot2", "simhelpers", "dplyr",
               "knitr"),                   # knitr: R engine quarto uses for code chunks
    system_pkgs = "quarto",                # quarto CLI to render the paper
    # LaTeX deps for the Quarto PDF; extend as the build reports missing .sty
    tex_pkgs = c("amsmath", "framed", "fvextra", "fancyvrb", "booktabs",
                 "caption", "etoolbox", "xcolor", "geometry", "hyperref",
                 "float", "pgf", "standalone",
                 # apaquarto (apa7.cls) + the packages it and its template load
                 "apa7", "scalerel", "endfloat", "substr", "xstring",
                 "threeparttable",
                 "threeparttablex", "multirow", "colortbl", "xpatch", "lineno",
                 "fontawesome5", "newtx", "fontaxes", "pbalance", "tcolorbox",
                 "tikzfill", "environ", "pdfcol", "mathspec", "xecjk",
                 "pdflscape", "babel-english"),
    # modsem: CRAN 1.0.24 (2026-11-02) post-dates the snapshot (which has 1.0.23);
    # this commit is source-identical to the CRAN 1.0.24 tarball (R/, src/, NAMESPACE)
    git_pkgs = list(
      list(package_name = "modsem",
           repo_url = "https://github.com/Kss2k/modsem",
           commit = "c72311915236f1d96bce4e087f45b43a3e08f7ac")
    ),
    project_path = ".",
    overwrite = TRUE
)

# rix bug workaround (per the official docs: "manually inspect and fix the file")
# https://docs.ropensci.org/rix/articles/remote-dependencies.html
# modsem's DESCRIPTION has `Depends: R (>= ...)`, which rix
# mis-parses as an R-package dependency. R lives at pkgs.R, not pkgs.rPackages.R
# current solution (28/05/26): strip the bad line so the nix build succeeds
nix <- readLines("default.nix")
nix <- nix[!grepl("^\\s+R\\s*$", nix)]
writeLines(nix, "default.nix")

