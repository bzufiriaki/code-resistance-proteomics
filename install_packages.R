# Install the R packages required by the scripts in this repository.
#
# NOTE: the scripts were developed with xgboost 1.7.x. xgboost >= 2.0 changed
# the R interface (e.g. `watchlist` -> `evals`) and may give different
# results, so a 1.7.x version is installed explicitly.

options(repos = c(CRAN = "https://cloud.r-project.org"))

cran_pkgs <- c(
  "caret", "cvTools", "dplyr", "tidyr", "ggplot2", "svglite", "future.apply",
  "dynamicTreeCut", "dendextend", "pheatmap", "viridis",
  "lmerTest", "ggnewscale", "tibble", "remotes"
)
missing <- setdiff(cran_pkgs, rownames(installed.packages()))
if (length(missing) > 0) install.packages(missing)

if (!requireNamespace("xgboost", quietly = TRUE) ||
    packageVersion("xgboost") >= "2.0.0") {
  remotes::install_version("xgboost", version = "1.7.8.1", upgrade = "never")
}
