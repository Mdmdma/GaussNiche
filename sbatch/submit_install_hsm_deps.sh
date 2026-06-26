#!/bin/bash
#SBATCH --job-name=gn-install-hsmdeps
#SBATCH --account=es_schin
#SBATCH --output=output/install-hsmdeps-%j.out
#SBATCH --error=output/install-hsmdeps-%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4G
#SBATCH --time=00:40:00

# Install the CRAN packages the downstream-HSM arm (run_5d_hsm.R / hsm_eval.R /
# hsm_plots.R) needs but the rocker container lacks, into the user R library
# (~/R/rocker-rstudio/4.5). Java-free by design: maxnet (pure-R Maxent) replaces
# dismo+rJava+maxent.jar, so NO JDK / javareconf is required. mgcv (GAM) ships
# with base R and is NOT installed here.
#
#   ranger    RF (regression forest -> OOB r.squared)
#   gbm       BRT (bernoulli)            [archived on CRAN; PPM binary usually OK]
#   maxnet    Maxent (pure R)
#   glmnet    maxnet backend (Fortran)
#   Metrics   convenience metric helpers (optional; hsm_eval.R is self-contained)
#   dunn.test downstream pairwise significance test (hsm_plots.R)

set -euo pipefail
cd "$HOME/GaussNiche"
mkdir -p output
source "$HOME/.config/euler/jupyterhub/config_r_studio"
SIF="/cluster/scratch/$USER/rocker_rstudio_4.5.sif"
[[ -f "$SIF" ]] || { echo "Missing $SIF"; exit 2; }
LIB="$HOME/R/rocker-rstudio/4.5"

# proxy passthrough so install.packages can reach the configured repo
export APPTAINERENV_http_proxy="${http_proxy:-${HTTP_PROXY:-http://proxy.ethz.ch:3128}}"
export APPTAINERENV_https_proxy="${https_proxy:-${HTTPS_PROXY:-http://proxy.ethz.ch:3128}}"
export APPTAINERENV_HTTP_PROXY="${APPTAINERENV_http_proxy}"
export APPTAINERENV_HTTPS_PROXY="${APPTAINERENV_https_proxy}"

apptainer exec \
  --bind "/cluster/home/$USER,/cluster/software,/cluster/scratch/$USER" \
  "$SIF" Rscript -e '
    .libPaths(c("'"$LIB"'", .libPaths()))
    r <- getOption("repos")
    if (is.null(r) || is.na(r["CRAN"]) || r["CRAN"] %in% c("@CRAN@", "")) {
      options(repos = c(CRAN = "https://cloud.r-project.org"))
    }
    pkgs <- c("ranger", "gbm", "maxnet", "glmnet", "Metrics", "dunn.test")
    need <- pkgs[!pkgs %in% rownames(installed.packages())]
    cat("Installing:", paste(need, collapse = ", "), "\n")
    if (length(need)) install.packages(need, lib = "'"$LIB"'")
    ## A package can install yet fail to load (truncated overlay, missing native
    ## dep): verify each loads, since a silently-missing backend becomes an
    ## all-NA column that masquerades as a sampler effect.
    bad <- character(0)
    for (p in pkgs) {
      okp <- requireNamespace(p, quietly = TRUE)
      cat(sprintf("  %-10s %s\n", p, if (okp) "OK" else "FAILED TO LOAD"))
      if (!okp) bad <- c(bad, p)
    }
    if (length(bad)) { cat("STILL MISSING/BROKEN:", paste(bad, collapse = ", "), "\n"); quit(status = 1) }
    cat("ALL HSM DEPS PRESENT + LOADABLE:", paste(pkgs, collapse = ", "), "\n")
  '
echo "[install] done"
