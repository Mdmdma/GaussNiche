#!/bin/bash
#SBATCH --job-name=gn-install-vigdeps
#SBATCH --account=es_schin
#SBATCH --output=output/install-vigdeps-%j.out
#SBATCH --error=output/install-vigdeps-%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4G
#SBATCH --time=00:40:00

# Install the CRAN packages the MCMC vignette needs but the container lacks, into
# the user R library (~/R/rocker-rstudio/4.5), so the paper figures can be re-knit.
#
# Deliberately NOT installed: gganimate / gifski / transformr. The vignette's
# animation chunks are eval=FALSE in every case (the GIFs are pre-rendered +
# committed) and gganimate is no longer loaded in `setup`, so they are never
# needed at knit time — and gifski is Rust-based and will not build here.

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
    ## Use the container-configured repo (rocker Rprofile.site, usually a binary
    ## PPM mirror); fall back to a plain CRAN mirror if unset.
    r <- getOption("repos")
    if (is.null(r) || is.na(r["CRAN"]) || r["CRAN"] %in% c("@CRAN@", "")) {
      options(repos = c(CRAN = "https://cloud.r-project.org"))
    }
    pkgs <- c("virtualspecies", "coda", "entropy", "showtext", "sysfonts")
    need <- pkgs[!pkgs %in% rownames(installed.packages())]
    cat("Installing:", paste(need, collapse = ", "), "\n")
    if (length(need)) install.packages(need, lib = "'"$LIB"'")
    miss <- pkgs[!pkgs %in% rownames(installed.packages())]
    if (length(miss)) { cat("STILL MISSING:", paste(miss, collapse = ", "), "\n"); quit(status = 1) }
    cat("ALL PRESENT:", paste(pkgs, collapse = ", "), "\n")
  '
echo "[install] done"
