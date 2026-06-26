#!/bin/bash
#SBATCH --job-name=gn-reinstall-test-usemcmc
#SBATCH --account=es_schin
#SBATCH --output=output/reinstall-test-%j.out
#SBATCH --error=output/reinstall-test-%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4G
#SBATCH --time=00:45:00

# Reinstall the edited USE.MCMC package (~/USE.MCMC source) into the container R
# library, recompiling src/mcmc_loop.cpp (the species.cutoff.threshold = 1 uniform
# branch is a C++ logic change, so the recompile is required), then run the
# MCMC-path testthat files that exercise the change. Clears any stale 00LOCK first.

set -euo pipefail
cd "$HOME/GaussNiche"
mkdir -p output
source "$HOME/.config/euler/jupyterhub/config_r_studio"
SIF="/cluster/scratch/$USER/rocker_rstudio_4.5.sif"
[[ -f "$SIF" ]] || { echo "Missing $SIF"; exit 2; }
LIB="$HOME/R/rocker-rstudio/4.5"

apptainer exec \
  --bind "/cluster/home/$USER,/cluster/software,/cluster/scratch/$USER" \
  "$SIF" bash -lc "
    set -e
    rm -rf '$LIB/00LOCK-USE.MCMC'
    R CMD INSTALL --library='$LIB' '$HOME/USE.MCMC'
    Rscript -e '
      Sys.setenv(NOT_CRAN = \"true\")  # else skip_on_cran() skips the MCMC test files
      .libPaths(c(\"$LIB\", .libPaths()))
      library(USE.MCMC)
      cat(\"\n== running MCMC-path tests against the freshly built package ==\n\")
      testthat::test_dir(
        file.path(\"$HOME\", \"USE.MCMC\", \"tests\", \"testthat\"),
        filter = \"mclustDensityFunction|paSamplingMcmc|mcmcSampling|precomputeMcmcEnvironment\",
        stop_on_failure = TRUE,
        reporter = \"summary\")
      cat(\"\n== all selected tests passed ==\n\")'
  "
