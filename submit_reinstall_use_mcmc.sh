#!/bin/bash
#SBATCH --job-name=gn-reinstall-usemcmc
#SBATCH --account=es_schin
#SBATCH --output=output/reinstall-%j.out
#SBATCH --error=output/reinstall-%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4G
#SBATCH --time=00:30:00

# Reinstall the edited USE.MCMC package (~/USE.MCMC source) into the container
# R library so changes (e.g. paSamplingMcmc's precomputed.env arg) take effect.
# Compiles src/mcmc_loop.cpp via RcppArmadillo using the container toolchain
# (no spack cc_wrap needed). Clears any stale 00LOCK first.

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
    rm -rf '$LIB/00LOCK-USE.MCMC'
    R CMD INSTALL --library='$LIB' '$HOME/USE.MCMC'
    Rscript -e '.libPaths(c(\"$LIB\", .libPaths()));
      cat(\"precomputed.env arg present:\",
          \"precomputed.env\" %in% names(formals(USE.MCMC::paSamplingMcmc)), \"\n\")'
  "
