#!/bin/bash
#SBATCH --job-name=gn-bench-cache
#SBATCH --account=es_schin
#SBATCH --output=output/bench-cache-%j.out
#SBATCH --error=output/bench-cache-%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4G
#SBATCH --time=01:00:00

# (1) document() + reinstall USE.MCMC (with the precomputeMcmcEnvironment cache)
#     into the user lib, then (2) run bench_cache.R to verify correctness and
#     measure the speedup. Compiles the Rcpp/RcppArmadillo inner loop, so it runs
#     inside the rocker apptainer container on a compute node (never the login node).
#
#   sbatch submit_bench_cache.sh                 # K=4 species, chain.length=8000
#   sbatch --export=ALL,K=8,CHAIN=12000 submit_bench_cache.sh

set -euo pipefail
cd "$HOME/GaussNiche"
mkdir -p output

# spack APPTAINERENV_* (LD paths, krb5 shim) so terra/sf/USE.MCMC load in-container
source "$HOME/.config/euler/jupyterhub/config_r_studio"

# keep BLAS/OMP single-threaded (deterministic timings, no thread thrash)
export APPTAINERENV_OMP_NUM_THREADS=1
export APPTAINERENV_OPENBLAS_NUM_THREADS=1
export APPTAINERENV_MKL_NUM_THREADS=1

SIF="/cluster/scratch/$USER/rocker_rstudio_4.5.sif"
[[ -f "$SIF" ]] || { echo "Missing $SIF"; exit 2; }

USER_LIB="$HOME/R/rocker-rstudio/4.5"
PKG_DIR="$HOME/USE.MCMC"
K="${K:-4}"
CHAIN="${CHAIN:-8000}"

echo "[1/3] roxygenise (regen man/NAMESPACE) + install USE.MCMC from $PKG_DIR -> $USER_LIB"
# --no-lock guards against a stale 00LOCK-* left by any earlier interrupted install.
rm -rf "$USER_LIB/00LOCK-USE.MCMC"
apptainer exec \
  --bind "/cluster/home/$USER,/cluster/software,/cluster/scratch/$USER" \
  "$SIF" bash -c "
    set -e
    Rscript -e '.libPaths(c(\"$USER_LIB\", .libPaths())); roxygen2::roxygenise(\"$PKG_DIR\")'
    R CMD INSTALL --no-multiarch --no-lock --with-keep.source -l \"$USER_LIB\" \"$PKG_DIR\"
  "

echo "[2/3] run precomputeMcmcEnvironment tests"
# NOT_CRAN=true so the skip_on_cran() guard does not skip the suite outside R CMD check.
export APPTAINERENV_NOT_CRAN=true
apptainer exec \
  --bind "/cluster/home/$USER,/cluster/software,/cluster/scratch/$USER" \
  "$SIF" \
  Rscript -e '.libPaths(c("'"$USER_LIB"'", .libPaths()));
    library(USE.MCMC); library(testthat);
    testthat::test_dir("'"$PKG_DIR"'/tests/testthat",
                       filter = "precomputeMcmcEnvironment", stop_on_failure = TRUE)'

echo "[3/3] benchmark (K=$K species, chain.length=$CHAIN)"
apptainer exec \
  --bind "/cluster/home/$USER,/cluster/software,/cluster/scratch/$USER" \
  "$SIF" \
  Rscript bench_cache.R "$K" "$CHAIN"

echo "[done]"
