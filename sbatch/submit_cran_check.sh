#!/bin/bash
#SBATCH --job-name=usemcmc-cran-check
#SBATCH --account=es_schin
#SBATCH --output=output/cran-check-%j.out
#SBATCH --error=output/cran-check-%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4G
#SBATCH --time=01:00:00

# Authoritative CRAN-readiness check: R CMD build + R CMD check --as-cran on the
# USE.MCMC source, inside the apptainer container. Runs in a scratch workdir (not
# $HOME) so artifacts don't hit the home quota. --no-manual skips the PDF manual
# (no full LaTeX in the container). If the full build fails (e.g. a vignette that
# needs internet), it retries without vignettes so the rest of the check still runs.

set -uo pipefail
cd "$HOME/GaussNiche"
mkdir -p output
source "$HOME/.config/euler/jupyterhub/config_r_studio"
source "$HOME/GaussNiche/sbatch/eth_proxy.sh"   # outbound internet via the ETH proxy
SIF="/cluster/scratch/$USER/rocker_rstudio_4.5.sif"
[[ -f "$SIF" ]] || { echo "Missing $SIF"; exit 2; }
LIB="$HOME/R/rocker-rstudio/4.5"
WORK="/cluster/scratch/$USER/usemcmc_crancheck"
rm -rf "$WORK"; mkdir -p "$WORK"

apptainer exec \
  --bind "/cluster/home/$USER,/cluster/software,/cluster/scratch/$USER" \
  "$SIF" bash -lc "
    export R_LIBS_USER='$LIB'
    export _R_CHECK_CRAN_INCOMING_=true
    export _R_CHECK_FORCE_SUGGESTS_=false
    cd '$WORK'
    echo '=== R CMD build (full) ==='
    R CMD build '$HOME/USE.MCMC' --no-manual 2>&1 | tee build.log
    TARBALL=\$(ls -t USE.MCMC_*.tar.gz 2>/dev/null | head -1)
    if [ -z \"\$TARBALL\" ]; then
      echo '=== full build failed; retrying --no-build-vignettes ==='
      R CMD build '$HOME/USE.MCMC' --no-build-vignettes --no-manual 2>&1 | tee -a build.log
      TARBALL=\$(ls -t USE.MCMC_*.tar.gz 2>/dev/null | head -1)
    fi
    [ -z \"\$TARBALL\" ] && { echo 'BUILD FAILED, no tarball'; exit 3; }
    echo \"=== R CMD check --as-cran \$TARBALL ===\"
    R CMD check --as-cran --no-manual \"\$TARBALL\" 2>&1 | tee check.log
    echo '=== 00check.log ==='
    cat USE.MCMC.Rcheck/00check.log 2>/dev/null
  "
