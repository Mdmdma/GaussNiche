#!/bin/bash
#SBATCH --job-name=usemcmc-build-verify
#SBATCH --account=es_schin
#SBATCH --output=output/build-verify-%j.out
#SBATCH --error=output/build-verify-%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4G
#SBATCH --time=00:20:00

# Verify the vignette-bloat fix: R CMD build the package and report the tarball
# size + whether the 668 MB climate cache / *_files render dirs / bad filenames
# leaked in. --no-build-vignettes keeps this fast and free of vignette Suggests
# (the offline-build correctness of the vignettes is verified separately). R_LIBS
# is set so the package's deps (Rcpp etc.) resolve during build.

set -uo pipefail
cd "$HOME/GaussNiche"; mkdir -p output
source "$HOME/.config/euler/jupyterhub/config_r_studio"
SIF="/cluster/scratch/$USER/rocker_rstudio_4.5.sif"
[[ -f "$SIF" ]] || { echo "Missing $SIF"; exit 2; }
LIB="$HOME/R/rocker-rstudio/4.5"
WORK="/cluster/scratch/$USER/usemcmc_buildverify"
rm -rf "$WORK"; mkdir -p "$WORK"

apptainer exec \
  --bind "/cluster/home/$USER,/cluster/software,/cluster/scratch/$USER" \
  "$SIF" bash -lc "
    export R_LIBS='$LIB'
    cd '$WORK'
    R CMD build '$HOME/USE.MCMC' --no-build-vignettes --no-manual 2>&1 | tail -n 20
    TARBALL=\$(ls -t USE.MCMC_*.tar.gz 2>/dev/null | head -1)
    [ -z \"\$TARBALL\" ] && { echo 'BUILD FAILED'; exit 3; }
    echo '=== tarball size ==='; ls -lh \"\$TARBALL\"
    echo '=== leaked bloat? (climate / *_files / non-portable names) ==='
    tar tzf \"\$TARBALL\" | grep -iE 'climate|_files|[[:space:]]' || echo 'NONE — clean'
    echo '=== gifs still shipped (expected, ~10 MB) ==='
    tar tzf \"\$TARBALL\" | grep -i '\.gif$' || echo '(none)'
  "
