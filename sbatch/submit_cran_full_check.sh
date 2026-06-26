#!/bin/bash
#SBATCH --job-name=usemcmc-full-check
#SBATCH --account=es_schin
#SBATCH --output=output/full-check-%j.out
#SBATCH --error=output/full-check-%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4G
#SBATCH --time=01:30:00

# Full R CMD check --as-cran WITH vignettes built (clears the 'no inst/doc'
# warnings the --no-vignettes run leaves). Installs the vignette Suggests through
# the ETH proxy, then builds (vignettes included) + checks. Does NOT re-run
# document() so it never touches the working tree (safe to run during commits).
# The vignette GIFs are .Rbuildignore'd, so file.exists() is FALSE in the build
# copy and the gif chunks skip; geodata download chunk is eval=FALSE; medium-scale
# ne_countries reads from the installed rnaturalearthdata (no build-time internet).

set -uo pipefail
cd "$HOME/GaussNiche"; mkdir -p output
source "$HOME/.config/euler/jupyterhub/config_r_studio"
module load eth_proxy 2>/dev/null || true
SIF="/cluster/scratch/$USER/rocker_rstudio_4.5.sif"
[[ -f "$SIF" ]] || { echo "Missing $SIF"; exit 2; }
LIB="$HOME/R/rocker-rstudio/4.5"
PKG="$HOME/USE.MCMC"
WORK="/cluster/scratch/$USER/usemcmc_fullcheck"
rm -rf "$WORK"; mkdir -p "$WORK"

apptainer exec \
  --bind "/cluster/home/$USER,/cluster/software,/cluster/scratch/$USER" \
  "$SIF" bash -lc "
    export http_proxy='http://proxy.ethz.ch:3128'  https_proxy='http://proxy.ethz.ch:3128'
    export HTTP_PROXY=\$http_proxy HTTPS_PROXY=\$https_proxy
    export R_LIBS='$LIB'
    export _R_CHECK_FORCE_SUGGESTS_=false
    echo '=== install vignette Suggests (via proxy) ==='
    Rscript -e '
      pkgs <- c(\"rnaturalearth\",\"rnaturalearthdata\",\"virtualspecies\",\"coda\",\"viridis\",\"entropy\")
      miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
      if (length(miss)) { cat(\"installing:\", paste(miss, collapse=\" \"), \"\n\")
        install.packages(miss, repos = \"https://cloud.r-project.org\", lib = \"$LIB\") }
      ok <- vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)
      cat(\"present:\", paste(pkgs[ok], collapse=\" \"), \"\n\")
      cat(\"still missing:\", paste(pkgs[!ok], collapse=\" \"), \"\n\")'
    cd '$WORK'
    echo '=== R CMD build (WITH vignettes) ==='
    R CMD build '$PKG' --no-manual 2>&1 | tail -n 25
    TARBALL=\$(ls -t USE.MCMC_*.tar.gz 2>/dev/null | head -1)
    [ -z \"\$TARBALL\" ] && { echo 'BUILD FAILED (vignettes?)'; exit 3; }
    echo '=== tarball size ==='; ls -lh \"\$TARBALL\"
    echo \"=== R CMD check --as-cran \$TARBALL ===\"
    R CMD check --as-cran --no-manual \"\$TARBALL\" 2>&1 | tail -n 5
    echo '=== 00check.log ==='; cat USE.MCMC.Rcheck/00check.log 2>/dev/null
    cat USE.MCMC.Rcheck/USE.MCMC-Ex.Rout.fail 2>/dev/null || true
  "
