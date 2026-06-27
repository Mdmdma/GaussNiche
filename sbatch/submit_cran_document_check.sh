#!/bin/bash
#SBATCH --job-name=usemcmc-doc-check
#SBATCH --account=es_schin
#SBATCH --output=output/doc-check-%j.out
#SBATCH --error=output/doc-check-%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4G
#SBATCH --time=00:45:00

# Regenerate man/*.Rd + NAMESPACE from the edited roxygen, then run the
# authoritative R CMD check --as-cran on the result. R_LIBS is set so the
# package's deps resolve during the check's internal install (fixes the earlier
# 'no package called Rcpp' artifact). Vignettes are skipped here (--no-vignettes
# + --no-build-vignettes) to avoid the vignette Suggests/internet path; the
# vignette fixes were verified separately. _R_CHECK_FORCE_SUGGESTS_=false so a
# missing optional Suggests reports as a NOTE rather than failing the check.

set -uo pipefail
cd "$HOME/GaussNiche"; mkdir -p output
source "$HOME/.config/euler/jupyterhub/config_r_studio"
source "$HOME/GaussNiche/sbatch/eth_proxy.sh"   # outbound internet via the ETH proxy
SIF="/cluster/scratch/$USER/rocker_rstudio_4.5.sif"
[[ -f "$SIF" ]] || { echo "Missing $SIF"; exit 2; }
LIB="$HOME/R/rocker-rstudio/4.5"
PKG="$HOME/USE.MCMC"
WORK="/cluster/scratch/$USER/usemcmc_doccheck"
rm -rf "$WORK"; mkdir -p "$WORK"

apptainer exec \
  --bind "/cluster/home/$USER,/cluster/software,/cluster/scratch/$USER" \
  "$SIF" bash -lc "
    export R_LIBS='$LIB'
    export _R_CHECK_FORCE_SUGGESTS_=false
    echo '=== regenerate docs (roxygen2) ==='
    Rscript -e '
      if (requireNamespace(\"devtools\", quietly = TRUE)) {
        devtools::document(\"$PKG\")
      } else {
        roxygen2::roxygenise(\"$PKG\")
      }
      cat(\"DOCUMENT DONE\n\")'
    cd '$WORK'
    echo '=== R CMD build (no vignettes) ==='
    R CMD build '$PKG' --no-build-vignettes --no-manual 2>&1 | tail -n 15
    TARBALL=\$(ls -t USE.MCMC_*.tar.gz 2>/dev/null | head -1)
    [ -z \"\$TARBALL\" ] && { echo 'BUILD FAILED'; exit 3; }
    echo \"=== R CMD check --as-cran \$TARBALL ===\"
    R CMD check --as-cran --no-manual --no-vignettes \"\$TARBALL\" 2>&1 | tail -n 5
    echo '=== 00check.log ==='
    cat USE.MCMC.Rcheck/00check.log 2>/dev/null
    echo '=== examples output (if any failures) ==='
    cat USE.MCMC.Rcheck/USE.MCMC-Ex.Rout.fail 2>/dev/null || echo '(no example failures)'
  "
