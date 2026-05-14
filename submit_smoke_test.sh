#!/bin/bash
#SBATCH --job-name=gn-smoke
#SBATCH --account=es_schin
#SBATCH --output=output/smoke-%j.out
#SBATCH --error=output/smoke-%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
#SBATCH --time=00:05:00

# Lightweight check that R packages from ~/R/rocker-rstudio/4.5 load cleanly
# inside the rocker apptainer container. Run BEFORE re-dispatching the sweep.

set -euo pipefail
cd "$HOME/GaussNiche"
mkdir -p output

# Sets APPTAINERENV_* (spack LD_LIBRARY_PATH, LD_PRELOAD, R_LD_LIBRARY_PATH)
# that apptainer translates into in-container env for Rscript.
source "$HOME/.config/euler/jupyterhub/config_r_studio"

SIF="/cluster/scratch/$USER/rocker_rstudio_4.5.sif"
[[ -f "$SIF" ]] || { echo "Missing $SIF — run apptainer pull first."; exit 2; }

apptainer exec \
  --bind "/cluster/home/$USER,/cluster/software,/cluster/scratch/$USER" \
  "$SIF" \
  Rscript imports.R
