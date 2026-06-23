#!/bin/bash
#SBATCH --job-name=gn-final-stack
#SBATCH --account=es_schin
#SBATCH --output=output/final-stack-%j.out
#SBATCH --error=output/final-stack-%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=8G
#SBATCH --time=00:20:00

# Lock the lean_natural 5-D environment: subset 12 layers from the cached
# candidate stack, run PCA, save the 5-D background. No downloads (uses cache).

set -euo pipefail
cd "$HOME/GaussNiche"
mkdir -p output
source "$HOME/.config/euler/jupyterhub/config_r_studio"

SIF="/cluster/scratch/$USER/rocker_rstudio_4.5.sif"
[[ -f "$SIF" ]] || { echo "Missing $SIF"; exit 2; }

apptainer exec \
  --bind "/cluster/home/$USER,/cluster/software,/cluster/scratch/$USER" \
  "$SIF" \
  Rscript build_final_stack.R
