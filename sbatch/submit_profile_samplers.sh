#!/bin/bash
#SBATCH --job-name=gn-profile-samplers
#SBATCH --account=es_schin
#SBATCH --output=output/profile-%j.out
#SBATCH --error=output/profile-%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=6G
#SBATCH --time=00:45:00

# Decompose the computational cost of uniform (USE paSampling) vs uniform+
# (paSamplingMcmc) into model-fitting vs sampling/chain, single-threaded, on a
# compute node. See profile_samplers.R for the breakdown.
set -euo pipefail
cd "$HOME/GaussNiche"; mkdir -p output
source "$HOME/.config/euler/jupyterhub/config_r_studio"
export APPTAINERENV_OMP_NUM_THREADS=1 APPTAINERENV_OPENBLAS_NUM_THREADS=1 APPTAINERENV_MKL_NUM_THREADS=1
SIF="/cluster/scratch/$USER/rocker_rstudio_4.5.sif"
[[ -f "$SIF" ]] || { echo "Missing $SIF"; exit 2; }
apptainer exec --bind "/cluster/home/$USER,/cluster/software,/cluster/scratch/$USER" \
  "$SIF" Rscript profile_samplers.R
