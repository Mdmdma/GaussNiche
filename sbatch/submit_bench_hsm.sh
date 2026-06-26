#!/bin/bash
#SBATCH --job-name=gn-bench-hsm
#SBATCH --account=es_schin
#SBATCH --output=output/benchhsm-%j.out
#SBATCH --error=output/benchhsm-%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=6G
#SBATCH --time=00:25:00

# Decompose the per-realisation cost of the 5-D HSM pipeline (sampler /
# hypervolume / SDM hook / background-prediction) to rank speedups before the
# settings sweep. Single-threaded timing on a compute node.
set -euo pipefail
cd "$HOME/GaussNiche"; mkdir -p output
source "$HOME/.config/euler/jupyterhub/config_r_studio"
export APPTAINERENV_OMP_NUM_THREADS=1 APPTAINERENV_OPENBLAS_NUM_THREADS=1 APPTAINERENV_MKL_NUM_THREADS=1
SIF="/cluster/scratch/$USER/rocker_rstudio_4.5.sif"
[[ -f "$SIF" ]] || { echo "Missing $SIF"; exit 2; }
apptainer exec --bind "/cluster/home/$USER,/cluster/software,/cluster/scratch/$USER" \
  "$SIF" Rscript bench_hsm.R
