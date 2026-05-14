#!/bin/bash
#SBATCH --job-name=gn-cutoff-sweep
#SBATCH --account=es_schin
#SBATCH --output=output/slurm-%A_%a.out
#SBATCH --error=output/slurm-%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH --mem-per-cpu=2G
#SBATCH --time=01:00:00
#SBATCH --array=0-8

# Sweep species.cutoff.threshold in {0.1, 0.2, ..., 0.9} via SLURM array.
# Each task runs 4_multi_species_comparison.R with one cutoff value.
#
# Override realisations and workers from the sbatch command line:
#   sbatch --export=ALL,N_REALIZATIONS=10,N_WORKERS=8 submit_cutoff_sweep.sh
# Defaults below are intentionally modest — bump for production runs.

set -euo pipefail

cd "$HOME/GaussNiche"
mkdir -p output

# Source config_r_studio for its APPTAINERENV_* exports (spack
# LD_LIBRARY_PATH, libkrb5support LD_PRELOAD, R_LD_LIBRARY_PATH with shim).
# Apptainer translates those into the in-container env for Rscript below.
# Host R cannot load these packages — see plan + skill for the GLIBCXX_3.4.32
# / dangling-shim analysis.
source "$HOME/.config/euler/jupyterhub/config_r_studio"

SIF="/cluster/scratch/$USER/rocker_rstudio_4.5.sif"
[[ -f "$SIF" ]] || { echo "Missing $SIF — run apptainer pull first."; exit 2; }

# Pin OMP/BLAS threads to 1 inside the container. Without this, each PSOCK
# worker auto-detects nproc and spawns ~N_WORKERS internal threads → N_WORKERS²
# threads competing for N_WORKERS cores → severe thrashing (caused the 1h
# TIMEOUT at N_REALIZATIONS=50, N_WORKERS=20). One thread per worker is what
# we want when furrr is the parallelism layer.
export APPTAINERENV_OMP_NUM_THREADS=1
export APPTAINERENV_OPENBLAS_NUM_THREADS=1
export APPTAINERENV_MKL_NUM_THREADS=1
export APPTAINERENV_OMP_THREAD_LIMIT=1

CUTOFFS=(0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9)
CUTOFF="${CUTOFFS[$SLURM_ARRAY_TASK_ID]}"
N_REALIZATIONS="${N_REALIZATIONS:-10}"
N_WORKERS="${N_WORKERS:-${SLURM_CPUS_PER_TASK:-4}}"

echo "[launch] array=${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}" \
     "cutoff=${CUTOFF} n_realizations=${N_REALIZATIONS} n_workers=${N_WORKERS}"

apptainer exec \
  --bind "/cluster/home/$USER,/cluster/software,/cluster/scratch/$USER" \
  "$SIF" \
  Rscript 4_multi_species_comparison.R "$CUTOFF" "$N_REALIZATIONS" "$N_WORKERS"
