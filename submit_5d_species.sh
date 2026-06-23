#!/bin/bash
#SBATCH --job-name=gn-5d-species
#SBATCH --account=es_schin
#SBATCH --output=output/5d-%j.out
#SBATCH --error=output/5d-%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=4G
#SBATCH --time=04:00:00

# Run the 4 virtual species in 5-D (random + mcmc samplers).
#   smoke validation:  sbatch --export=ALL,MODE=smoke --cpus-per-task=2 --time=00:20:00 submit_5d_species.sh
#   full run:          sbatch submit_5d_species.sh        (MODE defaults to "full")
#   override depth:    sbatch --export=ALL,N_REALIZATIONS=50 submit_5d_species.sh

set -euo pipefail
cd "$HOME/GaussNiche"
mkdir -p output
source "$HOME/.config/euler/jupyterhub/config_r_studio"

# Pin BLAS/OMP threads to 1 so each furrr worker stays single-threaded (avoids
# N_workers^2 thread thrashing — see submit_cutoff_sweep.sh).
export APPTAINERENV_OMP_NUM_THREADS=1
export APPTAINERENV_OPENBLAS_NUM_THREADS=1
export APPTAINERENV_MKL_NUM_THREADS=1
export APPTAINERENV_OMP_THREAD_LIMIT=1

SIF="/cluster/scratch/$USER/rocker_rstudio_4.5.sif"
[[ -f "$SIF" ]] || { echo "Missing $SIF"; exit 2; }

MODE="${MODE:-full}"
NW="${SLURM_CPUS_PER_TASK:-4}"
echo "[launch] mode=$MODE workers=$NW"

apptainer exec \
  --bind "/cluster/home/$USER,/cluster/software,/cluster/scratch/$USER" \
  "$SIF" \
  Rscript 5_highdim_species.R "$MODE" "$NW"
