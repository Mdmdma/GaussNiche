#!/bin/bash
#SBATCH --job-name=gn-run
#SBATCH --account=es_schin
#SBATCH --output=output/gnrun-%j.out
#SBATCH --error=output/gnrun-%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=6G
#SBATCH --time=12:00:00
#
# Generic SLURM wrapper behind `.claude/bin/gnrun` (login-node branch): runs an
# arbitrary GaussNiche R script inside the rocker apptainer container with the
# standard binds + single-threaded BLAS pins. `gnrun` passes the script/mode via
# --export; override resources on the sbatch command line, e.g.
#   sbatch --cpus-per-task=32 --time=06:00:00 \
#          --export=ALL,R_SCRIPT=tune_5d_hsm.R,MODE=full sbatch/submit_gnrun.sh
set -euo pipefail
cd "$HOME/GaussNiche"; mkdir -p output
source "$HOME/.config/euler/jupyterhub/config_r_studio"

# Pin BLAS/OMP threads to 1 so each furrr/PSOCK worker stays single-threaded.
export APPTAINERENV_OMP_NUM_THREADS=1 APPTAINERENV_OPENBLAS_NUM_THREADS=1 \
       APPTAINERENV_MKL_NUM_THREADS=1 APPTAINERENV_OMP_THREAD_LIMIT=1

R_SCRIPT="${R_SCRIPT:?set R_SCRIPT=<script.R> via --export}"
MODE="${MODE:-full}"
NW="${WORKERS:-${SLURM_CPUS_PER_TASK:-4}}"
SIF="/cluster/scratch/$USER/rocker_rstudio_4.5.sif"
[[ -f "$SIF" ]] || { echo "Missing $SIF"; exit 2; }

echo "[gnrun/sbatch] script=$R_SCRIPT mode=$MODE workers=$NW cpus=${SLURM_CPUS_PER_TASK:-?}"
apptainer exec --bind "/cluster/home/$USER,/cluster/software,/cluster/scratch/$USER" \
  "$SIF" Rscript "$R_SCRIPT" "$MODE" "$NW"
