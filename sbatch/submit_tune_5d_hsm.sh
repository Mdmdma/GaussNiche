#!/bin/bash
#SBATCH --job-name=gn-tune-hsm
#SBATCH --account=es_schin
#SBATCH --output=output/tunehsm-%j.out
#SBATCH --error=output/tunehsm-%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem-per-cpu=3G
#SBATCH --time=02:00:00

# uniform+ settings sweep (env_cutoff x species_cutoff) for downstream-HSM
# performance. ONE node, ONE job: the (cell x species) jobs are fanned across the
# node's cores in a single future pool (tune_5d_hsm.R), each job's realisations
# serial inside its worker. The default 4x4 grid = 16 cells x 4 species = 64 jobs,
# so 64 cores runs the whole grid in one wave (~the time of a single job). No job
# arrays -- trivially reproducible: just `sbatch sbatch/submit_tune_5d_hsm.sh`.
# Speedups A1-A4 baked in (uniform+ only / no hypervolume / pc5 / 25 realisations).
#   smoke (parallel, validates the compute_hypervolume flag + grid plumbing):
#     sbatch --export=ALL,MODE=smoke --cpus-per-task=2 --time=00:30:00 sbatch/submit_tune_5d_hsm.sh
#   full:
#     sbatch sbatch/submit_tune_5d_hsm.sh
#   custom grid:
#     sbatch --export=ALL,ENV_CUTOFFS=0.001,0.01,SPECIES_CUTOFFS=0.9,0.6 sbatch/submit_tune_5d_hsm.sh

set -euo pipefail
cd "$HOME/GaussNiche"; mkdir -p output
source "$HOME/.config/euler/jupyterhub/config_r_studio"
export APPTAINERENV_OMP_NUM_THREADS=1
export APPTAINERENV_OPENBLAS_NUM_THREADS=1
export APPTAINERENV_MKL_NUM_THREADS=1
export APPTAINERENV_OMP_THREAD_LIMIT=1
for v in ENV_CUTOFFS SPECIES_CUTOFFS N_REALIZATIONS MAX_PRES MIN_PRES MERGE_PREV; do
  if [[ -n "${!v:-}" ]]; then export "APPTAINERENV_$v=${!v}"; fi
done
SIF="/cluster/scratch/$USER/rocker_rstudio_4.5.sif"
[[ -f "$SIF" ]] || { echo "Missing $SIF"; exit 2; }
MODE="${MODE:-full}"; NW="${SLURM_CPUS_PER_TASK:-4}"
echo "[launch] mode=$MODE workers=$NW env=${ENV_CUTOFFS:-default} species=${SPECIES_CUTOFFS:-default}"
apptainer exec --bind "/cluster/home/$USER,/cluster/software,/cluster/scratch/$USER" \
  "$SIF" Rscript tune_5d_hsm.R "$MODE" "$NW"
