#!/bin/bash
#SBATCH --job-name=gn-5d-hsm
#SBATCH --account=es_schin
#SBATCH --output=output/5dhsm-%j.out
#SBATCH --error=output/5dhsm-%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=6G
#SBATCH --time=08:00:00

# Downstream-HSM arm of the 5-D experiment: the 4 species in 5-D (RND / buffer-out
# / uniform+) over R Bernoulli realisations, fitting the 5 USE model families
# (GLM/GAM/RF/BRT/Maxent, Java-free) to each realisation on two predictor sets
# (5 PCs and 12 raw layers). Built on run_5d_experiment.R's engine via the
# sdm_hook; reuses each PA draw single-pass so pa_mcmc is not re-paid.
#
# Install deps first:  sbatch sbatch/submit_install_hsm_deps.sh
#   smoke (REQUIRED before full -- the only parallel test of furrr worker pkg loading):
#     sbatch --export=ALL,MODE=smoke --cpus-per-task=2 --time=00:40:00 sbatch/submit_5d_hsm.sh
#   full:
#     sbatch sbatch/submit_5d_hsm.sh        (MODE defaults to "full", R=50)
#   tune:
#     sbatch --export=ALL,N_REALIZATIONS=50,MAX_PRES=300,MIN_PRES=12 sbatch/submit_5d_hsm.sh

set -euo pipefail
cd "$HOME/GaussNiche"
mkdir -p output
source "$HOME/.config/euler/jupyterhub/config_r_studio"

# Pin BLAS/OMP threads to 1 so each furrr worker (ranger/gbm/maxnet) stays
# single-threaded and does not oversubscribe the allocation.
export APPTAINERENV_OMP_NUM_THREADS=1
export APPTAINERENV_OPENBLAS_NUM_THREADS=1
export APPTAINERENV_MKL_NUM_THREADS=1
export APPTAINERENV_OMP_THREAD_LIMIT=1

# Forward optional overrides into the container (only when set).
for v in N_REALIZATIONS MAX_PRES SPECIES_CUTOFF BUFFER_KM MIN_PRES; do
  if [[ -n "${!v:-}" ]]; then export "APPTAINERENV_$v=${!v}"; fi
done

SIF="/cluster/scratch/$USER/rocker_rstudio_4.5.sif"
[[ -f "$SIF" ]] || { echo "Missing $SIF"; exit 2; }

MODE="${MODE:-full}"
NW="${SLURM_CPUS_PER_TASK:-4}"
echo "[launch] mode=$MODE workers=$NW n_realizations=${N_REALIZATIONS:-default} max_pres=${MAX_PRES:-default} min_pres=${MIN_PRES:-default}"

apptainer exec \
  --bind "/cluster/home/$USER,/cluster/software,/cluster/scratch/$USER" \
  "$SIF" \
  Rscript run_5d_hsm.R "$MODE" "$NW"
