#!/bin/bash
#SBATCH --job-name=gn-5d-experiment
#SBATCH --account=es_schin
#SBATCH --output=output/5dexp-%j.out
#SBATCH --error=output/5dexp-%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=6G
#SBATCH --time=12:00:00

# 4 virtual species in 5-D (RND / buffer-out / uniform+) over R Bernoulli
# realisations — the 5-D arm of the paper §2.3 experiment. USE (paSampling) is
# 2-D-only and is absent here; buffer-out is the non-random geographic baseline.
#   smoke:  sbatch --export=ALL,MODE=smoke --cpus-per-task=2 --time=00:30:00 sbatch/submit_5d_experiment.sh
#   full:   sbatch sbatch/submit_5d_experiment.sh        (MODE defaults to "full", R=50)
#   tune:   sbatch --export=ALL,SPECIES_CUTOFF=0.7,BUFFER_KM=50,N_REALIZATIONS=50 sbatch/submit_5d_experiment.sh

set -euo pipefail
cd "$HOME/GaussNiche"
mkdir -p output
source "$HOME/.config/euler/jupyterhub/config_r_studio"

# Pin BLAS/OMP threads to 1 so each furrr worker stays single-threaded.
export APPTAINERENV_OMP_NUM_THREADS=1
export APPTAINERENV_OPENBLAS_NUM_THREADS=1
export APPTAINERENV_MKL_NUM_THREADS=1
export APPTAINERENV_OMP_THREAD_LIMIT=1

# Forward optional overrides into the container (only when set).
if [[ -n "${N_REALIZATIONS:-}" ]]; then export APPTAINERENV_N_REALIZATIONS="$N_REALIZATIONS"; fi
if [[ -n "${MAX_PRES:-}"       ]]; then export APPTAINERENV_MAX_PRES="$MAX_PRES"; fi
if [[ -n "${SPECIES_CUTOFF:-}" ]]; then export APPTAINERENV_SPECIES_CUTOFF="$SPECIES_CUTOFF"; fi
if [[ -n "${BUFFER_KM:-}"      ]]; then export APPTAINERENV_BUFFER_KM="$BUFFER_KM"; fi

SIF="/cluster/scratch/$USER/rocker_rstudio_4.5.sif"
[[ -f "$SIF" ]] || { echo "Missing $SIF"; exit 2; }

MODE="${MODE:-full}"
NW="${SLURM_CPUS_PER_TASK:-4}"
echo "[launch] mode=$MODE workers=$NW n_realizations=${N_REALIZATIONS:-default} cutoff=${SPECIES_CUTOFF:-default} buffer_km=${BUFFER_KM:-default}"

apptainer exec \
  --bind "/cluster/home/$USER,/cluster/software,/cluster/scratch/$USER" \
  "$SIF" \
  Rscript run_5d_experiment.R "$MODE" "$NW"
