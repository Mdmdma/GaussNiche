#!/bin/bash
#SBATCH --job-name=gn-2d-experiment
#SBATCH --account=es_schin
#SBATCH --output=output/2d-%j.out
#SBATCH --error=output/2d-%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=4G
#SBATCH --time=04:00:00

# Run the 4 virtual species in 2-D (random + uniform + mcmc samplers) — the
# paper §2.3 experiment — over R Bernoulli realisations, and emit the boxplots.
#   smoke validation:  sbatch --export=ALL,MODE=smoke --cpus-per-task=2 --time=00:20:00 sbatch/submit_2d_experiment.sh
#   full run (R=50):   sbatch sbatch/submit_2d_experiment.sh        (MODE defaults to "full")
#   override depth:    sbatch --export=ALL,N_REALIZATIONS=50 sbatch/submit_2d_experiment.sh
#   tune MCMC cutoff:  sbatch --export=ALL,SPECIES_CUTOFF=0.75 sbatch/submit_2d_experiment.sh

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

# Forward optional run overrides into the container (only when set, so the R
# driver's Sys.getenv() defaults still apply when an override is absent). Use
# if-blocks, not `[[ ]] && export`, which would abort under `set -e` when false.
if [[ -n "${N_REALIZATIONS:-}" ]]; then export APPTAINERENV_N_REALIZATIONS="$N_REALIZATIONS"; fi
if [[ -n "${MAX_PRES:-}"       ]]; then export APPTAINERENV_MAX_PRES="$MAX_PRES"; fi
if [[ -n "${SPECIES_CUTOFF:-}" ]]; then export APPTAINERENV_SPECIES_CUTOFF="$SPECIES_CUTOFF"; fi

SIF="/cluster/scratch/$USER/rocker_rstudio_4.5.sif"
[[ -f "$SIF" ]] || { echo "Missing $SIF"; exit 2; }

MODE="${MODE:-full}"
NW="${SLURM_CPUS_PER_TASK:-4}"
echo "[launch] mode=$MODE workers=$NW n_realizations=${N_REALIZATIONS:-default} species_cutoff=${SPECIES_CUTOFF:-default}"

apptainer exec \
  --bind "/cluster/home/$USER,/cluster/software,/cluster/scratch/$USER" \
  "$SIF" \
  Rscript run_2d_experiment.R "$MODE" "$NW"
