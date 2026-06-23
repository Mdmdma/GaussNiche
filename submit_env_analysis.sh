#!/bin/bash
#SBATCH --job-name=gn-env-analysis
#SBATCH --account=es_schin
#SBATCH --output=output/env-analysis-%j.out
#SBATCH --error=output/env-analysis-%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=8G
#SBATCH --time=08:00:00

# Phase 1 of the 5d-niche branch, end to end:
#   (0) ensure the rocker apptainer SIF exists (re-pull if scratch purged it)
#   (1) build_env_stack.R   -> download + harmonise the candidate stack (cached)
#   (2) analyze_env_stack.R -> per-dataset metrics, correlations, PCA dimensionality
# Heavy/downloading work; runs on a compute node. The build is cache-guarded, so
# re-submitting after an analysis fix re-uses the downloaded stack instantly.

set -euo pipefail
cd "$HOME/GaussNiche"
mkdir -p output

# config_r_studio: spack APPTAINERENV_* (LD paths, krb5 shim) + eth_proxy (HTTP(S)_PROXY)
source "$HOME/.config/euler/jupyterhub/config_r_studio"

# proxy must reach inside the container (geodata downloads + apptainer pull)
export APPTAINERENV_http_proxy="${http_proxy:-${HTTP_PROXY:-http://proxy.ethz.ch:3128}}"
export APPTAINERENV_https_proxy="${https_proxy:-${HTTPS_PROXY:-http://proxy.ethz.ch:3128}}"
export APPTAINERENV_HTTP_PROXY="${APPTAINERENV_http_proxy}"
export APPTAINERENV_HTTPS_PROXY="${APPTAINERENV_https_proxy}"

SIF="/cluster/scratch/$USER/rocker_rstudio_4.5.sif"
if [[ ! -f "$SIF" ]]; then
  echo "[setup] SIF missing — pulling docker://rocker/rstudio:4.5 ..."
  export APPTAINER_TMPDIR="/cluster/scratch/$USER/.apptainer_tmp"
  export APPTAINER_CACHEDIR="/cluster/scratch/$USER/.apptainer_cache"
  mkdir -p "$APPTAINER_TMPDIR" "$APPTAINER_CACHEDIR"
  apptainer pull "$SIF" docker://rocker/rstudio:4.5
  echo "[setup] pulled $SIF"
fi

BIND="--bind /cluster/home/$USER,/cluster/software,/cluster/scratch/$USER"

echo "[run] build_env_stack.R"
apptainer exec $BIND "$SIF" Rscript build_env_stack.R

echo "[run] analyze_env_stack.R"
apptainer exec $BIND "$SIF" Rscript analyze_env_stack.R

echo "[done] diagnostics in /cluster/scratch/$USER/GaussNiche/env5d/diagnostics"
