#!/bin/bash
#SBATCH --job-name=gn-build-env
#SBATCH --account=es_schin
#SBATCH --output=output/build-env-%j.out
#SBATCH --error=output/build-env-%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=8G
#SBATCH --time=06:00:00

# Phase 1 of the 5d-niche branch: download + harmonise the higher-dimensional
# environmental stack (soil + terrain + vegetation via geodata) onto the
# WorldClim 10 arc-min Europe grid, run PCA, and report how many PCs reach 80%.
# One-time build; outputs cached to /cluster/scratch/$USER/GaussNiche/env5d.
# geodata caches downloads, so a timeout mid-run can simply be resubmitted.

set -euo pipefail
cd "$HOME/GaussNiche"
mkdir -p output

# config_r_studio exports the spack APPTAINERENV_* (LD paths, krb5 shim) AND
# loads eth_proxy (sets http_proxy/https_proxy) — needed because this job, unlike
# the others, DOWNLOADS data from the internet.
source "$HOME/.config/euler/jupyterhub/config_r_studio"

# Ensure the proxy reaches INSIDE the container (geodata downloads over HTTPS).
export APPTAINERENV_http_proxy="${http_proxy:-http://proxy.ethz.ch:3128}"
export APPTAINERENV_https_proxy="${https_proxy:-http://proxy.ethz.ch:3128}"
export APPTAINERENV_HTTP_PROXY="${APPTAINERENV_http_proxy}"
export APPTAINERENV_HTTPS_PROXY="${APPTAINERENV_https_proxy}"

SIF="/cluster/scratch/$USER/rocker_rstudio_4.5.sif"
[[ -f "$SIF" ]] || { echo "Missing $SIF — run apptainer pull first."; exit 2; }

# Optional overrides forwarded into the container:
#   ENV_OUT_DIR  build to a sandbox dir (e.g. to bit-identity-check a change)
#   FORCE_REBUILD=1  ignore the cached stack + per-block caches and rebuild all
for v in ENV_OUT_DIR FORCE_REBUILD; do
  if [[ -n "${!v:-}" ]]; then export "APPTAINERENV_$v=${!v}"; fi
done

apptainer exec \
  --bind "/cluster/home/$USER,/cluster/software,/cluster/scratch/$USER" \
  "$SIF" \
  Rscript build_env_stack.R
