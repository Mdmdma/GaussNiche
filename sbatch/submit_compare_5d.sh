#!/bin/bash
#SBATCH --job-name=gn-5d-compare
#SBATCH --account=es_schin
#SBATCH --output=output/5d-compare-%j.out
#SBATCH --error=output/5d-compare-%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=8G
#SBATCH --time=00:20:00

# Build the cross-species/sampler comparison PDF + combined metrics CSV from the
# results saved by 5_highdim_species.R.  MODE defaults to "full".

set -euo pipefail
cd "$HOME/GaussNiche"
mkdir -p output
source "$HOME/.config/euler/jupyterhub/config_r_studio"
SIF="/cluster/scratch/$USER/rocker_rstudio_4.5.sif"
[[ -f "$SIF" ]] || { echo "Missing $SIF"; exit 2; }
MODE="${MODE:-full}"

# proxy passthrough so pip (pypdf install) works inside the container
export APPTAINERENV_http_proxy="${http_proxy:-${HTTP_PROXY:-http://proxy.ethz.ch:3128}}"
export APPTAINERENV_https_proxy="${https_proxy:-${HTTPS_PROXY:-http://proxy.ethz.ch:3128}}"
export APPTAINERENV_HTTP_PROXY="${APPTAINERENV_http_proxy}"
export APPTAINERENV_HTTPS_PROXY="${APPTAINERENV_https_proxy}"

BIND="--bind /cluster/home/$USER,/cluster/software,/cluster/scratch/$USER"
apptainer exec $BIND "$SIF" Rscript 6_compare_5d.R "$MODE"

# --- inject clickable PDF bookmarks (no gs/pdftk/qpdf in the container; use
#     python/pypdf, installed to ~/.local on first run). Non-fatal on failure. ---
OUT="/cluster/scratch/$USER/GaussNiche/results5d"
export APPTAINERENV_GN_PDF="$OUT/report_${MODE}.pdf"
export APPTAINERENV_GN_TSV="$OUT/report_${MODE}_bookmarks.tsv"
export APPTAINERENV_GN_HOME="$HOME/GaussNiche"
if [[ -f "$OUT/report_${MODE}.pdf" && -f "$OUT/report_${MODE}_bookmarks.tsv" ]]; then
  apptainer exec $BIND "$SIF" bash -lc '
    export PYTHONPATH="$GN_HOME/pylib${PYTHONPATH:+:$PYTHONPATH}"
    if python3 -c "import pypdf" 2>/dev/null; then
      python3 "$GN_HOME/add_bookmarks.py" "$GN_PDF" "$GN_TSV" "$GN_PDF.bm" \
        && mv -f "$GN_PDF.bm" "$GN_PDF" && echo "[bookmarks] injected via vendored pypdf"
    else
      echo "[bookmarks] vendored pypdf missing in $GN_HOME/pylib (run: python3 vendor_pypdf.py); divider pages still navigate"
    fi
  ' || echo "[bookmarks] step failed (non-fatal)"
fi
