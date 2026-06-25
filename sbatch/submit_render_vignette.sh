#!/bin/bash
#SBATCH --job-name=mcmc-vignette-figs
#SBATCH --account=es_schin
#SBATCH --output=output/render-vignette-%j.out
#SBATCH --error=output/render-vignette-%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8G
#SBATCH --time=04:00:00

# Re-knit the USE.MCMC MCMC vignette on a COMPUTE NODE to regenerate the paper's
# figures as vector PDFs (the vignette defaults to dev="cairo_pdf", so every
# figure is vector), then copy them into the manuscript's graphics/ via
# `make figures`.
#
# Why this exists: the figures live in the vignette (see the publication-figures
# skill in ../USE.MCMC/.claude/skills/). The vignette chunks were restyled
# (colourblind-safe palette, Lato, no baked-in title, dev="cairo_pdf"); re-knitting
# runs real MCMC (optimRes, paSamplingMcmc, the 4-chain convergence ensemble, the
# gif chunks) -- minutes of compute, NOT a login-node task.
#
# Knitting is done figures-only (run_pandoc = FALSE) so there is no pandoc
# dependency; self_contained = FALSE keeps the figure-html/ directory on disk.

set -euo pipefail
cd "$HOME/GaussNiche"            # logs land in GaussNiche/output, like the other jobs
mkdir -p output
source "$HOME/.config/euler/jupyterhub/config_r_studio"
SIF="/cluster/scratch/$USER/rocker_rstudio_4.5.sif"
[[ -f "$SIF" ]] || { echo "Missing $SIF"; exit 2; }

USE_MCMC="$HOME/USE.MCMC"
PAPER="$HOME/markov-chain-sampler-paper"
VIGNETTE="vignettes/insights-on-MCMC-pseudo-absence-sampling-vignette.Rmd"
FIGDIR="$USE_MCMC/vignettes/insights-on-MCMC-pseudo-absence-sampling-vignette_files/figure-html"

[[ -d "$USE_MCMC" ]] || { echo "Missing $USE_MCMC"; exit 2; }
[[ -d "$PAPER"    ]] || { echo "Missing $PAPER (the manuscript repo must be a sibling)"; exit 2; }

# proxy passthrough (in case any chunk fetches data)
export APPTAINERENV_http_proxy="${http_proxy:-${HTTP_PROXY:-http://proxy.ethz.ch:3128}}"
export APPTAINERENV_https_proxy="${https_proxy:-${HTTPS_PROXY:-http://proxy.ethz.ch:3128}}"
export APPTAINERENV_HTTP_PROXY="${APPTAINERENV_http_proxy}"
export APPTAINERENV_HTTPS_PROXY="${APPTAINERENV_https_proxy}"

BIND="--bind /cluster/home/$USER,/cluster/software,/cluster/scratch/$USER"

echo "[render] knitting $VIGNETTE (figures only) ..."
apptainer exec $BIND "$SIF" Rscript -e '
  setwd("'"$USE_MCMC"'")
  rmarkdown::render(
    "'"$VIGNETTE"'",
    output_format = rmarkdown::html_vignette(self_contained = FALSE),
    run_pandoc = FALSE, quiet = FALSE)
'

# Verify the paper's referenced vector PDFs were produced before touching it.
# (trace_plots.png is static / not knit-produced -- see the publication-figures skill.)
missing=0
for f in generate-combined-plot-1.pdf gelman-plot-1.pdf autocorreltation-plot-1.pdf \
         plot-in-geo-1.pdf plot-pipeline-scheme-1.pdf plot-methods-overview-1.pdf \
         trace_plots-1.pdf; do
  if [[ -f "$FIGDIR/$f" ]]; then echo "  ok   $f"; else echo "  MISS $f"; missing=1; fi
done
if [[ "$missing" -ne 0 ]]; then
  echo "[render] one or more vector PDFs missing -- check the knit log above"
  echo "         (showtext/cairo unavailable, or a chunk errored before the figure)."
  exit 3
fi

echo "[render] copying figures into $PAPER/graphics via 'make figures' ..."
make -C "$PAPER" figures

echo "[done] vector-PDF figures regenerated and copied. trace_plots.png left static."
echo "       To compile: 'latexmk -xelatex main.tex' in $PAPER (then land changes on main)."
