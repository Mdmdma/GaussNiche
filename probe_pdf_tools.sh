#!/bin/bash
#SBATCH --job-name=gn-probe
#SBATCH --account=es_schin
#SBATCH --output=output/probe-%j.out
#SBATCH --error=output/probe-%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
#SBATCH --time=00:05:00

# Probe the container for a PDF-bookmark-capable tool (ghostscript is absent).
set -euo pipefail
cd "$HOME/GaussNiche"; mkdir -p output
source "$HOME/.config/euler/jupyterhub/config_r_studio"
SIF="/cluster/scratch/$USER/rocker_rstudio_4.5.sif"
apptainer exec --bind "/cluster/home/$USER,/cluster/software,/cluster/scratch/$USER" "$SIF" bash -lc '
  echo "== binaries =="
  for b in gs pdftk qpdf mutool cpdf pdftoppm python3 perl; do printf "%-10s " "$b"; command -v "$b" || echo "(none)"; done
  echo "== python pypdf =="
  python3 -c "import pypdf; print(\"pypdf\", pypdf.__version__)" 2>&1 | head -1
  echo "== R pkgs =="
  Rscript -e ".libPaths(c(\"~/R/rocker-rstudio/4.5\", .libPaths())); for (p in c(\"qpdf\",\"staplr\",\"pdftools\")) cat(p, requireNamespace(p, quietly=TRUE), \"\n\")"
'
