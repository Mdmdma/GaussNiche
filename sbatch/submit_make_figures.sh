#!/bin/bash
#SBATCH --job-name=gn-figures
#SBATCH --account=es_schin
#SBATCH --output=output/figures-%j.out
#SBATCH --error=output/figures-%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=6G
#SBATCH --time=01:00:00
#
# PLOT STEP of the parallel-jobs pipeline: renders ALL Part-B figures from the data
# the 4 compute jobs saved (no recompute). run_all_figures.sh submits this with
# --dependency=afterok on the 4 compute jobs, so it runs once they all finish. The
# 4 compute jobs themselves carry NO dependency and run fully in parallel.
#   MODE (default full) forwarded to each step.
set -euo pipefail
cd "$HOME/GaussNiche"; mkdir -p output
source "$HOME/.config/euler/jupyterhub/config_r_studio"
export APPTAINERENV_OMP_NUM_THREADS=1 APPTAINERENV_OPENBLAS_NUM_THREADS=1 APPTAINERENV_MKL_NUM_THREADS=1
SIF="/cluster/scratch/$USER/rocker_rstudio_4.5.sif"
[[ -f "$SIF" ]] || { echo "Missing $SIF"; exit 2; }

MODE="${MODE:-full}"
export APPTAINERENV_MODE="$MODE"
HSM_CSV="/cluster/scratch/$USER/GaussNiche/results5d_hsm/hsm_metrics_5d_${MODE}.csv"
BIND="/cluster/home/$USER,/cluster/software,/cluster/scratch/$USER"
# Each figure family is independent -> a failure in one must not abort the rest.
run() { echo "[figures] $*"; apptainer exec --bind "$BIND" "$SIF" "$@"; }

# 5-D boxplots + PC-matrix + geo-grid (re-render from the saved exp_*_full.rds)
run Rscript run_5d_experiment.R figures         || echo "[figures] 5d figures step FAILED (continuing)"
# downstream-HSM full report + cross-species aggregates (read the HSM CSV / rds)
run Rscript hsm_report.R "$HSM_CSV"             || echo "[figures] hsm_report FAILED (continuing)"
run Rscript hsm_aggregate_report.R "$HSM_CSV"   || echo "[figures] hsm_aggregate_report FAILED (continuing)"
# 2-D boxplots + HSM violins + Dunn summary + ablation heatmaps (value + Δ-vs-RND)
run Rscript make_figures.R "$MODE"              || echo "[figures] make_figures FAILED (continuing)"
echo "[figures] done (mode=$MODE)."
