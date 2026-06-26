# hsm_aggregate_report.R --- cross-species aggregate figures + tables
# =============================================================================
# Reads a tidy HSM metrics CSV (from run_5d_hsm.R) and produces the cross-species
# aggregate views (algorithm dimension collapsed), so the four virtual species
# can be compared in one place:
#   1. agg_overlay_<ps>.{pdf,png}     all 4 species in ONE plot, each a distinct
#                                     shade + central symbol, samplers on x,
#                                     faceted by metric (median +/- IQR, HSMs averaged).
#   2. agg_species_row_<ps>.{pdf,png} 4 species in a row (columns), HSMs averaged:
#                                     sampler boxplots, faceted metric x species.
#   3. agg_table_species_sampler.csv  per predictor_set x species x sampler medians.
#      agg_table_overall.csv          per predictor_set x sampler (pooled species).
#      agg_table_delta_vs_rnd.csv     sampler - RND deltas per metric/species.
#
# Pure post-processing (no model fitting) -> reads the CSV, renders, exits.
# Usage:  Rscript hsm_aggregate_report.R [hsm_csv] [out_dir]
# =============================================================================

suppressPackageStartupMessages({
  .libPaths(c("~/R/rocker-rstudio/4.5", .libPaths()))
  library(ggplot2)
})
source("experiment_plots.R")   # sampler colours/labels reused by hsm_plots.R
source("hsm_plots.R")

args   <- commandArgs(trailingOnly = TRUE)
user   <- Sys.getenv("USER"); scratch <- file.path("/cluster/scratch", user)
in_csv <- if (length(args) >= 1) args[1] else
  file.path(scratch, "GaussNiche", "results5d_hsm", "hsm_metrics_5d_full.csv")
out_dir <- if (length(args) >= 2) args[2] else dirname(in_csv)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
cat("input :", in_csv, "\noutput:", out_dir, "\n")

hsm <- read.csv(in_csv, stringsAsFactors = FALSE)
pred_sets <- unique(hsm$predictor_set)
fig_metrics <- c("cor_truth", "rmse_truth", "auc", "tss", "boyce")  # 2 truth + 3 discrimination

# --- 1 + 2. figures per predictor_set ----------------------------------------
for (ps in pred_sets) {
  ov <- plot_hsm_species_overlay(hsm, predictor_set = ps, metrics = fig_metrics)
  save_hsm_figure(ov, file.path(out_dir, paste0("agg_overlay_", ps, ".pdf")),
                  width = 12, height = 4.2)
  rw <- plot_hsm_species_row(hsm, predictor_set = ps, metrics = fig_metrics)
  save_hsm_figure(rw, file.path(out_dir, paste0("agg_species_row_", ps, ".pdf")),
                  width = 10, height = 8)
  cat("figures written for predictor_set =", ps, "\n")
}

# --- 3. tables ----------------------------------------------------------------
tbl_metrics <- c("auc", "tss", "boyce", "sensitivity", "specificity",
                 "cor", "rmse", "cor_truth", "rmse_truth", "r2")

t_species <- hsm_aggregate_table(hsm, metrics = tbl_metrics,
                                 by = c("predictor_set", "species", "sampler"))
write.csv(t_species, file.path(out_dir, "agg_table_species_sampler.csv"), row.names = FALSE)

t_overall <- hsm_aggregate_table(hsm, metrics = tbl_metrics,
                                 by = c("predictor_set", "sampler"))
write.csv(t_overall, file.path(out_dir, "agg_table_overall.csv"), row.names = FALSE)

t_delta <- hsm_delta_table(hsm, metrics = fig_metrics, ref = "random",
                           by = c("predictor_set", "species"))
write.csv(t_delta, file.path(out_dir, "agg_table_delta_vs_rnd.csv"), row.names = FALSE)

cat("\n================= OVERALL (pooled over species + algos + realisations) =================\n")
print(t_overall, row.names = FALSE)
cat("\n================= DELTA vs RND (sampler - RND; >0 means better than random, except RMSE) =================\n")
print(t_delta, row.names = FALSE)
cat("\nWrote: agg_overlay_*, agg_species_row_*, agg_table_{species_sampler,overall,delta_vs_rnd}.csv to\n  ", out_dir, "\n")
