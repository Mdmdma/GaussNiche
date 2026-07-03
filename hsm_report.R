# hsm_report.R --- scrollable full per-species HSM report
# =============================================================================
# Builds one multi-page PDF from a tidy HSM metrics CSV (one with ALL THREE
# samplers -- i.e. the baseline results5d_hsm/hsm_metrics_5d_full.csv, NOT a
# uniform+-only tuning run):
#   p1        overview tables: median of every metric per sampler (per predictor set)
#   p2-3      WIN MATRIX per predictor set -- where uniform+ beats both naive
#             samplers (RND + buffer-out), green/amber/red + signed margin
#   p4..      per (predictor_set x species) all-metrics violin ("flute") grids
# Plus hsm_win_matrix.csv (the full table behind the matrices).
#
# Usage:  Rscript hsm_report.R [hsm_csv] [out_pdf]
# =============================================================================

suppressPackageStartupMessages({
  .libPaths(c("~/R/rocker-rstudio/4.5", .libPaths()))
  library(ggplot2)
})
source("experiment_plots.R"); source("hsm_plots.R")
for (pk in c("gridExtra", "grid")) if (!requireNamespace(pk, quietly = TRUE)) stop("need ", pk)

args   <- commandArgs(trailingOnly = TRUE)
user   <- Sys.getenv("USER"); scratch <- file.path("/cluster/scratch", user)
in_csv <- if (length(args) >= 1) args[1] else
  file.path(scratch, "GaussNiche", "results5d_hsm", "hsm_metrics_5d_full.csv")
out_pdf <- if (length(args) >= 2) args[2] else
  file.path(dirname(in_csv), "hsm_full_report.pdf")
cat("input :", in_csv, "\noutput:", out_pdf, "\n")

hsm <- read_hsm_metrics(in_csv)   # CSV, or committed summary_5d_hsm_<mode>.rds$hsm fallback
ok  <- hsm[hsm$status == "ok", , drop = FALSE]
relab <- c(random = "RND", buffer = "buffer-out", mcmc = "uniform+")
sp_order <- { u <- unique(ok[c("species", "species_label")]); u$species_label[order(u$species)] }
pred_sets <- intersect(c("pc5", "raw12"), unique(ok$predictor_set))

# --- overview median table per predictor set (metrics x samplers) ------------
mk_tab <- function(ps) {
  d <- ok[ok$predictor_set == ps, , drop = FALSE]
  ag <- aggregate(d[HSM_ALL_METRICS], by = list(sampler = d$sampler),
                  FUN = function(v) round(median(v, na.rm = TRUE), 3))
  ag <- ag[match(c("random", "buffer", "mcmc"), ag$sampler), , drop = FALSE]
  m <- t(as.matrix(ag[HSM_ALL_METRICS])); colnames(m) <- unname(relab[ag$sampler])
  data.frame(metric = unname(HSM_METRIC_LABELS[rownames(m)]), m, check.names = FALSE,
             stringsAsFactors = FALSE)
}
tt <- gridExtra::ttheme_minimal(base_size = 8,
        core = list(fg_params = list(hjust = 0, x = 0.05)),
        colhead = list(fg_params = list(fontface = "bold")))

# --- assemble the multi-page PDF ---------------------------------------------
grDevices::cairo_pdf(out_pdf, width = 12, height = 8.5, onefile = TRUE)

# page 1: overview tables
hdr <- function(t, s = 13) grid::textGrob(t, x = 0.02, hjust = 0,
        gp = grid::gpar(fontface = "bold", fontsize = s))
tabs <- lapply(pred_sets, function(ps) gridExtra::arrangeGrob(
  hdr(sprintf("Median per metric x sampler — %s predictors", ps), 11),
  gridExtra::tableGrob(mk_tab(ps), rows = NULL, theme = tt), ncol = 1,
  heights = grid::unit(c(0.4, 0.22 * (length(HSM_ALL_METRICS) + 1) + 0.3), "in")))
grid::grid.draw(gridExtra::arrangeGrob(
  hdr("GaussNiche 5-D downstream HSM — full report (RND / buffer-out / uniform+)", 15),
  do.call(gridExtra::arrangeGrob, c(tabs, list(ncol = length(tabs)))),
  hdr("Internal metrics = held-out test vs sampled PA;  *_truth = predicted map vs known suitability.", 9),
  ncol = 1, heights = grid::unit(c(0.5, 7.0, 0.4), "in")))

# pages 2-3: win matrices
for (ps in pred_sets)
  print(plot_hsm_win_matrix(hsm, predictor_set = ps))

# pages 4+: per (predictor_set x species) all-metrics violins
for (ps in pred_sets) for (sp in unique(ok$species[order(match(ok$species_label, sp_order))])) {
  p <- tryCatch(plot_hsm_species_metrics(hsm, species = sp, predictor_set = ps),
                error = function(e) NULL)
  if (!is.null(p)) print(p)
}
grDevices::dev.off()

# --- the full win-matrix table (CSV) -----------------------------------------
wm <- hsm_win_matrix(hsm)
write.csv(wm, file.path(dirname(out_pdf), "hsm_win_matrix.csv"), row.names = FALSE)

# --- console summary: where uniform+ beats both naive samplers ---------------
cat("\n=== uniform+ vs naive (RND & buffer-out): cells where uniform+ BEATS BOTH ===\n")
beats <- wm[wm$category == "beats both", c("predictor_set", "species", "metric", "delta_vs_best_naive")]
print(beats[order(beats$predictor_set, beats$species), ], row.names = FALSE)
cat(sprintf("\nbeats both: %d / %d   beats one: %d   beats none: %d\n",
            sum(wm$category == "beats both"), nrow(wm),
            sum(wm$category == "beats one"), sum(wm$category == "beats none")))
cat("\nWrote report:", out_pdf, "\n")
