# make_figures.R --- the separate plot step for the parallel-jobs pipeline
# =============================================================================
# The 4 compute jobs (run_2d_experiment, run_5d_experiment, run_5d_hsm,
# tune_5d_hsm) are compute + save ONLY and carry NO cross-job dependency, so they
# run fully in parallel. This script (invoked by submit_make_figures.sh, which
# run_all_figures.sh submits after the 4 jobs finish) renders the figures that
# used to be inline, from the saved data:
#   * 2-D §2.3 boxplots             <- results2d/metrics_2d_<mode>.csv
#   * downstream-HSM violins         <- results5d_hsm/hsm_metrics_5d_<mode>.csv
#   * HSM Dunn significance summary  <- results5d_hsm/hsm_dunn_5d_<mode>.csv
#   * uniform+ ablation heatmaps     <- results5d_tune/summary_tune_<mode>.rds ($tune)
#       value + Δ-vs-RND. The Δ view merges run_5d_hsm's RND/buffer medians HERE
#       -- that cross-job merge is why tune could drop its HSM dependency.
# The 5-D boxplots + PC-matrix + geo-grid come from `Rscript run_5d_experiment.R
# figures`, and the full HSM report + aggregates from hsm_report.R /
# hsm_aggregate_report.R; submit_make_figures.sh runs those alongside this.
#
# Usage:  Rscript make_figures.R [mode]     (mode default "full")
# Pure post-processing (ggplot only) -> minutes; reads scratch, writes scratch.
# =============================================================================

suppressPackageStartupMessages({
  .libPaths(c("~/R/rocker-rstudio/4.5", .libPaths()))
  library(ggplot2)
})
source("experiment_plots.R"); source("hsm_plots.R")

mode <- { a <- commandArgs(trailingOnly = TRUE); if (length(a) >= 1) a[1] else "full" }
user <- Sys.getenv("USER"); scratch <- file.path("/cluster/scratch", user)
gn   <- function(...) file.path(scratch, "GaussNiche", ...)
dir_2d   <- Sys.getenv("GN_OUT_2D",   gn("results2d"))
dir_hsm  <- Sys.getenv("GN_OUT_HSM",  gn("results5d_hsm"))
dir_tune <- Sys.getenv("GN_OUT_TUNE", gn("results5d_tune"))

# --- 2-D §2.3 boxplots (mirrors the old run_2d_experiment figure block) -------
fig_2d_boxplots <- function() {
  csv <- file.path(dir_2d, paste0("metrics_2d_", mode, ".csv"))
  if (!file.exists(csv)) { message("[make_figures] skip 2d: ", csv, " absent"); return(invisible()) }
  metrics <- read.csv(csv, stringsAsFactors = FALSE)
  n_real  <- if ("realization" %in% names(metrics)) suppressWarnings(max(metrics$realization, na.rm = TRUE)) else NA
  jit     <- is.finite(n_real) && n_real <= 30
  np      <- experiment_panel_count(metrics)
  p <- plot_experiment_boxplots(metrics, jitter = jit,
         title = "Sampler comparison across virtual species")
  save_experiment_figure(p, file.path(dir_2d, paste0("boxplots_metrics_2d_", mode, ".pdf")),
                         width = 7.0, panel_count = np)
  for (mt in c("overlap", "coverage", "trueabs")) {
    marg <- if (mt == "trueabs") "prop_true_abs" else mt
    ph <- if (mt == "coverage") 1.7 * length(grep("^rel_cov_", names(metrics))) else 3.2
    pm <- plot_experiment_metric(metrics, metric = marg, jitter = jit)
    save_experiment_figure(pm, file.path(dir_2d, paste0("box_", mt, "_2d_", mode, ".pdf")),
                           width = 7.0, height = ph)
  }
  message("[make_figures] 2-D boxplots -> ", dir_2d)
}

# --- downstream-HSM violins + Dunn summary (mirrors run_5d_hsm figure block) ---
fig_hsm_violins <- function() {
  csv <- file.path(dir_hsm, paste0("hsm_metrics_5d_", mode, ".csv"))
  if (!file.exists(csv)) { message("[make_figures] skip hsm violins: ", csv, " absent"); return(invisible()) }
  hsm    <- read.csv(csv, stringsAsFactors = FALSE)
  sp_lab <- { u <- unique(hsm[c("species", "species_label")]); setNames(u$species_label, u$species) }
  for (ps in unique(hsm$predictor_set)) {
    comb <- tryCatch(plot_hsm_violin(hsm, predictor_set = ps,
              title = sprintf("Downstream HSM accuracy across species (%s predictors)", ps)),
              error = function(e) { message("violin ", ps, " skipped: ", conditionMessage(e)); NULL })
    if (!is.null(comb)) {
      w <- 3 + 2.4 * length(unique(hsm$species))
      save_hsm_figure(comb, file.path(dir_hsm, paste0("hsm_violin_", ps, "_", mode, ".pdf")),
                      width = min(20, w), height = 8)
    }
    for (nm in unique(hsm$species)) {
      p <- tryCatch(plot_hsm_violin(hsm, predictor_set = ps, species = nm,
             title = sprintf("%s — HSM accuracy (%s predictors)", sp_lab[[nm]], ps)),
             error = function(e) NULL)
      if (!is.null(p))
        save_hsm_figure(p, file.path(dir_hsm, paste0("hsm_violin_", ps, "_", nm, "_", mode, ".pdf")),
                        width = 9, height = 7)
    }
  }
  dunn_csv <- file.path(dir_hsm, paste0("hsm_dunn_5d_", mode, ".csv"))
  dunn <- if (file.exists(dunn_csv)) read.csv(dunn_csv, stringsAsFactors = FALSE)
          else tryCatch(hsm_dunn_test(hsm), error = function(e) NULL)
  if (!is.null(dunn) && nrow(dunn)) {
    ds <- tryCatch(plot_hsm_dunn_summary(dunn), error = function(e) NULL)
    if (!is.null(ds))
      save_hsm_figure(ds, file.path(dir_hsm, paste0("hsm_dunn_summary_5d_", mode, ".pdf")),
                      width = 11, height = 3.5)
  }
  message("[make_figures] HSM violins + Dunn summary -> ", dir_hsm)
}

# --- uniform+ ablation heatmaps: value + Δ-vs-RND ----------------------------
# The Δ view needs run_5d_hsm's RND/buffer medians (pc5); reading them HERE is
# the cross-job merge that used to live in tune_5d_hsm.R (mirrors that block).
fig_ablation_heatmaps <- function() {
  rds <- file.path(dir_tune, paste0("summary_tune_", mode, ".rds"))
  if (!file.exists(rds)) { message("[make_figures] skip ablation: ", rds, " absent"); return(invisible()) }
  s <- readRDS(rds); tune <- s$tune
  if (is.null(tune)) { message("[make_figures] tune rds has no $tune"); return(invisible()) }
  env_lv <- s$env_cutoffs; sp_lv <- s$species_cutoffs   # full intended axes (degenerate cols too)
  metr <- c("cor_truth", "rmse_truth", "auc", "tss", "boyce")
  ref  <- list()
  base_csv <- file.path(dir_hsm, paste0("hsm_metrics_5d_", mode, ".csv"))
  if (file.exists(base_csv)) {
    b <- read.csv(base_csv, stringsAsFactors = FALSE)
    b <- b[b$status == "ok" & b$predictor_set == "pc5" & b$sampler %in% c("random", "buffer"), ]
    for (m in metr)
      ref[[m]] <- c(random = stats::median(b[b$sampler == "random", m], na.rm = TRUE),
                    buffer = stats::median(b[b$sampler == "buffer", m], na.rm = TRUE))
  } else message("[make_figures] baseline ", base_csv, " absent -> ablation value heatmaps only")
  nice <- c(cor_truth = "Truth recovery (r)", rmse_truth = "Truth error (RMSE)",
            auc = "Discrimination (AUC)", tss = "Discrimination (TSS)", boyce = "Calibration (CBI)")
  for (m in metr) {
    hb <- (m != "rmse_truth")
    rr <- if (length(ref)) unname(ref[[m]]["random"]) else NULL
    for (kind in c("value", "delta")) {
      if (kind == "delta" && is.null(rr)) next
      p <- tryCatch(plot_tune_heatmap(tune, m,
             ref = if (kind == "delta") rr else NULL,
             baseline = if (length(ref)) ref[[m]] else NULL, higher_better = hb,
             env_levels = env_lv, sp_levels = sp_lv,
             title = sprintf("%s%s", if (m %in% names(nice)) nice[[m]] else m,
                             if (kind == "delta") "  —  Δ vs RND" else "")),
             error = function(e) { message("heatmap ", m, " ", kind, " skipped: ", conditionMessage(e)); NULL })
      if (!is.null(p))
        save_hsm_figure(p, file.path(dir_tune, sprintf("tune_heat_%s_%s_%s.pdf", m, kind, mode)),
                        width = 9.5, height = 6.4)
    }
  }
  message("[make_figures] ablation heatmaps -> ", dir_tune)
}

# --- per-species ablation heatmaps -------------------------------------------
# Same cor_truth + auc VALUE heatmaps as the aggregate pair (which pools all four
# species), but split by species -> for the appendix. tune_heat_<m>_value_<spN>_<mode>.pdf
fig_ablation_per_species <- function() {
  rds <- file.path(dir_tune, paste0("summary_tune_", mode, ".rds"))
  if (!file.exists(rds)) { message("[make_figures] skip per-species ablation: ", rds, " absent"); return(invisible()) }
  s <- readRDS(rds); tune <- s$tune
  if (is.null(tune) || !"species" %in% names(tune)) { message("[make_figures] tune rds lacks $tune/species"); return(invisible()) }
  env_lv <- s$env_cutoffs; sp_lv <- s$species_cutoffs   # full intended axes (degenerate cols too)
  sp_key <- c(sp1_generalist_common = "sp1", sp2_specialist_common = "sp2",
              sp3_generalist_rare = "sp3", sp4_specialist_rare = "sp4")
  for (sp in unique(tune$species)) {
    d   <- tune[tune$species == sp, , drop = FALSE]
    key <- if (sp %in% names(sp_key)) sp_key[[sp]] else sp
    for (m in c("cor_truth", "auc")) {           # match the Results aggregate pair
      p <- tryCatch(plot_tune_heatmap(d, m, ref = NULL, higher_better = (m != "rmse_truth"),
             env_levels = env_lv, sp_levels = sp_lv,
             title = sprintf("uniform+ %s — %s", m, sp)),
             error = function(e) { message("per-species heatmap ", m, " ", sp, " skipped: ", conditionMessage(e)); NULL })
      if (!is.null(p))
        save_hsm_figure(p, file.path(dir_tune, sprintf("tune_heat_%s_value_%s_%s.pdf", m, key, mode)),
                        width = 9.5, height = 6.4)
    }
  }
  message("[make_figures] per-species ablation heatmaps -> ", dir_tune)
}

# --- full-page per-species ablation composite (paper appendix) ---------------
# ALL 8 per-species heatmaps on ONE portrait page: rows = the 4 virtual species,
# columns = {truth-recovery (cor_truth), discrimination (AUC)}. Each panel is
# viridis-scaled on its OWN range (species differ by an order of magnitude in
# cor_truth), carries a bold best-cell ring, and states that species/metric's fixed
# RND (uniform-random) + buffer-out baselines. Human-readable axis titles (no code
# dot-names). One cairo-PDF: tune_heat_per_species_page_<mode>.pdf
fig_ablation_per_species_page <- function() {
  rds <- file.path(dir_tune, paste0("summary_tune_", mode, ".rds"))
  if (!file.exists(rds)) { message("[make_figures] skip per-species page: ", rds, " absent"); return(invisible()) }
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    message("[make_figures] patchwork absent -> skip per-species page"); return(invisible()) }
  s <- readRDS(rds); tune <- s$tune
  if (is.null(tune) || !"species" %in% names(tune)) { message("[make_figures] no $tune/species -> skip page"); return(invisible()) }
  env_lv <- s$env_cutoffs; sp_lv <- s$species_cutoffs

  # per-species / per-metric fixed baselines (RND + buffer-out, pc5 predictors)
  base_csv <- file.path(dir_hsm, paste0("hsm_metrics_5d_", mode, ".csv"))
  bl <- if (file.exists(base_csv)) {
    b <- read.csv(base_csv, stringsAsFactors = FALSE)
    b[b$status == "ok" & b$predictor_set == "pc5" & b$sampler %in% c("random", "buffer"), ]
  } else NULL
  base_txt <- function(sp, m) {
    if (is.null(bl)) return(NULL)
    r  <- stats::median(bl[bl$species == sp & bl$sampler == "random", m], na.rm = TRUE)
    bu <- stats::median(bl[bl$species == sp & bl$sampler == "buffer", m], na.rm = TRUE)
    sprintf("RND %.2f  ·  buffer-out %.2f", r, bu)
  }

  sp_lab   <- c(sp1_generalist_common = "Generalist · common", sp2_specialist_common = "Specialist · common",
                sp3_generalist_rare = "Generalist · rare",     sp4_specialist_rare = "Specialist · rare")
  metr_lab <- c(cor_truth = "Truth recovery (r)", auc = "Discrimination (AUC)")
  sp_order <- names(sp_lab)[names(sp_lab) %in% unique(tune$species)]

  panel <- function(sp, m, show_x, show_y) {
    d  <- tune[tune$species == sp, , drop = FALSE]
    hb <- (m != "rmse_truth")
    agg <- aggregate(d[[m]], by = list(env = d$env_cutoff, sp = d$species_cutoff),
                     FUN = function(v) stats::median(v, na.rm = TRUE))
    names(agg)[3] <- "val"
    full <- merge(expand.grid(env = env_lv, sp = sp_lv), agg, by = c("env", "sp"), all.x = TRUE)
    full$env <- factor(full$env, levels = env_lv)
    full$sp  <- factor(full$sp,  levels = rev(sp_lv))
    full$deg <- !is.finite(full$val)
    full$lab <- ifelse(full$deg, "", sprintf("%.2f", full$val))
    win_i <- if (hb) which.max(full$val) else which.min(full$val); if (!length(win_i)) win_i <- 1L
    win <- full[win_i, , drop = FALSE]
    deg <- full[full$deg, , drop = FALSE]
    if (nrow(deg)) { deg$x0 <- as.integer(deg$env) - .45; deg$x1 <- as.integer(deg$env) + .45
                     deg$y0 <- as.integer(deg$sp)  - .45; deg$y1 <- as.integer(deg$sp)  + .45 }
    bt  <- base_txt(sp, m)
    sub <- if (!is.null(bt)) paste0(metr_lab[[m]], "   ·   baselines: ", bt) else metr_lab[[m]]
    g <- ggplot(full, aes(.data$env, .data$sp, fill = .data$val)) +
      geom_tile(colour = "white", linewidth = 0.3) +
      geom_tile(data = win, aes(.data$env, .data$sp), fill = NA,
                colour = "black", linewidth = 1.1, inherit.aes = FALSE) +
      geom_text(aes(label = .data$lab), size = 1.75) +
      scale_fill_viridis_c(option = "viridis", direction = if (hb) 1 else -1, guide = "none") +
      labs(title = sp_lab[[sp]], subtitle = sub,
           x = if (show_x) "Environmental cutoff (percentile)" else NULL,
           y = if (show_y) "Species cutoff (quantile)" else NULL) +
      theme_minimal(base_size = 8) +
      theme(panel.grid = element_blank(),
            plot.title = element_text(face = "bold", size = 9),
            plot.subtitle = element_text(size = 6.6, colour = "grey25"),
            axis.title = element_text(size = 7.2),
            axis.text = element_text(size = 5.2),
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
    if (nrow(deg)) g <- g +
      geom_segment(data = deg, aes(x = .data$x0, xend = .data$x1, y = .data$y0, yend = .data$y1),
                   colour = "red", linewidth = 0.5, inherit.aes = FALSE) +
      geom_segment(data = deg, aes(x = .data$x0, xend = .data$x1, y = .data$y1, yend = .data$y0),
                   colour = "red", linewidth = 0.5, inherit.aes = FALSE)
    g
  }

  panels <- list()
  for (i in seq_along(sp_order)) {
    last <- (i == length(sp_order))
    panels[[length(panels) + 1L]] <- panel(sp_order[i], "cor_truth", show_x = last, show_y = TRUE)
    panels[[length(panels) + 1L]] <- panel(sp_order[i], "auc",       show_x = last, show_y = FALSE)
  }
  # No master title/subtitle: that explainer lives in the LaTeX \caption. Sized to
  # leave caption room when included at \textwidth on the appendix page.
  pg <- patchwork::wrap_plots(panels, ncol = 2)
  out <- file.path(dir_tune, paste0("tune_heat_per_species_page_", mode, ".pdf"))
  ggsave(out, pg, device = grDevices::cairo_pdf, width = 8.3, height = 9.3)
  message("[make_figures] per-species composite page -> ", out)
  invisible(pg)
}

# Only render when run directly (Rscript make_figures.R); sourcing the file for a
# single-figure preview just loads the fig_* definitions without re-running all.
if (sys.nframe() == 0L) {
  fig_2d_boxplots()
  fig_hsm_violins()
  fig_ablation_heatmaps()
  fig_ablation_per_species()
  fig_ablation_per_species_page()
}
message("[make_figures] done (", mode, ")")
