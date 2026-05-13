# Multi-species evaluation wrapper for virtualSpecies_fn.R
#
# Runs an arbitrary catalogue of virtual species through all pseudo-absence
# samplers and produces, for every figure that the single-species script
# (2_testing_wrapper_function.R) renders, an aggregate side-by-side view across
# the species in the catalogue. Layout adjusts dynamically to the catalogue
# size (1, 2, 4, … species).

setwd("/home/mat/Desktop/semesterarbeit10/GaussNiche/")

library(terra)
# library(geodata)
library(USE.MCMC)
library(mvtnorm)
library(MASS)
library(ggplot2)
library(viridis)
library(hypervolume)
library(patchwork)
library(sf)
library(tictoc)
library(future)
library(furrr)
source("virtualSpecies_fn.R")

# ──────────────────────────────────────────────────────────────────────────────
# 1. ENVIRONMENTAL DATA AND PCA
# ──────────────────────────────────────────────────────────────────────────────
envData <- rast(USE.MCMC::Worldclim_tmp, type = "xyz")
rpc     <- rastPCA(envData, stand = TRUE)

dt <- na.omit(as.data.frame(rpc$PCs[[c("PC1", "PC2")]], xy = TRUE))

var_exp <- rpc$pca$sdev^2 / sum(rpc$pca$sdev^2)
pc1_lab <- paste0("PC1 (", round(var_exp[1] * 100, 1), "% var)")
pc2_lab <- paste0("PC2 (", round(var_exp[2] * 100, 1), "% var)")

cat("Background cells (n):", nrow(dt), "\n")
cat("PC1 range:", round(range(dt$PC1), 2), "\n")
cat("PC2 range:", round(range(dt$PC2), 2), "\n")

# ──────────────────────────────────────────────────────────────────────────────
# PDF DEVICE
#
# Open a wide PDF up front so multi-panel aggregates have room to breathe.
# Filename includes the launch timestamp so successive runs never overwrite a
# previous report.  Closed at script end (and on error) via on.exit().
# ──────────────────────────────────────────────────────────────────────────────
pdf_path <- sprintf("plots_multi_species_%s.pdf",
                    format(Sys.time(), "%Y%m%d-%H%M%S"))
pdf(pdf_path, width = 14, height = 10)
on.exit(try(dev.off(), silent = TRUE), add = TRUE)
cat("Writing aggregate plots to: ", pdf_path, "\n", sep = "")

# ──────────────────────────────────────────────────────────────────────────────
# 2. BACKGROUND KDE (figure #0 — shared across all species)
# ──────────────────────────────────────────────────────────────────────────────
kde            <- MASS::kde2d(dt$PC1, dt$PC2, n = 200)
kde_df         <- expand.grid(PC1 = kde$x, PC2 = kde$y)
kde_df$density <- as.vector(kde$z)
kde_df$density <- kde_df$density / max(kde_df$density)

p_background <- ggplot(kde_df, aes(x = PC1, y = PC2)) +
  geom_contour_filled(aes(z = density),
                      breaks = c(0, 0.05, 0.25, 0.5, 0.75, 1),
                      alpha  = 0.80) +
  scale_fill_viridis_d(option = "mako", name = "Background\ndensity [0–1]",
                       direction = -1) +
  coord_equal() +
  labs(title    = "Available environmental space",
       subtitle = "Kernel density of background cells in PC1 x PC2",
       x = pc1_lab, y = pc2_lab) +
  theme_classic(base_size = 13)

print(p_background)

# ──────────────────────────────────────────────────────────────────────────────
# 3. ONE-TIME COMPUTATIONS (fixed across all species for comparability)
# ──────────────────────────────────────────────────────────────────────────────
bw_background <- compute_bandwidth(dt)
cat("Background bandwidth (PC1, PC2):", round(bw_background, 5), "\n")

dt_sf  <- sf::st_as_sf(dt, coords = c("PC1", "PC2"))
res_opt <- USE.MCMC::optimRes(
  sdf      = dt_sf,
  grid.res = 1:10,
  perc.thr = 20,
  showOpt  = FALSE,
  cr       = 5
)
grid_res_opt <- as.numeric(res_opt$Opt_res)
cat("Optimal USE grid resolution:", grid_res_opt, "\n")

# ──────────────────────────────────────────────────────────────────────────────
# 4. SPECIES CATALOGUE (user-editable)
#
# Drop or add entries here to control which species are evaluated. Layouts
# downstream adjust automatically to the length of this list.
# ──────────────────────────────────────────────────────────────────────────────
species_catalogue <- list(
  sp1_generalist_common = list(label = "Generalist · common",
                               mu = c(0, 0), sigma_PC1 = 3,   sigma_PC2 = 2),
  sp2_specialist_common = list(label = "Specialist · common",
                               mu = c(0, 0), sigma_PC1 = 1,   sigma_PC2 = 0.5),
  sp3_generalist_rare   = list(label = "Generalist · rare",
                               mu = c(4, 1), sigma_PC1 = 3,   sigma_PC2 = 2),
  sp4_specialist_rare   = list(label = "Specialist · rare",
                               mu = c(4, 1), sigma_PC1 = 1,   sigma_PC2 = 0.5)
)

N_REALIZATIONS <- 1    # smoke test; restore to 10 for production

# ──────────────────────────────────────────────────────────────────────────────
# 5. SHARED virtualSpecies ARGUMENTS
# ──────────────────────────────────────────────────────────────────────────────
shared_args <- list(
  dt                       = dt,
  envData                  = envData,
  rho                      = 0,
  max_pres                 = 1000,
  bgk_prev                 = 1,
  pa_samplers              = list(random  = pa_random,
                                  uniform = pa_uniform,
                                  mcmc    = pa_mcmc),
  n_realizations           = N_REALIZATIONS,
  seed_base                = 42,
  seed_pseudo_base         = 123,
  pc1_lab                  = pc1_lab,
  pc2_lab                  = pc2_lab,
  kde_df                   = kde_df,
  bw                       = bw_background,
  pa_env_rast              = envData,
  verbose                  = TRUE,
  parallel                 = TRUE,
  n_workers                = NULL,
  grid.res                 = grid_res_opt,
  thres                    = 0.75,
  chain.length             = 20000,
  burnIn                   = 1000,
  species.cutoff.threshold = 0.1
)

# ──────────────────────────────────────────────────────────────────────────────
# 6. RUN ALL SPECIES
# ──────────────────────────────────────────────────────────────────────────────
tic("All species")
species_results <- lapply(names(species_catalogue), function(nm) {
  cat("\n", strrep("█", 70), "\n",
      "▶ Running species: ", nm, "  (",
      species_catalogue[[nm]]$label, ")\n",
      strrep("█", 70), "\n", sep = "")
  spec <- species_catalogue[[nm]]
  do.call(virtualSpecies, c(shared_args, list(
    mu        = spec$mu,
    sigma_PC1 = spec$sigma_PC1,
    sigma_PC2 = spec$sigma_PC2
  )))
})
names(species_results) <- names(species_catalogue)
toc()

sampler_names  <- names(shared_args$pa_samplers)
species_labels <- vapply(species_catalogue, `[[`, character(1), "label")

# ──────────────────────────────────────────────────────────────────────────────
# 7. AGGREGATE-PLOT HELPERS
# ──────────────────────────────────────────────────────────────────────────────

# Pull a nested ggplot from the virtualSpecies return value by following a
# character path, e.g. c("response_curves","plots","p_rc_PC1").
.get_plot <- function(x, path) Reduce(`[[`, path, init = x)

# Each helper:
#   * replaces each panel's title with a short species (or species · sampler)
#     label,
#   * strips the original subtitle (long niche equations / legends overflow
#     into the next panel when plots are wrapped to <2-inch widths),
#   * collects legends via patchwork (otherwise every panel renders its own
#     "Background density [0–1]" legend and crushes the actual plots).

aggregate_row <- function(species_results, plot_path, title) {
  plots <- lapply(names(species_results), function(nm) {
    .get_plot(species_results[[nm]], plot_path) +
      labs(title = species_catalogue[[nm]]$label, subtitle = NULL) +
      theme(plot.title = element_text(size = 11))
  })
  n <- length(plots)
  # Near-square grid for >=4 species; single row for fewer.
  if (n >= 4L) {
    nc <- ceiling(sqrt(n)); nr <- ceiling(n / nc)
  } else {
    nc <- n; nr <- 1L
  }
  patchwork::wrap_plots(plots, nrow = nr, ncol = nc) +
    patchwork::plot_layout(guides = "collect") +
    patchwork::plot_annotation(
      title = title,
      theme = theme(plot.title = element_text(face = "bold"))
    )
}

aggregate_grid <- function(species_results, sampler_names, plot_path, title,
                           transpose = FALSE, strip_inner_axes = FALSE) {
  # transpose = TRUE swaps the grid so rows = samplers, cols = species. For
  # 4 species × 3 samplers in a wide PDF (14×10) this lets each panel grow
  # because there are more horizontal cells available.
  # strip_inner_axes = TRUE hides x-axis on every row except the bottom, and
  # y-axis on every col except the leftmost — useful when every panel shares
  # the same axis labels (e.g. "PC1 (43.3% var)") so duplicating them only
  # eats panel area.
  if (transpose) {
    nr <- length(sampler_names);            nc <- length(species_results)
    grid_iter <- expand.grid(s = sampler_names, sp = names(species_results),
                             stringsAsFactors = FALSE)
    grid_iter <- grid_iter[order(match(grid_iter$s, sampler_names),
                                 match(grid_iter$sp, names(species_results))), ]
  } else {
    nr <- length(species_results);          nc <- length(sampler_names)
    grid_iter <- expand.grid(sp = names(species_results), s = sampler_names,
                             stringsAsFactors = FALSE)
    grid_iter <- grid_iter[order(match(grid_iter$sp, names(species_results)),
                                 match(grid_iter$s,  sampler_names)), ]
  }

  plots <- vector("list", nrow(grid_iter))
  for (k in seq_len(nrow(grid_iter))) {
    sp_nm <- grid_iter$sp[k]
    s_nm  <- grid_iter$s[k]
    p <- .get_plot(species_results[[sp_nm]],
                   c("samplers", s_nm, plot_path)) +
      labs(title = paste(species_catalogue[[sp_nm]]$label, "·", s_nm),
           subtitle = NULL) +
      theme(plot.title = element_text(size = 10))

    if (strip_inner_axes) {
      row <- ((k - 1L) %/% nc) + 1L
      col <- ((k - 1L) %%  nc) + 1L
      if (col > 1L)
        p <- p + theme(axis.title.y = element_blank(),
                       axis.text.y  = element_blank(),
                       axis.ticks.y = element_blank())
      if (row < nr)
        p <- p + theme(axis.title.x = element_blank(),
                       axis.text.x  = element_blank(),
                       axis.ticks.x = element_blank())
    }
    plots[[k]] <- p
  }

  patchwork::wrap_plots(plots, nrow = nr, ncol = nc, byrow = TRUE) +
    patchwork::plot_layout(guides = "collect") +
    patchwork::plot_annotation(
      title = title,
      theme = theme(plot.title = element_text(face = "bold"))
    )
}

# ──────────────────────────────────────────────────────────────────────────────
# 8. AGGREGATE GGPLOT FIGURES  (#1–11)
# ──────────────────────────────────────────────────────────────────────────────

# #1–4 Response curves (per species)
p_agg_rc_PC1       <- aggregate_row(species_results,
                                    c("response_curves","plots","p_rc_PC1"),
                                    "Response curve · PC1")
p_agg_rc_PC2       <- aggregate_row(species_results,
                                    c("response_curves","plots","p_rc_PC2"),
                                    "Response curve · PC2")
p_agg_rc_logit_PC1 <- aggregate_row(species_results,
                                    c("response_curves","plots","p_rc_logit_PC1"),
                                    "Logit response curve · PC1")
p_agg_rc_logit_PC2 <- aggregate_row(species_results,
                                    c("response_curves","plots","p_rc_logit_PC2"),
                                    "Logit response curve · PC2")

# #5 Niche plot (per species)
p_agg_niche <- aggregate_row(species_results, c("plots","p_niche"),
                             "Niche vs available environment")

# #6–8 PA / bias plots (per species × sampler)
p_agg_pa      <- aggregate_grid(species_results, sampler_names,
                                c("plots","p_pa"),
                                "Pseudo-absences in PC space",
                                transpose = TRUE, strip_inner_axes = TRUE)
p_agg_biasPC1 <- aggregate_grid(species_results, sampler_names,
                                c("plots","p_bias_PC1"),
                                "Bias · PC1")
p_agg_biasPC2 <- aggregate_grid(species_results, sampler_names,
                                c("plots","p_bias_PC2"),
                                "Bias · PC2")

# #9–11 Within-species boxplots (per species × sampler)
p_agg_box_overlap_local <- aggregate_grid(species_results, sampler_names,
                                          c("plots","p_box_overlap"),
                                          "Within-species overlap boxplots")
p_agg_box_cov_local     <- aggregate_grid(species_results, sampler_names,
                                          c("plots","p_box_coverage"),
                                          "Within-species coverage boxplots")
p_agg_box_trueabs_local <- aggregate_grid(species_results, sampler_names,
                                          c("plots","p_box_trueabs"),
                                          "Within-species true-absence boxplots")

print(p_agg_rc_PC1)
print(p_agg_rc_PC2)
print(p_agg_rc_logit_PC1)
print(p_agg_rc_logit_PC2)
print(p_agg_niche)
print(p_agg_pa)
print(p_agg_biasPC1)
print(p_agg_biasPC2)
print(p_agg_box_overlap_local)
print(p_agg_box_cov_local)
print(p_agg_box_trueabs_local)

# ──────────────────────────────────────────────────────────────────────────────
# 9. CROSS-SAMPLER METRICS TABLE + BOXPLOTS  (#13–16)
# ──────────────────────────────────────────────────────────────────────────────
metrics_all <- do.call(rbind, lapply(names(species_results), function(sp_nm) {
  sp <- species_results[[sp_nm]]
  do.call(rbind, lapply(names(sp$samplers), function(s_nm) {
    m <- sp$samplers[[s_nm]]$metrics
    cbind(species       = sp_nm,
          species_label = species_catalogue[[sp_nm]]$label,
          sampler       = s_nm,
          m,
          stringsAsFactors = FALSE)
  }))
}))
metrics_all$species       <- factor(metrics_all$species,
                                    levels = names(species_catalogue))
metrics_all$species_label <- factor(metrics_all$species_label,
                                    levels = species_labels)
metrics_all$sampler       <- factor(metrics_all$sampler,
                                    levels = sampler_names)
metrics_ok <- subset(metrics_all, status == "ok")

plot_metric <- function(df, metric, ylab, title, subtitle = NULL) {
  ggplot(df, aes(x = sampler, y = .data[[metric]], fill = sampler)) +
    geom_boxplot(alpha = 0.7, width = 0.55, outlier.shape = 21) +
    geom_jitter(width = 0.1, alpha = 0.4, size = 0.8) +
    scale_fill_brewer(palette = "Set1") +
    facet_wrap(~ species_label, nrow = 1) +
    labs(title = title, subtitle = subtitle, x = NULL, y = ylab) +
    theme_classic(base_size = 13) +
    theme(legend.position = "none",
          strip.background = element_rect(fill = "grey95", colour = NA))
}

p_xs_overlap <- plot_metric(
  metrics_ok, "overlap", "Overlap (intersection volume)",
  "Hypervolume overlap across samplers")
p_xs_covPC1  <- plot_metric(
  metrics_ok, "rel_cov_PC1", "Coverage PC1",
  "PC1 coverage across samplers")
p_xs_covPC2  <- plot_metric(
  metrics_ok, "rel_cov_PC2", "Coverage PC2",
  "PC2 coverage across samplers")
p_xs_trueabs <- plot_metric(
  metrics_ok, "prop_true_abs", "Proportion",
  "Proportion of pseudo-absences on true-absence cells",
  "Higher = sampler avoids true-presence cells more effectively")

p_xs_all <- (p_xs_overlap / p_xs_covPC1 / p_xs_covPC2 / p_xs_trueabs) +
  patchwork::plot_annotation(
    title = "Sampler comparison across species",
    theme = theme(plot.title = element_text(face = "bold"))
  )

print(p_xs_overlap)
print(p_xs_covPC1)
print(p_xs_covPC2)
print(p_xs_trueabs)
print(p_xs_all)

# ──────────────────────────────────────────────────────────────────────────────
# 10. GEOGRAPHIC AGGREGATES  (#17–20, base-R)
# ──────────────────────────────────────────────────────────────────────────────
plot_geo_row <- function(species_results, slot, palette, title) {
  n <- length(species_results)
  op <- par(mfrow = c(1, n),
            mar  = c(2, 2, 3, 1),
            oma  = c(0, 0, 2.5, 0))
  on.exit(par(op))
  for (nm in names(species_results)) {
    plot(species_results[[nm]][[slot]],
         col  = palette,
         main = species_catalogue[[nm]]$label)
  }
  mtext(title, outer = TRUE, line = 0.5, font = 2, cex = 1.1)
}

plot_geo_grid_pa <- function(species_results, sampler_names,
                             rast_slot, pa_colour_by_sampler, title) {
  n_sp <- length(species_results)
  n_s  <- length(sampler_names)
  op <- par(mfrow = c(n_sp, n_s),
            mar  = c(2, 2, 3, 1),
            oma  = c(0, 0, 2.5, 0))
  on.exit(par(op))
  for (nm in names(species_results)) {
    sp <- species_results[[nm]]
    rast_pal <- if (rast_slot == "suit_rast")
      hcl.colors(100, "YlOrRd", rev = TRUE) else c("grey85", "firebrick")
    for (s in sampler_names) {
      plot(sp[[rast_slot]],
           col  = rast_pal,
           main = sprintf("%s · %s", species_catalogue[[nm]]$label, s))
      pv <- sp$samplers[[s]]$pseudo_vect
      if (!is.null(pv))
        points(pv, col = pa_colour_by_sampler[[s]], cex = 0.4, pch = 16)
    }
  }
  mtext(title, outer = TRUE, line = 0.5, font = 2, cex = 1.1)
}

pa_cols <- list(random = "steelblue", uniform = "darkorange", mcmc = "purple")

plot_geo_row(species_results, "suit_rast",
             hcl.colors(100, "YlOrRd", rev = TRUE),
             "Suitability (geographic) — per species")
plot_geo_row(species_results, "pa_rast",
             c("grey85", "firebrick"),
             "Presences (geographic) — per species")
plot_geo_grid_pa(species_results, sampler_names, "suit_rast", pa_cols,
                 "Pseudo-absences on suitability raster (species × sampler)")
plot_geo_grid_pa(species_results, sampler_names, "pa_rast", pa_cols,
                 "Pseudo-absences on presence raster (species × sampler)")

# ──────────────────────────────────────────────────────────────────────────────
# 11. STATUS / PREVALENCE SUMMARY
# ──────────────────────────────────────────────────────────────────────────────
cat("\n── Realisation status by species × sampler ─────────────────────\n")
status_tbl <- with(metrics_all,
                   table(species = species_label, sampler = sampler,
                         status = status))
print(status_tbl)

cat("\n── Prevalence (mean presence per background cell) ──────────────\n")
prevalence_summary <- sapply(species_results, function(sp)
  round(mean(sp$background$pa), 3))
print(prevalence_summary)

cat("\n── Niche equations per species ─────────────────────────────────\n")
for (nm in names(species_results)) {
  cat("• ", species_catalogue[[nm]]$label, " (", nm, ")\n",
      "   PC1  : ", species_results[[nm]]$niche$equations$PC1,   "\n",
      "   PC2  : ", species_results[[nm]]$niche$equations$PC2,   "\n",
      "   2-D  : ", species_results[[nm]]$niche$equations$two_d, "\n",
      "   Logit: ", species_results[[nm]]$niche$equations$logit, "\n",
      sep = "")
}

dev.off()
cat("\nDone. Aggregate plots saved to: ", pdf_path, "\n", sep = "")
