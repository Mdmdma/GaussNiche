# virtualSpecies_fn.R
# Function wrapping steps 3-8 of the virtual species exercise
# Source this script, then call virtualSpecies() from the main script

library(mvtnorm)
library(MASS)
library(ggplot2)
library(hypervolume)
library(patchwork)
library(terra)

# =============================================================================
# PSEUDO-ABSENCE SAMPLERS
# Each sampler takes (background, N_pa, ...) and returns a data.frame
# with at minimum columns x, y, PC1, PC2 — same structure as dt
# New strategies plug in here without touching the main function
# =============================================================================

# random: draw uniformly from the full background
pa_random <- function(background, N_pa, seed = 123, ...) {
  set.seed(seed)
  background[sample(nrow(background), size = N_pa, replace = FALSE), ]
}

# placeholder for future strategies, e.g.:
# pa_buffer_out <- function(background, N_pa, pres, buffer_dist, ...) { ... }
# pa_uniform    <- function(background, N_pa, grid_res, ...) { ... }

# =============================================================================
# MAIN FUNCTION
# =============================================================================

virtualSpecies <- function(
    dt,                        # background data.frame: x, y, PC1, PC2
    envData,                   # original SpatRaster (needed for CRS only)
    mu,                        # numeric[2]: niche optimum in PC space
    sigma_PC1,                 # niche breadth along PC1
    sigma_PC2,                 # niche breadth along PC2
    rho       = 0,             # correlation between PC1 and PC2 axes [-1, 1]
    bgk_prev  = 1,             # presence:absence ratio (1 = 1:1, 2 = 1:2)
    pa_sampler = pa_random,    # function(background, N_pa, ...) → data.frame
    seed_pa   = 42,            # seed for Bernoulli presence sampling
    seed_pseudo = 123,         # seed forwarded to pa_sampler
    pc1_lab   = "PC1",        # axis label (pass pre-formatted string)
    pc2_lab   = "PC2",
    kde_df    = NULL,          # pre-computed background KDE data.frame
    ...                        # extra args forwarded to pa_sampler
) {
  
  # --- 3. BUILD NICHE (bivariate Gaussian) ----------------------------------
  
  Sigma <- matrix(
    c(sigma_PC1^2,
      rho * sigma_PC1 * sigma_PC2,
      rho * sigma_PC1 * sigma_PC2,
      sigma_PC2^2),
    nrow = 2
  )
  
  # normalise by theoretical peak: suitability = 1.0 exactly at mu
  peak_density <- mvtnorm::dmvnorm(
    x     = matrix(mu, nrow = 1),
    mean  = mu,
    sigma = Sigma
  )
  
  suitability <- mvtnorm::dmvnorm(
    x     = as.matrix(dt[, c("PC1", "PC2")]),
    mean  = mu,
    sigma = Sigma
  ) / peak_density
  
  cat("Suitability range:", round(range(suitability), 4), "\n")
  
  # --- 4. BERNOULLI SAMPLING ------------------------------------------------
  
  set.seed(seed_pa)
  pa      <- rbinom(n = length(suitability), size = 1, prob = suitability)
  dt$suit <- suitability
  dt$pa   <- pa
  
  cat("Prevalence:", round(mean(pa), 3), "\n")
  cat("N presences:", sum(pa), "of", nrow(dt), "cells\n")
  
  # --- 5. CONTOUR GRID (for overlay plots) ----------------------------------
  
  grid_pc <- expand.grid(
    PC1 = seq(min(dt$PC1), max(dt$PC1), length.out = 300),
    PC2 = seq(min(dt$PC2), max(dt$PC2), length.out = 300)
  )
  grid_pc$suit <- mvtnorm::dmvnorm(
    as.matrix(grid_pc), mean = mu, sigma = Sigma
  ) / peak_density
  
  # compute KDE if not supplied (avoids recomputing when called iteratively)
  if (is.null(kde_df)) {
    kde        <- MASS::kde2d(dt$PC1, dt$PC2, n = 200)
    kde_df     <- expand.grid(PC1 = kde$x, PC2 = kde$y)
    kde_df$density <- as.vector(kde$z)
  }
  
  # --- 6. PSEUDO-ABSENCE SAMPLING -------------------------------------------
  
  N_pa   <- round(sum(pa) / bgk_prev)
  pseudo <- pa_sampler(background = dt, N_pa = N_pa, seed = seed_pseudo, ...)
  
  cat("PA sampler used:", deparse(substitute(pa_sampler)), "\n")
  cat("N pseudo-absences sampled:", nrow(pseudo), "\n")
  
  # --- 7. CLASS OVERLAP (hypervolume intersection) --------------------------
  
  hyp_pres <- hypervolume(subset(dt, pa == 1)[, c("PC1", "PC2")], verbose = FALSE)
  hyp_pa   <- hypervolume(pseudo[, c("PC1", "PC2")],               verbose = FALSE)
  
  ovrlp <- get_volume(
    hypervolume_set(hyp_pres, hyp_pa, check.memory = FALSE, verbose = FALSE)
  )[[3]]   # [[1]] pres volume, [[2]] abs volume, [[3]] intersection
  
  cat("Class overlap (hypervolume intersection):", round(ovrlp, 4), "\n")
  
  # --- 8. SAMPLING BIAS (range coverage on each PC axis) -------------------
  
  bg_range_PC1 <- diff(range(dt$PC1))
  bg_range_PC2 <- diff(range(dt$PC2))
  pa_range_PC1 <- diff(range(pseudo$PC1))
  pa_range_PC2 <- diff(range(pseudo$PC2))
  rel_cov_PC1  <- round(pa_range_PC1 / bg_range_PC1, 3)
  rel_cov_PC2  <- round(pa_range_PC2 / bg_range_PC2, 3)
  
  cat("PC1 — background range:", round(bg_range_PC1, 2),
      "| PA range:", round(pa_range_PC1, 2),
      "| relative coverage:", rel_cov_PC1, "\n")
  cat("PC2 — background range:", round(bg_range_PC2, 2),
      "| PA range:", round(pa_range_PC2, 2),
      "| relative coverage:", rel_cov_PC2, "\n")
  
  # --- PLOTS ----------------------------------------------------------------
  
  # niche vs background
  p_niche <- ggplot() +
    geom_contour_filled(data = kde_df, aes(x = PC1, y = PC2, z = density),
                        alpha = 0.75) +
    scale_fill_viridis_d(option = "mako", name = "Background\ndensity",
                         direction = -1) +
    geom_contour(data = grid_pc, aes(x = PC1, y = PC2, z = suit),
                 colour = "firebrick", linewidth = 0.7,
                 breaks = c(0.1, 0.25, 0.5, 0.75, 0.95)) +
    geom_point(data = subset(dt, pa == 1), aes(x = PC1, y = PC2),
               colour = "gray40", alpha = 0.25, size = 0.6) +
    annotate("point", x = mu[1], y = mu[2],
             colour = "firebrick", size = 4, shape = 3, stroke = 2) +
    coord_equal() +
    labs(
      title    = "Niche vs available environment",
      subtitle = "Red contours = suitability (0.1, 0.25, 0.50, 0.75, 0.95) | + = optimum | dots = presences",
      x = pc1_lab, y = pc2_lab
    ) +
    theme_classic(base_size = 13)
  
  # presences vs pseudo-absences
  p_pa <- ggplot() +
    geom_contour_filled(data = kde_df, aes(x = PC1, y = PC2, z = density),
                        alpha = 0.75) +
    scale_fill_viridis_d(option = "mako", name = "Background\ndensity",
                         direction = -1) +
    geom_point(data = pseudo, aes(x = PC1, y = PC2),
               colour = "steelblue", alpha = 0.35, size = 0.7) +
    geom_point(data = subset(dt, pa == 1), aes(x = PC1, y = PC2),
               colour = "firebrick", alpha = 0.35, size = 0.7) +
    geom_contour(data = grid_pc, aes(x = PC1, y = PC2, z = suit),
                 colour = "firebrick", linewidth = 0.5, linetype = "dashed",
                 breaks = c(0.5, 0.95)) +
    annotate("text", x = min(dt$PC1), y = max(dt$PC2),
             label = paste0("Hypervolume overlap = ", round(ovrlp, 3)),
             hjust = 0, vjust = 1, size = 4, colour = "white") +
    coord_equal() +
    labs(
      title    = "Presences vs pseudo-absences in PC space",
      subtitle = "Red = presences | Blue = pseudo-absences | Dashed = suitability 0.50, 0.95",
      x = pc1_lab, y = pc2_lab
    ) +
    theme_classic(base_size = 13)
  
  # sampling bias: marginal densities on each PC axis
  bias_df <- rbind(
    data.frame(PC1 = dt$PC1,                  PC2 = dt$PC2,                  group = "Background"),
    data.frame(PC1 = pseudo$PC1,              PC2 = pseudo$PC2,              group = "Pseudo-absences"),
    data.frame(PC1 = subset(dt, pa == 1)$PC1, PC2 = subset(dt, pa == 1)$PC2, group = "Presences")
  )
  bias_df$group <- factor(bias_df$group,
                          levels = c("Background", "Pseudo-absences", "Presences"))
  
  col_vals  <- c("Background" = "grey50", "Pseudo-absences" = "steelblue",
                 "Presences"  = "firebrick")
  fill_vals <- c("Background" = "grey70", "Pseudo-absences" = "steelblue",
                 "Presences"  = "firebrick")
  
  p_bias_PC1 <- ggplot(bias_df, aes(x = PC1, colour = group, fill = group)) +
    geom_density(alpha = 0.25, linewidth = 0.8) +
    scale_colour_manual(values = col_vals) +
    scale_fill_manual(values   = fill_vals) +
    annotate("text", x = Inf, y = Inf,
             label = paste0("Coverage = ", rel_cov_PC1),
             hjust = 1.1, vjust = 1.5, size = 4) +
    labs(title = "Sampling bias — PC1", x = pc1_lab, y = "Density",
         colour = NULL, fill = NULL) +
    theme_classic(base_size = 13) +
    theme(legend.position = "bottom")
  
  p_bias_PC2 <- ggplot(bias_df, aes(x = PC2, colour = group, fill = group)) +
    geom_density(alpha = 0.25, linewidth = 0.8) +
    scale_colour_manual(values = col_vals) +
    scale_fill_manual(values   = fill_vals) +
    annotate("text", x = Inf, y = Inf,
             label = paste0("Coverage = ", rel_cov_PC2),
             hjust = 1.1, vjust = 1.5, size = 4) +
    labs(title = "Sampling bias — PC2", x = pc2_lab, y = "Density",
         colour = NULL, fill = NULL) +
    theme_classic(base_size = 13) +
    theme(legend.position = "bottom")
  
  # --- RASTERS --------------------------------------------------------------
  
  suit_rast   <- rast(dt[, c("x", "y", "suit")], type = "xyz", crs = crs(envData))
  pa_rast     <- rast(dt[, c("x", "y", "pa")],   type = "xyz", crs = crs(envData))
  pseudo_vect <- vect(pseudo, geom = c("x", "y"), crs = crs(envData))
  
  # --- COMBINED DATASET (presences + pseudo-absences) ----------------------
  # pa column: 1 = presence, 0 = pseudo-absence
  
  pres_df        <- subset(dt, pa == 1)[, c("x", "y", "PC1", "PC2", "suit")]
  pres_df$pa     <- 1L
  pseudo_df      <- pseudo[, c("x", "y", "PC1", "PC2", "suit")]
  pseudo_df$pa   <- 0L
  dataset        <- rbind(pres_df, pseudo_df)
  
  # --- OUTPUT LIST ----------------------------------------------------------
  
  list(
    # spatial outputs
    suit_rast   = suit_rast,
    pa_rast     = pa_rast,
    pseudo_vect = pseudo_vect,
    # tabular outputs
    background  = dt,           # full background with suit and pa columns
    dataset     = dataset,      # presences + pseudo-absences only
    # metrics
    metrics = list(
      N_pres       = sum(pa),
      N_pseudo     = nrow(pseudo),
      prevalence   = round(mean(pa), 3),
      overlap      = round(ovrlp, 4),
      bg_range_PC1 = round(bg_range_PC1, 2),
      bg_range_PC2 = round(bg_range_PC2, 2),
      pa_range_PC1 = round(pa_range_PC1, 2),
      pa_range_PC2 = round(pa_range_PC2, 2),
      rel_cov_PC1  = rel_cov_PC1,
      rel_cov_PC2  = rel_cov_PC2
    ),
    # diagnostic plots
    plots = list(
      p_niche     = p_niche,
      p_pa        = p_pa,
      p_bias_PC1  = p_bias_PC1,
      p_bias_PC2  = p_bias_PC2
    )
  )
}