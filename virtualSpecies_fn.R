# virtualSpecies_fn.R  ─ v2
# =============================================================================
# Bivariate-Gaussian virtual-species simulation with:
#   1. Marginal response curves + log-suitability / logit equations
#   2. n_realizations Bernoulli draws for confidence intervals on metrics
#   3. Pluggable pseudo-absence samplers stored in a list-of-lists
#   4. Uniform-in-E-space sampler via USE::paSampling
#   5. Proportion of pseudo-absences falling on true-absence cells
#   6. hypervolume_gaussian with a fixed background bandwidth
#
# Source this file, then call virtualSpecies() from your main script.
# =============================================================================

suppressPackageStartupMessages({
  library(mvtnorm)
  library(MASS)
  library(ggplot2)
  library(hypervolume)
  library(patchwork)
  library(terra)
  library(USE)
})

# =============================================================================
# PSEUDO-ABSENCE SAMPLERS
# Shared interface:
#   function(background, N_pa, pres = NULL, seed = 123, ...)
#   'background'  : data.frame with columns x, y, PC1, PC2, suit
#   'N_pa'        : integer, number of pseudo-absences requested
#   'pres'        : data.frame of presence rows from background (needed by
#                   samplers that use presence information, e.g. pa_uniform)
#   '...'         : extra args forwarded from virtualSpecies()
#   Returns a data.frame with the same columns as background.
# =============================================================================

#' Random pseudo-absence sampler
#' Draws N_pa rows uniformly from the full background (no weighting).
pa_random <- function(background, N_pa, pres = NULL, seed = 123, ...) {
  set.seed(seed)
  background[sample(nrow(background), size = N_pa, replace = FALSE), ]
}


#' Uniform-in-E-space pseudo-absence sampler (USE package)
#'
#' Samples pseudo-absences so they are spread uniformly across the
#' environmental (PC) space, with a kernel-density filter that excludes the
#' region of E-space most likely to be suitable for the species.
#'
#' @param background  Full background data.frame (x, y, PC1, PC2, suit, …)
#' @param N_pa        Number of pseudo-absences to return
#' @param pres        Presence rows from background (required)
#' @param seed        RNG seed
#' @param env.rast    SpatRaster of the ORIGINAL environmental variables (e.g.
#'                    envData / WorldClim), NOT the PC rasters. USE::paSampling
#'                    always runs rastPCA() internally — if you feed it PC
#'                    rasters, it re-runs PCA on already-orthogonal variables,
#'                    producing a rotated E-space where the presence KDE filter
#'                    no longer excludes the correct region, making uniform and
#'                    random sampling indistinguishable. Pass envData; the
#'                    function retrieves our analysis PC1/PC2 values afterwards
#'                    by joining on geographic coordinates.
#' @param grid.res    Grid resolution for USE::paSampling (default 5).
#'                    Pre-compute with USE::optimRes() and pass via '...'.
#' @param thres       Kernel-density threshold for excluding suitable cells
#'                    (default 0.75; lower = more exclusion).
#' @param ...         Ignored additional arguments (for interface compatibility)
pa_uniform <- function(background, N_pa, pres = NULL, seed = 123,
                       env.rast = NULL, grid.res = 5, thres = 0.75, ...) {
  
  if (!requireNamespace("sf", quietly = TRUE))
    stop("Package 'sf' is required for pa_uniform. ",
         "Install with: install.packages('sf')")
  
  if (is.null(pres))
    stop("pa_uniform requires pres (presence rows from the background).")
  
  if (is.null(env.rast))
    stop("pa_uniform requires 'env.rast' (the original environmental SpatRaster, ",
         "e.g. envData). Pass pa_env_rast = envData in virtualSpecies(). ",
         "Do NOT pass rpc$PCs: USE runs rastPCA() internally and rePCA-ing ",
         "PC scores produces a rotated E-space that breaks the exclusion filter.")
  
  set.seed(seed)
  
  # Build sf object of presences in geographic space so USE can project them
  pres_sf <- sf::st_as_sf(
    pres[, c("x", "y"), drop = FALSE],
    coords = c("x", "y"),
    crs    = terra::crs(env.rast)
  )
  
  # Approximate n.tr: target N_pa total across all non-excluded grid cells.
  # Use a generous ceiling; we subsample down to N_pa afterwards.
  # n.tr = PA per grid cell. Over-request by 3x (unknown excluded fraction),
  # then subsample to N_pa. Floor of 5 ensures non-empty cells are covered.
  n_tr_approx <- max(5L, ceiling(N_pa * 3L / max(grid.res^2L, 10L)))
  
  pa_result <- tryCatch(
    USE::paSampling(
      env.rast  = env.rast,
      pres      = pres_sf,
      thres     = thres,
      H         = NULL,
      grid.res  = grid.res,
      n.tr      = n_tr_approx,
      prev      = 1,
      sub.ts    = FALSE,
      plot_proc = FALSE,
      verbose   = FALSE
    ),
    error = function(e) {
      warning("USE::paSampling failed ('", conditionMessage(e),
              "'). Falling back to pa_random.")
      NULL
    }
  )
  
  if (is.null(pa_result) || nrow(pa_result) == 0L) {
    warning("pa_uniform returned 0 pseudo-absences; falling back to pa_random.")
    return(pa_random(background, N_pa, pres, seed))
  }
  
  # When sub.ts = FALSE, USE returns the sf object directly (no $obs.tr).
  # The sf GEOMETRY is in PC / E-space; geographic coordinates are stored
  # as plain attribute columns 'x' and 'y' — retrieve them with st_drop_geometry().
  pa_attrs <- sf::st_drop_geometry(pa_result)
  
  # Match back to background rows by geographic coordinates (4 dp = ~11 m
  # precision, safe for any raster resolution >= 1 arc-second)
  bg_key <- paste(round(background$x, 4L), round(background$y, 4L))
  pa_key <- paste(round(pa_attrs$x,   4L), round(pa_attrs$y,   4L))
  
  idx       <- match(pa_key, bg_key)
  idx_valid <- unique(idx[!is.na(idx)])
  pa_bg     <- background[idx_valid, ]
  if (length(idx_valid) == 0L)
    warning("pa_uniform: NO geographic coordinates matched the background. ",
            "Check that env.rast covers the same extent as background and ",
            "that you passed envData (not PC rasters) as pa_env_rast.")
  if (length(idx_valid) < N_pa * 0.5)
    warning("pa_uniform: only ", length(idx_valid), " of ", N_pa,
            " requested pseudo-absences matched background coordinates. ",
            "Remainder will be padded with random draws.")
  
  # Adjust to exactly N_pa ─────────────────────────────────────────────────
  if (nrow(pa_bg) >= N_pa) {
    pa_bg <- pa_bg[sample(nrow(pa_bg), N_pa, replace = FALSE), ]
  } else {
    # Pad with random draws from non-presence background cells
    pres_key  <- paste(round(pres$x, 4L), round(pres$y, 4L))
    all_key   <- paste(round(background$x, 4L), round(background$y, 4L))
    non_pres  <- which(!all_key %in% pres_key)
    remaining <- setdiff(non_pres, idx_valid)
    extra_n   <- N_pa - nrow(pa_bg)
    if (length(remaining) >= extra_n) {
      pa_bg <- rbind(pa_bg, background[sample(remaining, extra_n), ])
    } else if (length(remaining) > 0L) {
      pa_bg <- rbind(pa_bg, background[remaining, ])
      warning("pa_uniform: obtained ", nrow(pa_bg),
              " pseudo-absences instead of ", N_pa, " (E-space constraint).")
    }
  }
  
  row.names(pa_bg) <- NULL
  pa_bg
}


# =============================================================================
# HELPER — one-time bandwidth estimation from the environmental background
# =============================================================================

#' Estimate KDE bandwidth from the full background E-space.
#' Call this ONCE and pass the result as bw = ... to virtualSpecies() so the
#' same bandwidth is used for every species and every sampler — ensuring
#' comparability of hypervolume estimates across the experiment.
#'
#' @param dt  Background data.frame with columns PC1, PC2
#' @return    Named numeric vector (length 2) — one bandwidth per dimension
compute_bandwidth <- function(dt) {
  hypervolume::estimate_bandwidth(dt[, c("PC1", "PC2")])
}


# =============================================================================
# MAIN FUNCTION
# =============================================================================

#' Simulate a virtual species and evaluate pseudo-absence samplers
#'
#' @param dt               Background data.frame: x, y, PC1, PC2
#' @param envData          Original SpatRaster (used for CRS and raster output)
#' @param mu               Numeric[2]: niche optimum in PC1, PC2
#' @param sigma_PC1        Niche breadth along PC1
#' @param sigma_PC2        Niche breadth along PC2
#' @param rho              PC1-PC2 correlation [-1, 1]; default 0
#' @param bgk_prev         Presence:pseudo-absence ratio (1 = 1:1)
#' @param pa_samplers      Named list of sampler functions; each must follow the
#'                         shared interface described above
#' @param n_realizations   Number of Bernoulli draws for CI computation
#' @param max_pres         Maximum number of presences used per realisation
#'                         for hypervolume computation and PA sampling.
#'                         The full Bernoulli draw is kept for prevalence
#'                         reporting and true-absence matching; only the
#'                         subsample is passed to hypervolume_gaussian() and
#'                         the PA sampler. Caps the dominant runtime cost.
#'                         Default 500.
#' @param seed_base        Base RNG seed; realization r uses seed_base + r
#' @param seed_pseudo_base Base seed for PA draws; realization r uses
#'                         seed_pseudo_base + r
#' @param pc1_lab          Axis label for PC1 (pass pre-formatted string)
#' @param pc2_lab          Axis label for PC2
#' @param kde_df           Pre-computed background KDE data.frame; computed
#'                         internally if NULL
#' @param bw               Pre-computed hypervolume bandwidth (from
#'                         compute_bandwidth()); computed internally if NULL.
#'                         Strongly recommended to pass externally so it is
#'                         fixed across multiple species.
#' @param pa_env_rast      SpatRaster forwarded to samplers as env.rast.
#'                         Must be the ORIGINAL environmental raster (envData),
#'                         NOT rpc$PCs. USE::paSampling always calls rastPCA()
#'                         internally; feeding it PC rasters re-rotates the E-space
#'                         and breaks the presence-exclusion filter, making
#'                         pa_uniform indistinguishable from pa_random.
#' @param verbose          Print progress messages
#' @param ...              Extra arguments forwarded to every sampler
#'                         (e.g. grid.res = 5, thres = 0.75 for pa_uniform)
#'
#' @return A named list — see "Output structure" section in comments below
virtualSpecies <- function(
    dt,
    envData,
    mu,
    sigma_PC1,
    sigma_PC2,
    rho              = 0,
    bgk_prev         = 1,
    pa_samplers      = list(random = pa_random),
    n_realizations   = 50L,
    max_pres         = 500L,
    seed_base        = 42L,
    seed_pseudo_base = 123L,
    pc1_lab          = "PC1",
    pc2_lab          = "PC2",
    kde_df           = NULL,
    bw               = NULL,
    pa_env_rast      = NULL,
    verbose          = TRUE,
    ...
) {
  
  # Default the PA environment raster to envData if not supplied
  if (is.null(pa_env_rast)) pa_env_rast <- envData
  
  # ── 1. BIVARIATE GAUSSIAN NICHE ────────────────────────────────────────────
  
  Sigma <- matrix(
    c(sigma_PC1^2,
      rho * sigma_PC1 * sigma_PC2,
      rho * sigma_PC1 * sigma_PC2,
      sigma_PC2^2),
    nrow = 2L
  )
  
  # Normalise so suitability = 1 exactly at the optimum mu
  peak_density <- mvtnorm::dmvnorm(
    x = matrix(mu, nrow = 1L), mean = mu, sigma = Sigma
  )
  suitability <- mvtnorm::dmvnorm(
    x     = as.matrix(dt[, c("PC1", "PC2")]),
    mean  = mu,
    sigma = Sigma
  ) / peak_density
  
  dt$suit <- suitability
  if (verbose) cat("Suitability range:", round(range(suitability), 4L), "\n")
  
  # ── 2. RESPONSE CURVES & EQUATIONS ─────────────────────────────────────────
  #
  # Marginal response along PC1 (PC2 fixed at optimum mu[2]):
  #   suit(PC1) = exp(-0.5 * (PC1 - mu1)^2 / sigma_PC1^2)
  #
  # Full 2-D log-suitability (quadratic form of Sigma^{-1}):
  #   log(suit) = -0.5 * [(PC-mu)^T Sigma^{-1} (PC-mu)]
  #
  # Logit of suitability:
  #   logit(suit) = log(suit) - log(1 - suit)
  #   This is non-linear in PC axes; its shape is shown in the logit plots.
  
  Sigma_inv <- solve(Sigma)
  a    <- Sigma_inv[1L, 1L]
  b    <- Sigma_inv[1L, 2L]   # off-diagonal = Sigma_inv[2,1]
  cc   <- Sigma_inv[2L, 2L]
  
  pc1_seq <- seq(min(dt$PC1), max(dt$PC1), length.out = 500L)
  pc2_seq <- seq(min(dt$PC2), max(dt$PC2), length.out = 500L)
  
  rc_PC1 <- data.frame(
    PC1  = pc1_seq,
    suit = exp(-0.5 * (pc1_seq - mu[1L])^2 / sigma_PC1^2)
  )
  rc_PC2 <- data.frame(
    PC2  = pc2_seq,
    suit = exp(-0.5 * (pc2_seq - mu[2L])^2 / sigma_PC2^2)
  )
  # logit is ±Inf at suit = 0 or 1; keep only finite values for plotting
  rc_PC1$logit_suit <- log(rc_PC1$suit / (1 - rc_PC1$suit))
  rc_PC2$logit_suit <- log(rc_PC2$suit / (1 - rc_PC2$suit))
  
  # Equation strings (stored and printed for record-keeping)
  eq_PC1 <- sprintf(
    "suit(PC1 | PC2=mu2) = exp(-0.5 * (PC1 - %.3f)^2 / %.4f)",
    mu[1L], sigma_PC1^2
  )
  eq_PC2 <- sprintf(
    "suit(PC2 | PC1=mu1) = exp(-0.5 * (PC2 - %.3f)^2 / %.4f)",
    mu[2L], sigma_PC2^2
  )
  eq_2d <- sprintf(
    paste0("log(suit) = -0.5 * [%.4f*(PC1-%.3f)^2",
           " + 2*%.4f*(PC1-%.3f)*(PC2-%.3f)",
           " + %.4f*(PC2-%.3f)^2]"),
    a, mu[1L], b, mu[1L], mu[2L], cc, mu[2L]
  )
  eq_logit <- paste(
    "logit(suit) = log(suit) - log(1-suit)",
    "[non-linear in PC axes; inspect logit response curve plots]"
  )
  
  if (verbose) {
    cat("── Equations ──────────────────────────────────────────────\n")
    cat(" PC1  :", eq_PC1,  "\n")
    cat(" PC2  :", eq_PC2,  "\n")
    cat(" 2-D  :", eq_2d,   "\n")
    cat(" Logit:", eq_logit, "\n")
  }
  
  # Response curve plots
  p_rc_PC1 <- ggplot(rc_PC1, aes(x = PC1)) +
    geom_line(aes(y = suit), colour = "firebrick", linewidth = 1) +
    geom_rug(data = dt, aes(x = PC1), alpha = 0.04, colour = "grey40",
             length = unit(0.02, "npc")) +
    labs(title    = "Marginal response curve — PC1",
         subtitle = eq_PC1,
         x = pc1_lab, y = "Suitability [0, 1]") +
    theme_classic(base_size = 13)
  
  p_rc_PC2 <- ggplot(rc_PC2, aes(x = PC2)) +
    geom_line(aes(y = suit), colour = "firebrick", linewidth = 1) +
    geom_rug(data = dt, aes(x = PC2), alpha = 0.04, colour = "grey40",
             length = unit(0.02, "npc")) +
    labs(title    = "Marginal response curve — PC2",
         subtitle = eq_PC2,
         x = pc2_lab, y = "Suitability [0, 1]") +
    theme_classic(base_size = 13)
  
  rc_PC1_fin <- rc_PC1[is.finite(rc_PC1$logit_suit), ]
  rc_PC2_fin <- rc_PC2[is.finite(rc_PC2$logit_suit), ]
  
  p_rc_logit_PC1 <- ggplot(rc_PC1_fin, aes(x = PC1, y = logit_suit)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
    geom_line(colour = "steelblue", linewidth = 1) +
    labs(title = "Logit(suit) — PC1", x = pc1_lab, y = "logit(suit)") +
    theme_classic(base_size = 13)
  
  p_rc_logit_PC2 <- ggplot(rc_PC2_fin, aes(x = PC2, y = logit_suit)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
    geom_line(colour = "steelblue", linewidth = 1) +
    labs(title = "Logit(suit) — PC2", x = pc2_lab, y = "logit(suit)") +
    theme_classic(base_size = 13)
  
  # ── 3. BACKGROUND KDE + NICHE CONTOUR GRID ─────────────────────────────────
  
  if (is.null(kde_df)) {
    kde_raw        <- MASS::kde2d(dt$PC1, dt$PC2, n = 200L)
    kde_df         <- expand.grid(PC1 = kde_raw$x, PC2 = kde_raw$y)
    kde_df$density <- as.vector(kde_raw$z)
  }
  # Normalise density to [0,1] so the colour scale is comparable across
  # species regardless of background sample size.
  kde_df$density <- kde_df$density / max(kde_df$density)
  
  grid_pc <- expand.grid(
    PC1 = seq(min(dt$PC1), max(dt$PC1), length.out = 300L),
    PC2 = seq(min(dt$PC2), max(dt$PC2), length.out = 300L)
  )
  grid_pc$suit <- mvtnorm::dmvnorm(
    as.matrix(grid_pc), mean = mu, sigma = Sigma
  ) / peak_density
  
  # ── 4. HYPERVOLUME BANDWIDTH (fixed across species and samplers) ────────────
  #
  # estimate_bandwidth() is called on the FULL BACKGROUND so that all
  # hypervolumes — regardless of species niche location or PA sampler — are
  # computed at the same KDE resolution. This is essential for comparability.
  
  if (is.null(bw)) {
    if (verbose) cat("\nEstimating hypervolume bandwidth from background…\n")
    bw <- hypervolume::estimate_bandwidth(dt[, c("PC1", "PC2")])
  }
  if (verbose) cat("Bandwidth (PC1, PC2):", round(bw, 5L), "\n")
  
  # ── 5. PRE-DRAW ALL BERNOULLI REALISATIONS ─────────────────────────────────
  #
  # Separating the Bernoulli draws from the PA sampling loop makes it clear
  # that the 50 realisations quantify the stochastic variation in which
  # cells happen to be "presences" from the same underlying suitability
  # surface — i.e., a source of uncertainty inherent to the simulation, not
  # to the PA strategy.
  
  if (verbose) cat("\nDrawing", n_realizations, "Bernoulli realisations…\n")
  
  pa_matrix <- vapply(
    seq_len(n_realizations),
    FUN = function(r) {
      set.seed(seed_base + r)
      rbinom(nrow(dt), size = 1L, prob = suitability)
    },
    FUN.VALUE = integer(nrow(dt))
  )  # result: nrow(dt) x n_realizations integer matrix
  
  # Realization 1 is used as the reference for representative plots / rasters
  ref_pa       <- pa_matrix[, 1L]
  dt$pa        <- ref_pa
  ref_pres_all <- dt[ref_pa == 1L, ]   # full set — used for prevalence/rasters
  
  # Subsample presences for plots, PA sampling, and hypervolumes.
  # The full ref_pres_all is kept separately so prevalence is unaffected.
  set.seed(seed_base)
  ref_pres <- if (nrow(ref_pres_all) > max_pres)
    ref_pres_all[sample(nrow(ref_pres_all), max_pres), ]
  else ref_pres_all
  
  if (verbose) {
    cat("Prevalence (ref. realisation 1):", round(mean(ref_pa), 3L), "\n")
    cat("N presences total:", nrow(ref_pres_all),
        "| subsampled to:", nrow(ref_pres), "\n")
  }
  
  # Fixed background ranges (denominator for relative coverage metric)
  bg_range_PC1 <- diff(range(dt$PC1))
  bg_range_PC2 <- diff(range(dt$PC2))
  
  # ── 6. SHARED NICHE PLOT (reference realisation) ───────────────────────────
  
  p_niche <- ggplot() +
    geom_contour_filled(data = kde_df,
                        aes(x = PC1, y = PC2, z = density),
                        breaks = c(0, 0.05, 0.25, 0.5, 0.75, 1),
                        alpha  = 0.80) +
    scale_fill_viridis_d(option = "mako", name = "Background\ndensity [0–1]",
                         direction = -1) +
    geom_contour(data = grid_pc, aes(x = PC1, y = PC2, z = suit),
                 colour = "firebrick", linewidth = 0.7,
                 breaks = c(0.1, 0.25, 0.5, 0.75, 0.95)) +
    geom_point(data = ref_pres, aes(x = PC1, y = PC2),
               colour = "gray40", alpha = 0.25, size = 0.6) +
    annotate("point", x = mu[1L], y = mu[2L],
             colour = "firebrick", size = 4, shape = 3, stroke = 2) +
    coord_equal() +
    labs(
      title    = "Niche vs available environment",
      subtitle = "Red contours = suitability (0.1, 0.25, 0.50, 0.75, 0.95) | + = optimum | dots = presences",
      x = pc1_lab, y = pc2_lab
    ) +
    theme_classic(base_size = 13)
  
  # ── 7. SAMPLER LOOP ────────────────────────────────────────────────────────
  #
  # Output structure:
  #   sampler_results[[sampler_name]] = list(
  #     pseudo_ref  : reference pseudo-absences (realisation 1)
  #     pseudo_vect : SpatVector of reference pseudo-absences
  #     dataset     : combined presence + pseudo-absence data.frame (ref. real.)
  #     metrics     : data.frame — one row per realisation, columns:
  #                     realization, N_pres, N_pseudo, overlap,
  #                     rel_cov_PC1, rel_cov_PC2, prop_true_abs
  #     plots       : list(p_pa, p_bias_PC1, p_bias_PC2, p_boxplot)
  #   )
  
  sampler_results <- vector("list", length(pa_samplers))
  names(sampler_results) <- names(pa_samplers)
  
  for (s_name in names(pa_samplers)) {
    
    sampler  <- pa_samplers[[s_name]]
    if (verbose) cat("\n══ Sampler:", s_name, "══\n")
    
    # -- Reference pseudo-absences for diagnostic plots (realisation 1) ------
    N_pa_ref   <- max(1L, round(nrow(ref_pres) / bgk_prev))
    pseudo_ref <- tryCatch(
      sampler(
        background = dt, N_pa = N_pa_ref, pres = ref_pres,
        seed = seed_pseudo_base, env.rast = pa_env_rast, ...
      ),
      error = function(e) {
        stop("Sampler '", s_name, "' failed on reference realisation: ", e$message)
      }
    )
    
    # -- Metrics loop across all n_realizations --------------------------------
    #
    # NOTE ON PERFORMANCE: each iteration calls hypervolume_gaussian() twice
    # (presences + pseudo-absences). Presences are capped at max_pres (default
    # 500) so N_pa is also capped at max_pres / bgk_prev. This keeps each
    # hypervolume call fast regardless of species prevalence.
    
    metrics_list <- vector("list", n_realizations)
    
    for (r in seq_len(n_realizations)) {
      
      pa_r         <- pa_matrix[, r]
      pres_r_all   <- dt[pa_r == 1L, ]   # full draw — used for prop_true_abs
      true_abs_r   <- dt[pa_r == 0L, ]
      
      # Subsample presences to max_pres to cap hypervolume runtime.
      # N_pa is derived from the subsampled count so the PA:presence
      # ratio (bgk_prev) is consistent across all realisations.
      set.seed(seed_base + r)
      pres_r <- if (nrow(pres_r_all) > max_pres)
        pres_r_all[sample(nrow(pres_r_all), max_pres), ]
      else pres_r_all
      
      N_pa_r <- max(1L, round(nrow(pres_r) / bgk_prev))
      
      if (nrow(pres_r) < 5L || N_pa_r < 5L) {
        if (verbose) message("  r=", r, ": too few points, skipping.")
        next
      }
      
      pseudo_r <- tryCatch(
        sampler(
          background = dt, N_pa = N_pa_r, pres = pres_r,
          seed = seed_pseudo_base + r, env.rast = pa_env_rast, ...
        ),
        error = function(e) {
          warning("Sampler '", s_name, "' failed at r=", r, ": ", e$message)
          NULL
        }
      )
      if (is.null(pseudo_r) || nrow(pseudo_r) < 5L) next
      
      # Hypervolume — Gaussian KDE with fixed background bandwidth
      hyp_pres_r <- tryCatch(
        hypervolume::hypervolume_gaussian(
          data                    = pres_r[, c("PC1", "PC2")],
          kde.bandwidth           = bw,
          sd.count                = 3,
          quantile.requested      = 0.95,     # use 0.99 if following Enrico
          quantile.requested.type = "probability",
          chunk.size              = 1000L,
          verbose                 = FALSE
        ),
        error = function(e) NULL
      )
      hyp_pa_r <- tryCatch(
        hypervolume::hypervolume_gaussian(
          data                    = pseudo_r[, c("PC1", "PC2")],
          kde.bandwidth           = bw,
          sd.count                = 3,
          quantile.requested      = 0.95,
          quantile.requested.type = "probability",
          chunk.size              = 1000L,
          verbose                 = FALSE
        ),
        error = function(e) NULL
      )
      
      if (is.null(hyp_pres_r) || is.null(hyp_pa_r)) {
        if (verbose) message("  r=", r, ": hypervolume failed, skipping.")
        next
      }
      
      hv_set <- hypervolume::hypervolume_set(
        hyp_pres_r, hyp_pa_r, check.memory = FALSE, verbose = FALSE
      )
      # get_volume() returns: [[1]] pres vol, [[2]] pa vol, [[3]] intersection
      ovrlp <- hypervolume::get_volume(hv_set)[[3L]]
      
      # Range coverage: how much of the background PC range do the PAs cover?
      rel_cov_PC1 <- diff(range(pseudo_r$PC1)) / bg_range_PC1
      rel_cov_PC2 <- diff(range(pseudo_r$PC2)) / bg_range_PC2
      
      # Proportion of pseudo-absences that fall on TRUE-absence cells.
      # Matches on rounded (x, y) geographic coordinates to avoid float issues.
      # This metric flags samplers that inadvertently draw from true-presence
      # cells, which would inflate observed class overlap.
      ta_key  <- paste(round(true_abs_r$x, 4L), round(true_abs_r$y, 4L))
      ps_key  <- paste(round(pseudo_r$x,   4L), round(pseudo_r$y,   4L))
      prop_ta <- mean(ps_key %in% ta_key)
      
      metrics_list[[r]] <- data.frame(
        realization   = r,
        N_pres_total  = nrow(pres_r_all),
        N_pres_used   = nrow(pres_r),
        N_pseudo      = nrow(pseudo_r),
        overlap       = round(ovrlp, 6L),
        rel_cov_PC1   = round(rel_cov_PC1, 4L),
        rel_cov_PC2   = round(rel_cov_PC2, 4L),
        prop_true_abs = round(prop_ta, 4L)
      )
      
      if (verbose)
        cat(sprintf("\r  [%s] realisation %d / %d — %d remaining   ",
                    s_name, r, n_realizations, n_realizations - r),
            sep = "")
    }  # end realisation loop
    if (verbose) cat("\n")   # move cursor past the \r counter
    
    metrics_df <- do.call(rbind, Filter(Negate(is.null), metrics_list))
    
    # -- Plots for this sampler -----------------------------------------------
    
    # Proportion of reference PAs on true-absence cells (ref. realisation)
    ref_ta_key  <- paste(round(dt[ref_pa == 0L, "x"], 4L),
                         round(dt[ref_pa == 0L, "y"], 4L))
    pref_key    <- paste(round(pseudo_ref$x, 4L), round(pseudo_ref$y, 4L))
    prop_ta_ref <- round(mean(pref_key %in% ref_ta_key), 3L)
    
    med_ovrlp <- if (!is.null(metrics_df) && nrow(metrics_df) > 0L)
      round(median(metrics_df$overlap, na.rm = TRUE), 3L)
    else NA_real_
    
    annotation_txt <- paste0(
      "Overlap (median) = ", med_ovrlp,
      "\nProp. true abs. = ", prop_ta_ref
    )
    
    p_pa <- ggplot() +
      geom_contour_filled(data = kde_df,
                          aes(x = PC1, y = PC2, z = density),
                          breaks = c(0, 0.05, 0.25, 0.5, 0.75, 1),
                          alpha  = 0.80) +
      scale_fill_viridis_d(option = "mako", name = "Background\ndensity [0–1]",
                           direction = -1) +
      geom_point(data = pseudo_ref, aes(x = PC1, y = PC2),
                 colour = "darkorange", alpha = 0.45, size = 0.7) +
      geom_point(data = ref_pres, aes(x = PC1, y = PC2),
                 colour = "firebrick", alpha = 0.35, size = 0.7) +
      geom_contour(data = grid_pc, aes(x = PC1, y = PC2, z = suit),
                   colour = "firebrick", linewidth = 0.5, linetype = "dashed",
                   breaks = c(0.5, 0.95)) +
      annotate("text",
               x = min(dt$PC1), y = max(dt$PC2),
               label = annotation_txt,
               hjust = 0, vjust = 1, size = 3.5, colour = "white") +
      coord_equal() +
      labs(
        title    = paste0("Presences vs pseudo-absences  [", s_name, "]"),
        subtitle = "Red = presences | Blue = PA | Dashed = suit 0.50, 0.95",
        x = pc1_lab, y = pc2_lab
      ) +
      theme_classic(base_size = 13)
    
    # Marginal density bias plots (reference realisation)
    bias_df <- rbind(
      data.frame(PC1 = dt$PC1,         PC2 = dt$PC2,         group = "Background"),
      data.frame(PC1 = pseudo_ref$PC1, PC2 = pseudo_ref$PC2, group = "Pseudo-absences"),
      data.frame(PC1 = ref_pres$PC1,   PC2 = ref_pres$PC2,   group = "Presences")
    )
    bias_df$group <- factor(
      bias_df$group,
      levels = c("Background", "Pseudo-absences", "Presences")
    )
    col_vals  <- c(Background = "grey50", "Pseudo-absences" = "darkorange",
                   Presences  = "firebrick")
    fill_vals <- c(Background = "grey70", "Pseudo-absences" = "darkorange",
                   Presences  = "firebrick")
    
    p_bias_PC1 <- ggplot(bias_df, aes(x = PC1, colour = group, fill = group)) +
      geom_density(alpha = 0.25, linewidth = 0.8) +
      scale_colour_manual(values = col_vals) +
      scale_fill_manual(values   = fill_vals) +
      annotate("text", x = Inf, y = Inf,
               label = paste0("Coverage = ",
                              round(diff(range(pseudo_ref$PC1)) / bg_range_PC1, 3L)),
               hjust = 1.1, vjust = 1.5, size = 4) +
      labs(title  = paste0("Sampling bias — PC1  [", s_name, "]"),
           x = pc1_lab, y = "Density", colour = NULL, fill = NULL) +
      theme_classic(base_size = 13) +
      theme(legend.position = "bottom")
    
    p_bias_PC2 <- ggplot(bias_df, aes(x = PC2, colour = group, fill = group)) +
      geom_density(alpha = 0.25, linewidth = 0.8) +
      scale_colour_manual(values = col_vals) +
      scale_fill_manual(values   = fill_vals) +
      annotate("text", x = Inf, y = Inf,
               label = paste0("Coverage = ",
                              round(diff(range(pseudo_ref$PC2)) / bg_range_PC2, 3L)),
               hjust = 1.1, vjust = 1.5, size = 4) +
      labs(title  = paste0("Sampling bias — PC2  [", s_name, "]"),
           x = pc2_lab, y = "Density", colour = NULL, fill = NULL) +
      theme_classic(base_size = 13) +
      theme(legend.position = "bottom")
    
    # Three separate boxplots — one per metric group
    box_theme <- theme_classic(base_size = 13) +
      theme(legend.position = "none", axis.text.x = element_text(size = 11))
    sub_txt <- paste0(n_realizations, " Bernoulli realisations  [", s_name, "]")
    
    if (!is.null(metrics_df) && nrow(metrics_df) > 0L) {
      
      # -- Boxplot 1: hypervolume overlap -----------------------------------
      p_box_overlap <- ggplot(metrics_df,
                              aes(x = factor(0), y = overlap)) +
        geom_boxplot(fill = "#66C2A5", alpha = 0.7,
                     outlier.shape = 21, outlier.size = 1.5, width = 0.4) +
        geom_jitter(width = 0.08, alpha = 0.45, size = 0.9) +
        labs(title    = paste0("Hypervolume overlap — ", sub_txt),
             subtitle = "Intersection volume (pres. ∩ pseudo-abs.)",
             x = NULL, y = "Overlap volume") +
        scale_x_discrete(labels = NULL) +
        box_theme
      
      # -- Boxplot 2: PC range coverage -------------------------------------
      cov_long <- data.frame(
        value = c(metrics_df$rel_cov_PC1, metrics_df$rel_cov_PC2),
        axis  = rep(c("PC1", "PC2"), each = nrow(metrics_df))
      )
      p_box_coverage <- ggplot(cov_long, aes(x = axis, y = value, fill = axis)) +
        geom_boxplot(alpha = 0.7, outlier.shape = 21, outlier.size = 1.5,
                     width = 0.4) +
        geom_jitter(width = 0.08, alpha = 0.45, size = 0.9) +
        scale_fill_manual(values = c(PC1 = "#8DA0CB", PC2 = "#FC8D62")) +
        scale_y_continuous(limits = c(0, 1)) +
        labs(title    = paste0("Range coverage — ", sub_txt),
             subtitle = "Fraction of background PC range covered by pseudo-absences",
             x = NULL, y = "Relative coverage [0, 1]") +
        box_theme
      
      # -- Boxplot 3: proportion true absences ------------------------------
      p_box_trueabs <- ggplot(metrics_df,
                              aes(x = factor(0), y = prop_true_abs)) +
        geom_boxplot(fill = "#E78AC3", alpha = 0.7,
                     outlier.shape = 21, outlier.size = 1.5, width = 0.4) +
        geom_jitter(width = 0.08, alpha = 0.45, size = 0.9) +
        scale_y_continuous(limits = c(0, 1)) +
        labs(title    = paste0("True absences — ", sub_txt),
             subtitle = "Proportion of pseudo-absences on true-absence cells",
             x = NULL, y = "Proportion [0, 1]") +
        scale_x_discrete(labels = NULL) +
        box_theme
      
    } else {
      empty_p <- ggplot() +
        annotate("text", x = 0.5, y = 0.5,
                 label = "No completed realisations\n(all skipped due to errors)") +
        theme_void()
      p_box_overlap  <- empty_p
      p_box_coverage <- empty_p
      p_box_trueabs  <- empty_p
    }
    
    # Combined presence + PA dataset (reference realisation)
    pres_out      <- ref_pres[, c("x", "y", "PC1", "PC2", "suit")]
    pres_out$pa   <- 1L
    pseudo_out    <- pseudo_ref[, c("x", "y", "PC1", "PC2", "suit")]
    pseudo_out$pa <- 0L
    
    sampler_results[[s_name]] <- list(
      pseudo_ref  = pseudo_ref,
      pseudo_vect = terra::vect(pseudo_ref, geom = c("x", "y"), crs = crs(envData)),
      dataset     = rbind(pres_out, pseudo_out),
      metrics     = metrics_df,
      plots       = list(
        p_pa           = p_pa,
        p_bias_PC1     = p_bias_PC1,
        p_bias_PC2     = p_bias_PC2,
        p_box_overlap  = p_box_overlap,
        p_box_coverage = p_box_coverage,
        p_box_trueabs  = p_box_trueabs
      )
    )
  }  # end sampler loop
  
  # ── 8. RASTERS (reference realisation) ─────────────────────────────────────
  
  suit_rast <- terra::rast(
    dt[, c("x", "y", "suit")],
    type = "xyz", crs = crs(envData)
  )
  pa_rast <- terra::rast(
    data.frame(x = dt$x, y = dt$y, pa = as.numeric(dt$pa)),
    type = "xyz", crs = crs(envData)
  )
  
  # ── 9. RETURN ───────────────────────────────────────────────────────────────
  
  list(
    # Niche definition and equations
    niche = list(
      mu        = mu,
      sigma_PC1 = sigma_PC1,
      sigma_PC2 = sigma_PC2,
      rho       = rho,
      Sigma     = Sigma,
      equations = list(
        PC1   = eq_PC1,
        PC2   = eq_PC2,
        two_d = eq_2d,
        logit = eq_logit
      )
    ),
    # Suitability surface and Bernoulli presence layer (reference realisation)
    suit_rast = suit_rast,
    pa_rast   = pa_rast,
    # Marginal response curves and logit plots
    response_curves = list(
      data  = list(PC1 = rc_PC1, PC2 = rc_PC2),
      plots = list(
        p_rc_PC1        = p_rc_PC1,
        p_rc_PC2        = p_rc_PC2,
        p_rc_logit_PC1  = p_rc_logit_PC1,
        p_rc_logit_PC2  = p_rc_logit_PC2
      )
    ),
    # Per-sampler results — list of lists
    samplers   = sampler_results,
    # Shared diagnostic plot
    plots      = list(p_niche = p_niche),
    # Bandwidth used (store so you can re-use across species)
    bw         = bw,
    background = dt
  )
}