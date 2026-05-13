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


#' MCMC pseudo-absence sampler (USE.MCMC package)
#'
#' Samples pseudo-absences via a Markov chain whose stationary distribution
#' is the environmental GMM density minus the species-presence GMM density.
#' Internally calls USE.MCMC::paSamplingMcmc(); see that function for the
#' full algorithm and the C++/R engine dispatch.
#'
#' @param background   Full background data.frame (x, y, PC1, PC2, suit, ...)
#' @param N_pa         Number of pseudo-absences to return
#' @param pres         Presence rows from background (required)
#' @param seed         RNG seed; forwarded to paSamplingMcmc() as seed.number
#' @param env.rast     SpatRaster of the ORIGINAL environmental variables
#'                     (envData), NOT rpc$PCs. paSamplingMcmc() always runs
#'                     rastPCA() internally — feeding it PC rasters re-rotates
#'                     an already-orthogonal space and breaks the
#'                     environmental/species GMM filters, making the sampler
#'                     degenerate to random sampling. Pass envData.
#' @param chain.length MCMC chain length (default 10000)
#' @param burnIn       Robbins-Monro adaptation steps (default 1000)
#' @param num.chains   Parallel chains (default 1)
#' @param num.cores    Cores for multi-chain parallelism (default 1)
#' @param engine       "auto" (default), "R", or "cpp". Leave "auto" unless
#'                     forcing the reference loop for comparison runs.
#' @param species.cutoff.threshold
#'                     Percentile of the species-presence GMM density used to
#'                     define the region the chain may visit. Forwarded to
#'                     USE.MCMC::paSamplingMcmc(); package default is 0.95.
#' @param ...          Ignored additional arguments (interface compatibility,
#'                     swallows pa_uniform's grid.res / thres)
pa_mcmc <- function(background, N_pa, pres = NULL, seed = 123,
                    env.rast = NULL,
                    chain.length = 10000,
                    burnIn = 1000,
                    num.chains = 1, num.cores = 1,
                    engine = "auto",
                    species.cutoff.threshold = 0.95, ...) {

  if (!requireNamespace("USE.MCMC", quietly = TRUE))
    stop("Package 'USE.MCMC' is required for pa_mcmc. ",
         "Install with: devtools::install('../USE.MCMC')")
  if (!requireNamespace("sf", quietly = TRUE))
    stop("Package 'sf' is required for pa_mcmc.")

  if (is.null(pres))
    stop("pa_mcmc requires pres (presence rows from the background).")
  if (is.null(env.rast))
    stop("pa_mcmc requires 'env.rast' (the original environmental SpatRaster, ",
         "e.g. envData). Pass pa_env_rast = envData in virtualSpecies(). ",
         "Do NOT pass rpc$PCs: USE.MCMC runs rastPCA() internally and ",
         "rePCA-ing PC scores breaks the GMM filters.")

  set.seed(seed)

  pres_sf <- sf::st_as_sf(
    pres[, c("x", "y"), drop = FALSE],
    coords = c("x", "y"),
    crs    = terra::crs(env.rast)
  )

  # Over-request to absorb paSamplingMcmc's internal distance-threshold filter
  # and the dedup step on the first PC dimension. Subsampled to N_pa below.
  n_target <- ceiling(N_pa * 1.2)

  pa_result <- tryCatch(
    USE.MCMC::paSamplingMcmc(
      env.data.raster = env.rast,
      pres            = pres_sf,
      n.samples       = n_target,
      chain.length    = chain.length,
      burnIn          = burnIn,
      num.chains      = num.chains,
      num.cores       = num.cores,
      seed.number     = seed,
      engine          = engine,
      species.cutoff.threshold = species.cutoff.threshold,
      verbose         = FALSE,
      plot_proc       = FALSE
    ),
    error = function(e) {
      warning("USE.MCMC::paSamplingMcmc failed ('", conditionMessage(e),
              "'). Falling back to pa_random.")
      NULL
    }
  )

  if (is.null(pa_result) || nrow(pa_result) == 0L) {
    warning("pa_mcmc returned 0 pseudo-absences; falling back to pa_random.")
    return(pa_random(background, N_pa, pres, seed))
  }

  # paSamplingMcmc returns an sf whose geometry IS the geographic (x, y)
  # location (line 100 of paSamplingMcmc.R passes coords = c("x", "y") to
  # st_as_sf, so x/y live in the geometry column, NOT as attributes).
  # Pull them back out with st_coordinates() rather than st_drop_geometry().
  pa_coords <- sf::st_coordinates(pa_result)
  pa_x <- pa_coords[, 1L]
  pa_y <- pa_coords[, 2L]

  # Match back to background rows by geographic coordinates (4 dp = ~11 m).
  bg_key <- paste(round(background$x, 4L), round(background$y, 4L))
  pa_key <- paste(round(pa_x,         4L), round(pa_y,         4L))
  idx       <- match(pa_key, bg_key)
  idx_valid <- unique(idx[!is.na(idx)])
  pa_bg     <- background[idx_valid, ]

  if (length(idx_valid) == 0L)
    warning("pa_mcmc: NO geographic coordinates matched the background. ",
            "Check that env.rast covers the same extent as background and ",
            "that you passed envData (not PC rasters) as pa_env_rast.")
  if (length(idx_valid) < N_pa * 0.5)
    warning("pa_mcmc: only ", length(idx_valid), " of ", N_pa,
            " requested pseudo-absences matched background coordinates. ",
            "Remainder will be padded with random draws.")

  # Adjust to exactly N_pa — same pad-or-subsample pattern as pa_uniform.
  if (nrow(pa_bg) >= N_pa) {
    pa_bg <- pa_bg[sample(nrow(pa_bg), N_pa, replace = FALSE), ]
  } else {
    pres_key  <- paste(round(pres$x, 4L), round(pres$y, 4L))
    all_key   <- paste(round(background$x, 4L), round(background$y, 4L))
    non_pres  <- which(!all_key %in% pres_key)
    remaining <- setdiff(non_pres, idx_valid)
    extra_n   <- N_pa - nrow(pa_bg)
    if (length(remaining) >= extra_n) {
      pa_bg <- rbind(pa_bg, background[sample(remaining, extra_n), ])
    } else if (length(remaining) > 0L) {
      pa_bg <- rbind(pa_bg, background[remaining, ])
      warning("pa_mcmc: obtained ", nrow(pa_bg),
              " pseudo-absences instead of ", N_pa, ".")
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
#' @param parallel         If TRUE, dispatch (sampler, realisation) tasks over
#'                         a future::multisession worker pool. Each task is
#'                         independent and assigned dynamically, so workers
#'                         pick up new tasks as they free up. Requires
#'                         'future' and 'furrr'. Default FALSE.
#' @param n_workers        Number of parallel workers when parallel = TRUE.
#'                         NULL (default) -> parallel::detectCores() - 1, so
#'                         one core stays free for system responsiveness.
#'                         When parallel = TRUE, pa_mcmc's num.cores is forced
#'                         to 1 inside each worker to avoid oversubscription.
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
    parallel         = FALSE,
    n_workers        = NULL,
    ...
) {

  # Default the PA environment raster to envData if not supplied
  if (is.null(pa_env_rast)) pa_env_rast <- envData

  # ── 0. PARALLEL BACKEND ────────────────────────────────────────────────────
  #
  # plan(multisession) so the same code path works on Linux, macOS and Windows
  # (multicore/fork is silently downgraded on Windows and inside RStudio).
  # on.exit restores the previous plan even if a sampler errors out.
  if (parallel) {
    if (!requireNamespace("future", quietly = TRUE) ||
        !requireNamespace("furrr",  quietly = TRUE))
      stop("parallel = TRUE requires packages 'future' and 'furrr'. ",
           "Install with: install.packages(c('future', 'furrr'))")
    if (is.null(n_workers))
      #n_workers <- max(1L, parallel::detectCores() - 1L)
      n_workers <- 4
    old_plan <- future::plan(future::multisession, workers = n_workers)
    on.exit(future::plan(old_plan), add = TRUE)
    if (verbose) cat("Parallel mode: ", n_workers, " workers\n", sep = "")
  }

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

  # ── Reference pseudo-absences (realisation 1) for diagnostic plots ─────────
  # One call per sampler, run sequentially before dispatch — cheap, and the
  # sequential path lets a failure on the reference realisation halt early.
  N_pa_ref    <- max(1L, round(nrow(ref_pres) / bgk_prev))
  dot_args    <- list(...)
  pseudo_refs <- lapply(names(pa_samplers), function(s_name) {
    if (verbose) cat("\n══ Sampler:", s_name, "══\n")
    sampler <- pa_samplers[[s_name]]
    tryCatch(
      do.call(sampler, c(
        list(background = dt, N_pa = N_pa_ref, pres = ref_pres,
             seed = seed_pseudo_base, env.rast = pa_env_rast),
        dot_args
      )),
      error = function(e) {
        stop("Sampler '", s_name, "' failed on reference realisation: ", e$message)
      }
    )
  })
  names(pseudo_refs) <- names(pa_samplers)

  # Pack the SpatRaster for multisession transport. SpatRaster is an Rcpp
  # external pointer that does NOT survive R's serialize/unserialize, so
  # without wrap() each worker would receive a SpatRaster with a null C++
  # pointer — terra::crs(env.rast) inside pa_uniform / pa_mcmc would error,
  # the outer tryCatch would swallow it, and every parallel realisation
  # would be mis-tagged "skip_sampler_null_or_short".
  pa_env_rast_packed <- if (parallel) terra::wrap(pa_env_rast) else NULL

  # Skip-row constructor: every realisation returns a row of the same shape so
  # downstream rbind never sees mismatched columns. status = "ok" rows carry
  # real metrics; status = "skip_*" rows carry NAs and the reason the
  # realisation was abandoned. Letting empty samplers produce 0-row metric
  # tables (or NULL) is what triggered the rbind error in the testing wrapper.
  .skip_row <- function(r, status,
                        N_pres_total = NA_integer_,
                        N_pres_used  = NA_integer_) {
    data.frame(
      realization   = r,
      status        = status,
      N_pres_total  = N_pres_total,
      N_pres_used   = N_pres_used,
      N_pseudo      = NA_integer_,
      overlap       = NA_real_,
      rel_cov_PC1   = NA_real_,
      rel_cov_PC2   = NA_real_,
      prop_true_abs = NA_real_,
      stringsAsFactors = FALSE
    )
  }

  # ── Closure: evaluate one (sampler, realisation) task ──────────────────────
  #
  # Always returns a one-row data.frame with the same columns (see .skip_row).
  # status == "ok" rows have real metrics; status == "skip_*" rows have NAs
  # and a reason code, so aggregations never see NULL. The closure captures
  # dt, pa_matrix, bw, pa_env_rast, bgk_prev, max_pres, seeds, bg ranges,
  # pa_samplers and dot_args. RNG: each task calls set.seed(seed_base + r)
  # explicitly, so output is deterministic regardless of execution order
  # (furrr's L'Ecuyer seed is overridden by these calls).
  eval_realization <- function(s_name, r) {
    # Pin RNG kind so results match between sequential and parallel runs.
    # furrr_options(seed = TRUE) switches workers to L'Ecuyer-CMRG; without
    # this reset, set.seed(seed_base + r) would seed a different RNG inside
    # workers than in the main session and produce different sample() output.
    RNGkind("Mersenne-Twister")

    sampler    <- pa_samplers[[s_name]]
    pa_r       <- pa_matrix[, r]
    pres_r_all <- dt[pa_r == 1L, ]   # full draw — used for prop_true_abs
    true_abs_r <- dt[pa_r == 0L, ]

    # Subsample presences to max_pres to cap hypervolume runtime. N_pa is
    # derived from the subsampled count so PA:presence ratio (bgk_prev)
    # stays consistent across realisations.
    set.seed(seed_base + r)
    pres_r <- if (nrow(pres_r_all) > max_pres)
      pres_r_all[sample(nrow(pres_r_all), max_pres), ]
    else pres_r_all

    N_pa_r <- max(1L, round(nrow(pres_r) / bgk_prev))
    if (nrow(pres_r) < 5L || N_pa_r < 5L)
      return(.skip_row(r, "skip_few_pres",
                       N_pres_total = nrow(pres_r_all),
                       N_pres_used  = nrow(pres_r)))

    # Under outer parallelism, force pa_mcmc to single-threaded so N workers
    # don't each spawn M chains and oversubscribe the CPU.
    task_dot_args <- dot_args
    if (parallel && s_name == "mcmc") task_dot_args$num.cores <- 1L

    # Unwrap the packed SpatRaster to a live one inside the worker. In
    # sequential mode pa_env_rast_packed is NULL and we use pa_env_rast as-is.
    env_rast_for_task <- if (parallel) terra::unwrap(pa_env_rast_packed) else pa_env_rast

    pseudo_r <- tryCatch(
      do.call(sampler, c(
        list(background = dt, N_pa = N_pa_r, pres = pres_r,
             seed = seed_pseudo_base + r, env.rast = env_rast_for_task),
        task_dot_args
      )),
      error = function(e) NULL
    )
    if (is.null(pseudo_r) || nrow(pseudo_r) < 5L)
      return(.skip_row(r, "skip_sampler_null_or_short",
                       N_pres_total = nrow(pres_r_all),
                       N_pres_used  = nrow(pres_r)))

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
    if (is.null(hyp_pres_r) || is.null(hyp_pa_r))
      return(.skip_row(r, "skip_hv_fail",
                       N_pres_total = nrow(pres_r_all),
                       N_pres_used  = nrow(pres_r)))

    hv_set <- hypervolume::hypervolume_set(
      hyp_pres_r, hyp_pa_r, check.memory = FALSE, verbose = FALSE
    )
    # get_volume() returns: [[1]] pres vol, [[2]] pa vol, [[3]] intersection
    ovrlp <- hypervolume::get_volume(hv_set)[[3L]]

    rel_cov_PC1 <- diff(range(pseudo_r$PC1)) / bg_range_PC1
    rel_cov_PC2 <- diff(range(pseudo_r$PC2)) / bg_range_PC2

    # Proportion of pseudo-absences on TRUE-absence cells. Flags samplers
    # that draw from true-presence cells, inflating observed class overlap.
    ta_key  <- paste(round(true_abs_r$x, 4L), round(true_abs_r$y, 4L))
    ps_key  <- paste(round(pseudo_r$x,   4L), round(pseudo_r$y,   4L))
    prop_ta <- mean(ps_key %in% ta_key)

    metrics_row <- data.frame(
      realization   = r,
      status        = "ok",
      N_pres_total  = nrow(pres_r_all),
      N_pres_used   = nrow(pres_r),
      N_pseudo      = nrow(pseudo_r),
      overlap       = round(ovrlp, 6L),
      rel_cov_PC1   = round(rel_cov_PC1, 4L),
      rel_cov_PC2   = round(rel_cov_PC2, 4L),
      prop_true_abs = round(prop_ta, 4L),
      stringsAsFactors = FALSE
    )

    # Force release of hypervolume / KDE working memory before returning.
    # Persistent multisession workers run several tasks back-to-back and R's
    # GC is lazy — without this, HV objects (100-500 MB each) drift up the
    # high-water mark per worker until RStudio OOMs. metrics_row holds no
    # references to anything below.
    rm(hyp_pres_r, hyp_pa_r, hv_set,
       pseudo_r, pres_r, pres_r_all, true_abs_r,
       ta_key, ps_key)
    gc(verbose = FALSE)

    metrics_row
  }

  # ── Flat task list: every (sampler, realisation) pair ──────────────────────
  tasks <- expand.grid(
    s_name = names(pa_samplers),
    r      = seq_len(n_realizations),
    stringsAsFactors = FALSE,
    KEEP.OUT.ATTRS   = FALSE
  )

  if (verbose)
    cat("\nEvaluating ", nrow(tasks), " tasks (",
        length(pa_samplers), " samplers x ", n_realizations,
        " realisations)…\n", sep = "")

  if (parallel) {
    results <- furrr::future_pmap(
      tasks, eval_realization,
      .options  = furrr::furrr_options(seed = TRUE, globals = TRUE),
      .progress = verbose
    )
  } else {
    results <- Map(eval_realization, tasks$s_name, tasks$r)
  }

  # Group results by sampler in the original sampler order
  results_by_sampler <- split(results, factor(tasks$s_name,
                                              levels = names(pa_samplers)))

  for (s_name in names(pa_samplers)) {

    pseudo_ref   <- pseudo_refs[[s_name]]
    metrics_list <- results_by_sampler[[s_name]]
    # Every task now returns a one-row data.frame (status = "ok" or "skip_*"),
    # so no Filter() / NULL handling is needed. rbind is shape-safe.
    metrics_df <- do.call(rbind, metrics_list)
    if (verbose) {
      tally <- table(factor(metrics_df$status,
                            levels = c("ok", "skip_few_pres",
                                       "skip_sampler_null_or_short",
                                       "skip_hv_fail")))
      cat(sprintf("  [%s] %s\n", s_name,
                  paste(sprintf("%s=%d", names(tally), as.integer(tally)),
                        collapse = "  ")))
    }
    # Plots / median etc. operate only on successful realisations.
    metrics_ok <- metrics_df[metrics_df$status == "ok", , drop = FALSE]

    # -- Plots for this sampler -----------------------------------------------

    # Proportion of reference PAs on true-absence cells (ref. realisation)
    ref_ta_key  <- paste(round(dt[ref_pa == 0L, "x"], 4L),
                         round(dt[ref_pa == 0L, "y"], 4L))
    pref_key    <- paste(round(pseudo_ref$x, 4L), round(pseudo_ref$y, 4L))
    prop_ta_ref <- round(mean(pref_key %in% ref_ta_key), 3L)

    med_ovrlp <- if (nrow(metrics_ok) > 0L)
      round(median(metrics_ok$overlap, na.rm = TRUE), 3L)
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

    if (nrow(metrics_ok) > 0L) {

      # -- Boxplot 1: hypervolume overlap -----------------------------------
      p_box_overlap <- ggplot(metrics_ok,
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
        value = c(metrics_ok$rel_cov_PC1, metrics_ok$rel_cov_PC2),
        axis  = rep(c("PC1", "PC2"), each = nrow(metrics_ok))
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
      p_box_trueabs <- ggplot(metrics_ok,
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

