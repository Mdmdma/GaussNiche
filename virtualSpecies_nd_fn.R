# virtualSpecies_nd_fn.R  —  k-dimensional GaussNiche (5d-niche branch)
# =============================================================================
# A dimension-general (k-D) sibling of virtualSpecies_fn.R. Same pipeline —
# multivariate-Gaussian niche on PCA axes -> Bernoulli presences -> pseudo-absence
# sampling -> hypervolume / coverage diagnostics -> geographic back-projection —
# but in an arbitrary number of PC dimensions instead of exactly two.
#
# What changes vs the 2-D module:
#   * mu is a length-k vector; niche breadth is the length-k vector `sigma`;
#     Sigma = diag(sigma) %*% Cor %*% diag(sigma) (Cor from scalar/`matrix` rho).
#   * suitability, Bernoulli draws, the fixed background bandwidth and
#     hypervolume_gaussian() all operate on the k PC columns.
#   * response curves are produced for every axis; the 2-D niche / PA scatter
#     plots are drawn on the PC1 x PC2 PROJECTION (suitability sliced at mu on the
#     remaining axes), since a k-D surface cannot be drawn directly.
#   * coverage is reported per axis (rel_cov_PC1 .. rel_cov_PCk).
#
# Samplers: pa_random (k-D trivially), pa_mcmc (k-D via USE.MCMC::paSamplingMcmc) and
# pa_nn (k-D via USE.MCMC::paSamplingNn). All three are dimension-agnostic — pass
# `dimensions = c("PC1",...,"PCk")`. The nearest-neighbour "uniform" sampler used to
# be 2-D by construction, but USE.MCMC::paSamplingNn is now dimension-general
# (uniform-box proposals + NN-remap, with a sqrt(d/2) correction on the support
# distance threshold); its acceptance rate falls with dimension, so it is compute-
# bound in high d (raise `n.candidates` / `n.tr`). Only the KDE backend
# USE.MCMC::paSampling remains 2-D. pa_nn is provided but not in the default sampler
# list; add it via `pa_samplers = list(random = pa_random, mcmc = pa_mcmc, nn = pa_nn)`.
#
# Source this file, then call virtualSpecies_nd().
# =============================================================================

suppressPackageStartupMessages({
  library(mvtnorm)
  library(MASS)
  library(ggplot2)
  library(hypervolume)
  library(patchwork)
  library(terra)
  library(USE.MCMC)
})

# =============================================================================
# PSEUDO-ABSENCE SAMPLERS  (shared interface)
#   function(background, N_pa, pres = NULL, seed = 123, env.rast = NULL, ...)
#   returns a data.frame with the same columns as `background`.
# =============================================================================

#' Random pseudo-absence sampler — draws N_pa rows uniformly from the background.
pa_random <- function(background, N_pa, pres = NULL, seed = 123, ...) {
  set.seed(seed)
  background[sample(nrow(background), size = N_pa, replace = FALSE), ]
}

#' MCMC pseudo-absence sampler (k-D) via USE.MCMC::paSamplingMcmc.
#'
#' paSamplingMcmc / mcmcSampling are dimension-agnostic; the only 2-D-ness in the
#' package is the default `dimensions = c("PC1","PC2")`. Pass the full PC name
#' vector to sample in k-D. Results are matched back to the background by
#' geographic coordinates (rounded to 4 dp), which is dimension-independent.
#'
#' @param env.rast  ORIGINAL environmental SpatRaster (e.g. the 12-layer
#'                  lean_natural stack), NOT the PC rasters — paSamplingMcmc runs
#'                  rastPCA() internally.
#' @param dimensions character vector of PC columns to sample in (default PC1..PC5).
#' @param precomputed.env optional bundle from USE.MCMC::precomputeMcmcEnvironment().
#'                  The PCA + environmental GMM fit are constant across species on
#'                  the same env.rast, so compute the bundle once and pass it here
#'                  to skip that work on every call. Must be built for the same
#'                  `dimensions` and `seed`.
pa_mcmc <- function(background, N_pa, pres = NULL, seed = 123,
                    env.rast = NULL,
                    dimensions = c("PC1", "PC2", "PC3", "PC4", "PC5"),
                    chain.length = 10000, burnIn = 1000,
                    num.chains = 1, num.cores = 1, engine = "auto",
                    species.cutoff.threshold = 0.95,
                    precomputed.env = NULL, ...) {
  if (!requireNamespace("USE.MCMC", quietly = TRUE))
    stop("Package 'USE.MCMC' is required for pa_mcmc.")
  if (!requireNamespace("sf", quietly = TRUE))
    stop("Package 'sf' is required for pa_mcmc.")
  if (is.null(pres))     stop("pa_mcmc requires pres (presence rows).")
  if (is.null(env.rast)) stop("pa_mcmc requires 'env.rast' (the original env SpatRaster, NOT PC rasters).")

  set.seed(seed)
  pres_sf <- sf::st_as_sf(pres[, c("x", "y"), drop = FALSE],
                          coords = c("x", "y"), crs = terra::crs(env.rast))
  n_target <- ceiling(N_pa * 1.2)

  pa_result <- tryCatch(
    USE.MCMC::paSamplingMcmc(
      env.data.raster = env.rast, pres = pres_sf,
      n.samples = n_target, chain.length = chain.length, burnIn = burnIn,
      dimensions = dimensions,
      num.chains = num.chains, num.cores = num.cores,
      seed.number = seed, engine = engine,
      species.cutoff.threshold = species.cutoff.threshold,
      precomputed.env = precomputed.env,
      verbose = FALSE, plot_proc = FALSE),
    error = function(e) {
      # A misconfiguration (e.g. an incompatible precomputed.env) must fail
      # loudly, not be silently downgraded to random sampling.
      if (inherits(e, "USE.MCMC_config_error")) stop(e)
      warning("paSamplingMcmc failed ('", conditionMessage(e), "'); falling back to pa_random.")
      NULL
    })

  if (is.null(pa_result) || nrow(pa_result) == 0L) {
    warning("pa_mcmc returned 0 pseudo-absences; falling back to pa_random.")
    return(pa_random(background, N_pa, pres, seed))
  }

  # geometry holds the geographic (x, y) location
  pa_coords <- sf::st_coordinates(pa_result)
  bg_key <- paste(round(background$x, 4L), round(background$y, 4L))
  pa_key <- paste(round(pa_coords[, 1L], 4L), round(pa_coords[, 2L], 4L))
  idx_valid <- unique(match(pa_key, bg_key))
  idx_valid <- idx_valid[!is.na(idx_valid)]
  pa_bg <- background[idx_valid, ]

  if (length(idx_valid) == 0L)
    warning("pa_mcmc: NO coordinates matched the background. Check env.rast extent.")

  # adjust to exactly N_pa (subsample or pad with random non-presence draws)
  if (nrow(pa_bg) >= N_pa) {
    pa_bg <- pa_bg[sample(nrow(pa_bg), N_pa, replace = FALSE), ]
  } else {
    pres_key <- paste(round(pres$x, 4L), round(pres$y, 4L))
    all_key  <- paste(round(background$x, 4L), round(background$y, 4L))
    remaining <- setdiff(which(!all_key %in% pres_key), idx_valid)
    extra_n   <- N_pa - nrow(pa_bg)
    if (length(remaining) >= extra_n) {
      pa_bg <- rbind(pa_bg, background[sample(remaining, extra_n), ])
    } else if (length(remaining) > 0L) {
      pa_bg <- rbind(pa_bg, background[remaining, ])
      warning("pa_mcmc: obtained ", nrow(pa_bg), " pseudo-absences instead of ", N_pa, ".")
    }
  }
  row.names(pa_bg) <- NULL
  pa_bg
}

#' Nearest-neighbour pseudo-absence sampler (k-D) via USE.MCMC::paSamplingNn.
#'
#' paSamplingNn is now dimension-general: it draws uniform proposals over the PC
#' bounding box, remaps each to its nearest realized environmental point (which is
#' uniform-over-support in any dimension), and rejects proposals whose remap distance
#' exceeds a sqrt(d/2)-corrected support threshold. Like pa_mcmc, results are matched
#' back to the background by geographic (x, y) rounded to 4 dp.
#'
#' Because the realized environment fills an exponentially smaller share of its
#' bounding box as the number of `dimensions` grows, the acceptance rate falls with
#' dimension: raise `n.candidates` (raw proposals per batch) and/or `n.tr` for k > ~4.
#'
#' @param env.rast ORIGINAL environmental SpatRaster (NOT the PC rasters);
#'                 paSamplingNn runs rastPCA() internally.
#' @param dimensions character vector of PC columns to sample in (default PC1..PC5).
#' @param n.candidates raw uniform proposals per batch (NULL -> grid.res^2 * n.tr).
#' @param dim.correction support-threshold dimension correction ("voronoi"/"simplex"/
#'                 "none" or a positive multiplier).
pa_nn <- function(background, N_pa, pres = NULL, seed = 123,
                  env.rast = NULL,
                  dimensions = c("PC1", "PC2", "PC3", "PC4", "PC5"),
                  grid.res = 25, n.tr = 5, thres = 0.75,
                  nn.based.presence.exclusion = TRUE,
                  data.based.distance.threshold = TRUE,
                  n.candidates = NULL, dim.correction = "voronoi",
                  precomputed.pca = NULL, ...) {
  if (!requireNamespace("USE.MCMC", quietly = TRUE))
    stop("Package 'USE.MCMC' is required for pa_nn.")
  if (!requireNamespace("sf", quietly = TRUE))
    stop("Package 'sf' is required for pa_nn.")
  if (is.null(pres))     stop("pa_nn requires pres (presence rows).")
  if (is.null(env.rast)) stop("pa_nn requires 'env.rast' (the original env SpatRaster, NOT PC rasters).")

  set.seed(seed)
  pres_sf <- sf::st_as_sf(pres[, c("x", "y"), drop = FALSE],
                          coords = c("x", "y"), crs = terra::crs(env.rast))
  n_target <- ceiling(N_pa * 1.2)

  pa_result <- tryCatch(
    USE.MCMC::paSamplingNn(
      env.rast = env.rast, pres = pres_sf,
      n.samples = n_target, dimensions = dimensions,
      grid.res = grid.res, n.tr = n.tr, thres = thres,
      nn.based.presence.exclusion = nn.based.presence.exclusion,
      data.based.distance.threshold = data.based.distance.threshold,
      n.candidates = n.candidates, dim.correction = dim.correction,
      precomputed.pca = precomputed.pca,
      verbose = FALSE, plot_proc = FALSE),
    error = function(e) {
      warning("paSamplingNn failed ('", conditionMessage(e), "'); falling back to pa_random.")
      NULL
    })

  if (is.null(pa_result) || nrow(pa_result) == 0L) {
    warning("pa_nn returned 0 pseudo-absences; falling back to pa_random.")
    return(pa_random(background, N_pa, pres, seed))
  }

  # geometry holds the geographic (x, y) location
  pa_coords <- sf::st_coordinates(pa_result)
  bg_key <- paste(round(background$x, 4L), round(background$y, 4L))
  pa_key <- paste(round(pa_coords[, 1L], 4L), round(pa_coords[, 2L], 4L))
  idx_valid <- unique(match(pa_key, bg_key))
  idx_valid <- idx_valid[!is.na(idx_valid)]
  pa_bg <- background[idx_valid, ]

  if (length(idx_valid) == 0L)
    warning("pa_nn: NO coordinates matched the background. Check env.rast extent.")

  # adjust to exactly N_pa (subsample or pad with random non-presence draws)
  if (nrow(pa_bg) >= N_pa) {
    pa_bg <- pa_bg[sample(nrow(pa_bg), N_pa, replace = FALSE), ]
  } else {
    pres_key <- paste(round(pres$x, 4L), round(pres$y, 4L))
    all_key  <- paste(round(background$x, 4L), round(background$y, 4L))
    remaining <- setdiff(which(!all_key %in% pres_key), idx_valid)
    extra_n   <- N_pa - nrow(pa_bg)
    if (length(remaining) >= extra_n) {
      pa_bg <- rbind(pa_bg, background[sample(remaining, extra_n), ])
    } else if (length(remaining) > 0L) {
      pa_bg <- rbind(pa_bg, background[remaining, ])
      warning("pa_nn: obtained ", nrow(pa_bg), " pseudo-absences instead of ", N_pa, ".")
    }
  }
  row.names(pa_bg) <- NULL
  pa_bg
}

# =============================================================================
# HELPERS
# =============================================================================

#' Detect PC columns in a background data.frame (PC1, PC2, ... in numeric order).
.detect_pcs <- function(dt) {
  pcs <- grep("^PC[0-9]+$", names(dt), value = TRUE)
  if (length(pcs) < 2L)
    stop("Need >=2 columns named PC1, PC2, ... in `dt`; found: ",
         paste(pcs, collapse = ", "))
  pcs[order(as.integer(sub("^PC", "", pcs)))]
}

#' One-time hypervolume bandwidth from the full background E-space (length-k).
#' Call once and pass as bw = ... so the same bandwidth is used for every species
#' and sampler — essential for comparability of hypervolume estimates.
compute_bandwidth_nd <- function(dt, pc_cols = NULL) {
  if (is.null(pc_cols)) pc_cols <- .detect_pcs(dt)
  hypervolume::estimate_bandwidth(dt[, pc_cols])
}

# Build a k x k covariance matrix from a sigma vector and a correlation spec.
.make_sigma <- function(sigma, rho) {
  k <- length(sigma)
  Cor <- if (is.matrix(rho)) rho else { m <- matrix(rho, k, k); diag(m) <- 1; m }
  D <- diag(sigma, nrow = k)
  D %*% Cor %*% D
}

# =============================================================================
# MAIN FUNCTION  (k-D)
# =============================================================================

#' Simulate a k-dimensional virtual species and evaluate pseudo-absence samplers.
#'
#' @param dt          Background data.frame with x, y and PC columns (PC1..PCk).
#' @param envData     Original SpatRaster (for CRS and raster back-projection).
#' @param mu          Numeric[k]: niche optimum in PC space.
#' @param sigma       Numeric[k]: niche breadth per axis.
#' @param rho         Scalar equicorrelation in [-1,1] OR a k x k correlation
#'                    matrix. Default 0 (diagonal Sigma).
#' @param pc_cols     PC column names to use (default: all PC* in dt, ordered).
#' @param bgk_prev    Presence:pseudo-absence ratio (1 = 1:1).
#' @param pa_samplers Named list of samplers. Default random + mcmc.
#' @param n_realizations Number of Bernoulli draws.
#' @param max_pres    Cap on presences used for hypervolume + PA sampling per
#'                    realisation (default 300; k-D hypervolumes are heavier).
#' @param seed_base / seed_pseudo_base  Layered RNG seeds (realisation r uses +r).
#' @param pc_labs     Optional character[k] of axis labels (default "PCk (xx% var)").
#' @param bw          Pre-computed length-k bandwidth (compute_bandwidth_nd()).
#' @param pa_env_rast SpatRaster forwarded to samplers as env.rast (ORIGINAL vars).
#' @param hv_quantile Probability quantile for hypervolume_gaussian (default 0.95).
#' @param verbose / parallel / n_workers  As in the 2-D module.
#' @param ...         Extra args forwarded to samplers (e.g. dimensions=,
#'                    chain.length=, burnIn= for pa_mcmc).
#' @return Named list: niche, parameters, suit_rast, pa_rast, response_curves,
#'         samplers, plots, bw, background.
virtualSpecies_nd <- function(
    dt, envData, mu, sigma,
    rho              = 0,
    pc_cols          = NULL,
    bgk_prev         = 1,
    pa_samplers      = list(random = pa_random, mcmc = pa_mcmc),
    n_realizations   = 20L,
    max_pres         = 300L,
    seed_base        = 42L,
    seed_pseudo_base = 123L,
    pc_labs          = NULL,
    bw               = NULL,
    pa_env_rast      = NULL,
    hv_quantile      = 0.95,
    verbose          = TRUE,
    parallel         = FALSE,
    n_workers        = NULL,
    ...
) {
  if (is.null(pa_env_rast)) pa_env_rast <- envData
  if (is.null(pc_cols))     pc_cols <- .detect_pcs(dt)
  k <- length(pc_cols)
  stopifnot(length(mu) == k, length(sigma) == k)
  if (is.null(pc_labs)) pc_labs <- pc_cols
  names(pc_labs) <- pc_cols

  # ---- 0. parallel backend ------------------------------------------------
  if (parallel) {
    if (!requireNamespace("future", quietly = TRUE) ||
        !requireNamespace("furrr", quietly = TRUE))
      stop("parallel = TRUE requires 'future' and 'furrr'.")
    if (is.null(n_workers)) n_workers <- 4L
    old_plan <- future::plan(future::multisession, workers = n_workers)
    on.exit(future::plan(old_plan), add = TRUE)
    if (verbose) cat("Parallel mode:", n_workers, "workers\n")
  }

  # ---- 1. k-D Gaussian niche ----------------------------------------------
  Sigma    <- .make_sigma(sigma, rho)
  peak     <- mvtnorm::dmvnorm(matrix(mu, nrow = 1L), mean = mu, sigma = Sigma)
  Xbg      <- as.matrix(dt[, pc_cols])
  suit     <- mvtnorm::dmvnorm(Xbg, mean = mu, sigma = Sigma) / peak
  dt$suit  <- suit
  if (verbose) cat("Suitability range:", round(range(suit), 4L),
                   "| k =", k, "dims\n")

  marg_sd <- sqrt(diag(Sigma))   # marginal sd per axis (exact marginal of MVN)

  # ---- 2. response curves (one per axis) ----------------------------------
  rc_data <- lapply(seq_len(k), function(i) {
    xs <- seq(min(dt[[pc_cols[i]]]), max(dt[[pc_cols[i]]]), length.out = 400L)
    data.frame(x = xs, suit = exp(-0.5 * (xs - mu[i])^2 / marg_sd[i]^2))
  })
  names(rc_data) <- pc_cols
  rc_plots <- lapply(seq_len(k), function(i) {
    ggplot(rc_data[[i]], aes(x = x, y = suit)) +
      geom_line(colour = "firebrick", linewidth = 1) +
      geom_rug(data = dt, aes(x = .data[[pc_cols[i]]]), inherit.aes = FALSE,
               alpha = 0.04, colour = "grey40", length = unit(0.02, "npc")) +
      labs(title = paste0("Marginal response — ", pc_cols[i]),
           subtitle = sprintf("mu=%.2f, sigma=%.2f", mu[i], sigma[i]),
           x = pc_labs[[pc_cols[i]]], y = "Suitability [0,1]") +
      theme_classic(base_size = 12)
  })
  names(rc_plots) <- paste0("p_rc_", pc_cols)

  # logit response curves per axis: logit(suit) = log(suit/(1-suit)), finite only
  rc_logit_plots <- lapply(seq_len(k), function(i) {
    d <- rc_data[[i]]; d$logit <- log(d$suit / (1 - d$suit))
    d <- d[is.finite(d$logit), , drop = FALSE]
    ggplot(d, aes(x = x, y = logit)) +
      geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
      geom_line(colour = "steelblue", linewidth = 1) +
      labs(title = paste0("Logit(suit) — ", pc_cols[i]),
           x = pc_labs[[pc_cols[i]]], y = "logit(suit)") +
      theme_classic(base_size = 12)
  })
  names(rc_logit_plots) <- paste0("p_rc_logit_", pc_cols)

  # equation string (generic quadratic form)
  eq_nd <- sprintf("log(suit) = -0.5 * (x-mu)^T Sigma^-1 (x-mu);  k=%d, mu=(%s), sigma=(%s)",
                   k, paste(round(mu, 2), collapse = ", "),
                   paste(round(sigma, 2), collapse = ", "))
  if (verbose) cat(eq_nd, "\n")

  # ---- 3. fixed background bandwidth (k-D) --------------------------------
  bw_source <- if (is.null(bw)) "computed" else "user-supplied"
  if (is.null(bw)) {
    if (verbose) cat("Estimating k-D hypervolume bandwidth from background...\n")
    bw <- hypervolume::estimate_bandwidth(dt[, pc_cols])
  }
  if (verbose) cat("Bandwidth:", paste(round(bw, 4L), collapse = ", "), "\n")

  # ---- 4. Bernoulli realisations ------------------------------------------
  if (verbose) cat("Drawing", n_realizations, "Bernoulli realisations...\n")
  pa_matrix <- vapply(seq_len(n_realizations), function(r) {
    set.seed(seed_base + r); rbinom(nrow(dt), 1L, suit)
  }, integer(nrow(dt)))

  ref_pa       <- pa_matrix[, 1L]
  dt$pa        <- ref_pa
  ref_pres_all <- dt[ref_pa == 1L, ]
  set.seed(seed_base)
  ref_pres <- if (nrow(ref_pres_all) > max_pres)
    ref_pres_all[sample(nrow(ref_pres_all), max_pres), ] else ref_pres_all
  if (verbose)
    cat("Prevalence (ref):", round(mean(ref_pa), 3L),
        "| N pres:", nrow(ref_pres_all), "-> used", nrow(ref_pres), "\n")

  bg_range <- vapply(pc_cols, function(pc) diff(range(dt[[pc]])), numeric(1))

  # ---- 5. PC1 x PC2 projection plots --------------------------------------
  px <- pc_cols[1]; py <- pc_cols[2]
  kde_raw <- MASS::kde2d(dt[[px]], dt[[py]], n = 150L)
  kde_df  <- expand.grid(X = kde_raw$x, Y = kde_raw$y)
  kde_df$density <- as.vector(kde_raw$z); kde_df$density <- kde_df$density / max(kde_df$density)

  # suitability sliced at mu on the non-displayed axes
  grid_xy <- expand.grid(
    X = seq(min(dt[[px]]), max(dt[[px]]), length.out = 200L),
    Y = seq(min(dt[[py]]), max(dt[[py]]), length.out = 200L))
  slice <- matrix(rep(mu, each = nrow(grid_xy)), nrow = nrow(grid_xy))
  colnames(slice) <- pc_cols
  slice[, px] <- grid_xy$X; slice[, py] <- grid_xy$Y
  grid_xy$suit <- mvtnorm::dmvnorm(slice[, pc_cols], mean = mu, sigma = Sigma) / peak

  p_niche <- ggplot() +
    geom_contour_filled(data = kde_df, aes(X, Y, z = density),
                        breaks = c(0, 0.05, 0.25, 0.5, 0.75, 1), alpha = 0.8) +
    scale_fill_viridis_d(option = "mako", name = "Background\ndensity", direction = -1) +
    geom_contour(data = grid_xy, aes(X, Y, z = suit), colour = "firebrick",
                 linewidth = 0.6, breaks = c(0.1, 0.25, 0.5, 0.75, 0.95)) +
    geom_point(data = ref_pres, aes(.data[[px]], .data[[py]]),
               colour = "grey30", alpha = 0.25, size = 0.5) +
    annotate("point", x = mu[1], y = mu[2], colour = "firebrick", shape = 3,
             size = 4, stroke = 2) +
    coord_equal() +
    labs(title = "Niche vs available environment (PC1 x PC2 projection)",
         subtitle = paste0("Suitability sliced at mu on PC3..PC", k,
                           "  |  + = optimum  |  dots = presences"),
         x = pc_labs[[px]], y = pc_labs[[py]]) +
    theme_classic(base_size = 12)

  # ---- 6. metrics-row builder (dynamic per-PC coverage columns) ------------
  relcov_names <- paste0("rel_cov_", pc_cols)
  make_row <- function(r, status, N_pres_total = NA_integer_, N_pres_used = NA_integer_,
                       N_pseudo = NA_integer_, overlap = NA_real_, prop_true_abs = NA_real_,
                       relcov = setNames(rep(NA_real_, k), relcov_names)) {
    base <- data.frame(realization = r, status = status,
                       N_pres_total = N_pres_total, N_pres_used = N_pres_used,
                       N_pseudo = N_pseudo, overlap = overlap, stringsAsFactors = FALSE)
    rc <- as.data.frame(as.list(relcov), stringsAsFactors = FALSE); names(rc) <- relcov_names
    cbind(base, rc, prop_true_abs = prop_true_abs)
  }

  # ---- 7. reference pseudo-absences per sampler ---------------------------
  N_pa_ref <- max(1L, round(nrow(ref_pres) / bgk_prev))
  dot_args <- list(...)
  pseudo_refs <- lapply(names(pa_samplers), function(s) {
    if (verbose) cat("== reference PA:", s, "==\n")
    do.call(pa_samplers[[s]], c(list(background = dt, N_pa = N_pa_ref, pres = ref_pres,
                                     seed = seed_pseudo_base, env.rast = pa_env_rast), dot_args))
  })
  names(pseudo_refs) <- names(pa_samplers)

  pa_env_packed <- if (parallel) terra::wrap(pa_env_rast) else NULL

  # ---- 8. one (sampler, realisation) task ---------------------------------
  eval_realization <- function(s_name, r) {
    RNGkind("Mersenne-Twister")
    pa_r       <- pa_matrix[, r]
    pres_r_all <- dt[pa_r == 1L, ]
    true_abs_r <- dt[pa_r == 0L, ]
    set.seed(seed_base + r)
    pres_r <- if (nrow(pres_r_all) > max_pres)
      pres_r_all[sample(nrow(pres_r_all), max_pres), ] else pres_r_all
    N_pa_r <- max(1L, round(nrow(pres_r) / bgk_prev))
    if (nrow(pres_r) < 5L || N_pa_r < 5L)
      return(make_row(r, "skip_few_pres", nrow(pres_r_all), nrow(pres_r)))

    task_dot <- dot_args
    if (parallel && s_name == "mcmc") task_dot$num.cores <- 1L
    env_for_task <- if (parallel) terra::unwrap(pa_env_packed) else pa_env_rast

    pseudo_r <- tryCatch(
      do.call(pa_samplers[[s_name]], c(list(background = dt, N_pa = N_pa_r, pres = pres_r,
              seed = seed_pseudo_base + r, env.rast = env_for_task), task_dot)),
      error = function(e) NULL)
    if (is.null(pseudo_r) || nrow(pseudo_r) < 5L)
      return(make_row(r, "skip_sampler_null_or_short", nrow(pres_r_all), nrow(pres_r)))

    hv_pres <- tryCatch(hypervolume::hypervolume_gaussian(
      pres_r[, pc_cols], kde.bandwidth = bw, sd.count = 3,
      quantile.requested = hv_quantile, quantile.requested.type = "probability",
      chunk.size = 1000L, verbose = FALSE), error = function(e) NULL)
    hv_pa <- tryCatch(hypervolume::hypervolume_gaussian(
      pseudo_r[, pc_cols], kde.bandwidth = bw, sd.count = 3,
      quantile.requested = hv_quantile, quantile.requested.type = "probability",
      chunk.size = 1000L, verbose = FALSE), error = function(e) NULL)
    if (is.null(hv_pres) || is.null(hv_pa))
      return(make_row(r, "skip_hv_fail", nrow(pres_r_all), nrow(pres_r)))

    hv_set <- hypervolume::hypervolume_set(hv_pres, hv_pa, check.memory = FALSE, verbose = FALSE)
    ovrlp  <- hypervolume::get_volume(hv_set)[[3L]]

    relcov <- setNames(vapply(pc_cols, function(pc)
      diff(range(pseudo_r[[pc]])) / bg_range[[pc]], numeric(1)), relcov_names)

    ta_key  <- paste(round(true_abs_r$x, 4L), round(true_abs_r$y, 4L))
    ps_key  <- paste(round(pseudo_r$x, 4L), round(pseudo_r$y, 4L))
    prop_ta <- mean(ps_key %in% ta_key)

    row <- make_row(r, "ok", nrow(pres_r_all), nrow(pres_r), nrow(pseudo_r),
                    round(ovrlp, 6L), round(prop_ta, 4L), round(relcov, 4L))
    rm(hv_pres, hv_pa, hv_set, pseudo_r, pres_r, pres_r_all, true_abs_r); gc(verbose = FALSE)
    row
  }

  # ---- 9. dispatch --------------------------------------------------------
  tasks <- expand.grid(s_name = names(pa_samplers), r = seq_len(n_realizations),
                       stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE)
  if (verbose) cat("Evaluating", nrow(tasks), "tasks (",
                   length(pa_samplers), "samplers x", n_realizations, "realisations)...\n")
  results <- if (parallel)
    furrr::future_pmap(tasks, eval_realization,
                       .options = furrr::furrr_options(seed = TRUE, globals = TRUE),
                       .progress = verbose)
  else Map(eval_realization, tasks$s_name, tasks$r)
  results_by_sampler <- split(results, factor(tasks$s_name, levels = names(pa_samplers)))

  # ---- 10. per-sampler assembly + plots -----------------------------------
  sampler_results <- vector("list", length(pa_samplers))
  names(sampler_results) <- names(pa_samplers)
  for (s_name in names(pa_samplers)) {
    metrics_df <- do.call(rbind, results_by_sampler[[s_name]])
    metrics_ok <- metrics_df[metrics_df$status == "ok", , drop = FALSE]
    pseudo_ref <- pseudo_refs[[s_name]]
    if (verbose) {
      tally <- table(factor(metrics_df$status,
                            levels = c("ok", "skip_few_pres",
                                       "skip_sampler_null_or_short", "skip_hv_fail")))
      cat(sprintf("  [%s] %s\n", s_name,
                  paste(sprintf("%s=%d", names(tally), as.integer(tally)), collapse = "  ")))
    }

    # reference-realisation scalar diagnostics for the PA-plot annotation
    med_ovrlp   <- if (nrow(metrics_ok) > 0L) round(median(metrics_ok$overlap, na.rm = TRUE), 3L) else NA_real_
    ref_ta_key  <- paste(round(dt$x[ref_pa == 0L], 4L), round(dt$y[ref_pa == 0L], 4L))
    pref_key    <- paste(round(pseudo_ref$x, 4L), round(pseudo_ref$y, 4L))
    prop_ta_ref <- round(mean(pref_key %in% ref_ta_key), 3L)
    annot_txt   <- paste0("Overlap (median) = ", med_ovrlp,
                          "\nProp. true abs. = ", prop_ta_ref)

    p_pa <- ggplot() +
      geom_contour_filled(data = kde_df, aes(X, Y, z = density),
                          breaks = c(0, 0.05, 0.25, 0.5, 0.75, 1), alpha = 0.8) +
      scale_fill_viridis_d(option = "mako", name = "Background\ndensity", direction = -1) +
      geom_point(data = pseudo_ref, aes(.data[[px]], .data[[py]]),
                 colour = "darkorange", alpha = 0.45, size = 0.6) +
      geom_point(data = ref_pres, aes(.data[[px]], .data[[py]]),
                 colour = "firebrick", alpha = 0.35, size = 0.6) +
      annotate("text", x = min(dt[[px]]), y = max(dt[[py]]), label = annot_txt,
               hjust = 0, vjust = 1, size = 3.3, colour = "grey15") +
      coord_equal() +
      labs(title = paste0("Presences vs pseudo-absences [", s_name, "]  (PC1 x PC2)"),
           subtitle = "Red = presences | Orange = pseudo-absences",
           x = pc_labs[[px]], y = pc_labs[[py]]) +
      theme_classic(base_size = 12)

    if (nrow(metrics_ok) > 0L) {
      p_box_overlap <- ggplot(metrics_ok, aes(factor(0), overlap)) +
        geom_boxplot(fill = "#66C2A5", alpha = 0.7, width = 0.4) +
        geom_jitter(width = 0.08, alpha = 0.45, size = 0.9) +
        labs(title = paste0("Hypervolume overlap [", s_name, "]"),
             x = NULL, y = "Intersection volume") +
        theme_classic(base_size = 12) + theme(axis.text.x = element_blank())
      cov_long <- do.call(rbind, lapply(pc_cols, function(pc)
        data.frame(axis = pc, value = metrics_ok[[paste0("rel_cov_", pc)]])))
      p_box_coverage <- ggplot(cov_long, aes(axis, value, fill = axis)) +
        geom_boxplot(alpha = 0.7, width = 0.5) + scale_y_continuous(limits = c(0, 1)) +
        labs(title = paste0("Per-axis coverage [", s_name, "]"),
             x = NULL, y = "rel. coverage [0,1]") +
        theme_classic(base_size = 12) + theme(legend.position = "none")
      p_box_trueabs <- ggplot(metrics_ok, aes(factor(0), prop_true_abs)) +
        geom_boxplot(fill = "#E78AC3", alpha = 0.7, width = 0.4) +
        scale_y_continuous(limits = c(0, 1)) +
        labs(title = paste0("True absences [", s_name, "]"),
             x = NULL, y = "Proportion") +
        theme_classic(base_size = 12) + theme(axis.text.x = element_blank())
    } else {
      empty_p <- ggplot() + annotate("text", x = 0.5, y = 0.5,
                 label = "No completed realisations") + theme_void()
      p_box_overlap <- p_box_coverage <- p_box_trueabs <- empty_p
    }

    # sampling-bias marginal densities per axis (reference realisation):
    # overlaid Background / Pseudo-absences / Presences density, annotated with
    # that axis's relative coverage.
    bias_df <- rbind(
      data.frame(grp = "Background",      dt[, pc_cols, drop = FALSE]),
      data.frame(grp = "Pseudo-absences", pseudo_ref[, pc_cols, drop = FALSE]),
      data.frame(grp = "Presences",       ref_pres[, pc_cols, drop = FALSE]))
    bias_df$grp <- factor(bias_df$grp,
                          levels = c("Background", "Pseudo-absences", "Presences"))
    bias_cols <- c(Background = "grey50", "Pseudo-absences" = "darkorange",
                   Presences = "firebrick")
    p_bias <- setNames(lapply(pc_cols, function(pc) {
      cov_pc <- round(diff(range(pseudo_ref[[pc]])) / bg_range[[pc]], 3L)
      ggplot(bias_df, aes(x = .data[[pc]], colour = grp, fill = grp)) +
        geom_density(alpha = 0.25, linewidth = 0.7) +
        scale_colour_manual(values = bias_cols) +
        scale_fill_manual(values = bias_cols) +
        annotate("text", x = Inf, y = Inf, label = paste0("Coverage = ", cov_pc),
                 hjust = 1.1, vjust = 1.5, size = 3.2) +
        labs(title = paste0("Bias — ", pc, " [", s_name, "]"),
             x = pc_labs[[pc]], y = "Density", colour = NULL, fill = NULL) +
        theme_classic(base_size = 11) + theme(legend.position = "bottom")
    }), pc_cols)

    pres_out <- ref_pres[, c("x", "y", pc_cols, "suit")]; pres_out$pa <- 1L
    pa_out   <- pseudo_ref[, c("x", "y", pc_cols, "suit")]; pa_out$pa <- 0L
    sampler_results[[s_name]] <- list(
      pseudo_ref  = pseudo_ref,
      pseudo_vect = terra::vect(pseudo_ref, geom = c("x", "y"), crs = terra::crs(envData)),
      dataset     = rbind(pres_out, pa_out),
      metrics     = metrics_df,
      plots       = list(p_pa = p_pa, p_box_overlap = p_box_overlap,
                         p_box_coverage = p_box_coverage, p_box_trueabs = p_box_trueabs,
                         p_bias = p_bias))
  }

  # ---- 11. geographic rasters (reference realisation) ---------------------
  suit_rast <- terra::rast(dt[, c("x", "y", "suit")], type = "xyz", crs = terra::crs(envData))
  pa_rast   <- terra::rast(data.frame(x = dt$x, y = dt$y, pa = as.numeric(dt$pa)),
                           type = "xyz", crs = terra::crs(envData))

  # ---- 12. parameter record -----------------------------------------------
  parameters <- list(
    k = k, pc_cols = pc_cols, mu = mu, sigma = sigma, rho = rho, Sigma = Sigma,
    bgk_prev = bgk_prev, n_realizations = n_realizations, max_pres = max_pres,
    seed_base = seed_base, seed_pseudo_base = seed_pseudo_base,
    bw = bw, bw_source = bw_source, hv_quantile = hv_quantile,
    sampler_names = names(pa_samplers), sampler_dots = dot_args,
    n_background = nrow(dt), n_pres_total = nrow(ref_pres_all), n_pres_ref = nrow(ref_pres),
    n_pa_ref = N_pa_ref, env_layer_names = names(envData),
    r_version = R.version.string, timestamp = Sys.time())

  list(
    niche = list(mu = mu, sigma = sigma, rho = rho, Sigma = Sigma, equation = eq_nd),
    parameters = parameters,
    suit_rast = suit_rast, pa_rast = pa_rast,
    response_curves = list(data = rc_data, plots = rc_plots, logit_plots = rc_logit_plots),
    samplers = sampler_results,
    plots = list(p_niche = p_niche),
    bw = bw, background = dt)
}

# =============================================================================
# REPLICATION-RECORD PAGE  (k-D port of summarize_parameters)
# =============================================================================

#' One-page replication record for a multi-species 5-D experiment.
#'
#' Reads each result's $parameters slot and renders three gridExtra tables on a
#' single page: (1) niche parameters per species — mu_PC1..PCk and sd_PC1..PCk
#' (generalised from the 2-D mu_PC1/PC2 + sigma_PC1/PC2), rho, prevalence,
#' n_pres_ref, n_pa_ref; (2) run controls incl. the fixed bandwidth per axis;
#' (3) environment & provenance. print() it inside an open PDF device.
#'
#' @param species_results Named list of virtualSpecies_nd() returns.
#' @param species_labels  Optional display labels (defaults to the names).
#' @param pca_var_exp     Optional numeric vector of PC variance fractions.
#' @param extra           Optional named list of free-form key/value rows.
summarize_parameters_nd <- function(species_results, species_labels = NULL,
                                    pca_var_exp = NULL, extra = list()) {
  if (!requireNamespace("gridExtra", quietly = TRUE))
    stop("summarize_parameters_nd() requires 'gridExtra'.")
  stopifnot(length(species_results) >= 1L)
  pl <- lapply(species_results, `[[`, "parameters")
  if (any(vapply(pl, is.null, logical(1))))
    stop("one or more results lack a $parameters slot.")
  sp_names <- names(species_results)
  if (is.null(species_labels)) species_labels <- sp_names
  p1  <- pl[[1L]]; k <- p1$k; pcc <- p1$pc_cols

  # Table 1: niche parameters per species (generalised to k axes)
  niche <- data.frame(species = sp_names, label = species_labels,
                      stringsAsFactors = FALSE)
  for (i in seq_len(k))
    niche[[paste0("mu_", pcc[i])]] <- round(vapply(pl, function(p) p$mu[i], numeric(1)), 2)
  for (i in seq_len(k))
    niche[[paste0("sd_", pcc[i])]] <- round(vapply(pl, function(p) p$sigma[i], numeric(1)), 2)
  niche$rho        <- vapply(pl, function(p) if (is.matrix(p$rho)) NA_real_ else as.numeric(p$rho), numeric(1))
  niche$prevalence <- round(vapply(species_results, function(s) mean(s$background$pa), numeric(1)), 4)
  niche$n_pres_ref <- vapply(pl, function(p) as.integer(p$n_pres_ref), integer(1))
  niche$n_pa_ref   <- vapply(pl, function(p) as.integer(p$n_pa_ref), integer(1))

  # Table 2: run controls (from species 1)
  run <- data.frame(
    parameter = c("k (dims)", "pc_cols", "n_realizations", "max_pres", "bgk_prev",
                  "seed_base", "seed_pseudo_base", "hv_quantile", "n_background",
                  "bw_source", "samplers",
                  paste0("bw (", paste(pcc, collapse = ","), ")")),
    value = c(p1$k, paste(pcc, collapse = ","), p1$n_realizations, p1$max_pres,
              p1$bgk_prev, p1$seed_base, p1$seed_pseudo_base, p1$hv_quantile,
              p1$n_background, p1$bw_source, paste(p1$sampler_names, collapse = ","),
              paste(round(p1$bw, 4), collapse = ", ")),
    stringsAsFactors = FALSE)

  # Table 3: environment & provenance
  pca_str <- if (!is.null(pca_var_exp))
    paste(sprintf("%s=%.1f%%", pcc, 100 * pca_var_exp[seq_len(k)]), collapse = ", ")
  else "(not provided)"
  env <- data.frame(
    field = c("env layers", "PCA variance (PC1..PCk)", "R version", "timestamp",
              names(extra)),
    value = c(paste(p1$env_layer_names, collapse = ", "), pca_str, p1$r_version,
              format(p1$timestamp, "%Y-%m-%d %H:%M:%S"),
              vapply(extra, function(x) paste(as.character(x), collapse = ", "),
                     character(1))),
    stringsAsFactors = FALSE)

  tt <- gridExtra::ttheme_minimal(
    base_size = 7,
    core    = list(fg_params = list(hjust = 0, x = 0.02)),
    colhead = list(fg_params = list(hjust = 0, x = 0.02, fontface = "bold")))
  g_niche <- gridExtra::tableGrob(niche, rows = NULL, theme = tt)
  g_run   <- gridExtra::tableGrob(run,   rows = NULL, theme = tt)
  g_env   <- gridExtra::tableGrob(env,   rows = NULL, theme = tt)
  hdr <- function(t, s = 11) grid::textGrob(
    t, x = 0.02, hjust = 0, gp = grid::gpar(fontface = "bold", fontsize = s))
  page <- gridExtra::arrangeGrob(
    hdr("GaussNiche 5-D run summary — replication record", 13),
    hdr("1. Niche parameters (per species)"), g_niche,
    hdr("2. Run controls"), g_run,
    hdr("3. Environment & provenance"), g_env,
    ncol = 1,
    heights = grid::unit(c(0.45, 0.3, 0.22 * (nrow(niche) + 1) + 0.2,
                           0.3, 0.22 * (nrow(run) + 1) + 0.2,
                           0.3, 0.22 * (nrow(env) + 1) + 0.2), "in"))
  patchwork::wrap_elements(page)
}
