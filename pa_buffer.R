# pa_buffer.R --- geographic "buffer-out" pseudo-absence sampler
# =============================================================================
# A spatial baseline benchmark (cf. the original USE paper, Da Re et al. 2023,
# which compares USE against random and "buffer-out" geographic samplers):
# pseudo-absences are drawn uniformly from the geographic background but ONLY
# from cells lying OUTSIDE a buffer of `buffer_km` around the presence points,
# keeping them away from the (sampled) species localities. Unlike the
# environmental samplers it needs no env.rast / PCA --- only the x,y coordinates
# --- so it works unchanged in 2-D and k-D. It is "better than random" in that it
# avoids false absences near presences, yet, being geographic, it still
# over-represents the commonest environments (the limitation USE / uniform+
# address).
#
# Interface (the GaussNiche sampler contract): same signature shape as
# pa_random / pa_mcmc; env.rast and any sampler-specific dots are swallowed.
# Distances are Euclidean in a metric CRS (default LAEA Europe, EPSG:3035, valid
# for the C/W-Europe study extent); pass crs_metric for other regions.
# Returns up to N_pa rows of `background` (fewer if the buffer leaves too few
# cells --- the caller then records a skip, exactly as for the other samplers; we
# deliberately do NOT pad with near-presence cells, which would defeat the
# method).
# =============================================================================

pa_buffer_out <- function(background, N_pa, pres = NULL, seed = 123,
                          env.rast = NULL, buffer_km = 50,
                          crs_ll = "EPSG:4326", crs_metric = "EPSG:3035", ...) {
  if (!requireNamespace("terra", quietly = TRUE))
    stop("pa_buffer_out requires 'terra'.")
  if (!requireNamespace("FNN", quietly = TRUE))
    stop("pa_buffer_out requires 'FNN'.")
  if (is.null(pres) || !nrow(pres))
    stop("pa_buffer_out requires non-empty 'pres' (presence rows from background).")
  set.seed(seed)

  # Project lon/lat -> metric CRS for fast Euclidean nearest-neighbour distances
  # (great-circle is unnecessary at a 50 km scale within Europe).
  to_m <- function(df) {
    v <- terra::vect(cbind(df$x, df$y), type = "points", crs = crs_ll)
    terra::crds(terra::project(v, crs_metric))
  }
  pr_m <- to_m(pres)
  bg_m <- to_m(background)
  if (anyNA(pr_m) || anyNA(bg_m))
    stop("pa_buffer_out: projection to ", crs_metric, " produced NA coordinates ",
         "(points outside the CRS domain?).")
  stopifnot(ncol(pr_m) == 2L, ncol(bg_m) == 2L)

  # Nearest-presence distance (m) per background cell; keep those beyond buffer.
  d_near  <- FNN::get.knnx(pr_m, bg_m, k = 1L)$nn.dist[, 1L]
  outside <- which(d_near > buffer_km * 1000)
  if (length(outside) == 0L) return(background[0, ])

  k   <- min(N_pa, length(outside))
  out <- background[sample(outside, k, replace = FALSE), ]
  row.names(out) <- NULL
  out
}
