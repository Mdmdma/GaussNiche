# 5_highdim_species.R  —  the 4 virtual species in 5-D (5d-niche branch)
# =============================================================================
# Reproduces the 2x2 factorial (generalist/specialist x common/rare) from
# 3_ virtualSpecies_cases.R, but in the locked 5-D 'lean_natural' E-space.
# mu/sigma are defined as MULTIPLES of each PC's background sd (read at runtime),
# so the design is reproducible from the data:
#   generalist sigma = 1.0 x sd     specialist sigma = 0.4 x sd
#   common     mu    = 0            rare        mu    = 0.8 x sd  (all axes)
#
# Usage (positional args):
#   Rscript 5_highdim_species.R [mode] [n_workers]
#     mode = "smoke" -> 1 species, 2 realisations, tiny chain (validation)
#     mode = "full"  -> all 4 species (default); N_REALIZATIONS env overrides
#
# Inputs : <scratch>/GaussNiche/env5d/{final_stack_lean_natural.tif, background_5d.rds}
# Outputs: <scratch>/GaussNiche/results5d/{<sp>_<mode>.rds, summary_<mode>.rds}
# Run on a compute node via submit_5d_species.sh.
# =============================================================================

suppressPackageStartupMessages({
  .libPaths(c("~/R/rocker-rstudio/4.5", .libPaths()))
  library(terra); library(future); library(furrr)
})
source("virtualSpecies_nd_fn.R")

args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args) >= 1) args[1] else "full"
n_workers <- if (length(args) >= 2) as.integer(args[2]) else
  as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4"))

if (mode == "smoke") {
  N_REAL <- 2L; MAXP <- 80L; CHAIN <- 2000L; BURN <- 500L
  SPECIES_SET <- "sp1"; PAR <- FALSE
} else {
  N_REAL <- as.integer(Sys.getenv("N_REALIZATIONS", "20"))
  MAXP <- 300L; CHAIN <- 20000L; BURN <- 1000L
  SPECIES_SET <- "all"; PAR <- TRUE
}
cat(sprintf("mode=%s  n_real=%d  max_pres=%d  chain=%d  parallel=%s  workers=%d\n",
            mode, N_REAL, MAXP, CHAIN, PAR, n_workers))

user    <- Sys.getenv("USER"); scratch <- file.path("/cluster/scratch", user)
env_dir <- file.path(scratch, "GaussNiche", "env5d")
out_dir <- file.path(scratch, "GaussNiche", "results5d")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

envData <- terra::rast(file.path(env_dir, "final_stack_lean_natural.tif"))
dt      <- readRDS(file.path(env_dir, "background_5d.rds"))
pc_cols <- paste0("PC", 1:5)
stopifnot(all(pc_cols %in% names(dt)))

pc_sd <- vapply(dt[pc_cols], stats::sd, numeric(1))
cat("PC background sd:", paste(round(pc_sd, 3), collapse = ", "), "\n")

# --- 4-species design --------------------------------------------------------
# Covariances: generalist sigma = 1.0 x per-PC sd, specialist = 0.5 x sd
# (diagonal Sigma, rho = 0).
sig_gen  <- 1.0 * pc_sd
sig_spec <- 0.5 * pc_sd

# COMMON optimum = mode of the 5-D background kernel density (the densest point
# in PC space). RARE optimum = user-specified peripheral point.
bg5 <- as.matrix(dt[, pc_cols])
mode_estimate <- function(X) {
  # point of maximum kernel density; ks::kde with a normal-scale bandwidth,
  # evaluated at the background points (argmax). Falls back to a kNN density
  # mode (smallest k-th-NN distance = densest) if ks is unavailable/errors.
  if (requireNamespace("ks", quietly = TRUE)) {
    est <- tryCatch(ks::kde(x = X, H = ks::Hns(X), eval.points = X)$estimate,
                    error = function(e) NULL)
    if (!is.null(est)) return(as.numeric(X[which.max(est), ]))
  }
  k <- min(200L, nrow(X) - 1L)
  as.numeric(X[which.min(FNN::knn.dist(X, k = k)[, k]), ])
}
mu_common <- mode_estimate(bg5);                          names(mu_common) <- pc_cols
mu_rare   <- c(PC1 = 1.55, PC2 = 1.13, PC3 = 1.4, PC4 = 0.8, PC5 = 1.2)[pc_cols]
cat("mu_common (KDE mode):", paste(round(mu_common, 3), collapse = ", "), "\n")
cat("mu_rare            :", paste(round(mu_rare, 3),   collapse = ", "), "\n")

SPECIES <- list(
  sp1_generalist_common = list(mu = mu_common, sigma = sig_gen),
  sp2_specialist_common = list(mu = mu_common, sigma = sig_spec),
  sp3_generalist_rare   = list(mu = mu_rare,   sigma = sig_gen),
  sp4_specialist_rare   = list(mu = mu_rare,   sigma = sig_spec))
if (SPECIES_SET == "sp1") SPECIES <- SPECIES["sp1_generalist_common"]

# --- fixed background bandwidth (once, for cross-species comparability) ------
bw <- compute_bandwidth_nd(dt, pc_cols)
cat("fixed bandwidth:", paste(round(bw, 4), collapse = ", "), "\n")

# --- precompute the MCMC environment ONCE and reuse across all species/realis. --
# precomputeMcmcEnvironment caches the PCA + environmental GMM fit/density +
# proposal covariance + distance threshold (constant across presence sets), so
# paSamplingMcmc skips that work on every call -> large speedup.
# CRITICAL for statistical soundness: the bundle captures an RNG state that
# paSamplingMcmc would restore on every call (collapsing the realisations to a
# shared random stream). We NEUTRALISE it -- keep the field name (so the bundle
# still validates) but set it to NULL so NO restore happens. Each pa_mcmc call
# then seeds its own chain with seed_pseudo_base + r -> independent realisations.
cat("\nPrecomputing MCMC environment (PCA + env GMM) once for reuse...\n")
env_bundle <- USE.MCMC::precomputeMcmcEnvironment(
  env.data.raster = envData, dimensions = pc_cols, seed.number = 123, verbose = TRUE)
env_bundle["rng_state"] <- list(NULL)   # keep name, NULL value -> restore skipped

if (mode == "smoke") {
  # soundness check: with the shared bundle, different seeds must give different
  # pseudo-absences (independent) and the same seed must reproduce (deterministic).
  cat("\n== MCMC RNG-independence soundness check ==\n")
  set.seed(1L); tp <- dt[sample(nrow(dt), 50L), ]
  mk <- function(sd) pa_mcmc(dt, 50L, tp, seed = sd, env.rast = envData,
                             dimensions = pc_cols, chain.length = 2000L, burnIn = 300L,
                             precomputed.env = env_bundle)[, pc_cols, drop = FALSE]
  A1 <- mk(11L); A2 <- mk(22L); A1b <- mk(11L)
  cat("  different seeds -> different draws (independent):",
      !isTRUE(all.equal(A1, A2)),  "\n")
  cat("  same seed -> identical draws (deterministic)    :",
      isTRUE(all.equal(A1, A1b)), "\n")
}

results <- list()
for (nm in names(SPECIES)) {
  cat("\n##########  ", nm, "  ##########\n", sep = "")
  sp <- SPECIES[[nm]]
  res <- virtualSpecies_nd(
    dt = dt, envData = envData, mu = sp$mu, sigma = sp$sigma, rho = 0,
    pc_cols = pc_cols, bgk_prev = 1,
    pa_samplers = list(random = pa_random, mcmc = pa_mcmc),
    n_realizations = N_REAL, max_pres = MAXP,
    bw = bw, pa_env_rast = envData,
    dimensions = pc_cols, chain.length = CHAIN, burnIn = BURN,
    environmental.cutof.percentile = 0.001,
    species.cutoff.threshold = 0.7,
    precomputed.env = env_bundle,
    verbose = TRUE, parallel = PAR, n_workers = n_workers)
  saveRDS(res, file.path(out_dir, paste0(nm, "_", mode, ".rds")))
  results[[nm]] <- list(metrics = lapply(res$samplers, `[[`, "metrics"),
                        prevalence = mean(res$background$pa))
}
saveRDS(results, file.path(out_dir, paste0("summary_", mode, ".rds")))

cat("\n================ SUMMARY ================\n")
for (nm in names(results)) {
  cat(sprintf("\n%s  (prevalence %.4f)\n", nm, results[[nm]]$prevalence))
  for (s in names(results[[nm]]$metrics)) {
    m <- results[[nm]]$metrics[[s]]; ok <- m[m$status == "ok", , drop = FALSE]
    cat(sprintf("  %-7s  ok=%d/%d  overlap(median)=%s\n", s, nrow(ok), nrow(m),
                if (nrow(ok)) sprintf("%.4f", median(ok$overlap, na.rm = TRUE)) else "NA"))
  }
}
cat("\nSaved to ", out_dir, "\n", sep = "")
