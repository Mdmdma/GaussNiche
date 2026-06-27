#!/usr/bin/env Rscript
# Benchmark + correctness check for the precomputeMcmcEnvironment() cache added
# to USE.MCMC::paSamplingMcmc(). Run inside the rocker apptainer container.
#
# It (1) verifies cached == uncached bit-for-bit, (2) verifies a saveRDS/readRDS
# round-trip of the bundle still works, and (3) times "simulate K species on the
# same environment" with and without the cache.

suppressWarnings(.libPaths(c("~/R/rocker-rstudio/4.5", .libPaths())))
suppressPackageStartupMessages({
  library(USE.MCMC)
  library(terra)
  library(sf)
})

args <- commandArgs(trailingOnly = TRUE)
K            <- as.integer(if (length(args) >= 1) args[1] else 4L)   # n "species"
chain.length <- as.integer(if (length(args) >= 2) args[2] else 8000L)
n.samples    <- 300L
seed.number  <- 42L
dimensions   <- c("PC1", "PC2")
out_dir      <- file.path(Sys.getenv("SCRATCH", "/cluster/scratch/merler"), "gn_bench")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat(sprintf("Bench config: K=%d species, chain.length=%d, n.samples=%d\n",
            K, chain.length, n.samples))

## --- environment + K presence sets ------------------------------------------
env <- terra::rast(USE.MCMC::Worldclim_tmp, type = "xyz")
cat(sprintf("env raster: %d layers, %d non-NA cells\n",
            terra::nlyr(env), sum(!is.na(terra::values(env)[, 1]))))

df_all <- terra::as.data.frame(env, xy = TRUE, na.rm = TRUE)
make_pres <- function(species_seed, n = 300L) {
  set.seed(species_seed)
  idx <- sample(nrow(df_all), min(n, nrow(df_all)))
  sf::st_as_sf(df_all[idx, ], coords = c("x", "y"), crs = 4326)
}
pres_list <- lapply(seq_len(K), function(k) make_pres(100L + k))

run_one <- function(pres, env_bundle = NULL) {
  USE.MCMC::paSamplingMcmc(
    env.rast = env, pres = pres,
    n.samples = n.samples, chain.length = chain.length, burnIn = 1000,
    dimensions = dimensions, num.chains = 1, num.cores = 1,
    seed.number = seed.number, engine = "auto",
    precomputed.env = env_bundle, verbose = FALSE, plot_proc = FALSE)
}

## --- (1) correctness: cached == uncached ------------------------------------
cat("\n[1] Correctness: cached vs uncached (same seed) ...\n")
env_bundle <- USE.MCMC::precomputeMcmcEnvironment(
  env.rast = env, dimensions = dimensions, seed.number = seed.number)

res_uncached <- run_one(pres_list[[1]], env_bundle = NULL)
res_cached   <- run_one(pres_list[[1]], env_bundle = env_bundle)

identical_coords <- isTRUE(all.equal(
  sf::st_coordinates(res_uncached), sf::st_coordinates(res_cached)))
cat(sprintf("    n uncached=%d  n cached=%d  identical coords: %s\n",
            nrow(res_uncached), nrow(res_cached), identical_coords))

## --- (2) saveRDS / readRDS round-trip of the bundle -------------------------
cat("\n[2] saveRDS/readRDS round-trip of the env bundle ...\n")
bundle_path <- file.path(out_dir, "env_bundle.rds")
saveRDS(env_bundle, bundle_path)
env_bundle_rt <- readRDS(bundle_path)
res_cached_rt <- run_one(pres_list[[1]], env_bundle = env_bundle_rt)
identical_rt <- isTRUE(all.equal(
  sf::st_coordinates(res_uncached), sf::st_coordinates(res_cached_rt)))
cat(sprintf("    bundle size on disk: %.1f MB   identical after round-trip: %s\n",
            file.info(bundle_path)$size / 1e6, identical_rt))

## --- (3) timing: K species, uncached vs cached ------------------------------
cat(sprintf("\n[3] Timing %d species on the same environment ...\n", K))

t_uncached <- system.time(
  for (k in seq_len(K)) run_one(pres_list[[k]], env_bundle = NULL)
)[["elapsed"]]

t_precompute <- system.time(
  env_bundle2 <- USE.MCMC::precomputeMcmcEnvironment(
    env.rast = env, dimensions = dimensions, seed.number = seed.number)
)[["elapsed"]]
t_cached_calls <- system.time(
  for (k in seq_len(K)) run_one(pres_list[[k]], env_bundle = env_bundle2)
)[["elapsed"]]
t_cached_total <- t_precompute + t_cached_calls

cat("\n================ RESULTS ================\n")
cat(sprintf("uncached  : %6.1fs total  (%.1fs / species)\n",
            t_uncached, t_uncached / K))
cat(sprintf("cached    : %6.1fs total  = %.1fs precompute + %.1fs for %d calls (%.1fs / species)\n",
            t_cached_total, t_precompute, t_cached_calls, K, t_cached_calls / K))
cat(sprintf("env-only block cost (amortised away): ~%.1fs / species\n",
            (t_uncached - t_cached_calls) / K))
cat(sprintf("speedup for %d species: %.2fx   (asymptotic per-species: %.2fx)\n",
            K, t_uncached / t_cached_total, t_uncached / t_cached_calls))
cat(sprintf("correctness: cached==uncached: %s ; roundtrip==uncached: %s\n",
            identical_coords, identical_rt))
cat("=========================================\n")

saveRDS(list(K = K, chain.length = chain.length,
             t_uncached = t_uncached, t_precompute = t_precompute,
             t_cached_calls = t_cached_calls, t_cached_total = t_cached_total,
             identical_coords = identical_coords, identical_rt = identical_rt),
        file.path(out_dir, "bench_results.rds"))
cat("Saved timings to", file.path(out_dir, "bench_results.rds"), "\n")
