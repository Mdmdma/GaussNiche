# bench_hsm.R --- decompose the per-realisation cost of the 5-D HSM pipeline
# =============================================================================
# Times each component of one (sampler, realisation) unit single-threaded so we
# can rank speedups before the settings sweep:
#   pa_random / pa_buffer_out / pa_mcmc(chain 20k/10k/5k)  -- the samplers
#   hypervolume_gaussian x2                                -- the diagnostics
#   sdm_hook (fit + predict-test + predict-background)     -- the downstream HSM
#   predict-over-background: full (~16.5k cells) vs 4k subsample
# Single species (sp1), a few realisations; means reported. Run on a compute node
# via sbatch/submit_bench_hsm.sh. NOT a login-node task (fits real models).
# =============================================================================

suppressPackageStartupMessages({
  .libPaths(c("~/R/rocker-rstudio/4.5", .libPaths()))
  library(terra); library(USE.MCMC); library(hypervolume); library(FNN)
})
source("virtualSpecies_nd_fn.R"); source("hsm_eval.R"); source("pa_buffer.R")

NB <- as.integer(Sys.getenv("BENCH_REALS", "6"))
et <- function(expr) as.numeric(system.time(expr)["elapsed"])

user <- Sys.getenv("USER"); scratch <- file.path("/cluster/scratch", user)
env_dir <- file.path(scratch, "GaussNiche", "env5d")
envData <- terra::rast(file.path(env_dir, "final_stack_lean_natural.tif"))
dt <- readRDS(file.path(env_dir, "background_5d.rds"))
pc_cols <- paste0("PC", 1:5); raw_cols <- names(envData)
rv <- as.data.frame(terra::extract(envData, as.matrix(dt[, c("x", "y")])))[, raw_cols, drop = FALSE]
dt <- cbind(dt, rv); rownames(dt) <- NULL
cat(sprintf("background cells: %d\n", nrow(dt)))

# sp1 generalist-common niche (matches run_5d_hsm)
pc_sd <- vapply(dt[pc_cols], stats::sd, numeric(1))
mu <- vapply(dt[pc_cols], function(v) v[which.max(table(round(v, 1)))], numeric(1))  # cheap mode proxy
mu <- setNames(as.numeric(colMeans(dt[pc_cols])), pc_cols)                            # central optimum (common)
sigma <- 1.0 * pc_sd
Sigma <- diag(sigma) %*% diag(length(sigma)) %*% diag(sigma)
peak <- mvtnorm::dmvnorm(matrix(mu, 1), mu, Sigma)
suit <- mvtnorm::dmvnorm(as.matrix(dt[pc_cols]), mu, Sigma) / peak
dt$suit <- suit
bw <- compute_bandwidth_nd(dt, pc_cols)

cat("precompute MCMC env...\n")
env_bundle <- USE.MCMC::precomputeMcmcEnvironment(env.data.raster = envData,
              dimensions = pc_cols, seed.number = 123, verbose = FALSE)
env_bundle["rng_state"] <- list(NULL)

predictor_sets <- list(pc5 = pc_cols, raw12 = raw_cols)
hook <- make_sdm_hook("sp1", predictor_sets, algorithms = HSM_ALGORITHMS,
                      test_frac = 0.3, min_pres = 12L, seed = 1100L, libpaths = .libPaths(),
                      source_file = normalizePath("hsm_eval.R"))

hv1 <- function(df) hypervolume::hypervolume_gaussian(df[, pc_cols], kde.bandwidth = bw,
         sd.count = 3, quantile.requested = 0.95, quantile.requested.type = "probability",
         chunk.size = 1000L, verbose = FALSE)
mcmc_call <- function(pres, N_pa, chain) pa_mcmc(dt, N_pa, pres = pres, seed = 123,
         env.rast = envData, dimensions = pc_cols, chain.length = chain, burnIn = 1000,
         species.cutoff.threshold = 0.7, environmental.cutof.percentile = 0.001,
         precomputed.env = env_bundle)

acc <- list()
add <- function(k, v) acc[[k]] <<- c(acc[[k]], v)
for (r in seq_len(NB)) {
  set.seed(42 + r); pa <- rbinom(nrow(dt), 1L, suit)
  pres_all <- dt[pa == 1L, ]; set.seed(42 + r)
  pres <- if (nrow(pres_all) > 300) pres_all[sample(nrow(pres_all), 300), ] else pres_all
  Npa <- max(1L, round(nrow(pres)))
  add("pa_random",  et(rnd <- pa_random(dt, Npa, pres, 124)))
  add("pa_buffer",  et(pa_buffer_out(dt, Npa, pres, 124, env.rast = envData, buffer_km = 50)))
  add("pa_mcmc_20k", et(mc <- mcmc_call(pres, Npa, 20000)))
  add("pa_mcmc_10k", et(mcmc_call(pres, Npa, 10000)))
  add("pa_mcmc_5k",  et(mcmc_call(pres, Npa, 5000)))
  add("hypervol_x2", et({ hv1(pres); hv1(mc) }))
  add("sdm_hook_all", et(hook(pres_r = pres, pseudo_r = mc, background = dt,
                              pc_cols = pc_cols, realization = r, sampler = "mcmc")))
  # predict-over-background: full vs 4k subsample, for one RF model
  tr <- rbind(cbind(pres[pc_cols], pa = 1L), cbind(mc[pc_cols], pa = 0L))
  m_rf <- ranger::ranger(stats::reformulate(pc_cols, "pa"), data = tr, num.trees = 500L, num.threads = 1L)
  sub <- dt[sample(nrow(dt), 4000L), ]
  add("predict_bg_full", et(predict(m_rf, data = dt, num.threads = 1L)))
  add("predict_bg_4k",   et(predict(m_rf, data = sub, num.threads = 1L)))
  cat(sprintf("  realisation %d done\n", r))
}

cat("\n================= PER-REALISATION COST (mean over", NB, "realisations, seconds) =================\n")
ord <- c("pa_random","pa_buffer","pa_mcmc_5k","pa_mcmc_10k","pa_mcmc_20k",
         "hypervol_x2","sdm_hook_all","predict_bg_full","predict_bg_4k")
for (k in ord) if (!is.null(acc[[k]]))
  cat(sprintf("  %-16s %7.3f  (sd %.3f)\n", k, mean(acc[[k]]), stats::sd(acc[[k]])))

cat("\nDecomposition note: pa_mcmc(chain) ~ fixed + k*chain.length;\n",
    "  fixed (densityMclust + PC-unwrap + presence-extract + remap) ~ 2*t(5k) - t(10k);\n",
    "  chain-linear part per 1k steps ~ (t(20k) - t(5k))/15.\n")
f5 <- mean(acc$pa_mcmc_5k); f10 <- mean(acc$pa_mcmc_10k); f20 <- mean(acc$pa_mcmc_20k)
cat(sprintf("  est fixed cost      ~ %.3f s\n", max(0, 2*f5 - f10)))
cat(sprintf("  est per-1k-steps    ~ %.3f s\n", (f20 - f5)/15))
