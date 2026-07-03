# pcmatrix_at_setting.R --- per-species PC-matrix figures (like the appendix
# pcmatrix-5d-sp* on pages 22-25) for uniform+ at a CHOSEN cutoff setting.
# Lower triangle = RND, upper triangle = uniform+ at (env, pres). One figure per
# species. Minimal compute: only the one-off reference PA draw per sampler is
# needed (no realisation loop, no hypervolume). Faithful to run_5d_experiment's
# species definitions + plot_pc_matrix.
#
# Usage:  Rscript pcmatrix_at_setting.R
#   env overrides: ENV_CUTOFF (default 0.45), SPECIES_CUTOFF (default 0.5),
#                  CHAIN (30000), MAX_PRES (300), GN_ENV_DIR, GN_OUT_DIR
# =============================================================================
suppressPackageStartupMessages({
  .libPaths(c("~/R/rocker-rstudio/4.5", .libPaths()))
  library(terra); library(USE.MCMC); library(ggplot2)
  library(hypervolume); library(FNN); library(future); library(furrr)
})
source("virtualSpecies_nd_fn.R")
source("experiment_plots.R")

ENVc  <- as.numeric(Sys.getenv("ENV_CUTOFF", "0.45"))
PRESc <- as.numeric(Sys.getenv("SPECIES_CUTOFF", "0.5"))
CHAIN <- as.integer(Sys.getenv("CHAIN", "30000")); BURN <- 1000L
MAXP  <- as.integer(Sys.getenv("MAX_PRES", "300"))
tag   <- sprintf("env%02d_pres%02d", round(ENVc * 100), round(PRESc * 100))
cat(sprintf("PC-matrix at uniform+ env=%.3f pres=%.3f chain=%d  tag=%s\n", ENVc, PRESc, CHAIN, tag))

user    <- Sys.getenv("USER"); scratch <- file.path("/cluster/scratch", user)
env_dir <- Sys.getenv("GN_ENV_DIR", file.path(scratch, "GaussNiche", "env5d"))
out_dir <- Sys.getenv("GN_OUT_DIR", file.path(scratch, "GaussNiche", "results5d_experiment"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# --- environment + species catalogue (mirrors run_5d_experiment) -------------
envData <- terra::rast(file.path(env_dir, "final_stack_lean_natural.tif"))
dt      <- readRDS(file.path(env_dir, "background_5d.rds"))
pc_cols <- paste0("PC", 1:5)
stopifnot(all(pc_cols %in% names(dt)))
pc_sd   <- vapply(dt[pc_cols], stats::sd, numeric(1))
sig_gen <- 1.0 * pc_sd; sig_spec <- 0.5 * pc_sd
bg5     <- as.matrix(dt[, pc_cols])
mode_estimate <- function(X) {
  if (requireNamespace("ks", quietly = TRUE)) {
    est <- tryCatch(ks::kde(x = X, H = ks::Hns(X), eval.points = X)$estimate, error = function(e) NULL)
    if (!is.null(est)) return(as.numeric(X[which.max(est), ]))
  }
  k <- min(200L, nrow(X) - 1L)
  as.numeric(X[which.min(FNN::knn.dist(X, k = k)[, k]), ])
}
mu_common <- mode_estimate(bg5); names(mu_common) <- pc_cols
mu_rare   <- c(PC1 = 1.55, PC2 = 1.13, PC3 = 1.4, PC4 = 0.8, PC5 = 1.2)[pc_cols]
SPECIES <- list(
  sp1_generalist_common = list(label = "Generalist · common", mu = mu_common, sigma = sig_gen),
  sp2_specialist_common = list(label = "Specialist · common", mu = mu_common, sigma = sig_spec),
  sp3_generalist_rare   = list(label = "Generalist · rare",   mu = mu_rare,   sigma = sig_gen),
  sp4_specialist_rare   = list(label = "Specialist · rare",   mu = mu_rare,   sigma = sig_spec))

bw <- compute_bandwidth_nd(dt, pc_cols)
env_bundle <- USE.MCMC::precomputeMcmcEnvironment(
  env.rast = envData, dimensions = pc_cols, seed.number = 123, verbose = FALSE)
env_bundle["rng_state"] <- list(NULL)

# --- one PC-matrix per species (RND vs uniform+ at env/pres) ------------------
for (nm in names(SPECIES)) {
  cat("\n##########  ", nm, "  ##########\n", sep = "")
  sp  <- SPECIES[[nm]]
  res <- virtualSpecies_nd(
    dt = dt, envData = envData, mu = sp$mu, sigma = sp$sigma, rho = 0,
    pc_cols = pc_cols, bgk_prev = 1,
    pa_samplers = list(random = pa_random, mcmc = pa_mcmc),
    n_realizations = 1, max_pres = MAXP, bw = bw, pa_env_rast = envData,
    dimensions = pc_cols, chain.length = CHAIN, burnIn = BURN,
    environmental.cutof.percentile = ENVc, species.cutoff.threshold = PRESc,
    precomputed.env = env_bundle,
    compute_hypervolume = FALSE, make_plots = FALSE, compute_reference = TRUE,
    verbose = FALSE, parallel = FALSE)
  pmx <- plot_pc_matrix(res, samplers = c("random", "mcmc"),
    title = sprintf("%s — niche & pseudo-absences (lower RND, upper uniform+ at env=%.2f, pres=%.2f)",
                    sp$label, ENVc, PRESc))
  f <- file.path(out_dir, sprintf("pcmatrix_5d_%s_%s.pdf", nm, tag))
  save_experiment_figure(pmx, f, width = 9.5, height = 9.5)
  cat("saved ", f, "\n", sep = "")
}
cat("\nDone.\n")
