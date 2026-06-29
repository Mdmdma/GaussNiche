# profile_samplers.R --- decompose the computational cost of `uniform` (USE,
# paSampling) vs `uniform+` (paSamplingMcmc) into MODEL-FITTING vs SAMPLING/CHAIN.
# =============================================================================
# Goal (for the paper's "Computational performance" results): show that for
# uniform+ the Markov chain itself is cheap, and the cost is dominated by FITTING
# the models (the environmental + presence Gaussian mixtures). The chain is fast
# precisely because, once fitted, the GMMs evaluate quickly (C++ inner loop).
#
# We profile, single-threaded:
#   uniform+ (paSamplingMcmc), 2-D and 5-D, decomposed into
#       PCA  /  env-GMM fit  /  env-GMM density eval  /  presence-GMM fit  /
#       Markov chain (per-step, from a chain-length regression)  /  remap+finalise
#   uniform  (paSampling, 2-D only), decomposed into
#       PCA  /  KDE bandwidth (Hpi)  /  KDE fit (kde)  /  grid sampling  /
#       optimRes grid-resolution search (the dominant cost of using USE properly)
#
# Run on a compute node via sbatch/submit_profile_samplers.sh. NOT a login-node
# task (fits real GMM/KDE models). Outputs a table to stdout + an rds/csv under
# <scratch>/GaussNiche/results_profile.
# =============================================================================

suppressPackageStartupMessages({
  .libPaths(c("~/R/rocker-rstudio/4.5", .libPaths()))
  library(terra); library(USE.MCMC); library(sf); library(mclust); library(ks); library(FNN)
})

REPS  <- as.integer(Sys.getenv("PROFILE_REPS", "3"))
et    <- function(expr) tryCatch(as.numeric(system.time(expr)["elapsed"]), error = function(e) NA_real_)
medet <- function(expr_fun, reps = REPS) stats::median(vapply(seq_len(reps), function(i) expr_fun(), numeric(1)), na.rm = TRUE)
say   <- function(...) cat(sprintf(...), "\n")

user    <- Sys.getenv("USER"); scratch <- file.path("/cluster/scratch", user)
out_dir <- file.path(scratch, "GaussNiche", "results_profile"); dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
res     <- list()  # accumulate named timings (seconds)

# Build a generalist-common virtual species in PC space and draw a capped
# presence set (matches the paper's max_pres = 500 convention).
make_pres <- function(dt, pc_cols, cap = 500L, seed = 1L) {
  mu    <- setNames(as.numeric(colMeans(dt[pc_cols])), pc_cols)
  sg    <- vapply(dt[pc_cols], stats::sd, numeric(1))
  Sig   <- diag(sg) %*% diag(length(sg)) %*% diag(sg)
  peak  <- mvtnorm::dmvnorm(matrix(mu, 1), mu, Sig)
  suit  <- mvtnorm::dmvnorm(as.matrix(dt[pc_cols]), mu, Sig) / peak
  set.seed(seed); pa <- stats::rbinom(nrow(dt), 1L, suit)
  pc <- dt[pa == 1L, , drop = FALSE]
  if (nrow(pc) > cap) { set.seed(seed); pc <- pc[sample(nrow(pc), cap), ] }
  pc
}

# Regression-based split of the Markov chain cost: t(chain) = fixed + k*chain.
# Returns per-1000-step cost and the chain cost at 10k steps. env is precomputed
# (passed as a bundle) so only chain.length varies between calls.
chain_regression <- function(env.rast, pres_sf, dims, env_bundle, lengths = c(5000, 10000, 20000, 40000)) {
  ts <- vapply(lengths, function(L) medet(function()
    et(USE.MCMC::paSamplingMcmc(env.rast = NULL, pres = pres_sf, dimensions = dims,
        chain.length = L, burnIn = 1000, n.samples = 300, seed.number = 42,
        species.cutoff.threshold = 0.95, environmental.cutof.percentile = 0.001,
        precomputed.env = env_bundle, verbose = FALSE)), reps = max(2, REPS - 1)), numeric(1))
  fit <- stats::lm(ts ~ lengths)
  k1000 <- unname(stats::coef(fit)[2]) * 1000
  list(lengths = lengths, times = ts, per_1000 = k1000,
       intercept = unname(stats::coef(fit)[1]), chain_10k = unname(stats::coef(fit)[2]) * 10000)
}

# =====================================================================
# PART 1 --- 2-D environment (Worldclim_tmp): uniform vs uniform+
# =====================================================================
say("\n================ PART 1: 2-D (Worldclim_tmp) ================")
env2  <- terra::rast(USE.MCMC::Worldclim_tmp, type = "xyz")
dims2 <- c("PC1", "PC2")
rpc2  <- rastPCA(env2, stand = TRUE)
dt2   <- terra::as.data.frame(rpc2$PCs[[dims2]], xy = TRUE, na.rm = TRUE)
pres_pc2 <- make_pres(dt2, dims2, cap = 500L, seed = 1L)
pres2 <- sf::st_as_sf(pres_pc2[, c("x", "y", dims2)], coords = c("x", "y"), crs = 4326)
say("2-D background cells: %d   presence points: %d", nrow(dt2), nrow(pres2))

## ---- uniform+ (paSamplingMcmc) 2-D decomposition ----
say("\n-- uniform+ (paSamplingMcmc), 2-D --")
res[["u+_2d_PCA"]]          <- medet(function() et(rastPCA(env2, stand = TRUE)))
env_bundle2                 <- USE.MCMC::precomputeMcmcEnvironment(env.rast = env2, dimensions = dims2, seed.number = 42, verbose = FALSE)
res[["u+_2d_precompute"]]   <- medet(function() et(USE.MCMC::precomputeMcmcEnvironment(env.rast = env2, dimensions = dims2, seed.number = 42, verbose = FALSE)))
# env-GMM fit (on the <=2000-row subsample the function uses) + density eval over full bg
env_pc_full2 <- env_bundle2$env.data.cleaned
sub2 <- env_pc_full2[sample(nrow(env_pc_full2), min(2000L, nrow(env_pc_full2))), , drop = FALSE]
res[["u+_2d_envGMM_fit"]]   <- medet(function() et(mclust::densityMclust(sub2, plot = FALSE, verbose = FALSE)))
env_model2 <- mclust::densityMclust(sub2, plot = FALSE, verbose = FALSE)
res[["u+_2d_envGMM_eval"]]  <- medet(function() et(stats::predict(env_model2, env_pc_full2)))
# presence-GMM fit
pres_pc_sc2 <- as.data.frame(terra::extract(terra::unwrap(env_bundle2$pcs_packed), pres2, ID = FALSE))[, dims2, drop = FALSE]
res[["u+_2d_presGMM_fit"]]  <- medet(function() et(mclust::densityMclust(pres_pc_sc2, plot = FALSE, verbose = FALSE)))
# Markov chain (per-step regression) + cold one-shot total
cr2 <- chain_regression(env2, pres2, dims2, env_bundle2)
res[["u+_2d_chain_per1k"]]  <- cr2$per_1000
res[["u+_2d_chain_10k"]]    <- cr2$chain_10k
res[["u+_2d_total_cold10k"]]<- medet(function() et(USE.MCMC::paSamplingMcmc(env.rast = env2, pres = pres2, dimensions = dims2,
                                 chain.length = 10000, burnIn = 1000, n.samples = 300, seed.number = 42,
                                 species.cutoff.threshold = 0.95, environmental.cutof.percentile = 0.001, verbose = FALSE)), reps = max(2, REPS - 1))
# pure-R vs C++ chain at 10k (warm env) to show the chain is cheap with the C++ inner loop
res[["u+_2d_total_warm10k_cpp"]] <- medet(function() et(USE.MCMC::paSamplingMcmc(env.rast = NULL, pres = pres2, dimensions = dims2,
                                 chain.length = 10000, burnIn = 1000, n.samples = 300, seed.number = 42,
                                 species.cutoff.threshold = 0.95, environmental.cutof.percentile = 0.001,
                                 precomputed.env = env_bundle2, engine = "cpp", verbose = FALSE)), reps = max(2, REPS - 1))
res[["u+_2d_total_warm10k_R"]]   <- medet(function() et(USE.MCMC::paSamplingMcmc(env.rast = NULL, pres = pres2, dimensions = dims2,
                                 chain.length = 10000, burnIn = 1000, n.samples = 300, seed.number = 42,
                                 species.cutoff.threshold = 0.95, environmental.cutof.percentile = 0.001,
                                 precomputed.env = env_bundle2, engine = "R", verbose = FALSE)), reps = max(2, REPS - 1))

## ---- uniform (USE, paSampling) 2-D decomposition ----
say("\n-- uniform (USE, paSampling), 2-D --")
res[["u_2d_PCA"]] <- res[["u+_2d_PCA"]]   # same rastPCA
# KDE model fit (bandwidth + kde) on presence PC scores -- USE's "model"
res[["u_2d_KDE_Hpi"]] <- medet(function() et(ks::Hpi(x = pres_pc_sc2)))
Hkde <- ks::Hpi(x = pres_pc_sc2)
res[["u_2d_KDE_fit"]] <- medet(function() et(ks::kde(pres_pc_sc2, eval.points = pres_pc_sc2, h = Hkde)))
# whole paSampling (one-shot, fixed grid.res; no optimRes)
n_tr <- max(5L, ceiling(nrow(pres2) * 3L / max(10L^2L, 10L)))
res[["u_2d_total_grid10"]] <- medet(function() et(suppressMessages(USE.MCMC::paSampling(env.rast = env2, pres = pres2,
                                 thres = 0.75, grid.res = 10, n.tr = n_tr, verbose = FALSE))), reps = max(2, REPS - 1))
# optimRes grid-resolution search (the slow step a USE user runs to choose grid.res)
env_pc_sf2 <- sf::st_as_sf(stats::na.omit(terra::as.data.frame(rpc2$PCs[[dims2]])), coords = dims2)
res[["u_2d_optimRes"]] <- et(suppressMessages(USE.MCMC::optimRes(sdf = env_pc_sf2, grid.res = seq(2, 30, by = 2), cr = 1, showOpt = FALSE)))

# =====================================================================
# PART 2 --- 5-D environment (env5d): uniform+ scaling (USE not applicable)
# =====================================================================
say("\n================ PART 2: 5-D (env5d lean_natural) ================")
env5_path <- file.path(scratch, "GaussNiche", "env5d", "final_stack_lean_natural.tif")
bg5_path  <- file.path(scratch, "GaussNiche", "env5d", "background_5d.rds")
if (file.exists(env5_path) && file.exists(bg5_path)) {
  env5  <- terra::rast(env5_path)
  dims5 <- paste0("PC", 1:5)
  dt5   <- readRDS(bg5_path)
  pres_pc5 <- make_pres(dt5, dims5, cap = 500L, seed = 1L)
  pres5 <- sf::st_as_sf(pres_pc5[, c("x", "y", dims5)], coords = c("x", "y"), crs = 4326)
  say("5-D background cells: %d   presence points: %d", nrow(dt5), nrow(pres5))

  env_bundle5 <- USE.MCMC::precomputeMcmcEnvironment(env.rast = env5, dimensions = dims5, seed.number = 42, verbose = FALSE)
  res[["u+_5d_precompute"]]  <- medet(function() et(USE.MCMC::precomputeMcmcEnvironment(env.rast = env5, dimensions = dims5, seed.number = 42, verbose = FALSE)), reps = max(2, REPS - 1))
  env_pc_full5 <- env_bundle5$env.data.cleaned
  sub5 <- env_pc_full5[sample(nrow(env_pc_full5), min(2000L, nrow(env_pc_full5))), , drop = FALSE]
  res[["u+_5d_envGMM_fit"]]  <- medet(function() et(mclust::densityMclust(sub5, plot = FALSE, verbose = FALSE)))
  env_model5 <- mclust::densityMclust(sub5, plot = FALSE, verbose = FALSE)
  res[["u+_5d_envGMM_eval"]] <- medet(function() et(stats::predict(env_model5, env_pc_full5)))
  pres_pc_sc5 <- as.data.frame(terra::extract(terra::unwrap(env_bundle5$pcs_packed), pres5, ID = FALSE))[, dims5, drop = FALSE]
  res[["u+_5d_presGMM_fit"]] <- medet(function() et(mclust::densityMclust(pres_pc_sc5, plot = FALSE, verbose = FALSE)))
  cr5 <- chain_regression(env5, pres5, dims5, env_bundle5)
  res[["u+_5d_chain_per1k"]] <- cr5$per_1000
  res[["u+_5d_chain_10k"]]   <- cr5$chain_10k
  res[["u+_5d_total_cold10k"]] <- medet(function() et(USE.MCMC::paSamplingMcmc(env.rast = env5, pres = pres5, dimensions = dims5,
                                   chain.length = 10000, burnIn = 1000, n.samples = 300, seed.number = 42,
                                   species.cutoff.threshold = 0.95, environmental.cutof.percentile = 0.001, verbose = FALSE)), reps = max(2, REPS - 1))
} else {
  say("env5d artifacts not found (%s); skipping 5-D part.", env5_path)
}

# =====================================================================
# REPORT
# =====================================================================
say("\n================ TIMINGS (median of %d, seconds) ================", REPS)
ord <- c(
  "u+_2d_PCA", "u+_2d_precompute", "u+_2d_envGMM_fit", "u+_2d_envGMM_eval", "u+_2d_presGMM_fit",
  "u+_2d_chain_per1k", "u+_2d_chain_10k", "u+_2d_total_warm10k_cpp", "u+_2d_total_warm10k_R", "u+_2d_total_cold10k",
  "u_2d_PCA", "u_2d_KDE_Hpi", "u_2d_KDE_fit", "u_2d_total_grid10", "u_2d_optimRes",
  "u+_5d_precompute", "u+_5d_envGMM_fit", "u+_5d_envGMM_eval", "u+_5d_presGMM_fit",
  "u+_5d_chain_per1k", "u+_5d_chain_10k", "u+_5d_total_cold10k")
for (k in ord) if (!is.null(res[[k]])) say("  %-26s %8.3f", k, res[[k]])

say("\n-- headline ratios --")
fit2  <- sum(c(res[["u+_2d_precompute"]], res[["u+_2d_presGMM_fit"]]), na.rm = TRUE)
if (!is.null(res[["u+_2d_chain_10k"]]))
  say("  2-D uniform+: model fitting %.2fs vs chain(10k) %.3fs  -> fitting is %.1fx the chain",
      fit2, res[["u+_2d_chain_10k"]], fit2 / res[["u+_2d_chain_10k"]])
if (!is.null(res[["u+_5d_chain_10k"]])) {
  fit5 <- sum(c(res[["u+_5d_precompute"]], res[["u+_5d_presGMM_fit"]]), na.rm = TRUE)
  say("  5-D uniform+: model fitting %.2fs vs chain(10k) %.3fs  -> fitting is %.1fx the chain",
      fit5, res[["u+_5d_chain_10k"]], fit5 / res[["u+_5d_chain_10k"]])
}
if (!is.null(res[["u+_2d_total_warm10k_R"]]) && !is.null(res[["u+_2d_total_warm10k_cpp"]]))
  say("  2-D chain engine: C++ %.3fs vs pure-R %.3fs", res[["u+_2d_total_warm10k_cpp"]], res[["u+_2d_total_warm10k_R"]])

saveRDS(res, file.path(out_dir, "profile_samplers.rds"))
write.csv(data.frame(step = names(res), seconds = unlist(res)), file.path(out_dir, "profile_samplers.csv"), row.names = FALSE)
say("\nsaved -> %s", file.path(out_dir, "profile_samplers.{rds,csv}"))
