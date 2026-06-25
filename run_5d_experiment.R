# run_5d_experiment.R --- the 4 virtual species in 5-D (paper §2.3, 5-D arm)
# =============================================================================
# Mirrors run_2d_experiment.R but in the locked 5-D 'lean_natural' E-space
# (5 PCs ~82.5% var), per the GaussNiche 5d-niche setup. Same 2x2 factorial
# (generalist/specialist x common/rare); mu/sigma as multiples of each PC's
# background sd (generalist 1.0*sd, specialist 0.5*sd; common = KDE mode of the
# 5-D background, rare = fixed peripheral point). Samplers: RND (pa_random),
# buffer-out (pa_buffer_out, 50 km), and uniform+ (pa_mcmc). USE (paSampling) is
# 2-D-only and therefore absent; buffer-out is the non-random geographic baseline.
#
# Usage:  Rscript run_5d_experiment.R [mode] [n_workers]
#   mode = smoke -> 1 species, 2 realisations, tiny chain (validation)
#   mode = full  -> 4 species, R from N_REALIZATIONS (default 50)
#   env overrides: N_REALIZATIONS, MAX_PRES, SPECIES_CUTOFF, BUFFER_KM
# Inputs : <scratch>/GaussNiche/env5d/{final_stack_lean_natural.tif, background_5d.rds}
# Outputs: <scratch>/GaussNiche/results5d_experiment/
# Run on a compute node via sbatch/submit_5d_experiment.sh.
# =============================================================================

suppressPackageStartupMessages({
  .libPaths(c("~/R/rocker-rstudio/4.5", .libPaths()))
  library(terra); library(USE.MCMC); library(ggplot2)
  library(hypervolume); library(FNN); library(future); library(furrr)
})
source("virtualSpecies_nd_fn.R")
source("experiment_plots.R")
source("pa_buffer.R")

args      <- commandArgs(trailingOnly = TRUE)
mode      <- if (length(args) >= 1) args[1] else "full"
n_workers <- if (length(args) >= 2) as.integer(args[2]) else
  as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4"))

if (mode == "smoke") {
  N_REAL <- 2L; MAXP <- 80L; CHAIN <- 2000L; BURN <- 500L
  SPECIES_SET <- "sp1"; PAR <- FALSE
} else {
  N_REAL <- as.integer(Sys.getenv("N_REALIZATIONS", "50"))
  MAXP   <- as.integer(Sys.getenv("MAX_PRES", "300"))
  CHAIN  <- 20000L; BURN <- 1000L
  SPECIES_SET <- "all"; PAR <- TRUE
}
CUTOFF    <- as.numeric(Sys.getenv("SPECIES_CUTOFF", "0.7"))
BUFFER_KM <- as.numeric(Sys.getenv("BUFFER_KM", "50"))
cat(sprintf("mode=%s  n_real=%d  max_pres=%d  chain=%d  cutoff=%.2f  buffer_km=%g  parallel=%s  workers=%d\n",
            mode, N_REAL, MAXP, CHAIN, CUTOFF, BUFFER_KM, PAR, n_workers))

user    <- Sys.getenv("USER"); scratch <- file.path("/cluster/scratch", user)
env_dir <- file.path(scratch, "GaussNiche", "env5d")
out_dir <- file.path(scratch, "GaussNiche", "results5d_experiment")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# --- 1. Locked 5-D environment ----------------------------------------------
envData <- terra::rast(file.path(env_dir, "final_stack_lean_natural.tif"))
dt      <- readRDS(file.path(env_dir, "background_5d.rds"))
pc_cols <- paste0("PC", 1:5)
stopifnot(all(pc_cols %in% names(dt)), all(c("x", "y") %in% names(dt)))
cat(sprintf("background cells: %d\n", nrow(dt)))

pc_sd    <- vapply(dt[pc_cols], stats::sd, numeric(1))
sig_gen  <- 1.0 * pc_sd
sig_spec <- 0.5 * pc_sd

# COMMON optimum = mode of the 5-D background KDE; RARE = fixed peripheral point.
bg5 <- as.matrix(dt[, pc_cols])
mode_estimate <- function(X) {
  if (requireNamespace("ks", quietly = TRUE)) {
    est <- tryCatch(ks::kde(x = X, H = ks::Hns(X), eval.points = X)$estimate,
                    error = function(e) NULL)
    if (!is.null(est)) return(as.numeric(X[which.max(est), ]))
  }
  k <- min(200L, nrow(X) - 1L)
  as.numeric(X[which.min(FNN::knn.dist(X, k = k)[, k]), ])
}
mu_common <- mode_estimate(bg5); names(mu_common) <- pc_cols
mu_rare   <- c(PC1 = 1.55, PC2 = 1.13, PC3 = 1.4, PC4 = 0.8, PC5 = 1.2)[pc_cols]
cat("mu_common:", paste(round(mu_common, 3), collapse = ", "), "\n")
cat("mu_rare  :", paste(round(mu_rare,   3), collapse = ", "), "\n")

SPECIES <- list(
  sp1_generalist_common = list(label = "Generalist · common", mu = mu_common, sigma = sig_gen),
  sp2_specialist_common = list(label = "Specialist · common", mu = mu_common, sigma = sig_spec),
  sp3_generalist_rare   = list(label = "Generalist · rare",   mu = mu_rare,   sigma = sig_gen),
  sp4_specialist_rare   = list(label = "Specialist · rare",   mu = mu_rare,   sigma = sig_spec))
if (SPECIES_SET == "sp1") SPECIES <- SPECIES["sp1_generalist_common"]

# --- 2. One-time computations (fixed across species) -------------------------
bw <- compute_bandwidth_nd(dt, pc_cols)
cat("fixed bandwidth:", paste(round(bw, 4), collapse = ", "), "\n")

cat("\nPrecomputing MCMC environment once (PCA + env GMM) for reuse...\n")
env_bundle <- USE.MCMC::precomputeMcmcEnvironment(
  env.data.raster = envData, dimensions = pc_cols, seed.number = 123, verbose = TRUE)
env_bundle["rng_state"] <- list(NULL)   # keep name, NULL value -> restore skipped

# --- 3. Run all species ------------------------------------------------------
species_results <- list()
for (nm in names(SPECIES)) {
  cat("\n##########  ", nm, "  ##########\n", sep = "")
  sp  <- SPECIES[[nm]]
  res <- virtualSpecies_nd(
    dt = dt, envData = envData, mu = sp$mu, sigma = sp$sigma, rho = 0,
    pc_cols = pc_cols, bgk_prev = 1,
    pa_samplers = list(random = pa_random, buffer = pa_buffer_out, mcmc = pa_mcmc),
    n_realizations = N_REAL, max_pres = MAXP,
    bw = bw, pa_env_rast = envData,
    dimensions = pc_cols, chain.length = CHAIN, burnIn = BURN,
    environmental.cutof.percentile = 0.001,
    species.cutoff.threshold = CUTOFF,
    precomputed.env = env_bundle, buffer_km = BUFFER_KM,
    verbose = TRUE, parallel = PAR, n_workers = n_workers)
  saveRDS(res, file.path(out_dir, paste0("exp_", nm, "_", mode, ".rds")))
  species_results[[nm]] <- res
}

# --- 4. Collect metrics, CSV + lightweight summary ---------------------------
species_labels <- vapply(SPECIES, `[[`, character(1), "label")
metrics <- collect_experiment_metrics(species_results, species_labels = species_labels)
write.csv(metrics, file.path(out_dir, paste0("metrics_5d_", mode, ".csv")), row.names = FALSE)

summary_bundle <- list(
  metrics = metrics, species_labels = species_labels, pc_cols = pc_cols, bw = bw,
  prevalence = vapply(species_results, function(r) mean(r$background$pa), numeric(1)),
  mode = mode, n_realizations = N_REAL, max_pres = MAXP,
  species_cutoff = CUTOFF, buffer_km = BUFFER_KM)
saveRDS(summary_bundle, file.path(out_dir, paste0("summary_5d_", mode, ".rds")))

# --- 5. Figures (vector PDF + PNG; capped height since 5-D has 7 metric rows) -
np <- experiment_panel_count(metrics)
p_box <- plot_experiment_boxplots(metrics, jitter = (N_REAL <= 30),
  title = "Sampler comparison across virtual species (5-D)")
for (ext in c("pdf", "png"))
  save_experiment_figure(p_box, file.path(out_dir, paste0("boxplots_metrics_5d_", mode, ".", ext)),
                         width = 7.0, height = min(9.5, 1.6 * np))
for (mt in c("overlap", "coverage", "trueabs")) {
  metric_arg <- if (mt == "trueabs") "prop_true_abs" else mt
  ph <- if (mt == "coverage") min(9.5, 1.5 * length(grep("^rel_cov_", names(metrics)))) else 3.2
  pm <- plot_experiment_metric(metrics, metric = metric_arg, jitter = (N_REAL <= 30))
  save_experiment_figure(pm, file.path(out_dir, paste0("box_", mt, "_5d_", mode, ".pdf")),
                         width = 7.0, height = ph)
}

# --- 6. Console summary -------------------------------------------------------
cat("\n================ SUMMARY ================\n")
for (nm in names(species_results)) {
  cat(sprintf("\n%s  (prevalence %.4f)\n", nm, mean(species_results[[nm]]$background$pa)))
  for (s in names(species_results[[nm]]$samplers)) {
    m  <- species_results[[nm]]$samplers[[s]]$metrics
    ok <- m[m$status == "ok", , drop = FALSE]
    cat(sprintf("  %-8s ok=%d/%d  overlap(med)=%s  trueabs(med)=%s\n", s, nrow(ok), nrow(m),
                if (nrow(ok)) sprintf("%.4f", median(ok$overlap, na.rm = TRUE)) else "NA",
                if (nrow(ok)) sprintf("%.3f", median(ok$prop_true_abs, na.rm = TRUE)) else "NA"))
  }
}
cat("\nSaved to ", out_dir, "\n", sep = "")
