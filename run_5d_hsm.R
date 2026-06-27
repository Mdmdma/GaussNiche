# run_5d_hsm.R --- downstream HSM evaluation arm of the 5-D experiment
# =============================================================================
# The "downstream" half of the USE-paper experiment (Da Re et al. 2023), in
# GaussNiche's locked 5-D 'lean_natural' E-space. Mirrors run_5d_experiment.R
# (same 4 species, same 3 samplers RND/buffer-out/uniform+, same R=50 Bernoulli
# realisations, IDENTICAL seeds) but additionally fits the five USE model
# families (GLM/GAM/RF/BRT/Maxent, Java-free) to each realisation's presence /
# pseudo-absence set via the sdm_hook, scoring predictive accuracy both
# internally (held-out 30% test) and against the KNOWN suitability (USE's
# rmse.distribution / cor.distribution). Two predictor sets, like USE's 5-var vs
# 19-var arms: the 5 PCs and the 12 raw lean_natural layers.
#
# The hook reuses the PA set already drawn for the diagnostics, so the expensive
# pa_mcmc draw is NOT re-paid; identical seeds mean the standard overlap/coverage
# this produces reproduce metrics_5d_full.csv bit-for-bit (a built-in check).
#
# Usage:  Rscript run_5d_hsm.R [mode] [n_workers]
#   mode = smoke -> 1 species, 2 realisations, all 5 algos, PARALLEL (2 workers)
#                   -- the parallel smoke is REQUIRED: it is the only thing that
#                   exercises furrr worker package loading before the full run.
#   mode = full  -> 4 species, R from N_REALIZATIONS (default 50)
#   env overrides: N_REALIZATIONS, MAX_PRES, SPECIES_CUTOFF, BUFFER_KM, MIN_PRES
# Inputs : <scratch>/GaussNiche/env5d/{final_stack_lean_natural.tif, background_5d.rds}
# Outputs: <scratch>/GaussNiche/results5d_hsm/
#            hsm_metrics_5d_<mode>.csv     tidy stats (sp x sampler x real x set x algo)
#            metrics_5d_<mode>.csv         standard diagnostics (non-perturbation check)
#            hsm_dunn_5d_<mode>.csv        Dunn pairwise (uniform+ vs others)
#            hsm_violin_*_<mode>.pdf/.png  per-(predictor_set[,species]) violins
#            summary_5d_hsm_<mode>.rds     lightweight bundle (no 90 MB result rds)
# Run on a compute node via sbatch/submit_5d_hsm.sh.
# =============================================================================

suppressPackageStartupMessages({
  .libPaths(c("~/R/rocker-rstudio/4.5", .libPaths()))
  library(terra); library(USE.MCMC); library(ggplot2)
  library(hypervolume); library(FNN); library(future); library(furrr)
})
source("virtualSpecies_nd_fn.R")
source("experiment_plots.R")
source("hsm_eval.R")
source("hsm_plots.R")
source("pa_buffer.R")
HSM_EVAL_PATH <- normalizePath("hsm_eval.R")   # for the furrr-worker source guard

args      <- commandArgs(trailingOnly = TRUE)
mode      <- if (length(args) >= 1) args[1] else "full"
n_workers <- if (length(args) >= 2) as.integer(args[2]) else
  as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4"))

if (mode == "smoke") {
  N_REAL <- 2L; MAXP <- 120L; CHAIN <- 2000L; BURN <- 500L
  SPECIES_SET <- "sp1"; PAR <- TRUE                       # parallel smoke (required)
  if (is.na(n_workers) || n_workers < 2L) n_workers <- 2L
} else {
  N_REAL <- as.integer(Sys.getenv("N_REALIZATIONS", "50"))
  MAXP   <- as.integer(Sys.getenv("MAX_PRES", "300"))
  CHAIN  <- 20000L; BURN <- 1000L
  SPECIES_SET <- "all"; PAR <- TRUE
}
CUTOFF    <- as.numeric(Sys.getenv("SPECIES_CUTOFF", "0.7"))
BUFFER_KM <- as.numeric(Sys.getenv("BUFFER_KM", "50"))
MIN_PRES  <- as.integer(Sys.getenv("MIN_PRES", "12"))
ALGOS     <- HSM_ALGORITHMS                               # glm, gam, rf, brt, maxent
cat(sprintf("mode=%s n_real=%d max_pres=%d chain=%d cutoff=%.2f buffer_km=%g min_pres=%d parallel=%s workers=%d\n",
            mode, N_REAL, MAXP, CHAIN, CUTOFF, BUFFER_KM, MIN_PRES, PAR, n_workers))

user    <- Sys.getenv("USER"); scratch <- file.path("/cluster/scratch", user)
env_dir <- file.path(scratch, "GaussNiche", "env5d")
out_dir <- file.path(scratch, "GaussNiche", "results5d_hsm")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# --- 1. Locked 5-D environment + raw-layer augmentation ----------------------
envData <- terra::rast(file.path(env_dir, "final_stack_lean_natural.tif"))
dt      <- readRDS(file.path(env_dir, "background_5d.rds"))
pc_cols <- paste0("PC", 1:5)
stopifnot(all(pc_cols %in% names(dt)), all(c("x", "y") %in% names(dt)))

raw_cols <- names(envData)                                 # the 12 lean_natural layers
raw_vals <- terra::extract(envData, as.matrix(dt[, c("x", "y")]))
raw_vals <- as.data.frame(raw_vals)[, raw_cols, drop = FALSE]
ok_raw   <- stats::complete.cases(raw_vals)
if (!all(ok_raw)) {
  cat(sprintf("dropping %d background cells with NA raw-layer values\n", sum(!ok_raw)))
  dt <- dt[ok_raw, , drop = FALSE]; raw_vals <- raw_vals[ok_raw, , drop = FALSE]
}
dt <- cbind(dt, raw_vals); rownames(dt) <- NULL
predictor_sets <- list(pc5 = pc_cols, raw12 = raw_cols)
cat(sprintf("background cells: %d   predictors: pc5=%d raw12=%d\n",
            nrow(dt), length(pc_cols), length(raw_cols)))

pc_sd    <- vapply(dt[pc_cols], stats::sd, numeric(1))
sig_gen  <- 1.0 * pc_sd
sig_spec <- 0.5 * pc_sd

# COMMON optimum = mode of the 5-D background KDE; RARE = fixed peripheral point
# (identical construction to run_5d_experiment.R so the species match).
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
  env.rast = envData, dimensions = pc_cols, seed.number = 123, verbose = TRUE)
env_bundle["rng_state"] <- list(NULL)

LIBS <- .libPaths()   # captured so the hook can restore it inside furrr workers

# --- 3. Run all species (with the SDM hook) ----------------------------------
species_hsm     <- list()
species_metrics <- list()
prevalence      <- numeric(0)
sp_idx <- 0L
for (nm in names(SPECIES)) {
  sp_idx <- sp_idx + 1L
  cat("\n##########  ", nm, "  ##########\n", sep = "")
  sp   <- SPECIES[[nm]]
  hook <- make_sdm_hook(species = nm, predictor_sets = predictor_sets,
                        algorithms = ALGOS, test_frac = 0.3, min_pres = MIN_PRES,
                        seed = 1000L + sp_idx * 100L, libpaths = LIBS,
                        source_file = HSM_EVAL_PATH)
  res <- virtualSpecies_nd(
    dt = dt, envData = envData, mu = sp$mu, sigma = sp$sigma, rho = 0,
    pc_cols = pc_cols, bgk_prev = 1,
    pa_samplers = list(random = pa_random, buffer = pa_buffer_out, mcmc = pa_mcmc),
    n_realizations = N_REAL, max_pres = MAXP,
    bw = bw, pa_env_rast = envData,
    sdm_hook = hook,
    dimensions = pc_cols, chain.length = CHAIN, burnIn = BURN,
    environmental.cutof.percentile = 0.001,
    species.cutoff.threshold = CUTOFF,
    precomputed.env = env_bundle, buffer_km = BUFFER_KM,
    verbose = TRUE, parallel = PAR, n_workers = n_workers)
  if (!is.null(res$hsm)) {
    res$hsm$species_label <- sp$label
    species_hsm[[nm]] <- res$hsm
  }
  m <- collect_experiment_metrics(setNames(list(res), nm),
                                  species_labels = setNames(sp$label, nm))
  species_metrics[[nm]] <- m
  prevalence[nm] <- mean(res$background$pa)
  rm(res); gc(verbose = FALSE)
}

# --- 4. Persist tidy stats (HSM + standard diagnostics) ----------------------
hsm <- do.call(rbind, species_hsm); rownames(hsm) <- NULL
write.csv(hsm, file.path(out_dir, paste0("hsm_metrics_5d_", mode, ".csv")), row.names = FALSE)
metrics <- do.call(rbind, species_metrics); rownames(metrics) <- NULL
write.csv(metrics, file.path(out_dir, paste0("metrics_5d_", mode, ".csv")), row.names = FALSE)

cat(sprintf("\nHSM rows: %d  (ok=%d, skip_few_pres=%d, fit_failed=%d)\n",
            nrow(hsm), sum(hsm$status == "ok"), sum(hsm$status == "skip_few_pres"),
            sum(hsm$status == "fit_failed")))

species_labels <- vapply(SPECIES, `[[`, character(1), "label")

# --- 5. Dunn pairwise test (uniform+ vs random / buffer-out) -----------------
dunn <- tryCatch(hsm_dunn_test(hsm), error = function(e) {
  cat("Dunn test skipped:", conditionMessage(e), "\n"); NULL })
if (!is.null(dunn))
  write.csv(dunn, file.path(out_dir, paste0("hsm_dunn_5d_", mode, ".csv")), row.names = FALSE)

# --- 6. Figures (per predictor_set: combined-across-species + per species) ----
for (ps in names(predictor_sets)) {
  comb <- tryCatch(plot_hsm_violin(hsm, predictor_set = ps, species_labels = species_labels,
            title = sprintf("Downstream HSM accuracy across species (%s predictors)", ps)),
            error = function(e) { cat("violin", ps, "skipped:", conditionMessage(e), "\n"); NULL })
  if (!is.null(comb)) {
    w <- 3 + 2.4 * length(unique(hsm$species))
    save_hsm_figure(comb, file.path(out_dir, paste0("hsm_violin_", ps, "_", mode, ".pdf")),
                    width = min(20, w), height = 8)
  }
  for (nm in names(species_hsm)) {
    p <- tryCatch(plot_hsm_violin(hsm, predictor_set = ps, species = nm,
           species_labels = species_labels,
           title = sprintf("%s — HSM accuracy (%s predictors)", species_labels[[nm]], ps)),
           error = function(e) NULL)
    if (!is.null(p))
      save_hsm_figure(p, file.path(out_dir, paste0("hsm_violin_", ps, "_", nm, "_", mode, ".pdf")),
                      width = 9, height = 7)
  }
}
if (!is.null(dunn)) {
  ds <- tryCatch(plot_hsm_dunn_summary(dunn), error = function(e) NULL)
  if (!is.null(ds))
    save_hsm_figure(ds, file.path(out_dir, paste0("hsm_dunn_summary_5d_", mode, ".pdf")),
                    width = 11, height = 3.5)
}

# --- 7. Lightweight summary bundle (NO 90 MB per-species result rds) ---------
summary_bundle <- list(
  hsm = hsm, metrics = metrics, dunn = dunn,
  medians = tryCatch(summarize_hsm(hsm), error = function(e) NULL),
  species_labels = species_labels, prevalence = prevalence,
  predictor_sets = predictor_sets, algorithms = ALGOS, pc_cols = pc_cols, bw = bw,
  mode = mode, n_realizations = N_REAL, max_pres = MAXP, min_pres = MIN_PRES,
  species_cutoff = CUTOFF, buffer_km = BUFFER_KM)
saveRDS(summary_bundle, file.path(out_dir, paste0("summary_5d_hsm_", mode, ".rds")))

# --- 8. Console summary -------------------------------------------------------
cat("\n================ HSM SUMMARY (median over realisations) ================\n")
med <- summary_bundle$medians
if (!is.null(med)) {
  for (ps in names(predictor_sets)) {
    cat("\n--- predictor set:", ps, "---\n")
    sub <- med[med$predictor_set == ps, , drop = FALSE]
    for (nm in unique(sub$species)) {
      cat(sprintf(" %s\n", nm))
      ss <- sub[sub$species == nm, , drop = FALSE]
      for (a in unique(ss$algorithm)) {
        aa <- ss[ss$algorithm == a, , drop = FALSE]
        cat(sprintf("  %-7s  ", a))
        cat(paste(sprintf("%s:AUC=%.3f/CORt=%.3f", aa$sampler, aa$auc, aa$cor_truth),
                  collapse = "  "), "\n")
      }
    }
  }
}
cat("\nSaved to ", out_dir, "\n", sep = "")
