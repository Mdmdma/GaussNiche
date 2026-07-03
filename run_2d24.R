# run_2d24.R --- reevaluate §2.3 + downstream HSM on PC1xPC2 of the 24-layer env
# =============================================================================
# A 2-D analysis in the LEADING TWO PCs of the full 24-layer candidate PCA
# (env5d_all24). Uses the dimension-agnostic nd engine with pc_cols = PC1,PC2 so
# ONE pass yields BOTH the §2.3 sampler diagnostics (overlap / per-axis coverage /
# prop-true-absence) AND the downstream-HSM metrics (via the sdm_hook). Because we
# are in 2-D, the 2-D-only USE sampler (pa_uniform) rejoins the comparison next to
# RND / buffer-out / uniform+ (it is absent from the 5-D arm). random + buffer-out
# are drawn once and feed both halves.
#
# Species = the paper's 2-D Table-1 design (FIXED PC-unit mu/sigma), reevaluated
# as-is on the 24-layer PC1xPC2 plane (per user: keep Table-1). PC scales differ
# from WorldClim-2D, so prevalences shift -- reported in the console summary.
#
# Usage:  Rscript run_2d24.R [mode] [n_workers]
#   mode=smoke -> 1 sp, 2 real, tiny chain, serial (validation)
#   mode=full  -> 4 sp, R from N_REALIZATIONS (default 10)
#   env overrides: N_REALIZATIONS, MAX_PRES, SPECIES_CUTOFF, GRID_RES, MIN_PRES, BUFFER_KM
#   path overrides: GN_ENV_DIR, GN_STACK_TIF, GN_BG_RDS, GN_OUT_DIR
# Inputs : <scratch>/GaussNiche/env5d_all24/{final_stack_all24.tif, background_all24.rds}
# Outputs: <scratch>/GaussNiche/results2d24/
# =============================================================================

suppressPackageStartupMessages({
  .libPaths(c("~/R/rocker-rstudio/4.5", .libPaths()))
  library(terra); library(USE.MCMC); library(ggplot2)
  library(sf); library(hypervolume); library(FNN); library(MASS)
  library(future); library(furrr)
})
source("virtualSpecies_fn.R")      # pa_uniform (2-D USE sampler)
source("virtualSpecies_nd_fn.R")   # nd engine + nd pa_random/pa_mcmc (override the 2-D ones)
source("experiment_plots.R")
source("hsm_eval.R")
source("hsm_plots.R")
source("pa_buffer.R")
HSM_EVAL_PATH <- normalizePath("hsm_eval.R")   # furrr-worker source guard

args      <- commandArgs(trailingOnly = TRUE)
mode      <- if (length(args) >= 1) args[1] else "full"
n_workers <- if (length(args) >= 2) as.integer(args[2]) else
  as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4"))

if (mode == "smoke") {
  N_REAL <- 2L; MAXP <- 120L; CHAIN <- 2000L; BURN <- 500L
  SPECIES_SET <- "sp1"; PAR <- TRUE
  if (is.na(n_workers) || n_workers < 2L) n_workers <- 2L
} else {
  N_REAL <- as.integer(Sys.getenv("N_REALIZATIONS", "10"))
  MAXP   <- as.integer(Sys.getenv("MAX_PRES", "300"))
  CHAIN  <- 20000L; BURN <- 1000L
  SPECIES_SET <- "all"; PAR <- TRUE
}
CUTOFF    <- as.numeric(Sys.getenv("SPECIES_CUTOFF", "0.1"))   # 2-D convention
GRID_RES  <- as.integer(Sys.getenv("GRID_RES", "10"))          # USE grid (fixed; not optimRes)
BUFFER_KM <- as.numeric(Sys.getenv("BUFFER_KM", "50"))
MIN_PRES  <- as.integer(Sys.getenv("MIN_PRES", "12"))
ALGOS     <- HSM_ALGORITHMS
cat(sprintf("mode=%s n_real=%d max_pres=%d chain=%d cutoff=%.2f grid.res=%d min_pres=%d parallel=%s workers=%d\n",
            mode, N_REAL, MAXP, CHAIN, CUTOFF, GRID_RES, MIN_PRES, PAR, n_workers))

user    <- Sys.getenv("USER"); scratch <- file.path("/cluster/scratch", user)
env_dir <- Sys.getenv("GN_ENV_DIR",  file.path(scratch, "GaussNiche", "env5d_all24"))
out_dir <- Sys.getenv("GN_OUT_DIR",  file.path(scratch, "GaussNiche", "results2d24"))
stack_f <- Sys.getenv("GN_STACK_TIF", "final_stack_all24.tif")
bg_f    <- Sys.getenv("GN_BG_RDS",    "background_all24.rds")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# --- 1. Env + PC1,PC2 background + raw-layer augmentation ---------------------
envData <- terra::rast(file.path(env_dir, stack_f))
dt      <- readRDS(file.path(env_dir, bg_f))
pc_cols <- c("PC1", "PC2")
stopifnot(all(pc_cols %in% names(dt)), all(c("x", "y") %in% names(dt)))
dt <- dt[, c("x", "y", pc_cols)]      # 2-D analysis: keep only the leading two PCs

raw_cols <- names(envData)            # the 24 candidate layers
raw_vals <- as.data.frame(terra::extract(envData, as.matrix(dt[, c("x", "y")])))[, raw_cols, drop = FALSE]
ok_raw   <- stats::complete.cases(raw_vals)
if (!all(ok_raw)) { cat(sprintf("dropping %d cells with NA raw values\n", sum(!ok_raw)))
  dt <- dt[ok_raw, ]; raw_vals <- raw_vals[ok_raw, ] }
dt <- cbind(dt, raw_vals); rownames(dt) <- NULL
predictor_sets <- list(pc2 = pc_cols, raw24 = raw_cols)
cat(sprintf("background cells: %d   predictors: pc2=%d raw24=%d\n",
            nrow(dt), length(pc_cols), length(raw_cols)))

# --- 2. Paper Table-1 species (fixed PC-unit mu/sigma), rho = 0 ---------------
SPECIES <- list(
  sp1_generalist_common = list(label = "Generalist · common", mu = c(PC1 = 0, PC2 = 0), sigma = c(PC1 = 3, PC2 = 2)),
  sp2_specialist_common = list(label = "Specialist · common", mu = c(PC1 = 0, PC2 = 0), sigma = c(PC1 = 1, PC2 = 0.5)),
  sp3_generalist_rare   = list(label = "Generalist · rare",   mu = c(PC1 = 4, PC2 = 1), sigma = c(PC1 = 3, PC2 = 2)),
  sp4_specialist_rare   = list(label = "Specialist · rare",   mu = c(PC1 = 4, PC2 = 1), sigma = c(PC1 = 1, PC2 = 0.5)))
if (SPECIES_SET == "sp1") SPECIES <- SPECIES["sp1_generalist_common"]

# --- 3. One-time computations ------------------------------------------------
bw <- compute_bandwidth_nd(dt, pc_cols)
cat("fixed bandwidth:", paste(round(bw, 4), collapse = ", "), "\n")
cat("\nPrecomputing MCMC environment once (PCA + env GMM) for reuse...\n")
env_bundle <- USE.MCMC::precomputeMcmcEnvironment(
  env.rast = envData, dimensions = pc_cols, seed.number = 123, verbose = TRUE)
env_bundle["rng_state"] <- list(NULL)
LIBS <- .libPaths()

# --- 4. Run all species (4 samplers incl USE, + the SDM hook) ----------------
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
    pa_samplers = list(random = pa_random, buffer = pa_buffer_out,
                       uniform = pa_uniform, mcmc = pa_mcmc),
    n_realizations = N_REAL, max_pres = MAXP,
    bw = bw, pa_env_rast = envData, sdm_hook = hook,
    dimensions = pc_cols, chain.length = CHAIN, burnIn = BURN,
    environmental.cutof.percentile = 0.001, species.cutoff.threshold = CUTOFF,
    grid.res = GRID_RES, thres = 0.75,                 # -> pa_uniform (USE)
    precomputed.env = env_bundle, buffer_km = BUFFER_KM,
    verbose = TRUE, parallel = PAR, n_workers = n_workers)
  if (!is.null(res$hsm)) { res$hsm$species_label <- sp$label; species_hsm[[nm]] <- res$hsm }
  species_metrics[[nm]] <- collect_experiment_metrics(setNames(list(res), nm),
                                                      species_labels = setNames(sp$label, nm))
  prevalence[nm] <- mean(res$background$pa)
  rm(res); gc(verbose = FALSE)
}

# --- 5. Persist tidy stats ---------------------------------------------------
hsm <- do.call(rbind, species_hsm); rownames(hsm) <- NULL
write.csv(hsm, file.path(out_dir, paste0("hsm_metrics_2d24_", mode, ".csv")), row.names = FALSE)
metrics <- do.call(rbind, species_metrics); rownames(metrics) <- NULL
write.csv(metrics, file.path(out_dir, paste0("metrics_2d24_", mode, ".csv")), row.names = FALSE)
cat(sprintf("\nHSM rows: %d  (ok=%d, skip_few_pres=%d, fit_failed=%d)\n",
            nrow(hsm), sum(hsm$status == "ok"), sum(hsm$status == "skip_few_pres"),
            sum(hsm$status == "fit_failed")))

species_labels <- vapply(SPECIES, `[[`, character(1), "label")

# --- 6. Dunn test + figures (both halves) ------------------------------------
dunn <- tryCatch(hsm_dunn_test(hsm), error = function(e) { cat("Dunn skipped:", conditionMessage(e), "\n"); NULL })
if (!is.null(dunn)) write.csv(dunn, file.path(out_dir, paste0("hsm_dunn_2d24_", mode, ".csv")), row.names = FALSE)

# §2.3 sampler-comparison boxplots (all 4 samplers incl USE)
np <- experiment_panel_count(metrics)
p_box <- plot_experiment_boxplots(metrics, jitter = (N_REAL <= 30),
  title = "Sampler comparison — PC1×PC2 of the 24-layer env (2-D)")
save_experiment_figure(p_box, file.path(out_dir, paste0("boxplots_metrics_2d24_", mode, ".pdf")),
                       width = 7.0, panel_count = np)

# HSM violins per predictor set + Dunn summary
for (ps in names(predictor_sets)) {
  comb <- tryCatch(plot_hsm_violin(hsm, predictor_set = ps, species_labels = species_labels,
            title = sprintf("Downstream HSM accuracy (2-D, %s predictors)", ps)),
            error = function(e) { cat("violin", ps, "skipped:", conditionMessage(e), "\n"); NULL })
  if (!is.null(comb))
    save_hsm_figure(comb, file.path(out_dir, paste0("hsm_violin_", ps, "_", mode, ".pdf")),
                    width = min(20, 3 + 2.4 * length(unique(hsm$species))), height = 8)
}
if (!is.null(dunn)) {
  ds <- tryCatch(plot_hsm_dunn_summary(dunn), error = function(e) NULL)
  if (!is.null(ds)) save_hsm_figure(ds, file.path(out_dir, paste0("hsm_dunn_summary_2d24_", mode, ".pdf")),
                                    width = 11, height = 3.5)
}

# --- 7. Lightweight summary + console ----------------------------------------
saveRDS(list(hsm = hsm, metrics = metrics, dunn = dunn,
             medians = tryCatch(summarize_hsm(hsm), error = function(e) NULL),
             species_labels = species_labels, prevalence = prevalence,
             predictor_sets = predictor_sets, pc_cols = pc_cols, bw = bw,
             mode = mode, n_realizations = N_REAL, max_pres = MAXP,
             species_cutoff = CUTOFF, grid_res = GRID_RES),
        file.path(out_dir, paste0("summary_2d24_", mode, ".rds")))

cat("\n================ SUMMARY ================\n")
for (nm in names(prevalence)) cat(sprintf("  %s  prevalence=%.4f\n", nm, prevalence[nm]))
cat("\nSaved to ", out_dir, "\n", sep = "")
