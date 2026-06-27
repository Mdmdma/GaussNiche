# tune_5d_hsm.R --- uniform+ settings sweep for downstream-HSM performance
# =============================================================================
# Sweeps a grid of the two uniform+ (pa_mcmc) tunables that control how far the
# pseudo-absences are pushed from the presences / the realised environment, to
# find a setting that closes uniform+'s truth-recovery gap (cor/RMSE-to-truth,
# CBI) without giving back its discrimination edge (AUC/TSS) vs the fixed
# RND / buffer-out baselines (results5d_hsm/hsm_metrics_5d_full.csv).
#
#   environmental.cutof.percentile  in ENV_CUTOFFS  (default 0.001,0.005,0.01,0.05)
#   species.cutoff.threshold        in SPECIES_CUTOFFS (default 0.9,0.75,0.6,0.5)
#
# Speedups (measured by bench_hsm.R) baked in so the 16-cell grid runs in minutes:
#   A1  uniform+ only (RND/buffer are knob-invariant -> reuse the baseline)
#   A2  compute_hypervolume = FALSE   (the ~70%-of-task 5-D hypervolume cost)
#   A3  pc5 predictors only           (raw12 told the same story, milder)
#   A4  25 realisations               (ranks knob cells; re-confirm winner at 50)
#
# Usage:  Rscript tune_5d_hsm.R [mode] [n_workers]
#   mode = smoke -> 1 species, 2x2 grid, 2 realisations (validation)
#   mode = full  -> 4 species, full grid, N_REALIZATIONS (default 25)
#   env overrides: ENV_CUTOFFS, SPECIES_CUTOFFS, N_REALIZATIONS, MAX_PRES, MIN_PRES
# Inputs : env5d/{final_stack_lean_natural.tif, background_5d.rds}
#          results5d_hsm/hsm_metrics_5d_full.csv (baseline RND/buffer reference)
# Outputs: <scratch>/GaussNiche/results5d_tune/
# Run via sbatch/submit_tune_5d_hsm.sh.
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
HSM_EVAL_PATH <- normalizePath("hsm_eval.R")

args      <- commandArgs(trailingOnly = TRUE)
mode      <- if (length(args) >= 1) args[1] else "full"
n_workers <- if (length(args) >= 2) as.integer(args[2]) else
  as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4"))

ENV_CUTOFFS <- as.numeric(strsplit(Sys.getenv("ENV_CUTOFFS", "0.001,0.005,0.01,0.05"), ",")[[1]])
SP_CUTOFFS  <- as.numeric(strsplit(Sys.getenv("SPECIES_CUTOFFS", "0.9,0.75,0.6,0.5"), ",")[[1]])
if (mode == "smoke") {
  N_REAL <- 2L; MAXP <- 120L; CHAIN <- 2000L; BURN <- 500L
  SPECIES_SET <- "sp1"; ENV_CUTOFFS <- c(0.001, 0.05); SP_CUTOFFS <- c(0.9, 0.5)
  if (is.na(n_workers) || n_workers < 2L) n_workers <- 2L
} else {
  N_REAL <- as.integer(Sys.getenv("N_REALIZATIONS", "25"))
  MAXP   <- as.integer(Sys.getenv("MAX_PRES", "300"))
  CHAIN  <- 20000L; BURN <- 1000L
  SPECIES_SET <- "all"
}
MIN_PRES <- as.integer(Sys.getenv("MIN_PRES", "12"))
ALGOS    <- HSM_ALGORITHMS
cat(sprintf("mode=%s n_real=%d max_pres=%d workers=%d\n  env_cutoffs=%s\n  species_cutoffs=%s\n",
            mode, N_REAL, MAXP, n_workers, paste(ENV_CUTOFFS, collapse=","), paste(SP_CUTOFFS, collapse=",")))

user    <- Sys.getenv("USER"); scratch <- file.path("/cluster/scratch", user)
env_dir <- file.path(scratch, "GaussNiche", "env5d")
out_dir <- file.path(scratch, "GaussNiche", "results5d_tune")
base_csv <- file.path(scratch, "GaussNiche", "results5d_hsm", "hsm_metrics_5d_full.csv")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# --- 1. Environment (pc5 only -> no raw-layer augmentation) -------------------
envData <- terra::rast(file.path(env_dir, "final_stack_lean_natural.tif"))
dt      <- readRDS(file.path(env_dir, "background_5d.rds"))
pc_cols <- paste0("PC", 1:5)
stopifnot(all(pc_cols %in% names(dt)))
predictor_sets <- list(pc5 = pc_cols)
cat(sprintf("background cells: %d\n", nrow(dt)))

pc_sd    <- vapply(dt[pc_cols], stats::sd, numeric(1))
sig_gen  <- 1.0 * pc_sd; sig_spec <- 0.5 * pc_sd
bg5 <- as.matrix(dt[, pc_cols])
mode_estimate <- function(X) {
  if (requireNamespace("ks", quietly = TRUE)) {
    est <- tryCatch(ks::kde(x = X, H = ks::Hns(X), eval.points = X)$estimate, error = function(e) NULL)
    if (!is.null(est)) return(as.numeric(X[which.max(est), ]))
  }
  k <- min(200L, nrow(X) - 1L); as.numeric(X[which.min(FNN::knn.dist(X, k = k)[, k]), ])
}
mu_common <- mode_estimate(bg5); names(mu_common) <- pc_cols
mu_rare   <- c(PC1 = 1.55, PC2 = 1.13, PC3 = 1.4, PC4 = 0.8, PC5 = 1.2)[pc_cols]
SPECIES <- list(
  sp1_generalist_common = list(label = "Generalist · common", mu = mu_common, sigma = sig_gen),
  sp2_specialist_common = list(label = "Specialist · common", mu = mu_common, sigma = sig_spec),
  sp3_generalist_rare   = list(label = "Generalist · rare",   mu = mu_rare,   sigma = sig_gen),
  sp4_specialist_rare   = list(label = "Specialist · rare",   mu = mu_rare,   sigma = sig_spec))
if (SPECIES_SET == "sp1") SPECIES <- SPECIES["sp1_generalist_common"]
species_labels <- vapply(SPECIES, `[[`, character(1), "label")

# --- 2. One-time computations ------------------------------------------------
bw <- compute_bandwidth_nd(dt, pc_cols)
cat("\nPrecomputing MCMC environment once...\n")
env_bundle <- USE.MCMC::precomputeMcmcEnvironment(
  env.rast = envData, dimensions = pc_cols, seed.number = 123, verbose = FALSE)
env_bundle["rng_state"] <- list(NULL)
# Save the (multi-MB) bundle to disk and pass its PATH to the sampler: furrr
# ships the tiny string, each worker loads it once (cached). See .gn_load_env_bundle.
bundle_path <- file.path(out_dir, "env_bundle.rds")
saveRDS(env_bundle, bundle_path)
LIBS <- .libPaths()

# One hook per species (knob-independent: it just scores whatever PA it is given;
# reused across grid cells so the train/test split is identical for a fair compare).
hooks <- list(); sp_i <- 0L
for (nm in names(SPECIES)) {
  sp_i <- sp_i + 1L
  hooks[[nm]] <- make_sdm_hook(species = nm, predictor_sets = predictor_sets,
                               algorithms = ALGOS, test_frac = 0.3, min_pres = MIN_PRES,
                               seed = 1000L + sp_i * 100L, libpaths = LIBS,
                               source_file = HSM_EVAL_PATH)
}

# --- 3. Sweep the grid -------------------------------------------------------
# WHY this shape (read before changing it):
#   * Each virtualSpecies_nd() call evaluates ONE (cell, species) over N_REAL
#     realisations; those calls are fully independent.
#   * So we run the WHOLE set of (cell x species) jobs in ONE future pool on ONE
#     node, each job's realisations serial INSIDE its worker (parallel = FALSE).
#     The default 4x4 grid is 16 cells x 4 species = 64 jobs -> a clean fit for a
#     64-core node: one wave, every core busy, ~the wall time of a single job.
#   * It stays a SINGLE node / SINGLE SLURM job (no job arrays) -- the easiest
#     thing for someone else to reproduce: `sbatch sbatch/submit_tune_5d_hsm.sh`.
#   * Bit-identical to the old per-cell loop: every result depends only on the
#     layered seeds (seed_base + r / seed_pseudo_base + r) set inside
#     eval_realization, never on which worker runs a job.
grid <- expand.grid(env_c = ENV_CUTOFFS, sp_c = SP_CUTOFFS, KEEP.OUT.ATTRS = FALSE)
jobs <- expand.grid(cell = seq_len(nrow(grid)), species = names(SPECIES),
                    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
cat(sprintf("\nGrid: %d cells x %d species = %d jobs x %d realisations = %d uniform+ tasks  (workers=%d)\n",
            nrow(grid), length(SPECIES), nrow(jobs), N_REAL, nrow(jobs) * N_REAL, n_workers))

# A SpatRaster cannot cross a future worker boundary (its C++ pointer is invalid
# in the worker), so wrap it ONCE here and unwrap it inside each job.
envData_packed <- terra::wrap(envData)

# One job = one (cell, species). Returns its $hsm rows tagged with the two swept
# knobs (or NULL if the cell is degenerate). The uniform+-only / no-hypervolume /
# no-plots / no-reference flags strip the call to just the downstream-HSM numbers.
run_job <- function(cell, species) {
  ed <- terra::unwrap(envData_packed)                   # local, valid SpatRaster
  env_c <- grid$env_c[cell]; sp_c <- grid$sp_c[cell]; sp <- SPECIES[[species]]
  res <- virtualSpecies_nd(
    dt = dt, envData = ed, mu = sp$mu, sigma = sp$sigma, rho = 0,
    pc_cols = pc_cols, bgk_prev = 1,
    pa_samplers = list(mcmc = pa_mcmc),                 # A1: uniform+ only
    n_realizations = N_REAL, max_pres = MAXP,
    bw = bw, pa_env_rast = ed, sdm_hook = hooks[[species]],
    compute_hypervolume = FALSE,                        # A2: skip the 5-D hypervolume
    make_plots = FALSE, compute_reference = FALSE,      # no plots, no reference draw
    dimensions = pc_cols, chain.length = CHAIN, burnIn = BURN,
    environmental.cutof.percentile = env_c,             # swept knob
    species.cutoff.threshold = sp_c,                    # swept knob
    precomputed.env = bundle_path,
    verbose = FALSE, parallel = FALSE)                  # realisations serial INSIDE the job
  h <- res$hsm
  if (!is.null(h)) { h$env_cutoff <- env_c; h$species_cutoff <- sp_c; h$species_label <- sp$label }
  h
}

future::plan(future::multisession, workers = n_workers)   # one pool for the whole grid
rows <- furrr::future_pmap(jobs, run_job,
                           .options = furrr::furrr_options(seed = TRUE, globals = TRUE))
future::plan(future::sequential)                          # tear the pool down

rows <- rows[!vapply(rows, is.null, logical(1))]
tune <- do.call(rbind, rows); rownames(tune) <- NULL
write.csv(tune, file.path(out_dir, paste0("tune_hsm_metrics_5d_", mode, ".csv")), row.names = FALSE)
cat(sprintf("\nTune rows: %d  (ok=%d, skip_few_pres=%d, fit_failed=%d)\n",
            nrow(tune), sum(tune$status == "ok"), sum(tune$status == "skip_few_pres"),
            sum(tune$status == "fit_failed")))

# --- 4. Baseline RND / buffer-out reference (pc5, knob-invariant) -------------
ref <- list()
if (file.exists(base_csv)) {
  b <- read.csv(base_csv, stringsAsFactors = FALSE)
  b <- b[b$status == "ok" & b$predictor_set == "pc5" & b$sampler %in% c("random", "buffer"), ]
  for (m in c("auc","tss","boyce","cor_truth","rmse_truth"))
    ref[[m]] <- c(random = stats::median(b[b$sampler=="random", m], na.rm=TRUE),
                  buffer = stats::median(b[b$sampler=="buffer", m], na.rm=TRUE))
  cat("\nBaseline pc5 reference (median over species+algos+reals):\n")
  print(round(do.call(rbind, ref), 3))
} else cat("\n(baseline", base_csv, "not found -> heatmaps without reference)\n")

# --- 5. Per-cell aggregate table + heatmaps ----------------------------------
metr <- c("cor_truth","rmse_truth","auc","tss","boyce")
to   <- tune[tune$status == "ok", ]
cell <- aggregate(to[metr], by = list(env_cutoff = to$env_cutoff, species_cutoff = to$species_cutoff),
                  FUN = function(v) round(stats::median(v, na.rm = TRUE), 4))
cell <- cell[order(cell$env_cutoff, -cell$species_cutoff), ]
write.csv(cell, file.path(out_dir, paste0("tune_cell_medians_", mode, ".csv")), row.names = FALSE)
cat("\n================= PER-CELL MEDIANS (pooled species+algos+reals) =================\n")
print(cell, row.names = FALSE)

for (m in metr) {
  hb <- (m != "rmse_truth")
  rr <- if (length(ref)) unname(ref[[m]]["random"]) else NULL
  for (kind in c("value", "delta")) {
    if (kind == "delta" && is.null(rr)) next
    p <- tryCatch(plot_tune_heatmap(tune, m,
           ref = if (kind == "delta") rr else NULL,
           baseline = if (length(ref)) ref[[m]] else NULL, higher_better = hb,
           title = sprintf("uniform+ %s%s", m, if (kind == "delta") " (Δ vs RND)" else "")),
           error = function(e) { cat("heatmap", m, kind, "skipped:", conditionMessage(e), "\n"); NULL })
    if (!is.null(p))
      save_hsm_figure(p, file.path(out_dir, sprintf("tune_heat_%s_%s_%s.pdf", m, kind, mode)),
                      width = 6.2, height = 4.6)
  }
}

saveRDS(list(tune = tune, cell = cell, ref = ref, env_cutoffs = ENV_CUTOFFS,
             species_cutoffs = SP_CUTOFFS, n_realizations = N_REAL, mode = mode),
        file.path(out_dir, paste0("summary_tune_", mode, ".rds")))
cat("\nSaved to ", out_dir, "\n", sep = "")
