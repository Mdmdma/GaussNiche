# run_2d_experiment.R  —  the 4 virtual species in 2-D (paper §2.3)
# =============================================================================
# Reproduces the "Experimental setting" of the markov-chain-sampler paper: the
# 2x2 factorial (breadth x position) of four virtual species in the leading two
# PCs of the WorldClim bioclimatic stack (USE.MCMC::Worldclim_tmp), benchmarking
# three pseudo-absence samplers — RND (pa_random), USE (pa_uniform), and uniform+
# (pa_mcmc) — over R independent Bernoulli realisations.
#
# Niche design (Table 1): rho = 0, prevalence 1:1, all four share the breadth x
# position grid:
#   sp1 generalist common  mu=(0,0)  sigma=(3,2)
#   sp2 specialist common  mu=(0,0)  sigma=(1,0.5)
#   sp3 generalist rare    mu=(4,1)  sigma=(3,2)
#   sp4 specialist rare    mu=(4,1)  sigma=(1,0.5)
#
# Usage (positional args), run on a compute node via sbatch/submit_2d_experiment.sh:
#   Rscript run_2d_experiment.R [mode] [n_workers]
#     mode = "smoke" -> 1 species, 2 realisations, tiny chain (validation)
#     mode = "full"  -> all 4 species, R from N_REALIZATIONS env (default 50)
#   env overrides: N_REALIZATIONS, MAX_PRES, SPECIES_CUTOFF
#
# Inputs : USE.MCMC::Worldclim_tmp (bundled; no scratch input needed)
# Outputs: <scratch>/GaussNiche/results2d/
#            <sp>_<mode>.rds                full per-species result objects
#            metrics_2d_<mode>.csv          tidy combined metrics (one row/realis.)
#            summary_2d_<mode>.rds           lightweight bundle for the vignette
#            boxplots_metrics_2d_<mode>.pdf/.png   the §2.3 figure (all metrics)
#            box_{overlap,coverage,trueabs}_2d_<mode>.pdf   per-metric figures
# =============================================================================

suppressPackageStartupMessages({
  .libPaths(c("~/R/rocker-rstudio/4.5", .libPaths()))
  library(terra); library(USE.MCMC); library(ggplot2)
  library(sf); library(hypervolume); library(MASS)
  library(future); library(furrr)
})
source("virtualSpecies_fn.R")
source("experiment_plots.R")

args      <- commandArgs(trailingOnly = TRUE)
mode      <- if (length(args) >= 1) args[1] else "full"
n_workers <- if (length(args) >= 2) as.integer(args[2]) else
  as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4"))

if (mode == "smoke") {
  N_REAL <- 2L; MAXP <- 80L; CHAIN <- 2000L; BURN <- 500L
  SPECIES_SET <- "sp1"; PAR <- FALSE; CALIBRATE_GRID <- FALSE
} else {
  N_REAL <- as.integer(Sys.getenv("N_REALIZATIONS", "50"))
  MAXP   <- as.integer(Sys.getenv("MAX_PRES", "1000"))
  CHAIN  <- 20000L; BURN <- 1000L
  SPECIES_SET <- "all"; PAR <- TRUE; CALIBRATE_GRID <- TRUE
}
CUTOFF <- as.numeric(Sys.getenv("SPECIES_CUTOFF", "0.1"))
cat(sprintf("mode=%s  n_real=%d  max_pres=%d  chain=%d  cutoff=%.2f  parallel=%s  workers=%d\n",
            mode, N_REAL, MAXP, CHAIN, CUTOFF, PAR, n_workers))

user    <- Sys.getenv("USER"); scratch <- file.path("/cluster/scratch", user)
out_dir <- file.path(scratch, "GaussNiche", "results2d")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# --- 1. Environmental data and PCA -------------------------------------------
envData <- terra::rast(USE.MCMC::Worldclim_tmp, type = "xyz")
rpc     <- rastPCA(envData, stand = TRUE)
dt      <- na.omit(as.data.frame(rpc$PCs[[c("PC1", "PC2")]], xy = TRUE))
var_exp <- rpc$pca$sdev^2 / sum(rpc$pca$sdev^2)
pc1_lab <- paste0("PC1 (", round(var_exp[1] * 100, 1), "% var)")
pc2_lab <- paste0("PC2 (", round(var_exp[2] * 100, 1), "% var)")
cat(sprintf("background cells: %d   PC1 var=%.1f%%   PC2 var=%.1f%%\n",
            nrow(dt), var_exp[1] * 100, var_exp[2] * 100))

# Background KDE (once; reused by every species for the niche plot).
kde        <- MASS::kde2d(dt$PC1, dt$PC2, n = 200)
kde_df     <- expand.grid(PC1 = kde$x, PC2 = kde$y)
kde_df$density <- as.vector(kde$z); kde_df$density <- kde_df$density / max(kde_df$density)

# --- 2. One-time computations (fixed across species for comparability) -------
bw_background <- compute_bandwidth(dt)
cat("background bandwidth (PC1, PC2):", paste(round(bw_background, 5), collapse = ", "), "\n")

# USE grid resolution: calibrate once with optimRes (cached to scratch), else the
# pa_uniform default. optimRes is the slow one-off step (minutes) — see CLAUDE.md.
grid_res_opt <- 5
if (CALIBRATE_GRID) {
  cache <- file.path(out_dir, "optimres_2d.rds")
  if (file.exists(cache)) {
    grid_res_opt <- readRDS(cache); cat("optimRes (cached) grid.res:", grid_res_opt, "\n")
  } else {
    grid_res_opt <- tryCatch({
      dt_sf  <- sf::st_as_sf(dt, coords = c("PC1", "PC2"))
      g <- as.numeric(USE.MCMC::optimRes(sdf = dt_sf, grid.res = 1:10,
                                         perc.thr = 20, showOpt = FALSE, cr = 5)$Opt_res)
      saveRDS(g, cache); g
    }, error = function(e) { cat("optimRes failed (", conditionMessage(e),
                                 "); using grid.res = 5\n"); 5 })
    cat("optimRes grid.res:", grid_res_opt, "\n")
  }
}

# --- 3. Species catalogue (Table 1) ------------------------------------------
species_catalogue <- list(
  sp1_generalist_common = list(label = "Generalist · common",
                               mu = c(0, 0), sigma_PC1 = 3, sigma_PC2 = 2),
  sp2_specialist_common = list(label = "Specialist · common",
                               mu = c(0, 0), sigma_PC1 = 1, sigma_PC2 = 0.5),
  sp3_generalist_rare   = list(label = "Generalist · rare",
                               mu = c(4, 1), sigma_PC1 = 3, sigma_PC2 = 2),
  sp4_specialist_rare   = list(label = "Specialist · rare",
                               mu = c(4, 1), sigma_PC1 = 1, sigma_PC2 = 0.5))
if (SPECIES_SET == "sp1") species_catalogue <- species_catalogue["sp1_generalist_common"]

# --- 4. Shared virtualSpecies() arguments ------------------------------------
shared_args <- list(
  dt = dt, envData = envData, rho = 0,
  max_pres = MAXP, bgk_prev = 1,
  pa_samplers = list(random = pa_random, uniform = pa_uniform, mcmc = pa_mcmc),
  n_realizations = N_REAL, seed_base = 42, seed_pseudo_base = 123,
  pc1_lab = pc1_lab, pc2_lab = pc2_lab, kde_df = kde_df,
  bw = bw_background, pa_env_rast = envData,
  verbose = TRUE, parallel = PAR,
  n_workers = if (is.na(n_workers)) NULL else n_workers,
  grid.res = grid_res_opt, thres = 0.75,
  chain.length = CHAIN, burnIn = BURN,
  species.cutoff.threshold = CUTOFF)

# --- 5. Run all species ------------------------------------------------------
species_results <- list()
for (nm in names(species_catalogue)) {
  cat("\n##########  ", nm, "  ##########\n", sep = "")
  spec <- species_catalogue[[nm]]
  res  <- do.call(virtualSpecies, c(shared_args, list(
    mu = spec$mu, sigma_PC1 = spec$sigma_PC1, sigma_PC2 = spec$sigma_PC2)))
  saveRDS(res, file.path(out_dir, paste0(nm, "_", mode, ".rds")))
  species_results[[nm]] <- res
}

# --- 6. Collect metrics, write CSV + lightweight summary ----------------------
species_labels <- vapply(species_catalogue, `[[`, character(1), "label")
metrics <- collect_experiment_metrics(species_results, species_labels = species_labels)
write.csv(metrics, file.path(out_dir, paste0("metrics_2d_", mode, ".csv")), row.names = FALSE)

summary_bundle <- list(
  metrics = metrics,
  catalogue = species_catalogue, species_labels = species_labels,
  var_exp = var_exp, grid_res_opt = grid_res_opt, bw = bw_background,
  prevalence = vapply(species_results, function(r) mean(r$background$pa), numeric(1)),
  mode = mode, n_realizations = N_REAL, max_pres = MAXP)
saveRDS(summary_bundle, file.path(out_dir, paste0("summary_2d_", mode, ".rds")))

# --- 7. The demanded figures (boxplots; vector PDF + PNG) --------------------
np <- experiment_panel_count(metrics)
p_box <- plot_experiment_boxplots(
  metrics, jitter = (N_REAL <= 30),
  title = "Sampler comparison across virtual species")
save_experiment_figure(p_box, file.path(out_dir, paste0("boxplots_metrics_2d_", mode, ".pdf")),
                       width = 7.0, panel_count = np)
save_experiment_figure(p_box, file.path(out_dir, paste0("boxplots_metrics_2d_", mode, ".png")),
                       width = 7.0, panel_count = np)
for (mt in c("overlap", "coverage", "trueabs")) {
  metric_arg <- if (mt == "trueabs") "prop_true_abs" else mt
  ph <- if (mt == "coverage") 1.7 * length(grep("^rel_cov_", names(metrics))) else 3.2
  pm <- plot_experiment_metric(metrics, metric = metric_arg, jitter = (N_REAL <= 30))
  save_experiment_figure(pm, file.path(out_dir, paste0("box_", mt, "_2d_", mode, ".pdf")),
                         width = 7.0, height = ph)
}

# --- 8. Console summary -------------------------------------------------------
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
