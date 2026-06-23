# build_final_stack.R  —  lock the lean_natural 5-D environment (5d-niche branch)
# =============================================================================
# Subsets the locked 12-layer 'lean_natural' set from the cached candidate stack,
# runs the standardised PCA, keeps PC1-5 as the analysis E-space, and reports the
# per-PC loadings (axis interpretation) + per-PC background distribution (needed
# to place the 4 virtual species' mu/sigma). No downloads — reuses the cache.
#
# Outputs (to <scratch>/GaussNiche/env5d/):
#   final_stack_lean_natural.tif  the 12 chosen env layers
#   final_pcs.tif                 PC-score rasters (all PCs)
#   final_rastPCA.rds             PCA model (princomp-style)
#   background_5d.rds             na.omit data.frame (x, y, PC1..PC5)  <- niche input
#   final_loadings.csv            variable loadings
#   final_pc_stats.csv            per-PC mean/sd/quantiles
# =============================================================================

suppressPackageStartupMessages({
  .libPaths(c("~/R/rocker-rstudio/4.5", .libPaths()))
  library(terra)
  library(USE.MCMC)
})

user     <- Sys.getenv("USER")
scratch  <- file.path("/cluster/scratch", user)
out_dir  <- file.path(scratch, "GaussNiche", "env5d")
cand_tif <- file.path(out_dir, "env_stack_candidate.tif")
stopifnot("candidate stack not found; run build_env_stack.R first" = file.exists(cand_tif))
stk <- terra::rast(cand_tif)

# locked 'lean_natural' set: 5 climate + 4 soil + 2 terrain + 1 vegetation
LEAN <- c("clim_wc2.1_10m_bio_3", "clim_wc2.1_10m_bio_4", "clim_wc2.1_10m_bio_9",
          "clim_wc2.1_10m_bio_14", "clim_wc2.1_10m_bio_15",
          "soil_phh2o", "soil_clay", "soil_bdod", "soil_soc",
          "ter_TRI", "ter_TPI", "veg_trees")
miss <- setdiff(LEAN, names(stk))
if (length(miss)) stop("missing layers in candidate stack: ", paste(miss, collapse = ", "))
env <- stk[[LEAN]]
cat("Final env stack: ", terra::nlyr(env), " layers\n  ",
    paste(names(env), collapse = ", "), "\n", sep = "")

# --- PCA ----------------------------------------------------------------------
rpc <- USE.MCMC::rastPCA(env, stand = TRUE)
ve  <- rpc$pca$sdev^2 / sum(rpc$pca$sdev^2)
cum <- cumsum(ve)
n80 <- which(cum >= 0.80)[1]
cat(sprintf("\nVariance: first 5 PCs = %.1f%%   (80%% reached at %d PCs)\n",
            100 * cum[5], n80))
print(data.frame(PC = paste0("PC", seq_along(ve)),
                 var_pct = round(100 * ve, 2), cum_pct = round(100 * cum, 2)),
      row.names = FALSE)

# --- loadings: how the 5 axes map to datasets ---------------------------------
L <- tryCatch(unclass(rpc$pca$loadings), error = function(e) rpc$pca$rotation)
cat("\nLoadings (PC1-PC5):\n"); print(round(L[, 1:5], 3))
cat("\nDominant variable per PC:\n")
for (k in 1:5) {
  lk <- abs(L[, k])
  cat(sprintf("  PC%d  %-22s |load|=%.2f  var=%.1f%%\n",
              k, rownames(L)[which.max(lk)], max(lk), 100 * ve[k]))
}

# --- 5-D background + per-PC distribution (for species mu/sigma) ---------------
K   <- 5L
pcs <- rpc$PCs[[paste0("PC", 1:K)]]
dt  <- na.omit(as.data.frame(pcs, xy = TRUE))
cat(sprintf("\nBackground: %d cells in 5-D PC space\n", nrow(dt)))
pcc <- paste0("PC", 1:K)
pc_stats <- data.frame(
  PC   = pcc,
  mean = sapply(dt[pcc], mean),
  sd   = sapply(dt[pcc], stats::sd),
  min  = sapply(dt[pcc], min),
  q05  = sapply(dt[pcc], stats::quantile, 0.05),
  q50  = sapply(dt[pcc], stats::median),
  q95  = sapply(dt[pcc], stats::quantile, 0.95),
  max  = sapply(dt[pcc], max),
  row.names = NULL)
cat("\nPer-PC background distribution (basis for common/rare mu and generalist/specialist sigma):\n")
ps <- pc_stats; isnum <- sapply(ps, is.numeric); ps[isnum] <- lapply(ps[isnum], round, 3)
print(ps, row.names = FALSE)

# --- save ---------------------------------------------------------------------
terra::writeRaster(env,     file.path(out_dir, "final_stack_lean_natural.tif"), overwrite = TRUE)
terra::writeRaster(rpc$PCs, file.path(out_dir, "final_pcs.tif"),                overwrite = TRUE)
saveRDS(rpc$pca, file.path(out_dir, "final_rastPCA.rds"))
saveRDS(dt,      file.path(out_dir, "background_5d.rds"))
write.csv(round(L, 4),        file.path(out_dir, "final_loadings.csv"))
write.csv(pc_stats, file.path(out_dir, "final_pc_stats.csv"), row.names = FALSE)
cat("\nSaved final 5-D environment to ", out_dir, "\n", sep = "")
