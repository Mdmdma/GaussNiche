# analyze_env_stack.R  —  Phase 1 (diagnostics) of the 5d-niche branch
# =============================================================================
# Reads the cached candidate stack (build_env_stack.R output) and computes
# standard data-science diagnostics to decide WHICH layers/blocks are needed to
# reach >=5 useful PCA dimensions:
#
#   1. Per-variable descriptive statistics (per dataset/block).
#   2. Full Pearson + Spearman correlation matrices; VIF (multicollinearity);
#      highly-correlated pairs.
#   3. Per-block diagnostics: intra-block effective dimensionality, mean intra
#      correlation, and mean |r| with the climate block (orthogonality = how much
#      genuinely NEW signal each dataset adds).
#   4. Block x block mean |r| matrix.
#   5. Full standardised PCA: eigenvalues, variance explained, and THREE
#      "how many useful PCs" criteria — cumulative 80%, Kaiser (eig>1), and the
#      broken-stick model. Loadings -> which block leads each PC.
#   6. Forward block inclusion: climate -> +soil -> +terrain -> +veg -> +anth,
#      reporting the marginal gain in useful PCs at each step.
#   7. A parsimonious recommendation (drop within-block redundancy) + its PCA.
#
# Outputs (to <scratch>/GaussNiche/env5d/diagnostics/): CSV tables, PNG plots,
# a metrics .rds, and a human-readable report.txt. Base R only (no new installs).
# =============================================================================

suppressPackageStartupMessages({
  .libPaths(c("~/R/rocker-rstudio/4.5", .libPaths()))
  library(terra)
})

user    <- Sys.getenv("USER")
scratch <- file.path("/cluster/scratch", user)
out_dir <- file.path(scratch, "GaussNiche", "env5d")
diag_dir <- file.path(out_dir, "diagnostics")
dir.create(diag_dir, recursive = TRUE, showWarnings = FALSE)
cand_tif <- file.path(out_dir, "env_stack_candidate.tif")
stopifnot("cached stack not found; run build_env_stack.R first" = file.exists(cand_tif))

# Tee all console output to report.txt as well.
report_con <- file(file.path(diag_dir, "report.txt"), open = "wt")
sink(report_con, split = TRUE)
on.exit({ sink(); close(report_con) }, add = TRUE)

hr <- function(t) cat("\n", strrep("=", 78), "\n", t, "\n", strrep("=", 78), "\n", sep = "")
BLOCK_ORDER <- c("clim", "soil", "ter", "veg", "anth")
BLOCK_NAME  <- c(clim = "Climate (WorldClim bioclim)", soil = "Soil (SoilGrids)",
                 ter = "Terrain (SRTM)", veg = "Vegetation (WorldCover)",
                 anth = "Anthropogenic (Footprint)")

# --- load ---------------------------------------------------------------------
stk   <- terra::rast(cand_tif)
vars  <- names(stk)
block <- sub("_.*", "", vars)
block <- factor(block, levels = intersect(BLOCK_ORDER, unique(block)))
p_all <- length(vars)

# per-layer % complete (fraction of non-NA cells)
pct_all <- 100 * (1 - as.numeric(terra::global(is.na(stk), "mean")[, 1]))
names(pct_all) <- vars

# value matrix on complete cases
df0  <- as.data.frame(stk, na.rm = FALSE)
keep <- stats::complete.cases(df0)
df   <- df0[keep, , drop = FALSE]
n_cc <- nrow(df)

# drop zero-variance variables (scale()/prcomp would divide by 0)
sdv <- sapply(df, stats::sd)
zv  <- names(df)[!is.finite(sdv) | sdv == 0]
if (length(zv)) {
  message("Dropping zero-variance variables: ", paste(zv, collapse = ", "))
  df <- df[, !(names(df) %in% zv), drop = FALSE]
}
# re-align variable/block bookkeeping to the surviving columns
vars  <- names(df)
bsrc  <- sub("_.*", "", vars)
block <- factor(bsrc, levels = intersect(BLOCK_ORDER, unique(bsrc)))
p_all <- length(vars)
pct_complete <- pct_all[vars]
hr("INPUT")
cat(sprintf("Stack: %d variables x %d cells (%d complete cases, %.1f%%)\n",
            p_all, nrow(df0), n_cc, 100 * n_cc / nrow(df0)))
cat("Blocks:\n"); print(table(block))

# helper metrics ---------------------------------------------------------------
skew <- function(x) { x <- x[is.finite(x)]; m <- mean(x); s <- stats::sd(x)
                      if (s == 0) return(0); mean((x - m)^3) / s^3 }
# number of PCs to reach a cumulative-variance threshold
n_pc_cum <- function(sdev, thr = 0.80) {
  ve <- sdev^2 / sum(sdev^2); which(cumsum(ve) >= thr)[1]
}
# broken-stick: contiguous count of PCs whose variance fraction exceeds the
# null expectation bs_k = (1/p) * sum_{i=k}^{p} 1/i
n_pc_brokenstick <- function(sdev) {
  ve <- sdev^2 / sum(sdev^2); p <- length(ve)
  bs <- sapply(seq_len(p), function(k) sum(1 / (k:p)) / p)
  ok <- ve > bs; w <- which(!ok)[1]; if (is.na(w)) p else w - 1L
}
pca_summary <- function(X) {                       # X: data.frame/matrix of vars
  pr <- stats::prcomp(X, center = TRUE, scale. = TRUE)
  ev <- pr$sdev^2
  list(prcomp = pr, eig = ev, ve = ev / sum(ev), cum = cumsum(ev) / sum(ev),
       n80 = n_pc_cum(pr$sdev, 0.80), n90 = n_pc_cum(pr$sdev, 0.90),
       kaiser = sum(ev > 1), bstick = n_pc_brokenstick(pr$sdev))
}

# =============================================================================
# 1. PER-VARIABLE DESCRIPTIVE STATISTICS
# =============================================================================
hr("1. PER-VARIABLE DESCRIPTIVE STATISTICS")
desc <- data.frame(
  variable = vars, block = as.character(block),
  pct_complete = round(pct_complete, 1),
  mean   = round(sapply(df, mean), 4),
  sd     = round(sapply(df, stats::sd), 4),
  min    = round(sapply(df, min), 4),
  median = round(sapply(df, stats::median), 4),
  max    = round(sapply(df, max), 4),
  skew   = round(sapply(df, skew), 3),
  row.names = NULL
)
print(desc)
write.csv(desc, file.path(diag_dir, "01_variable_stats.csv"), row.names = FALSE)

# =============================================================================
# 2. CORRELATION + MULTICOLLINEARITY
# =============================================================================
hr("2. CORRELATION & MULTICOLLINEARITY")
R_pear <- stats::cor(df, method = "pearson")
R_spear <- stats::cor(df, method = "spearman")
write.csv(round(R_pear, 4),  file.path(diag_dir, "02_corr_pearson.csv"))
write.csv(round(R_spear, 4), file.path(diag_dir, "02_corr_spearman.csv"))

# VIF = diagonal of the inverse correlation matrix (guard against singularity)
vif <- tryCatch({ d <- diag(solve(R_pear)); setNames(round(d, 2), vars) },
                error = function(e) setNames(rep(NA_real_, p_all), vars))
cat("\nVariance Inflation Factor (>5 = notable, >10 = severe collinearity):\n")
print(vif)

# highly correlated pairs
hi <- which(abs(R_pear) > 0.8 & upper.tri(R_pear), arr.ind = TRUE)
if (nrow(hi)) {
  hipairs <- data.frame(v1 = vars[hi[, 1]], v2 = vars[hi[, 2]],
                        r = round(R_pear[hi], 3))
  hipairs <- hipairs[order(-abs(hipairs$r)), ]
  cat("\nHighly correlated pairs (|r| > 0.8):\n"); print(hipairs, row.names = FALSE)
  write.csv(hipairs, file.path(diag_dir, "02_high_corr_pairs.csv"), row.names = FALSE)
} else cat("\nNo variable pairs with |r| > 0.8.\n")

# =============================================================================
# 3. PER-BLOCK DIAGNOSTICS
# =============================================================================
hr("3. PER-BLOCK DIAGNOSTICS (effective dim & orthogonality to climate)")
clim_vars <- vars[block == "clim"]
block_rows <- lapply(levels(block), function(b) {
  bv <- vars[block == b]
  intra <- if (length(bv) >= 2) mean(abs(R_pear[bv, bv][upper.tri(R_pear[bv, bv])])) else NA
  effdim <- if (length(bv) >= 2) pca_summary(df[, bv, drop = FALSE])$n80 else 1L
  # mean |r| of this block's vars vs the climate block (lower = more orthogonal)
  vsclim <- if (b == "clim") NA else
    mean(abs(R_pear[bv, clim_vars, drop = FALSE]))
  data.frame(block = b, label = BLOCK_NAME[[b]], n_vars = length(bv),
             eff_dim_80 = effdim, mean_intra_absr = round(intra, 3),
             mean_absr_vs_climate = round(vsclim, 3))
})
block_tbl <- do.call(rbind, block_rows)
print(block_tbl, row.names = FALSE)
cat("\n(eff_dim_80 = PCs to explain 80% WITHIN the block; lower mean_absr_vs_climate",
    "\n = more orthogonal to climate = adds genuinely new PCA dimensions.)\n")
write.csv(block_tbl, file.path(diag_dir, "03_block_diagnostics.csv"), row.names = FALSE)

# block x block mean |r|
bxb <- outer(levels(block), levels(block), Vectorize(function(a, b) {
  va <- vars[block == a]; vb <- vars[block == b]
  m <- abs(R_pear[va, vb, drop = FALSE]); if (a == b) m <- m[upper.tri(m)]
  if (length(m) == 0) NA else mean(m)
}))
dimnames(bxb) <- list(levels(block), levels(block))
cat("\nBlock x block mean |r|:\n"); print(round(bxb, 3))
write.csv(round(bxb, 3), file.path(diag_dir, "03_block_x_block_absr.csv"))

# =============================================================================
# 4. FULL PCA + "HOW MANY USEFUL PCs"
# =============================================================================
hr("4. FULL PCA (all candidate variables)")
full <- pca_summary(df)
pca_tbl <- data.frame(PC = paste0("PC", seq_along(full$eig)),
                      eigenvalue = round(full$eig, 3),
                      var_pct = round(100 * full$ve, 2),
                      cum_pct = round(100 * full$cum, 2))
print(pca_tbl, row.names = FALSE)
write.csv(pca_tbl, file.path(diag_dir, "04_pca_variance.csv"), row.names = FALSE)
cat(sprintf("\nUseful-PC criteria:\n  cumulative 80%%   : %d PCs\n  cumulative 90%%   : %d PCs\n  Kaiser (eig>1)   : %d PCs\n  broken-stick     : %d PCs\n",
            full$n80, full$n90, full$kaiser, full$bstick))

# loadings: dominant block/variable per PC
load_mat <- full$prcomp$rotation
n_show <- min(ncol(load_mat), max(6L, full$n80))
cat("\nDominant variable per PC (max |loading|):\n")
dom <- sapply(seq_len(n_show), function(k) {
  lk <- abs(load_mat[, k]); i <- which.max(lk)
  sprintf("  PC%-2d  %-16s [%s]  |load|=%.2f  var=%.1f%%",
          k, vars[i], block[i], lk[i], 100 * full$ve[k])
})
cat(paste(dom, collapse = "\n"), "\n")
write.csv(round(load_mat, 4), file.path(diag_dir, "04_pca_loadings.csv"))

# =============================================================================
# 5. FORWARD BLOCK INCLUSION (marginal contribution of each dataset)
# =============================================================================
hr("5. FORWARD BLOCK INCLUSION (which datasets are needed for >=5 useful PCs)")
present <- levels(block)
fwd <- list(); cumv <- character(0)
for (b in present) {
  cumv <- c(cumv, vars[block == b])
  s <- pca_summary(df[, cumv, drop = FALSE])
  fwd[[b]] <- data.frame(added_block = b, n_vars = length(cumv),
                         n80 = s$n80, n90 = s$n90,
                         kaiser = s$kaiser, broken_stick = s$bstick)
}
fwd_tbl <- do.call(rbind, fwd)
print(fwd_tbl, row.names = FALSE)
write.csv(fwd_tbl, file.path(diag_dir, "05_forward_inclusion.csv"), row.names = FALSE)
# first cumulative set reaching >=5 useful (broken-stick) PCs
reach <- fwd_tbl$added_block[which(fwd_tbl$broken_stick >= 5)[1]]
if (!is.na(reach)) {
  cat(sprintf("\n>>> >=5 broken-stick PCs first reached after adding block: '%s'\n", reach))
} else {
  cat("\n>>> Candidate pool does NOT reach 5 broken-stick PCs by that strict criterion.\n")
}

# --- 5b. curated candidate configurations -> PCs for ~80% --------------------
hr("5b. CANDIDATE CONFIGURATIONS (curated subsets vs the '5 PCs ~ 80%' goal)")
tryCatch({
  pick <- function(nm) vars[vars %in% nm]
  cfg <- list(
    all_24           = vars,
    drop_aspect_tpi  = setdiff(vars, c("ter_eastness", "ter_northness", "ter_TPI", "ter_roughness")),
    lean_5axis       = c(vars[block == "clim"],
                         pick(paste0("soil_", c("phh2o", "clay", "bdod", "soc"))),
                         "ter_TRI", "veg_trees"),
    lean_plus_anth   = c(vars[block == "clim"],
                         pick(paste0("soil_", c("phh2o", "clay", "bdod", "soc"))),
                         "ter_TRI", "veg_trees", "veg_grassland", "anth_footprint"),
    # natural-only 5-axis variants (no anthropogenic footprint): use terrain
    # POSITION (TPI) as the 5th orthogonal axis instead of human footprint
    lean_natural     = c(vars[block == "clim"],
                         pick(paste0("soil_", c("phh2o", "clay", "bdod", "soc"))),
                         "ter_TRI", "ter_TPI", "veg_trees"),
    lean_natural_veg = c(vars[block == "clim"],
                         pick(paste0("soil_", c("phh2o", "clay", "bdod", "soc"))),
                         "ter_TRI", "ter_TPI", "veg_trees", "veg_grassland")
  )
  cfg_tbl <- do.call(rbind, lapply(names(cfg), function(nm) {
    v <- intersect(cfg[[nm]], vars)
    s <- pca_summary(df[, v, drop = FALSE])
    data.frame(config = nm, n_vars = length(v),
               cum_5pc_pct = round(100 * s$cum[min(5, length(v))], 1),
               n80 = s$n80, n90 = s$n90, kaiser = s$kaiser, broken_stick = s$bstick)
  }))
  cat("cum_5pc_pct = % variance captured by the FIRST 5 PCs of that subset.\n\n")
  print(cfg_tbl, row.names = FALSE)
  write.csv(cfg_tbl, file.path(diag_dir, "05b_candidate_configs.csv"), row.names = FALSE)
}, error = function(e) message("5b configs failed: ", conditionMessage(e)))

# =============================================================================
# 6. PARSIMONIOUS RECOMMENDATION (drop within-block redundancy)
# =============================================================================
hr("6. PARSIMONIOUS RECOMMENDATION")
# greedily keep variables in block order, dropping any whose |r| with an
# already-kept variable exceeds 0.90 (near-duplicate signal).
kept <- character(0)
for (v in vars[order(block)]) {
  if (length(kept) == 0 || max(abs(R_pear[v, kept])) <= 0.90) kept <- c(kept, v)
}
dropped <- setdiff(vars, kept)
cat("Kept (", length(kept), "):\n", paste(" ", kept, collapse = "\n"), "\n", sep = "")
cat("\nDropped as redundant (|r|>0.90 with a kept var):\n",
    if (length(dropped)) paste(" ", dropped, collapse = "\n") else "  (none)", "\n", sep = "")
red <- pca_summary(df[, kept, drop = FALSE])
cat(sprintf("\nReduced-set PCA: %d vars -> 80%%=%d PCs, Kaiser=%d, broken-stick=%d\n",
            length(kept), red$n80, red$kaiser, red$bstick))
write.csv(data.frame(variable = vars, block = as.character(block),
                     kept = vars %in% kept, vif = vif[vars]),
          file.path(diag_dir, "06_recommended_layers.csv"), row.names = FALSE)

# =============================================================================
# 7. PLOTS (base graphics — no extra packages)
# =============================================================================
# scree + broken-stick + Kaiser
png(file.path(diag_dir, "scree.png"), width = 1100, height = 700, res = 130)
ve <- full$ve; p <- length(ve)
bs <- sapply(seq_len(p), function(k) sum(1 / (k:p)) / p)
plot(seq_len(p), 100 * ve, type = "b", pch = 19, xlab = "Principal component",
     ylab = "Variance explained (%)", main = "Scree: observed vs broken-stick")
lines(seq_len(p), 100 * bs, type = "b", pch = 1, lty = 2, col = "red")
abline(h = 100 / p, col = "grey60", lty = 3)  # Kaiser threshold (eig=1) as % var
legend("topright", c("observed", "broken-stick null", "Kaiser (eig=1)"),
       col = c("black", "red", "grey60"), lty = c(1, 2, 3), pch = c(19, 1, NA), bty = "n")
dev.off()

# cumulative variance
png(file.path(diag_dir, "cum_variance.png"), width = 1100, height = 700, res = 130)
plot(seq_len(p), 100 * full$cum, type = "b", pch = 19, ylim = c(0, 100),
     xlab = "Number of PCs", ylab = "Cumulative variance (%)",
     main = "Cumulative variance explained")
abline(h = 80, col = "blue", lty = 2); abline(v = full$n80, col = "blue", lty = 3)
dev.off()

# correlation heatmap (ordered by block)
png(file.path(diag_dir, "corr_heatmap.png"), width = 1200, height = 1100, res = 120)
ord <- order(block, vars); Rord <- R_pear[ord, ord]
par(mar = c(9, 9, 3, 2))
image(seq_len(p_all), seq_len(p_all), t(Rord)[, p_all:1], zlim = c(-1, 1),
      col = hcl.colors(41, "Blue-Red 3"), axes = FALSE, xlab = "", ylab = "",
      main = "Pearson correlation (ordered by block)")
axis(1, at = seq_len(p_all), labels = colnames(Rord), las = 2, cex.axis = 0.55)
axis(2, at = seq_len(p_all), labels = rev(colnames(Rord)), las = 2, cex.axis = 0.55)
box()
dev.off()

# save everything machine-readable
saveRDS(list(desc = desc, corr_pearson = R_pear, corr_spearman = R_spear, vif = vif,
             block_tbl = block_tbl, block_x_block = bxb, pca = pca_tbl,
             loadings = load_mat, forward = fwd_tbl, kept = kept, dropped = dropped,
             criteria = c(n80 = full$n80, n90 = full$n90, kaiser = full$kaiser,
                          broken_stick = full$bstick)),
        file.path(diag_dir, "metrics.rds"))

hr("DONE")
cat("Diagnostics written to: ", diag_dir, "\n", sep = "")
cat("Key files: report.txt, 0*_*.csv, scree.png, cum_variance.png, corr_heatmap.png, metrics.rds\n")
