# Reproducible example for virtualSpecies_fn.R


# Uncomment if the wd is not the repo. Leave commented in commits for transportability
setwd("/home/mat/Desktop/semesterarbeit10/GaussNiche/")

library(terra)
# library(geodata)
library(USE.MCMC)
library(mvtnorm)
library(MASS)
library(ggplot2)
library(viridis)
library(hypervolume)
library(patchwork)
library(sf)                          # needed by pa_uniform
library(tictoc)
library(future)                      # parallel backend
library(furrr)                       # furrr::future_pmap dispatch
source("virtualSpecies_fn.R")

#1. ENVIRONMENTAL DATA AND PCA ----
 envData <- rast(USE.MCMC::Worldclim_tmp, type = "xyz")
rpc     <- rastPCA(envData, stand = TRUE)

# Full background: every non-NA cell projected into PC space
dt <- na.omit(as.data.frame(rpc$PCs[[c("PC1", "PC2")]], xy = TRUE))

# Axis labels with variance explained
var_exp <- rpc$pca$sdev^2 / sum(rpc$pca$sdev^2)
cat(var_exp)
pc1_lab <- paste0("PC1 (", round(var_exp[1] * 100, 1), "% var)")
pc2_lab <- paste0("PC2 (", round(var_exp[2] * 100, 1), "% var)")

cat("Background cells (n):", nrow(dt), "\n")
cat("PC1 range:", round(range(dt$PC1), 2), "\n")
cat("PC2 range:", round(range(dt$PC2), 2), "\n")

# 2. BACKGROUND KDE (pre-computed once, reused across all species)----
kde        <- MASS::kde2d(dt$PC1, dt$PC2, n = 200)
kde_df     <- expand.grid(PC1 = kde$x, PC2 = kde$y)
kde_df$density <- as.vector(kde$z)
kde_df$density <- kde_df$density / max(kde_df$density)  # normalise to [0,1]

p_background <- ggplot(kde_df, aes(x = PC1, y = PC2)) +
  geom_contour_filled(aes(z = density),
                      breaks = c(0, 0.05, 0.25, 0.5, 0.75, 1),
                      alpha  = 0.80) +
  scale_fill_viridis_d(option = "mako", name = "Background\ndensity [0–1]",
                       direction = -1) +
  coord_equal() +
  labs(
    title    = "Available environmental space",
    subtitle = "Kernel density of background cells in PC1 x PC2",
    x = pc1_lab, y = pc2_lab
  ) +
  theme_classic(base_size = 13)

p_background

# 3. ONE-TIME COMPUTATIONS (fix across all species for comparability)----

## 3a. Hypervolume bandwidth from the FULL BACKGROUND ----
#     This is fixed for every species and every sampler so hypervolume
#     estimates are comparable across the experiment.
bw_background <- compute_bandwidth(dt)
cat("Background bandwidth (PC1, PC2):", round(bw_background, 5), "\n")

## 3b. Optimal USE grid resolution (pre-computed once for the E-space)----
#     optimRes needs the background projected into PC space as an sf object
#     (using USE's internal PC axes, which it derives from envData).
dt_sf  <- sf::st_as_sf(dt, coords = c("PC1", "PC2"))
res_opt <- USE.MCMC::optimRes(
  sdf      = dt_sf,
  grid.res = 1:10,
  perc.thr = 20,
  showOpt  = FALSE,
  cr       = 5
)
grid_res_opt <- as.numeric(res_opt$Opt_res)
cat("Optimal USE grid resolution:", grid_res_opt, "\n")

# 4. VIRTUAL SPECIES — run with both samplers----

# WHEN WILL THE TWO SAMPLERS VISIBLY DIFFER?
#   pa_random draws uniformly from the full background (already approximately
#   uniform in E-space since dt = all raster cells).
#   pa_uniform adds an exclusion step: it removes the high-suitability region
#   from the sampling pool (controlled by thres). The difference is therefore
#   largest when:
#     (a) the niche is BROAD relative to the PC space (mu near centre, large sigma)
#     (b) thres is LOW (e.g. 0.25) so more of the suitable E-space is excluded
#   With mu = c(4,1) and narrow sigmas, the suitable region is small (~3-5%
#   of E-space) and the two samplers look similar.  Try mu = c(0,0) and
#   sigma_PC1 = 3, sigma_PC2 = 2 with thres = 0.25 to see a clear divergence.
#
# NOTE ON RUNTIME:
#   n_realizations = 50 means ~100 hypervolume_gaussian() calls per sampler.
#   Reduce to n_realizations = 10 for quick exploration.
# debugonce(virtualSpecies)
tic()
sp1 <- virtualSpecies(
  dt               = dt,
  envData          = envData,
  mu               = c(4, 1),
  sigma_PC1        = 3,
  sigma_PC2        = 2,
  rho              = 0,
  max_pres         = 1000, # keep 1000 to do not oversample the env space
  bgk_prev         = 1,
  pa_samplers      = list(
    random  = pa_random,
    uniform = pa_uniform,
    mcmc    = pa_mcmc
  ),
  n_realizations   = 10,          # set to 10 for rapid exploration
  seed_base        = 42,
  seed_pseudo_base = 123,
  pc1_lab          = pc1_lab,
  pc2_lab          = pc2_lab,
  kde_df           = kde_df,
  bw               = bw_background,   # fixed bandwidth from section 3a
  pa_env_rast      = envData,          # MUST be original vars, NOT rpc$PCs
  verbose          = TRUE,
  parallel         = TRUE,           # dispatch (sampler, realisation) tasks
  n_workers        = NULL,           # NULL -> detectCores() - 1
  # extra args forwarded to pa_uniform:
  grid.res         = grid_res_opt,    # from section 3b
  thres            = 0.75,
  # extra args forwarded to pa_mcmc:
  chain.length     = 20000,
  burnIn           = 1000,
  species.cutoff.threshold = 0.1   # USE.MCMC species GMM percentile
)
toc()

# 5. INSPECT NICHE EQUATIONS----
cat("\n── Species niche equations ────────────────────────────────────\n")
cat("PC1  :", sp1$niche$equations$PC1,   "\n")
cat("PC2  :", sp1$niche$equations$PC2,   "\n")
cat("2-D  :", sp1$niche$equations$two_d, "\n")
cat("Logit:", sp1$niche$equations$logit, "\n")

# 6. RESPONSE CURVE PLOTS----
# Suitability response curves (marginal, other axis fixed at optimum)
sp1$response_curves$plots$p_rc_PC1 + sp1$response_curves$plots$p_rc_PC2

# Logit of suitability (shows the non-linear transformation clearly)
sp1$response_curves$plots$p_rc_logit_PC1 + sp1$response_curves$plots$p_rc_logit_PC2

# Combined 2x2 panel
# (sp1$response_curves$plots$p_rc_PC1 + sp1$response_curves$plots$p_rc_PC2) /
#   (sp1$response_curves$plots$p_rc_logit_PC1 + sp1$response_curves$plots$p_rc_logit_PC2)

# 7. NICHE PLOT (shared across samplers)----
sp1$plots$p_niche

# 8. PER-SAMPLER DIAGNOSTIC PLOTS----

## Random sampler----
sp1$samplers$random$plots$p_pa
sp1$samplers$random$plots$p_bias_PC1 + sp1$samplers$random$plots$p_bias_PC2
sp1$samplers$random$plots$p_box_overlap
sp1$samplers$random$plots$p_box_coverage
sp1$samplers$random$plots$p_box_trueabs

## Uniform sampler----
sp1$samplers$uniform$plots$p_pa
sp1$samplers$uniform$plots$p_bias_PC1 + sp1$samplers$uniform$plots$p_bias_PC2
sp1$samplers$uniform$plots$p_box_overlap
sp1$samplers$uniform$plots$p_box_coverage
sp1$samplers$uniform$plots$p_box_trueabs

## MCMC sampler----
sp1$samplers$mcmc$plots$p_pa
sp1$samplers$mcmc$plots$p_bias_PC1 + sp1$samplers$mcmc$plots$p_bias_PC2
sp1$samplers$mcmc$plots$p_box_overlap
sp1$samplers$mcmc$plots$p_box_coverage
sp1$samplers$mcmc$plots$p_box_trueabs

# Side-by-side niche plots across samplers
sp1$samplers$random$plots$p_pa + sp1$samplers$uniform$plots$p_pa + sp1$samplers$mcmc$plots$p_pa

# 9. INSPECT METRICS TABLES----
# Per-sampler task status tally — spot samplers that returned no successful
# realisations (skip_few_pres / skip_sampler_null_or_short / skip_hv_fail).
sapply(sp1$samplers, function(s) table(s$metrics$status))

# Summary statistics across 50 Bernoulli realisations — random sampler
summary(sp1$samplers$random$metrics)

# Summary — uniform sampler
summary(sp1$samplers$uniform$metrics)

# rbind is shape-safe now: every metrics table has the same columns, including
# rows tagged with status = "skip_*" for realisations that were abandoned.
metrics_combined <- do.call(rbind, lapply(names(sp1$samplers), function(nm) {
  cbind(sampler = nm, sp1$samplers[[nm]]$metrics)
}))
metrics_combined$sampler <- factor(metrics_combined$sampler,
                                   levels = names(sp1$samplers))

# Only successful realisations contribute to cross-sampler comparisons.
metrics_ok <- subset(metrics_combined, status == "ok")

# Cross-sampler boxplot of overlap (ggplot2)
ggplot(metrics_ok, aes(x = sampler, y = overlap, fill = sampler)) +
  geom_boxplot(alpha = 0.7, width = 0.45, outlier.shape = 21) +
  geom_jitter(width = 0.1, alpha = 0.4, size = 0.9) +
  scale_fill_brewer(palette = "Set1") +
  scale_x_discrete(drop = FALSE) +
  labs(title    = "Hypervolume overlap across samplers",
       subtitle = paste0(sum(sp1$samplers$random$metrics$status == "ok"),
                         " successful Bernoulli realisations (random sampler)"),
       x = NULL, y = "Overlap (intersection volume)") +
  theme_classic(base_size = 13) +
  theme(legend.position = "none")

# Cross-sampler comparison of ranges
ggplot(metrics_ok, aes(x = sampler, y = rel_cov_PC1, fill = sampler)) +
  geom_boxplot(alpha = 0.7, width = 0.45, outlier.shape = 21) +
  geom_jitter(width = 0.1, alpha = 0.4, size = 0.9) +
  scale_fill_brewer(palette = "Set1") +
  scale_x_discrete(drop = FALSE) +
  # labs(title    = "Proportion of pseudo-absences on true-absence cells",
       # subtitle = "Higher = sampler avoids true-presence cells more effectively",
       # x = NULL, y = "coverage PC1") +
  theme_classic(base_size = 13) +
  theme(legend.position = "none")

ggplot(metrics_ok, aes(x = sampler, y = rel_cov_PC2, fill = sampler)) +
  geom_boxplot(alpha = 0.7, width = 0.45, outlier.shape = 21) +
  geom_jitter(width = 0.1, alpha = 0.4, size = 0.9) +
  scale_fill_brewer(palette = "Set1") +
  scale_x_discrete(drop = FALSE) +
  # labs(title    = "Proportion of pseudo-absences on true-absence cells",
  # subtitle = "Higher = sampler avoids true-presence cells more effectively",
  # x = NULL, y = "coverage PC1") +
  theme_classic(base_size = 13) +
  theme(legend.position = "none")

# Cross-sampler comparison of prop_true_abs
ggplot(metrics_ok, aes(x = sampler, y = prop_true_abs, fill = sampler)) +
  geom_boxplot(alpha = 0.7, width = 0.45, outlier.shape = 21) +
  geom_jitter(width = 0.1, alpha = 0.4, size = 0.9) +
  scale_fill_brewer(palette = "Set1") +
  scale_x_discrete(drop = FALSE) +
  labs(title    = "Proportion of pseudo-absences on true-absence cells",
       subtitle = "Higher = sampler avoids true-presence cells more effectively",
       x = NULL, y = "Proportion") +
  theme_classic(base_size = 13) +
  theme(legend.position = "none")

# 10. GEOGRAPHIC MAPS (reference realisation) ----
plot(sp1$suit_rast,
     col  = hcl.colors(100, "YlOrRd", rev = TRUE),
     main = "Suitability (geographic space)")

plot(sp1$pa_rast,
     col  = c("grey85", "firebrick"),
     main = paste0("Presences — realisation 1"))

# Geographic pseudo-absences: random
plot(sp1$suit_rast, col = hcl.colors(100, "YlOrRd", rev = TRUE),
     main = "Random pseudo-absences (geographic space)")
points(sp1$samplers$random$pseudo_vect, col = "steelblue", cex = 0.4, pch = 16)

# Geographic pseudo-absences: uniform
plot(sp1$suit_rast, col = hcl.colors(100, "YlOrRd", rev = TRUE),
     main = "Uniform pseudo-absences (geographic space)")
points(sp1$samplers$uniform$pseudo_vect, col = "darkorange", cex = 0.4, pch = 16)

# Geographic pseudo-absences: mcmc
plot(sp1$suit_rast, col = hcl.colors(100, "YlOrRd", rev = TRUE),
     main = "MCMC pseudo-absences (geographic space)")
points(sp1$samplers$mcmc$pseudo_vect, col = "purple", cex = 0.4, pch = 16)

plot(sp1$pa_rast,
     col  = c("grey85", "firebrick"),
     main = paste0("Presences — realisation 1"))
points(sp1$samplers$random$pseudo_vect, col = "steelblue", cex = 0.8, pch = 16)

plot(sp1$pa_rast,
     col  = c("grey85", "firebrick"),
     main = paste0("Presences — realisation 1"))
points(sp1$samplers$uniform$pseudo_vect, col = "darkorange", cex = 0.8, pch = 16)

#one of the issue here might be that we are using all the presences!!!

# 11. ACCESS COMBINED DATASETS FOR MODELLING----

# dataset for each sampler: pa = 1 (presence) or 0 (pseudo-absence)
head(sp1$samplers$random$dataset)
head(sp1$samplers$uniform$dataset)

# Bandwidth used (reuse for other species to keep hypervolumes comparable)
sp1$bw


