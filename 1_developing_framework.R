# Virtual species in PC-environmental space
# pseudo-absence sampling + class overlap + sampling bias diagnostics

setwd("/home/dared/Downloads/")

library(terra)
library(geodata)
library(USE)
library(mvtnorm)
library(MASS)
library(ggplot2)
library(viridis)
library(hypervolume)  
library(patchwork)

# 1. ENVIRONMENTAL DATA AND PCA ----
envData <- rast(USE::Worldclim_tmp, type = "xyz")
rpc     <- rastPCA(envData, stand = TRUE)

# full background: every non-NA cell in the raster
dt <- na.omit(as.data.frame(rpc$PCs[[c("PC1", "PC2")]], xy = TRUE))

# axis labels with variance explained
var_exp <- rpc$pca$sdev^2 / sum(rpc$pca$sdev^2)
pc1_lab <- paste0("PC1 (", round(var_exp[1] * 100, 1), "% var)")
pc2_lab <- paste0("PC2 (", round(var_exp[2] * 100, 1), "% var)")

cat("Background cells (n):", nrow(dt), "\n")
cat("PC1 range:", round(range(dt$PC1), 2), "\n")
cat("PC2 range:", round(range(dt$PC2), 2), "\n")

# 2. VISUALISE BACKGROUND IN ENVIRONMENTAL SPACE ----
kde        <- MASS::kde2d(dt$PC1, dt$PC2, n = 200)
kde_df     <- expand.grid(PC1 = kde$x, PC2 = kde$y)
kde_df$density <- as.vector(kde$z)

p_background <- ggplot(kde_df, aes(x = PC1, y = PC2)) +
  geom_contour_filled(aes(z = density), alpha = 0.9) +
  scale_fill_viridis_d(option = "mako", name = "Background\ndensity", direction = -1) +
  coord_equal() +
  labs(
    title    = "Available environmental space",
    subtitle = "Kernel density of background cells in PC1 × PC2",
    x = pc1_lab, y = pc2_lab
  ) +
  theme_classic(base_size = 13)

p_background

# 3. VIRTUAL SPECIES — DEFINE THE NICHE (bivariate Gaussian) ----
# mu    = niche optimum in PC space
# Sigma = niche breadth: large variance → generalist, small → specialist
# rho   = tolerance correlation between PC1 and PC2
mu        <- c(2.0, -1.0)    # <-- TUNE: niche optimum; look at background first
sigma_PC1 <- 1.5              # <-- TUNE: breadth along PC1
sigma_PC2 <- 1.0              # <-- TUNE: breadth along PC2
rho       <- 0.2              # <-- TUNE: axis correlation [-1, 1]

Sigma <- matrix(
  c(sigma_PC1^2,
    rho * sigma_PC1 * sigma_PC2,
    rho * sigma_PC1 * sigma_PC2,
    sigma_PC2^2),
  nrow = 2
)

# normalise by the theoretical peak at mu: suitability = 1.0 at optimum always,
# regardless of how much background falls near mu
peak_density <- mvtnorm::dmvnorm(x = matrix(mu, nrow = 1), mean = mu, sigma = Sigma)

suitability <- mvtnorm::dmvnorm(
  x     = as.matrix(dt[, c("PC1", "PC2")]),
  mean  = mu,
  sigma = Sigma
) / peak_density

cat("Suitability range:", round(range(suitability), 4), "\n")

# 4. BERNOULLI SAMPLING → presence / absence ----
set.seed(42)
pa      <- rbinom(n = length(suitability), size = 1, prob = suitability)
dt$suit <- suitability
dt$pa   <- pa

cat("Prevalence:", round(mean(pa), 3), "\n")
cat("N presences:", sum(pa), "of", nrow(dt), "cells\n")

# 5. OVERLAY: niche vs background in environmental space ----
grid_pc      <- expand.grid(
  PC1 = seq(min(dt$PC1), max(dt$PC1), length.out = 300),
  PC2 = seq(min(dt$PC2), max(dt$PC2), length.out = 300)
)
grid_pc$suit <- mvtnorm::dmvnorm(as.matrix(grid_pc), mean = mu, sigma = Sigma) / peak_density

p_niche <- ggplot() +
  geom_contour_filled(data = kde_df, aes(x = PC1, y = PC2, z = density), alpha = 0.75) +
  scale_fill_viridis_d(option = "mako", name = "Background\ndensity", direction = -1) +
  geom_contour(data = grid_pc, aes(x = PC1, y = PC2, z = suit),
               colour = "firebrick", linewidth = 0.7,
               breaks = c(0.1, 0.25, 0.5, 0.75, 0.95)) +
  geom_point(data = subset(dt, pa == 1), aes(x = PC1, y = PC2),
             colour = "gray40", alpha = 0.25, size = 0.6) +
  annotate("point", x = mu[1], y = mu[2],
           colour = "firebrick", size = 4, shape = 3, stroke = 2) +
  coord_equal() +
  labs(
    title    = "Niche vs available environment",
    subtitle = "Red contours = suitability (0.1, 0.25, 0.50, 0.75, 0.95) | + = optimum | dots = presences",
    x = pc1_lab, y = pc2_lab
  ) +
  theme_classic(base_size = 13)

p_niche

# 6. RANDOM PSEUDO-ABSENCE SAMPLING IN ENVIRONMENTAL SPACE ----
# draw uniformly from the full background (all cells, including cells where pa = 1)
# this is the "random" strategy from 5_ClassOverlap.R / 7_SampleBias.R
# N is matched to presences at a 1:1 ratio (bgk_prev = 1)

bgk_prev <- 1                           # <-- TUNE: 1 = 1:1; 2 = 1:2 pres:abs
N_pa     <- round(sum(pa) / bgk_prev)

set.seed(123)
pseudo <- dt[sample(nrow(dt), size = N_pa, replace = FALSE), ]

cat("N pseudo-absences sampled:", nrow(pseudo), "\n")

# 7. CLASS OVERLAP IN PC SPACE ----
# from 5_ClassOverlap.R: build hypervolumes on PC scores of presences and
# pseudo-absences separately; [[3]] of get_volume() is the intersection volume
# NOTE: hypervolume is stochastic and slow on large N — subsample if needed

hyp_pres <- hypervolume(subset(dt, pa == 1)[, c("PC1", "PC2")], verbose = FALSE)
hyp_pa   <- hypervolume(pseudo[, c("PC1", "PC2")],               verbose = FALSE)

ovrlp <- get_volume(
  hypervolume_set(hyp_pres, hyp_pa, check.memory = FALSE, verbose = FALSE)
)[[3]]   # [[1]] pres, [[2]] abs, [[3]] intersection

cat("Class overlap (hypervolume intersection):", round(ovrlp, 4), "\n")

# 8. SAMPLING BIAS IN PC SPACE ----
# from 7_SampleBias.R: total range of PC1 and PC2 covered by pseudo-absences
# relative coverage = pseudo-absence range / background range (1 = perfect coverage)

bg_range_PC1  <- diff(range(dt$PC1))
bg_range_PC2  <- diff(range(dt$PC2))
pa_range_PC1  <- diff(range(pseudo$PC1))
pa_range_PC2  <- diff(range(pseudo$PC2))
rel_cov_PC1   <- round(pa_range_PC1 / bg_range_PC1, 3)
rel_cov_PC2   <- round(pa_range_PC2 / bg_range_PC2, 3)

cat("PC1 — background range:", round(bg_range_PC1, 2),
    "| PA range:", round(pa_range_PC1, 2),
    "| relative coverage:", rel_cov_PC1, "\n")
cat("PC2 — background range:", round(bg_range_PC2, 2),
    "| PA range:", round(pa_range_PC2, 2),
    "| relative coverage:", rel_cov_PC2, "\n")

# 9. VISUALISE: presences + pseudo-absences in environmental space ----
p_pa <- ggplot() +
  geom_contour_filled(data = kde_df, aes(x = PC1, y = PC2, z = density), alpha = 0.75) +
  scale_fill_viridis_d(option = "mako", name = "Background\ndensity", direction = -1) +
  geom_point(data = pseudo, aes(x = PC1, y = PC2),
             colour = "steelblue", alpha = 0.35, size = 0.7) +
  geom_point(data = subset(dt, pa == 1), aes(x = PC1, y = PC2),
             colour = "firebrick", alpha = 0.35, size = 0.7) +
  geom_contour(data = grid_pc, aes(x = PC1, y = PC2, z = suit),
               colour = "firebrick", linewidth = 0.5, linetype = "dashed",
               breaks = c(0.5, 0.95)) +
  annotate("text", x = min(dt$PC1), y = max(dt$PC2),
           label = paste0("Hypervolume overlap = ", round(ovrlp, 3)),
           hjust = 0, vjust = 1, size = 4, colour = "white") +
  coord_equal() +
  labs(
    title    = "Presences vs pseudo-absences in PC space",
    subtitle = "Red = presences | Blue = pseudo-absences | Dashed = suitability 0.50, 0.95",
    x = pc1_lab, y = pc2_lab
  ) +
  theme_classic(base_size = 13)

p_pa

# sampling bias: marginal density of Background, Presences, Pseudo-absences on each PC axis
bias_df <- rbind(
  data.frame(PC1 = dt$PC1,                       PC2 = dt$PC2,                       group = "Background"),
  data.frame(PC1 = pseudo$PC1,                   PC2 = pseudo$PC2,                   group = "Pseudo-absences"),
  data.frame(PC1 = subset(dt, pa == 1)$PC1,      PC2 = subset(dt, pa == 1)$PC2,      group = "Presences")
)
bias_df$group <- factor(bias_df$group, levels = c("Background", "Pseudo-absences", "Presences"))

col_vals  <- c("Background" = "grey50",   "Pseudo-absences" = "steelblue", "Presences" = "firebrick")
fill_vals <- c("Background" = "grey70",   "Pseudo-absences" = "steelblue", "Presences" = "firebrick")

p_bias_PC1 <- ggplot(bias_df, aes(x = PC1, colour = group, fill = group)) +
  geom_density(alpha = 0.25, linewidth = 0.8) +
  scale_colour_manual(values = col_vals)  +
  scale_fill_manual(values   = fill_vals) +
  annotate("text", x = Inf, y = Inf,
           label = paste0("Coverage = ", rel_cov_PC1),
           hjust = 1.1, vjust = 1.5, size = 4) +
  labs(title = "Sampling bias — PC1", x = pc1_lab, y = "Density",
       colour = NULL, fill = NULL) +
  theme_classic(base_size = 13) +
  theme(legend.position = "bottom")

p_bias_PC2 <- ggplot(bias_df, aes(x = PC2, colour = group, fill = group)) +
  geom_density(alpha = 0.25, linewidth = 0.8) +
  scale_colour_manual(values = col_vals)  +
  scale_fill_manual(values   = fill_vals) +
  annotate("text", x = Inf, y = Inf,
           label = paste0("Coverage = ", rel_cov_PC2),
           hjust = 1.1, vjust = 1.5, size = 4) +
  labs(title = "Sampling bias — PC2", x = pc2_lab, y = "Density",
       colour = NULL, fill = NULL) +
  theme_classic(base_size = 13) +
  theme(legend.position = "bottom")

p_bias_PC1 + p_bias_PC2

# 10. BACK-PROJECT TO GEOGRAPHIC SPACE ----
suit_rast   <- rast(dt[, c("x", "y", "suit")], type = "xyz", crs = crs(envData))
pa_rast     <- rast(dt[, c("x", "y", "pa")],   type = "xyz", crs = crs(envData))
pseudo_vect <- vect(pseudo, geom = c("x", "y"), crs = crs(envData))

plot(suit_rast,
     col  = hcl.colors(100, "YlOrRd", rev = TRUE),
     main = "Suitability (geographic space)")

plot(pa_rast,
     col  = c("grey85", "firebrick"),
     main = paste0("Presences — prevalence = ", round(mean(pa), 3)))

plot(suit_rast,
     col  = hcl.colors(100, "YlOrRd", rev = TRUE),
     main = "Pseudo-absences in geographic space")
points(pseudo_vect, col = "steelblue", cex = 0.4, pch = 16)
