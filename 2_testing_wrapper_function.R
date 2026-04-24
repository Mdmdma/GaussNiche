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

# testing functions
source("virtualSpecies_fn.R")

sp1 <- virtualSpecies(
  dt          = dt,
  envData     = envData,
  mu          = c(2.0, -1.0),
  sigma_PC1   = 1.5,
  sigma_PC2   = 1.0,
  rho         = 0.2,
  bgk_prev    = 1,
  pa_sampler  = pa_random,    # <-- swap this for future strategies
  pc1_lab     = pc1_lab,
  pc2_lab     = pc2_lab,
  kde_df      = kde_df        # pass pre-computed KDE from section 2
)

# inspect metrics
sp1$metrics

# print plots
sp1$plots$p_niche
sp1$plots$p_pa
sp1$plots$p_bias_PC1 + sp1$plots$p_bias_PC2

# geographic maps
plot(sp1$suit_rast, col = hcl.colors(100, "YlOrRd", rev = TRUE),
     main = "Suitability (geographic space)")
plot(sp1$pa_rast, col = c("grey85", "firebrick"),
     main = paste0("Presences — prevalence = ", sp1$metrics$prevalence))
plot(sp1$suit_rast, col = hcl.colors(100, "YlOrRd", rev = TRUE),
     main = "Pseudo-absences in geographic space")
points(sp1$pseudo_vect, col = "steelblue", cex = 0.4, pch = 16)

