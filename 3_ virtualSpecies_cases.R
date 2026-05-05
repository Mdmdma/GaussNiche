# =============================================================================
# Four virtual species arranged as a 2x2 factorial:
#
#                     | Common environment  | Rare environment
#                     | (mu near PC origin) | (mu at PC periphery)
#   ------------------|---------------------|---------------------
#   Generalist        |       sp1           |       sp3
#   (broad niche)     |                     |
#   ------------------|---------------------|---------------------
#   Specialist        |       sp2           |       sp4
#   (narrow niche)    |                     |
#
# "Common environment"  = mu near (0, 0), the high-density region of E-space
# "Rare environment"    = mu near (4, 1), the low-density periphery
# "Generalist"          = sigma_PC1=3, sigma_PC2=2  (broad tolerance)
# "Specialist"          = sigma_PC1=1, sigma_PC2=0.5 (narrow tolerance)
#
# The two axes stress-test the pseudo-absence samplers differently:
#   - Breadth controls how much E-space the USE exclusion filter removes
#   - Position controls the overlap between background density and the niche
# =============================================================================

# Shared call arguments — only mu and sigmas change across species
shared_args <- list(
  dt               = dt,
  envData          = envData,
  rho              = 0,
  bgk_prev         = 1,
  pa_samplers      = list(random = pa_random, uniform = pa_uniform),
  n_realizations   = 5,
  max_pres         = 1000,
  seed_base        = 42,
  seed_pseudo_base = 123,
  pc1_lab          = pc1_lab,
  pc2_lab          = pc2_lab,
  kde_df           = kde_df,
  bw               = bw_background,
  pa_env_rast      = envData,
  verbose          = TRUE,
  grid.res         = grid_res_opt,
  thres            = 0.75
)

# =============================================================================
# sp1 — Generalist, common environment
# Broad niche centred in the densest region of E-space.
# Most background cells are potential presences; USE exclusion removes a large
# fraction → strongest contrast between random and uniform PA sampling.
# =============================================================================
sp1 <- do.call(virtualSpecies, c(shared_args, list(
  mu        = c(0, 0),
  sigma_PC1 = 3,
  sigma_PC2 = 2
)))

# =============================================================================
# sp2 — Specialist, common environment
# Narrow niche in a common environment: many background cells surround the niche
# but the suitable region is small. High prevalence relative to niche volume.
# =============================================================================
sp2 <- do.call(virtualSpecies, c(shared_args, list(
  mu        = c(0, 0),
  sigma_PC1 = 1,
  sigma_PC2 = 0.5
)))

# =============================================================================
# sp3 — Generalist, rare environment
# Broad niche at the periphery of E-space. Few background cells near the
# optimum, so prevalence is low despite wide tolerance. The niche extends
# into dense background regions, creating an asymmetric suitability gradient.
# =============================================================================
sp3 <- do.call(virtualSpecies, c(shared_args, list(
  mu        = c(4, 1),
  sigma_PC1 = 3,
  sigma_PC2 = 2
)))

# =============================================================================
# sp4 — Specialist, rare environment
# Narrow niche at the periphery: low prevalence, small suitable volume, and
# few background cells near mu. Hardest case for both samplers — small
# presence set and sparse background in the relevant E-space region.
# =============================================================================
sp4 <- do.call(virtualSpecies, c(shared_args, list(
  mu        = c(4, 1),
  sigma_PC1 = 1,
  sigma_PC2 = 0.5
)))

# =============================================================================
# COLLECT ALL SPECIES FOR COMPARISON
# =============================================================================

species_list <- list(
  sp1_generalist_common    = sp1,
  sp2_specialist_common    = sp2,
  sp3_generalist_rare      = sp3,
  sp4_specialist_rare      = sp4
)

# Quick prevalence summary
prevalence_summary <- sapply(species_list, function(sp)
  round(mean(sp$background$pa), 3))
print(prevalence_summary)