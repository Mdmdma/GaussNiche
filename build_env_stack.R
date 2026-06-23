# build_env_stack.R  —  Phase 1 (build) of the 5d-niche branch
# =============================================================================
# Download + harmonise a BROAD candidate pool of environmental layers onto the
# WorldClim baseline grid (USE.MCMC::Worldclim_tmp: 10 arc-min, Central & Western
# Europe, 6 bioclim layers = bio1/3/9 temp + bio12/13/15 precip), then cache the
# stack to scratch. All sources via the already-installed `geodata` package.
#
# This script ONLY builds + caches (it does no PCA). The diagnostics that decide
# which layers/blocks to keep live in analyze_env_stack.R, which reads the cache.
# Re-runs are cheap: if the cached .tif exists it is loaded and the build skipped
# (set env FORCE_REBUILD=1 to force). geodata also caches each download.
#
# Candidate blocks (name prefix = block id used by the analysis):
#   clim_  6 WorldClim bioclim (baseline, kept)
#   soil_  SoilGrids 2.0  geodata::soil_world()   pH, clay, sand, bulk dens, SOC, N
#   ter_   SRTM elevation geodata::elevation_global(res=10) + terra::terrain()
#   veg_   ESA WorldCover geodata::landcover()    trees/grass/shrubs/bare/cropland
#   anth_  Human Footprint geodata::footprint()   (anthropogenic candidate)
#
# Output: <scratch>/GaussNiche/env5d/env_stack_candidate.tif
#
# Run on a COMPUTE NODE (never login). See sbatch/submit_env_analysis.sh.
# =============================================================================

suppressPackageStartupMessages({
  .libPaths(c("~/R/rocker-rstudio/4.5", .libPaths()))
  library(terra)
  library(geodata)
  library(USE.MCMC)
})
cat("terra", as.character(packageVersion("terra")),
    "| geodata", as.character(packageVersion("geodata")), "\n")

# --- network: compute nodes reach the internet only via the ETH proxy --------
if (Sys.getenv("https_proxy") == "" && Sys.getenv("HTTPS_PROXY") == "") {
  Sys.setenv(http_proxy  = "http://proxy.ethz.ch:3128",
             https_proxy = "http://proxy.ethz.ch:3128")
  message("Set ETH proxy fallback (http://proxy.ethz.ch:3128)")
}

# --- paths: bulk data + outputs on scratch, NEVER $HOME ----------------------
user         <- Sys.getenv("USER")
scratch      <- file.path("/cluster/scratch", user)
path_geodata <- file.path(scratch, "GaussNiche", "geodata")
out_dir      <- file.path(scratch, "GaussNiche", "env5d")
dir.create(path_geodata, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir,      recursive = TRUE, showWarnings = FALSE)
cand_tif <- file.path(out_dir, "env_stack_candidate.tif")

FORCE <- identical(Sys.getenv("FORCE_REBUILD"), "1")
if (file.exists(cand_tif) && !FORCE) {
  r <- terra::rast(cand_tif)
  message("Cached stack exists, skipping rebuild: ", cand_tif)
  message(sprintf("  %d layers: %s", terra::nlyr(r), paste(names(r), collapse = ", ")))
  quit(save = "no", status = 0)
}
message("geodata cache: ", path_geodata)
message("outputs:       ", out_dir)

# --- candidate pool (edit here to tune what gets built) ----------------------
SOIL_VARS  <- c("phh2o", "clay", "sand", "bdod", "soc", "nitrogen")  # 'cec' is NOT valid
SOIL_DEPTH <- 5L                                                     # topsoil 0-5 cm
VEG_VARS   <- c("trees", "grassland", "shrubs", "bare", "cropland")  # continuous fractions
INCLUDE_FOOTPRINT <- TRUE

# =============================================================================
# 1. BASELINE = target grid, extent, mask
# =============================================================================
envData <- terra::rast(USE.MCMC::Worldclim_tmp, type = "xyz")
if (is.na(terra::crs(envData)) || terra::crs(envData) == "")
  terra::crs(envData) <- "EPSG:4326"
names(envData) <- paste0("clim_", names(envData))
tmpl   <- envData[[1]]
eu_ext <- terra::ext(envData)
message(sprintf("Baseline: %d layers (%s); %d cells; res %.4f deg",
                terra::nlyr(envData), paste(names(envData), collapse = ","),
                terra::ncell(envData), terra::res(envData)[1]))

# crop -> (aggregate) -> snap onto template -> mask to baseline footprint
harmonise <- function(r, agg = NULL, fun = "mean") {
  r <- terra::crop(r, eu_ext)
  if (!is.null(agg)) r <- terra::aggregate(r, fact = agg, fun = fun, na.rm = TRUE)
  r <- terra::resample(r, tmpl, method = "bilinear")
  terra::mask(r, tmpl)
}

# =============================================================================
# 2. SOIL  (SoilGrids 2.0 via geodata; 30 arc-sec -> /20 to 10 arc-min)
# =============================================================================
message("\n== SOIL ==")
soil <- terra::rast(lapply(SOIL_VARS, function(v) {
  message("  soil_world: ", v, " (depth ", SOIL_DEPTH, " cm)")
  r <- geodata::soil_world(var = v, depth = SOIL_DEPTH, stat = "mean", path = path_geodata)
  harmonise(r, agg = 20)
}))
names(soil) <- paste0("soil_", SOIL_VARS)

# =============================================================================
# 3. TERRAIN  (SRTM elevation via geodata at res=10; derive on the global grid)
# =============================================================================
message("\n== TERRAIN ==")
elev <- geodata::elevation_global(res = 10, path = path_geodata)
# Aspect is circular -> split into eastness/northness before PCA. Elevation is
# kept as a candidate even though it is lapse-rate-coupled to temperature: the
# analysis will SHOW that correlation rather than us assuming it.
slp  <- terra::terrain(elev, v = "slope", unit = "degrees")
tri  <- terra::terrain(elev, v = "TRI")
tpi  <- terra::terrain(elev, v = "TPI")
rgh  <- terra::terrain(elev, v = "roughness")
asp  <- terra::terrain(elev, v = "aspect", unit = "radians")
east <- sin(asp)
north <- cos(asp)
terr <- c(elev, slp, tri, tpi, rgh, east, north)
names(terr) <- c("ter_elevation", "ter_slope", "ter_TRI", "ter_TPI",
                 "ter_roughness", "ter_eastness", "ter_northness")
terr <- harmonise(terr)   # res=10 already: crop + snap + mask, no aggregate

# =============================================================================
# 4. VEGETATION  (ESA WorldCover fractions via geodata; 30 arc-sec -> /20)
# =============================================================================
message("\n== VEGETATION ==")
veg <- terra::rast(lapply(VEG_VARS, function(v) {
  message("  landcover: ", v)
  r <- geodata::landcover(var = v, path = path_geodata)
  harmonise(r, agg = 20)
}))
names(veg) <- paste0("veg_", VEG_VARS)

# =============================================================================
# 5. ANTHROPOGENIC  (Human Footprint via geodata; 30 arc-sec -> /20)
# =============================================================================
anth <- NULL
if (INCLUDE_FOOTPRINT) {
  message("\n== FOOTPRINT ==")
  fp <- geodata::footprint(year = 2009, path = path_geodata)
  anth <- harmonise(fp, agg = 20)
  names(anth) <- "anth_footprint"
}

# =============================================================================
# 6. COMBINE + CACHE
# =============================================================================
message("\n== COMBINE ==")
allEnv <- c(envData, soil, terr, veg)
if (!is.null(anth)) allEnv <- c(allEnv, anth)
allEnv <- terra::mask(allEnv, tmpl)

n_base <- terra::global(!is.na(tmpl), "sum")[1, 1]
n_keep <- terra::global(!any(is.na(allEnv)), "sum")[1, 1]
message(sprintf("Cells: baseline=%d, complete-across-all-layers=%d (%.1f%%)",
                n_base, n_keep, 100 * n_keep / n_base))
message(sprintf("Candidate stack: %d layers: %s",
                terra::nlyr(allEnv), paste(names(allEnv), collapse = ", ")))

terra::writeRaster(allEnv, cand_tif, overwrite = TRUE)
message("Wrote ", cand_tif)
