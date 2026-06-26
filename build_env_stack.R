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

# Threaded GDAL warp/resample -- value-neutral (changes scheduling, not numbers),
# a free speedup on the per-block harmonise(). No furrr layer in this script, so a
# multithreaded BLAS is also fine here -- do NOT pin OMP/OPENBLAS=1 in its sbatch.
Sys.setenv(GDAL_NUM_THREADS = "ALL_CPUS")

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
# ENV_OUT_DIR lets you build to a SANDBOX (e.g. to bit-identity-check a change)
# without clobbering the live env5d the ablation reads. Default = live env5d.
out_dir      <- Sys.getenv("ENV_OUT_DIR", file.path(scratch, "GaussNiche", "env5d"))
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

# Per-block cache. A block (soil / terrain / veg / anth) is re-harmonised only
# when its layer set (`layer_names`, the cache key) changed, or FORCE_REBUILD=1.
# So editing ONE block's variables re-runs only that block on the next build, not
# the whole 17-layer harmonise -- the common case when tweaking the environment
# slightly. Bit-identical: a cached block tif == a fresh harmonise (harmonise is
# deterministic given the geodata download cache), and the clean (no-cache) build
# is unchanged, so this never alters env_stack_candidate.tif.
build_block <- function(name, layer_names, builder) {
  bt <- file.path(out_dir, paste0("block_", name, ".tif"))
  mf <- file.path(out_dir, paste0("block_", name, ".names"))
  if (!FORCE && file.exists(bt) && file.exists(mf) &&
      identical(readLines(mf, warn = FALSE), layer_names)) {
    message("  [cache hit] block '", name, "': ", paste(layer_names, collapse = ", "))
    r <- terra::rast(bt); names(r) <- layer_names; return(r)
  }
  r <- builder(); names(r) <- layer_names
  terra::writeRaster(r, bt, overwrite = TRUE); writeLines(layer_names, mf)
  message("  [rebuilt] block '", name, "'")
  r
}

# =============================================================================
# 2. SOIL  (SoilGrids 2.0 via geodata; 30 arc-sec -> /20 to 10 arc-min)
# =============================================================================
message("\n== SOIL ==")
soil <- build_block("soil", paste0("soil_", SOIL_VARS), function()
  terra::rast(lapply(SOIL_VARS, function(v) {
    message("  soil_world: ", v, " (depth ", SOIL_DEPTH, " cm)")
    harmonise(geodata::soil_world(var = v, depth = SOIL_DEPTH, stat = "mean",
                                  path = path_geodata), agg = 20)
  })))

# =============================================================================
# 3. TERRAIN  (SRTM elevation via geodata at res=10; derive on the global grid)
# =============================================================================
message("\n== TERRAIN ==")
TER_NAMES <- c("ter_elevation", "ter_slope", "ter_TRI", "ter_TPI",
               "ter_roughness", "ter_eastness", "ter_northness")
terr <- build_block("terr", TER_NAMES, function() {
  elev <- geodata::elevation_global(res = 10, path = path_geodata)
  # #6a: crop to a BUFFERED European window BEFORE deriving terrain. slope/TRI/
  # TPI/roughness use a 3x3 neighbourhood, so a generous halo (1 deg ~ 6 cells at
  # 10 arc-min, >> the 1-cell focal radius) keeps the European cells bit-identical
  # to deriving on the global grid; the final harmonise() mask drops the halo.
  # Avoids running all 5 terrain passes over the whole globe (the biggest waste).
  elev <- terra::crop(elev, terra::ext(terra::xmin(eu_ext) - 1, terra::xmax(eu_ext) + 1,
                                       terra::ymin(eu_ext) - 1, terra::ymax(eu_ext) + 1))
  # Aspect is circular -> split into eastness/northness before PCA. Elevation is
  # kept even though it is lapse-rate-coupled to temperature: the analysis SHOWS
  # that correlation rather than us assuming it.
  asp <- terra::terrain(elev, v = "aspect", unit = "radians")
  terr0 <- c(elev,
             terra::terrain(elev, v = "slope", unit = "degrees"),
             terra::terrain(elev, v = "TRI"),
             terra::terrain(elev, v = "TPI"),
             terra::terrain(elev, v = "roughness"),
             sin(asp), cos(asp))
  harmonise(terr0)   # res=10 already: crop + snap + mask, no aggregate
})

# =============================================================================
# 4. VEGETATION  (ESA WorldCover fractions via geodata; 30 arc-sec -> /20)
# =============================================================================
message("\n== VEGETATION ==")
veg <- build_block("veg", paste0("veg_", VEG_VARS), function()
  terra::rast(lapply(VEG_VARS, function(v) {
    message("  landcover: ", v)
    harmonise(geodata::landcover(var = v, path = path_geodata), agg = 20)
  })))

# =============================================================================
# 5. ANTHROPOGENIC  (Human Footprint via geodata; 30 arc-sec -> /20)
# =============================================================================
anth <- NULL
if (INCLUDE_FOOTPRINT) {
  message("\n== FOOTPRINT ==")
  anth <- build_block("anth", "anth_footprint", function()
    harmonise(geodata::footprint(year = 2009, path = path_geodata), agg = 20))
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
