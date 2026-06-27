#!/bin/bash
#SBATCH --job-name=usemcmc-gated-check
#SBATCH --account=es_schin
#SBATCH --output=output/gated-check-%j.out
#SBATCH --error=output/gated-check-%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4G
#SBATCH --time=01:30:00

# Verify the vignette chunks that are eval=FALSE / gated for CRAN compatibility
# actually RUN (they are never exercised by R CMD check). Each gated chunk's machinery
# is run standalone in the apptainer env with the ETH proxy on and the optional packages
# installed: the geodata WorldClim download (USE chunks 34/46), the optimRes call
# (USE 112), rnaturalearth::ne_countries (NN 41/73), and a gganimate animation +
# anim_save built from the exact animate-sample-chain chunk pattern (MCMC). Each block
# is tryCatch'd and reported PASS/FAIL; the job never aborts on a single failure.

set -uo pipefail
cd "$HOME/GaussNiche"; mkdir -p output
source "$HOME/.config/euler/jupyterhub/config_r_studio"
module load eth_proxy 2>/dev/null || true
SIF="/cluster/scratch/$USER/rocker_rstudio_4.5.sif"
[[ -f "$SIF" ]] || { echo "Missing $SIF"; exit 2; }
LIB="$HOME/R/rocker-rstudio/4.5"
WORK="/cluster/scratch/$USER/usemcmc_gated"; rm -rf "$WORK"; mkdir -p "$WORK"

apptainer exec \
  --bind "/cluster/home/$USER,/cluster/software,/cluster/scratch/$USER" \
  "$SIF" bash -lc "
    export http_proxy='http://proxy.ethz.ch:3128' https_proxy='http://proxy.ethz.ch:3128'
    export HTTP_PROXY=\$http_proxy HTTPS_PROXY=\$https_proxy
    export R_LIBS='$LIB'
    echo '=== install optional packages used by the gated chunks (best-effort) ==='
    Rscript -e 'p <- c(\"gganimate\",\"gifski\",\"transformr\",\"rnaturalearthdata\")
      m <- p[!vapply(p, requireNamespace, logical(1), quietly=TRUE)]
      if (length(m)) try(install.packages(m, repos=\"https://cloud.r-project.org\", lib=\"$LIB\"))
      for (q in p) cat(sprintf(\"  %-18s %s\n\", q, if (requireNamespace(q, quietly=TRUE)) \"OK\" else \"MISSING\"))'
    cd '$WORK'
    Rscript -e '
      .libPaths(c(\"$LIB\", .libPaths())); library(USE.MCMC)
      ok <- function(tag, expr) {
        r <- tryCatch({ force(expr); \"PASS\" },
                      error = function(e) paste0(\"FAIL: \", conditionMessage(e)))
        cat(sprintf(\"[%s] %s\n\", r, tag))
      }
      env <- terra::rast(USE.MCMC::Worldclim_tmp, type = \"xyz\")

      ## USE chunks 34/46 : geodata WorldClim download + names() + crop()
      ok(\"USE geodata::worldclim_global(res=10) + names + crop\", {
        wc <- geodata::worldclim_global(var = \"bio\", res = 10, path = tempdir())
        names(wc) <- paste0(\"bio\", seq_len(terra::nlyr(wc)))
        terra::crop(wc, terra::ext(-12, 25, 36, 60))
      })

      ## USE chunk 112 : optimRes on the PC-space sf (myRes is hardcoded in the build)
      ok(\"USE optimRes(grid.res=1:10, perc.thr=20, cr=5)\", {
        rpc <- USE.MCMC::rastPCA(env, stand = TRUE)
        dt  <- na.omit(as.data.frame(rpc\$PCs[[c(\"PC1\",\"PC2\")]], xy = TRUE))
        dt  <- sf::st_as_sf(dt, coords = c(\"PC1\",\"PC2\"))
        USE.MCMC::optimRes(sdf = dt, grid.res = 1:10, perc.thr = 20, cr = 5, showOpt = FALSE)
      })

      ## NN chunks 41/73 : rnaturalearth::ne_countries(scale=medium) (needs rnaturalearthdata)
      ok(\"NN rnaturalearth::ne_countries(scale=medium)\", {
        if (!requireNamespace(\"rnaturalearthdata\", quietly = TRUE))
          stop(\"rnaturalearthdata not installed\")
        cs <- rnaturalearth::ne_countries(scale = \"medium\", returnclass = \"sf\")
        cs[cs\$admin == \"Italy\", ]
      })

      ## MCMC animate-* chunks : exact gganimate pattern from animate-sample-chain
      ok(\"MCMC gganimate animate() + anim_save() (animation toolchain)\", {
        if (!requireNamespace(\"gganimate\", quietly = TRUE) ||
            !requireNamespace(\"gifski\", quietly = TRUE)) stop(\"gganimate/gifski not installed\")
        library(ggplot2); library(gganimate)
        sp <- data.frame(x = cumsum(rnorm(30)), y = cumsum(rnorm(30)),
                         density = 1, step = 1:30)
        a <- ggplot(sp) + geom_point(aes(x = y, y = x, color = step)) +
          geom_path(aes(x = y, y = x, color = step)) + theme_minimal() +
          transition_reveal(step) + ease_aes(\"linear\")
        g <- gganimate::animate(a, renderer = gifski::gifski_renderer(),
                                width = 300, height = 250, fps = 5, duration = 2,
                                nframes = 10)
        gganimate::anim_save(file.path(tempdir(), \"smoke.gif\"), g)
        stopifnot(file.exists(file.path(tempdir(), \"smoke.gif\")))
      })

      ## NN chunk 423 : naive distance threshold (trivial base R)
      ok(\"NN naive distance threshold (max/2 + subset)\", {
        step.x <- 0.4; step.y <- 0.3
        msp <- data.frame(distance = runif(20, 0, 1))
        thr <- max(step.y, step.x) / 2
        msp[msp\$distance < thr, , drop = FALSE]
      })
      cat(\"=== gated-chunk verification done ===\n\")'
  "
