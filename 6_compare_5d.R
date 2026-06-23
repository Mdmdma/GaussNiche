# 6_compare_5d.R  —  shareable, self-contained report (5d-niche branch)
# =============================================================================
# Renders report_<mode>.pdf: a framed, self-contained document for collaborators.
# Structure:
#   FRONT MATTER  title/framing -> how-to-read -> contents -> run configuration
#   CROSS-SPECIES sampler comparison (overlap / coverage / true-absence)
#   PER SPECIES   divider (params + plain-language description) ->
#                 pairwise PC matrices (environment+niche, and pseudo-absences
#                 per sampler) -> niche+PA projection -> response/logit curves ->
#                 per-sampler sampling-bias + within-species boxplots
#   GEOGRAPHIC    suitability / presences / pseudo-absences back-projected to maps
# Every figure carries a one-line explainer. Report-only: reads the cached .rds.
#
# Usage:  Rscript 6_compare_5d.R [mode]      (mode = "full" default, or "smoke")
# =============================================================================

suppressPackageStartupMessages({
  .libPaths(c("~/R/rocker-rstudio/4.5", .libPaths()))
  library(terra); library(ggplot2); library(patchwork); library(grid)
})
source("virtualSpecies_nd_fn.R")   # summarize_parameters_nd()

mode    <- if (length(commandArgs(TRUE)) >= 1) commandArgs(TRUE)[1] else "full"
user    <- Sys.getenv("USER"); scratch <- file.path("/cluster/scratch", user)
env_dir <- file.path(scratch, "GaussNiche", "env5d")
out_dir <- file.path(scratch, "GaussNiche", "results5d")

files <- list.files(out_dir, pattern = paste0("_", mode, "\\.rds$"), full.names = TRUE)
files <- files[!grepl("^summary_", basename(files))]
stopifnot("no result files found; run 5_highdim_species.R first" = length(files) > 0)
sp_names  <- sub(paste0("_", mode, ".rds$"), "", basename(files))
res_list  <- setNames(lapply(files, readRDS), sp_names)
samplers  <- names(res_list[[1]]$samplers)
pa_cols   <- c(random = "steelblue", mcmc = "purple", uniform = "darkorange")
var_exp   <- tryCatch({ pm <- readRDS(file.path(env_dir, "final_rastPCA.rds"))
                        pm$sdev^2 / sum(pm$sdev^2) }, error = function(e) NULL)

# plain-language species descriptions (2x2 generalist/specialist x common/rare)
sp_meta <- list(
  sp1_generalist_common = list(label = "Generalist · common",
    desc = "Broad niche (sigma = 1.0 x sd on every axis) centred on the MODE of the background density. The easy case: many suitable cells, high prevalence."),
  sp2_specialist_common = list(label = "Specialist · common",
    desc = "Narrow niche (sigma = 0.5 x sd) at the density mode: few but well-defined presences embedded in a common environment."),
  sp3_generalist_rare   = list(label = "Generalist · rare",
    desc = "Broad niche at a peripheral optimum (low-density E-space). The niche spills back toward denser regions, giving an asymmetric suitability gradient."),
  sp4_specialist_rare   = list(label = "Specialist · rare",
    desc = "Narrow niche at the periphery: lowest prevalence and smallest suitable volume — the hardest case for pseudo-absence sampling."))
sp_label <- function(nm) if (!is.null(sp_meta[[nm]])) sp_meta[[nm]]$label else nm

# --- combined metrics + console summary --------------------------------------
comb <- do.call(rbind, lapply(sp_names, function(nm)
  do.call(rbind, lapply(samplers, function(s) {
    m <- res_list[[nm]]$samplers[[s]]$metrics
    if (is.null(m) || !nrow(m)) return(NULL)
    cbind(species = nm, sampler = s, m)
  }))))
ok <- comb[comb$status == "ok", , drop = FALSE]
relcov_cols <- grep("^rel_cov_", names(comb), value = TRUE)
write.csv(comb, file.path(out_dir, paste0("metrics_combined_", mode, ".csv")), row.names = FALSE)

cat("\n==================== SUMMARY ====================\n")
for (nm in sp_names) {
  cat(sprintf("\n%s  (prevalence %.4f)\n", nm, mean(res_list[[nm]]$background$pa)))
  for (s in samplers) {
    o <- ok[ok$species == nm & ok$sampler == s, , drop = FALSE]
    cat(sprintf("  %-7s ok=%2d  overlap(med)=%-8s\n", s, nrow(o),
                if (nrow(o)) sprintf("%.3f", median(o$overlap, na.rm = TRUE)) else "NA"))
  }
}
cat("\n", paste(capture.output(print(table(species = comb$species, sampler = comb$sampler,
                                            status = comb$status))), collapse = "\n"), "\n")

# =============================================================================
# helpers
# =============================================================================
# page counter + bookmark targets (for the clickable PDF outline added after render)
PAGE <- 0L; BM <- list()
mark <- function(title) BM[[length(BM) + 1L]] <<- list(title = title, page = PAGE)
pg <- function(expr) tryCatch({ print(expr); PAGE <<- PAGE + 1L },
        error = function(e) message("page failed: ", conditionMessage(e)))
cov_long <- function(df) do.call(rbind, lapply(relcov_cols, function(cc)
  data.frame(species = df$species, sampler = df$sampler,
             axis = sub("rel_cov_", "", cc), value = df[[cc]])))
clean <- function(p, title = NULL) {        # strip subtitle / shrink title for composites
  p <- p + labs(subtitle = NULL) + theme(plot.title = element_text(size = 9))
  if (!is.null(title)) p <- p + labs(title = title)
  p
}
# unified per-axis sampling-bias page: background vs presences vs BOTH samplers'
# pseudo-absences overlaid per axis, sharing one vertical legend (more room for plots).
unified_bias <- function(r, title, samplers, pa_cols, max_pts = 1500L) {
  pcs <- grep("^PC[0-9]+$", names(r$background), value = TRUE)
  pcs <- pcs[order(as.integer(sub("PC", "", pcs)))]
  bg  <- r$background
  sub_ <- function(d) if (nrow(d) > max_pts) d[sample(nrow(d), max_pts), ] else d
  pres <- sub_(bg[bg$pa == 1L, , drop = FALSE])
  pas  <- setNames(lapply(samplers, function(s) sub_(r$samplers[[s]]$pseudo_ref)), samplers)
  cols <- c(Background = "grey55", Presences = "firebrick")
  for (s in samplers) cols[s] <- pa_cols[[s]]
  levs <- c("Background", "Presences", samplers)
  plots <- lapply(pcs, function(pc) {
    dl <- rbind(data.frame(val = bg[[pc]],   group = "Background"),
                data.frame(val = pres[[pc]], group = "Presences"),
                do.call(rbind, lapply(samplers, function(s)
                  data.frame(val = pas[[s]][[pc]], group = s))))
    dl$group <- factor(dl$group, levels = levs)
    ggplot(dl, aes(val, colour = group, fill = group)) +
      geom_density(alpha = .15, linewidth = .6) +
      scale_colour_manual(values = cols, name = NULL, drop = FALSE) +
      scale_fill_manual(values = cols, name = NULL, drop = FALSE) +
      labs(title = pc, x = pc, y = "Density") +
      theme_classic(base_size = 10) +
      theme(legend.position = "right", legend.direction = "vertical")
  })
  wrap_plots(plots, ncol = 3) + plot_layout(guides = "collect") +
    plot_annotation(title = title,
      caption = "Per-axis density: background vs presences vs each sampler's pseudo-absences. Closer overlap of a sampler with the presence curve = worse separation on that axis.",
      theme = theme(plot.title = element_text(face = "bold")))
}
# unified within-species metrics page: both samplers side by side.
unified_box <- function(title, oks) {
  p_o <- ggplot(oks, aes(sampler, overlap, fill = sampler)) +
    geom_boxplot(alpha = .7, width = .5) + geom_jitter(width = .1, alpha = .4, size = .8) +
    labs(title = "Hypervolume overlap (lower better)", x = NULL, y = "overlap") +
    theme_classic(base_size = 10) + theme(legend.position = "none")
  p_t <- ggplot(oks, aes(sampler, prop_true_abs, fill = sampler)) +
    geom_boxplot(alpha = .7, width = .5) + ylim(0, 1) +
    labs(title = "True-absence fraction (higher better)", x = NULL, y = "proportion") +
    theme_classic(base_size = 10) + theme(legend.position = "none")
  p_c <- ggplot(cov_long(oks), aes(sampler, value, fill = sampler)) +
    geom_boxplot(alpha = .7, width = .6) + ylim(0, 1) + facet_wrap(~axis, nrow = 1) +
    labs(title = "Per-axis coverage of pseudo-absences (higher better)", x = NULL, y = "rel. coverage") +
    theme_classic(base_size = 9) + theme(legend.position = "none")
  ((p_o | p_t) / p_c) + plot_layout(heights = c(1, 1.05)) +
    plot_annotation(title = title,
      caption = "Both samplers side by side across the independent Bernoulli realisations.",
      theme = theme(plot.title = element_text(face = "bold")))
}
# all-species response curves on one figure: per axis, every species' marginal
# (suitability or logit) overlaid -> compare niche position & breadth at a glance.
all_species_response <- function(res_list, sp_names, sp_label, logit = FALSE) {
  pcs <- grep("^PC[0-9]+$", names(res_list[[1]]$response_curves$data), value = TRUE)
  pcs <- pcs[order(as.integer(sub("PC", "", pcs)))]
  plots <- lapply(pcs, function(pc) {
    dl <- do.call(rbind, lapply(sp_names, function(nm) {
      d <- res_list[[nm]]$response_curves$data[[pc]]
      yv <- if (logit) log(d$suit / (1 - d$suit)) else d$suit
      data.frame(x = d$x, y = yv, species = sp_label(nm))
    }))
    dl$species <- factor(dl$species, levels = vapply(sp_names, sp_label, character(1)))
    if (logit) dl <- dl[is.finite(dl$y), ]
    ggplot(dl, aes(x, y, colour = species)) + geom_line(linewidth = .7) +
      scale_colour_brewer(palette = "Dark2", name = NULL, drop = FALSE) +
      labs(title = pc, x = pc, y = if (logit) "logit(suit)" else "suitability") +
      theme_classic(base_size = 10) + theme(legend.position = "right")
  })
  ttl <- if (logit) "Logit response curves — all species" else "Suitability response curves — all species"
  wrap_plots(plots, ncol = 3) + plot_layout(guides = "collect") +
    plot_annotation(title = ttl,
      caption = "Marginal niche shape along each PC for every species (other axes at that species' optimum). Compare where each peak sits (position) and how wide it is (breadth).",
      theme = theme(plot.title = element_text(face = "bold")))
}
# full-page text block (title + wrapped paragraphs). Drawn straight to the device.
text_page <- function(title, body, subtitle = NULL, wrap = 108, body_size = 10.5, bookmark = NULL) {
  grid.newpage()
  y <- unit(1, "npc") - unit(9, "mm")
  grid.text(title, x = unit(8, "mm"), y = y, just = c("left", "top"),
            gp = gpar(fontface = "bold", fontsize = 17))
  y <- y - unit(10, "mm")
  if (!is.null(subtitle)) {
    grid.text(subtitle, x = unit(8, "mm"), y = y, just = c("left", "top"),
              gp = gpar(fontface = "italic", fontsize = 11, col = "grey35"))
    y <- y - unit(8, "mm")
  }
  y <- y - unit(2, "mm")
  for (para in body) {
    for (ln in strwrap(para, width = wrap)) {
      grid.text(ln, x = unit(8, "mm"), y = y, just = c("left", "top"),
                gp = gpar(fontsize = body_size))
      y <- y - unit(5.2, "mm")
    }
    y <- y - unit(2.6, "mm")
  }
  PAGE <<- PAGE + 1L
  if (!is.null(bookmark)) mark(bookmark)
}
# section-divider page (large centred title + optional sub-lines)
divider_page <- function(big, lines = NULL, bookmark = NULL) {
  grid.newpage()
  grid.rect(gp = gpar(fill = "grey94", col = NA))
  grid.text(big, x = 0.5, y = 0.62, just = "centre", gp = gpar(fontface = "bold", fontsize = 24))
  if (!is.null(lines)) {
    yy <- 0.50
    for (ln in lines) { grid.text(ln, x = 0.5, y = yy, just = "centre",
                                  gp = gpar(fontsize = 12, col = "grey25")); yy <- yy - 0.045 }
  }
  PAGE <<- PAGE + 1L
  if (!is.null(bookmark)) mark(bookmark)
}

# Pairwise PC scatter-matrix for ONE species. Lower triangle = PCj (x) vs PCi (y)
# panels with the background KDE (mako) + presences + niche suitability contours;
# diagonal = per-axis marginal density (bg vs presences). If `pa` is supplied,
# that sampler's pseudo-absences are overlaid (and presences are dimmed).
pc_matrix_figure <- function(r, title, subtitle = NULL,
                             pa_lower = NULL, col_lower = "steelblue",
                             pa_upper = NULL, col_upper = "purple",
                             max_pts = 1200L, ng = 70L) {
  pcs <- grep("^PC[0-9]+$", names(r$background), value = TRUE)
  pcs <- pcs[order(as.integer(sub("PC", "", pcs)))]
  k   <- length(pcs); bg <- r$background
  sub_ <- function(d) if (!is.null(d) && nrow(d) > max_pts) d[sample(nrow(d), max_pts), ] else d
  pres <- sub_(bg[bg$pa == 1L, , drop = FALSE])
  pa_lower <- sub_(pa_lower); pa_upper <- sub_(pa_upper)
  has_pa <- !is.null(pa_lower) || !is.null(pa_upper)
  mu   <- r$niche$mu; Sigma <- r$niche$Sigma
  peak <- mvtnorm::dmvnorm(matrix(mu, nrow = 1L), mean = mu, sigma = Sigma)
  rng  <- setNames(lapply(pcs, function(p) range(bg[[p]])), pcs)
  pres_alpha <- if (has_pa) 0.07 else 0.16

  # cell axes: ALWAYS lower-rank PC on x, higher-rank PC on y (so an upper-triangle
  # cell and its lower-triangle mirror share the same orientation -> comparable).
  scatter <- function(i, j, pa, pa_col) {
    xc <- pcs[min(i, j)]; yc <- pcs[max(i, j)]
    g  <- expand.grid(X = seq(rng[[xc]][1], rng[[xc]][2], length.out = ng),
                      Y = seq(rng[[yc]][1], rng[[yc]][2], length.out = ng))
    sl <- matrix(rep(mu, each = nrow(g)), nrow = nrow(g)); colnames(sl) <- pcs
    sl[, xc] <- g$X; sl[, yc] <- g$Y
    g$suit <- mvtnorm::dmvnorm(sl[, pcs], mean = mu, sigma = Sigma) / peak
    kd  <- MASS::kde2d(bg[[xc]], bg[[yc]], n = ng, lims = c(rng[[xc]], rng[[yc]]))
    kdf <- expand.grid(X = kd$x, Y = kd$y); kdf$d <- as.vector(kd$z); kdf$d <- kdf$d / max(kdf$d)
    p <- ggplot() +
      geom_contour_filled(data = kdf, aes(X, Y, z = d),
                          breaks = c(0, .05, .25, .5, .75, 1), alpha = .9) +
      # neutral grey background so the coloured points (incl. blue random PAs,
      # which were camouflaged on the teal mako fill) all stand out.
      scale_fill_grey(start = .92, end = .45, guide = "none") +
      # presences first (small, bright, no ring) as the reference layer ...
      geom_point(data = pres, aes(.data[[xc]], .data[[yc]]), colour = "#FFC400",
                 size = .4, alpha = .6)
    # ... then the pseudo-absences drawn OVER the presences.
    if (!is.null(pa))
      p <- p + geom_point(data = pa, aes(.data[[xc]], .data[[yc]]), colour = pa_col,
                          alpha = .55, size = .4)
    p + geom_contour(data = g, aes(X, Y, z = suit), colour = "firebrick",
                     linewidth = .4, breaks = c(0.1, 0.5, 0.95)) +
      labs(x = xc, y = yc) + theme_classic(base_size = 7) +
      theme(legend.position = "none", axis.title = element_text(size = 6),
            plot.margin = margin(1, 1, 1, 1))
  }
  diagonal <- function(i) {
    p <- pcs[i]
    g <- ggplot() +
      geom_density(data = bg,   aes(.data[[p]]), fill = "grey70", colour = "grey40", alpha = .5) +
      geom_density(data = pres, aes(.data[[p]]), fill = "firebrick", colour = "firebrick", alpha = .3)
    if (!is.null(pa_lower)) g <- g + geom_density(data = pa_lower, aes(.data[[p]]), colour = col_lower, fill = NA, linewidth = .5)
    if (!is.null(pa_upper)) g <- g + geom_density(data = pa_upper, aes(.data[[p]]), colour = col_upper, fill = NA, linewidth = .5)
    g + labs(title = p, x = NULL, y = NULL) + theme_classic(base_size = 7) +
      theme(legend.position = "none", plot.title = element_text(size = 8, face = "bold"),
            axis.text = element_blank(), axis.ticks = element_blank(),
            plot.margin = margin(1, 1, 1, 1))
  }
  # lower triangle -> pa_lower (one sampler); upper triangle -> pa_upper (the other)
  cells <- vector("list", k * k)
  for (i in seq_len(k)) for (j in seq_len(k))
    cells[[(i - 1L) * k + j]] <-
      if (i == j) diagonal(i)
      else if (i > j) scatter(i, j, pa_lower, col_lower)
      else if (has_pa) scatter(i, j, pa_upper, col_upper)
      else patchwork::plot_spacer()
  if (is.null(subtitle))
    subtitle <- "background KDE (grey) + presences (gold) + niche suitability 0.1/0.5/0.95 (red); diagonal = marginal density"
  wrap_plots(cells, nrow = k, ncol = k) +
    plot_annotation(title = title, subtitle = subtitle,
                    theme = theme(plot.title = element_text(face = "bold")))
}

# =============================================================================
# PDF REPORT
# =============================================================================
pdf(file.path(out_dir, paste0("report_", mode, ".pdf")), width = 13, height = 9)

# ---- FRONT MATTER -----------------------------------------------------------
text_page(
  "GaussNiche — pseudo-absence sampling in a 5-dimensional environment",
  subtitle = sprintf("Virtual-species experiment | %d species x %d samplers x %d realisations | run: %s",
                     length(sp_names), length(samplers),
                     max(comb$realization, na.rm = TRUE), mode),
  body = c(
    "WHAT THIS IS. GaussNiche simulates 'virtual species' as multivariate-Gaussian niches in an environmental space, then asks a methodological question: how well do different PSEUDO-ABSENCE sampling strategies separate a species' presences from its absences? Because the species are simulated, we know the ground truth exactly, so we can score each sampler fairly.",
    "THE PIPELINE. (1) Stack environmental rasters and reduce them to principal components (PCA). (2) Define a species as a Gaussian 'suitability' bump in PC space, centred on an optimum mu with breadth sigma. (3) Draw presences: each background cell becomes a presence with probability equal to its suitability (a Bernoulli trial), repeated over many realisations to capture stochastic variation. (4) Draw pseudo-absences with each sampler. (5) Score the sampler with hypervolume class-overlap and coverage diagnostics. (6) Back-project to maps.",
    "WHAT IS NEW HERE. The original framework worked in 2 dimensions (PC1, PC2). This is the 5-DIMENSIONAL extension: the environment is enriched (soil + climate + terrain + vegetation) so that 5 PCs are needed to explain >80% of the variance, and the whole pipeline - niche, presences, hypervolumes, samplers - runs in 5-D. Two samplers are compared: RANDOM and MCMC (see next page).",
    "THE TAKEAWAY TO LOOK FOR. A good pseudo-absence sampler should draw absences from environments UNLIKE the species' presences (low class overlap) while still spanning the whole available environment (high coverage). The MCMC sampler is designed to do this; random sampling is the naive baseline. The figures that follow quantify the difference, per species."),
  wrap = 116)
mark("Overview")

text_page(
  "How to read this report",
  body = c(
    "THE ENVIRONMENT (5 axes). The PC axes come from a 12-layer stack over Central & Western Europe at 10 arc-min. Roughly: PC1 = soil fertility/texture, PC2 = temperature continentality, PC3 = terrain ruggedness / forest cover, PC4 = precipitation seasonality, PC5 = topographic position. The first 5 PCs explain ~82% of the environmental variance.",
    "THE 4 SPECIES (a 2x2 design). Two niche BREADTHS x two niche POSITIONS:  Generalist (broad, sigma = 1.0 x sd) vs Specialist (narrow, sigma = 0.5 x sd);  Common (optimum on the densest environment) vs Rare (optimum at a sparse, peripheral environment). This spans easy-to-hard cases for the samplers.",
    "THE 2 SAMPLERS.  RANDOM = pseudo-absences drawn uniformly from all background cells (the naive baseline).  MCMC = a Markov chain (USE.MCMC) whose target is 'environmental density minus presence density', i.e. it spreads pseudo-absences uniformly across the environment while avoiding presence-like conditions. MCMC is the high-dimensional generalisation of the 2-D 'uniform-in-E-space' grid sampler (which does not scale beyond ~2 axes).",
    "THE METRICS.  PREVALENCE = fraction of background cells that are presences.  HYPERVOLUME OVERLAP = volume shared by the presence and pseudo-absence point-clouds in 5-D (LOWER = better class separation).  PER-AXIS COVERAGE = fraction of each PC's range spanned by the pseudo-absences (HIGHER = absences span the full gradient, near 1 is ideal).  PROP. TRUE-ABSENCE = fraction of pseudo-absences landing on genuine absence cells (HIGHER = sampler avoids true-presence cells).",
    "STRUCTURE & NAVIGATION. After this front matter and a cross-species comparison, each species gets its own SECTION (a grey divider page starts each one) so you can jump straight to a species - the PDF also has a clickable bookmark per section. Within a section: a pairwise PC matrix (niche + each sampler's pseudo-absences, one sampler per triangle), suitability and logit response curves, a per-axis sampling-bias comparison, and a side-by-side metrics page. A geographic section closes the report. Every figure has a one-line explainer at the top."),
  wrap = 116)

text_page(
  "Contents",
  body = c(
    "1.  Run configuration / replication record",
    "2.  Cross-species sampler comparison  (overlap, coverage, true-absence)",
    "3.  Niche shapes across species  (suitability & logit response curves, all 4)",
    paste0("4.  Per-species sections:"),
    paste0("       ", sprintf("4.%d  %s  (%s)", seq_along(sp_names),
                              vapply(sp_names, sp_label, character(1)), sp_names)),
    "5.  Geographic projection  (suitability, presences, pseudo-absences on maps)",
    "",
    "Each per-species section contains: a pairwise PC matrix (niche + pseudo-absences, one sampler per triangle); suitability & logit response curves; a unified per-axis sampling-bias comparison (both samplers); and a side-by-side metrics page (overlap, coverage, true-absence)."),
  wrap = 116)

pg(summarize_parameters_nd(res_list, pca_var_exp = var_exp,
                           extra = list(mode = mode, n_species = length(sp_names))))
mark("Run configuration")

# ---- CROSS-SPECIES COMPARISON -----------------------------------------------
divider_page("Cross-species sampler comparison",
             c("How random vs MCMC pseudo-absences compare, across all 4 species"),
             bookmark = "Cross-species comparison")
text_page(
  "Cross-species comparison — how to read the next pages",
  body = c(
    "Each of the following panels shows one diagnostic with the two samplers side by side, faceted by species. Each box summarises the spread across the independent Bernoulli realisations.",
    "OVERLAP (lower is better): the hypervolume shared by presences and pseudo-absences. A sampler that draws absences from presence-like environments inflates this and hurts model discrimination.",
    "PER-AXIS COVERAGE (higher is better, max 1): how much of each PC's gradient the pseudo-absences span. Random tracks background density and tends to under-cover the minor axes; MCMC spreads more uniformly.",
    "PROPORTION ON TRUE-ABSENCE CELLS (higher is better): fraction of pseudo-absences on cells that were genuine absences in the Bernoulli draw.",
    "WATCH FOR: MCMC usually shows lower overlap and higher coverage than random. The generalist-rare species is the instructive exception - a broad niche at the periphery fills so much of E-space that even a good sampler cannot avoid overlap."),
  wrap = 116)

if (nrow(ok) > 0) {
  p_ov <- ggplot(ok, aes(sampler, overlap, fill = sampler)) +
    geom_boxplot(alpha = .7, width = .5) + geom_jitter(width = .1, alpha = .4, size = .7) +
    facet_wrap(~species, nrow = 1, scales = "free_y") +
    labs(title = "Hypervolume overlap (presence vs pseudo-absence) — LOWER is better",
         x = NULL, y = "overlap volume") +
    theme_classic(base_size = 11) + theme(legend.position = "none")
  p_cv <- ggplot(cov_long(ok), aes(sampler, value, fill = sampler)) +
    geom_boxplot(alpha = .7, width = .6) + ylim(0, 1) + facet_grid(axis ~ species) +
    labs(title = "Per-axis coverage of pseudo-absences (PC1..PC5) — HIGHER is better",
         x = NULL, y = "rel. coverage") +
    theme_classic(base_size = 10) + theme(legend.position = "none")
  p_ta <- ggplot(ok, aes(sampler, prop_true_abs, fill = sampler)) +
    geom_boxplot(alpha = .7, width = .5) + ylim(0, 1) + facet_wrap(~species, nrow = 1) +
    labs(title = "Pseudo-absences on true-absence cells — HIGHER is better",
         x = NULL, y = "proportion") +
    theme_classic(base_size = 11) + theme(legend.position = "none")
  pg(p_ov); pg(p_ta); pg(p_cv)
}

# ---- NICHE SHAPES ACROSS SPECIES (all-species response curves up front) ------
text_page("Niche shapes across species (response curves)",
  body = c(
    "These two pages overlay ALL species' marginal response curves per PC axis, so you can compare niche POSITION (where each peak sits on an axis) and BREADTH (how wide) across the 2x2 generalist/specialist x common/rare design at a glance.",
    "Suitability curves (next page) show the niche bump directly; logit curves (the page after) show the non-linear transform a logistic SDM would fit. Each species also has these curves within its own section."),
  bookmark = "Niche shapes (all species)", wrap = 116)
asr <- all_species_response(res_list, sp_names, sp_label, logit = FALSE)
pg(asr)
tryCatch(ggsave(file.path(out_dir, "preview_response_all.png"), asr,
                width = 12, height = 7, dpi = 110), error = function(e) NULL)
pg(all_species_response(res_list, sp_names, sp_label, logit = TRUE))

# ---- PER-SPECIES SECTIONS ---------------------------------------------------
for (si in seq_along(sp_names)) {
  nm <- sp_names[si]; r <- res_list[[nm]]
  prm <- r$parameters
  prev <- mean(r$background$pa)
  ov <- sapply(samplers, function(s) { o <- ok[ok$species == nm & ok$sampler == s, ]
                if (nrow(o)) round(median(o$overlap, na.rm = TRUE), 1) else NA })
  divider_page(
    sprintf("Species %d / %d   —   %s", si, length(sp_names), sp_label(nm)),
    c(if (!is.null(sp_meta[[nm]])) sp_meta[[nm]]$desc else "",
      sprintf("mu = (%s)", paste(round(prm$mu, 2), collapse = ", ")),
      sprintf("sigma = (%s)", paste(round(prm$sigma, 2), collapse = ", ")),
      sprintf("prevalence = %.4f   |   median overlap:  %s",
              prev, paste(sprintf("%s=%s", samplers, ov), collapse = "   "))),
    bookmark = sprintf("Species %d — %s", si, sp_label(nm)))

  # pairwise PC matrix: niche + presences + pseudo-absences (sampler 1 LOWER, sampler 2 UPPER).
  # This single matrix subsumes the old environment-only matrix and the PC1xPC2 projection.
  s1 <- samplers[1]; s2 <- samplers[min(2L, length(samplers))]
  pa_mx <- pc_matrix_figure(
    r, paste0(sp_label(nm), " — niche & pseudo-absences across all PC pairs"),
    subtitle = sprintf("LOWER triangle = %s (%s)  |  UPPER triangle = %s (%s)  |  presences gold  |  niche red 0.1/0.5/0.95",
                       s1, pa_cols[[s1]], s2, pa_cols[[s2]]),
    pa_lower = r$samplers[[s1]]$pseudo_ref, col_lower = pa_cols[[s1]],
    pa_upper = r$samplers[[s2]]$pseudo_ref, col_upper = pa_cols[[s2]])
  pg(pa_mx)
  if (si == 1L) tryCatch(ggsave(file.path(out_dir, "preview_pa_matrix.png"), pa_mx,
                                width = 11, height = 10, dpi = 110), error = function(e) NULL)
  # response curves (suitability + logit)
  pg(wrap_plots(lapply(r$response_curves$plots, clean), ncol = 3) +
       plot_annotation(title = paste0(sp_label(nm), " — suitability response curves"),
                       caption = "Marginal suitability along each PC (other axes fixed at the optimum). Rug = background."))
  pg(wrap_plots(lapply(r$response_curves$logit_plots, clean), ncol = 3) +
       plot_annotation(title = paste0(sp_label(nm), " — logit response curves"),
                       caption = "logit(suitability) per axis; the non-linear transform a logistic SDM would see."))
  # unified sampling-bias: both samplers overlaid per axis (one vertical legend)
  ub <- unified_bias(r, paste0(sp_label(nm), " — sampling bias per axis (both samplers)"),
                     samplers, pa_cols)
  pg(ub)
  if (si == 1L) tryCatch(ggsave(file.path(out_dir, "preview_bias.png"), ub,
                                width = 12, height = 7, dpi = 110), error = function(e) NULL)
  # unified within-species metrics: both samplers side by side
  pg(unified_box(paste0(sp_label(nm), " — metrics across realisations (both samplers)"),
                 ok[ok$species == nm, , drop = FALSE]))
}

# ---- GEOGRAPHIC -------------------------------------------------------------
divider_page("Geographic projection",
             c("Suitability, presences and pseudo-absences back-projected to the map"),
             bookmark = "Geographic projection")
text_page(
  "Geographic maps — how to read",
  body = c(
    "The niche lives in PC space, but every cell maps back to a location. These maps show, per species: the suitability surface (yellow->red), the reference presences (red), and where each sampler placed its pseudo-absences (points).",
    "The PC-space patterns above translate here: a sampler with good coverage scatters pseudo-absences broadly across the map; one that tracks the species clusters them near the presences."),
  wrap = 116)
n_sp <- length(sp_names); n_s <- length(samplers); geo_crs <- "EPSG:4326"
suit_r <- function(r) terra::rast(r$background[, c("x", "y", "suit")], type = "xyz", crs = geo_crs)
pa_r   <- function(r) terra::rast(data.frame(x = r$background$x, y = r$background$y,
                                  pa = as.numeric(r$background$pa)), type = "xyz", crs = geo_crs)
tryCatch({
  op <- par(mfrow = c(1, n_sp), mar = c(2, 2, 3, 1), oma = c(0, 0, 2.5, 0))
  for (nm in sp_names) plot(suit_r(res_list[[nm]]), col = hcl.colors(100, "YlOrRd", rev = TRUE),
                            main = sp_label(nm))
  mtext("Suitability (geographic) — per species", outer = TRUE, line = .5, font = 2); par(op); PAGE <<- PAGE + 1L
}, error = function(e) message("geo suit failed: ", conditionMessage(e)))
tryCatch({
  op <- par(mfrow = c(1, n_sp), mar = c(2, 2, 3, 1), oma = c(0, 0, 2.5, 0))
  for (nm in sp_names) plot(pa_r(res_list[[nm]]), col = c("grey85", "firebrick"), main = sp_label(nm))
  mtext("Presences (geographic) — per species", outer = TRUE, line = .5, font = 2); par(op); PAGE <<- PAGE + 1L
}, error = function(e) message("geo pa failed: ", conditionMessage(e)))
tryCatch({
  op <- par(mfrow = c(n_sp, n_s), mar = c(2, 2, 3, 1), oma = c(0, 0, 2.5, 0))
  for (nm in sp_names) { sr <- suit_r(res_list[[nm]])
    for (s in samplers) {
      plot(sr, col = hcl.colors(100, "YlOrRd", rev = TRUE), main = sprintf("%s · %s", sp_label(nm), s))
      pts <- res_list[[nm]]$samplers[[s]]$pseudo_ref[, c("x", "y"), drop = FALSE]
      if (!is.null(pts) && nrow(pts)) points(pts$x, pts$y, col = pa_cols[[s]], cex = .35, pch = 16)
    } }
  mtext("Pseudo-absences on the suitability map (species x sampler)", outer = TRUE, line = .5, font = 2); par(op); PAGE <<- PAGE + 1L
}, error = function(e) message("geo PA grid failed: ", conditionMessage(e)))

dev.off()
cat(sprintf("\nRendered %d pages; %d bookmark targets.\n", PAGE, length(BM)))

# Write the bookmark targets (page<TAB>title) for the outline injection step in
# sbatch/submit_compare_5d.sh (python/pypdf), which runs after this script.
bm_tsv <- file.path(out_dir, paste0("report_", mode, "_bookmarks.tsv"))
con <- file(bm_tsv, "wt")
for (b in BM) cat(sprintf("%d\t%s\n", b$page, gsub("[\t\n]", " ", b$title)), file = con)
close(con)

cat("\nWrote framed report_", mode, ".pdf (+ bookmarks.tsv) to ", out_dir, "\n", sep = "")
