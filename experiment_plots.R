# experiment_plots.R --- dimension-agnostic boxplots of the sampler-comparison
# =============================================================================
# Diagnostics for the GaussNiche virtual-species experiment (paper §2.3,
# "Experimental setting"). Consumes results from either engine:
#   - virtualSpecies()      (2-D, virtualSpecies_fn.R)
#   - virtualSpecies_nd()   (k-D, virtualSpecies_nd_fn.R)
# Both store per-realisation diagnostics in result$samplers[[name]]$metrics with
# the SAME column conventions:
#   overlap         hypervolume intersection (presence vs pseudo-absence clouds)
#   prop_true_abs   fraction of pseudo-absences on true-absence (y=0) cells
#   rel_cov_<PC>    range coverage = PA PC-range / background PC-range, ONE per axis
# The number and names of the rel_cov_* columns are discovered at runtime, so the
# same plotting code produces the demanded boxplots whether the experiment ran in
# 2, 5, or k dimensions — the coverage rows simply scale with the axis count.
#
# Public API:
#   collect_experiment_metrics(x, species_labels)   -> tidy wide metrics df
#   experiment_metrics_long(metrics, ...)            -> tidy long df (one panel col)
#   plot_experiment_boxplots(x, ...)                 -> the §2.3 figure (ggplot)
#   plot_experiment_metric(x, metric, ...)           -> one metric in isolation
#   save_experiment_figure(plot, file, ...)          -> vector PDF / 600-dpi PNG
#
# Modules return ggplot objects; saving is the caller's job (repo convention).
# =============================================================================

suppressPackageStartupMessages(library(ggplot2))

# --- Sampler identity ---------------------------------------------------------
# Maps the pa_samplers list KEYS used in the drivers to the paper's notation, and
# to a colourblind-safe (Okabe-Ito) fill. 'mcmc' is the paper's "uniform+" method
# (Markov-chain sampler, §2.2). Override both for other sampler sets.
EXPERIMENT_SAMPLER_LABELS <- c(random = "RND", buffer = "buffer-out",
                               uniform = "USE", mcmc = "uniform+", nn = "NN")
EXPERIMENT_SAMPLER_COLORS <- c(random = "#E69F00", buffer = "#D55E00",
                               uniform = "#56B4E9", mcmc = "#009E73", nn = "#CC79A7")

# --- internal helpers ---------------------------------------------------------

# Coverage columns, ordered by their trailing PC index (PC1, PC2, ..., PC10).
.experiment_coverage_cols <- function(df) {
  cc <- grep("^rel_cov_", names(df), value = TRUE)
  if (!length(cc)) return(character(0))
  idx <- suppressWarnings(as.integer(gsub("[^0-9]", "", sub("^rel_cov_", "", cc))))
  cc[order(idx, na.last = TRUE)]
}

# Resolve sampler factor order + display labels from the present keys: paper order
# first (per EXPERIMENT_SAMPLER_LABELS), then any unknown samplers appended.
.experiment_sampler_levels <- function(sampler_keys, sampler_labels) {
  present <- intersect(names(sampler_labels), unique(sampler_keys))
  extra   <- setdiff(unique(sampler_keys), present)
  keys    <- c(present, extra)
  labs    <- ifelse(keys %in% names(sampler_labels),
                    unname(sampler_labels[keys]), keys)
  list(keys = keys, labels = labs)
}

# --- collect: results object(s) -> one tidy wide metrics data.frame -----------
# `x` may be:
#   (a) a data.frame of metrics (must carry a 'species' column)  -> passthrough
#   (b) a single result object (has $samplers)                   -> one species
#   (c) a named list of result objects (each with $samplers)     -> many species
#   (d) a named list of {metrics = list(<sampler> = df), ...}     -> the lightweight
#       summary shape saved by the drivers                       -> many species
# `species_labels` is an optional named character vector (key -> pretty label).
collect_experiment_metrics <- function(x, species_labels = NULL) {
  if (is.data.frame(x)) {
    if (is.null(x$species)) stop("metrics data.frame must have a 'species' column")
    if (is.null(x$species_label)) x$species_label <- as.character(x$species)
    return(x)
  }
  if (!is.list(x)) stop("`x` must be a results list or a metrics data.frame")
  if (!is.null(x$samplers)) x <- list(species = x)   # wrap a single result

  sp_names <- names(x)
  if (is.null(sp_names)) sp_names <- paste0("sp", seq_along(x))

  rows <- list()
  for (i in seq_along(x)) {
    nm  <- sp_names[i]; el <- x[[i]]
    lab <- if (!is.null(species_labels) && nm %in% names(species_labels))
      unname(species_labels[[nm]]) else nm
    sms <- if (!is.null(el$samplers)) lapply(el$samplers, `[[`, "metrics")
           else if (!is.null(el$metrics)) el$metrics
           else stop(sprintf("element '%s' has neither $samplers nor $metrics", nm))
    for (s in names(sms)) {
      m <- sms[[s]]
      if (is.null(m) || !nrow(m)) next
      rows[[length(rows) + 1L]] <- cbind(
        species = nm, species_label = lab, sampler = s, m,
        stringsAsFactors = FALSE)
    }
  }
  if (!length(rows)) stop("no metrics found in `x`")
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

# --- reshape: wide metrics -> long (one row per species/sampler/realisation/panel)
# `panel` is an ordered factor: overlap, then one "Range coverage (PCk)" per axis,
# then proportion of true absences. `ok_only` keeps only status == "ok" rows.
experiment_metrics_long <- function(metrics, ok_only = TRUE,
                                    sampler_labels = EXPERIMENT_SAMPLER_LABELS) {
  metrics <- collect_experiment_metrics(metrics)
  if (ok_only && "status" %in% names(metrics))
    metrics <- metrics[metrics$status == "ok", , drop = FALSE]
  if (!nrow(metrics)) stop("no successful (status == 'ok') realisations to plot")

  cov_cols <- .experiment_coverage_cols(metrics)
  pcs      <- sub("^rel_cov_", "", cov_cols)
  cov_lab  <- paste0("Range coverage (", pcs, ")")
  panel_levels <- c("Hypervolume overlap", cov_lab, "Proportion true absences")

  reali <- if ("realization" %in% names(metrics)) metrics$realization else NA_integer_
  mk <- function(value, panel) data.frame(
    species       = metrics$species,
    species_label = metrics$species_label,
    sampler       = metrics$sampler,
    realization   = reali,
    panel         = panel,
    value         = value,
    stringsAsFactors = FALSE)

  parts <- list()
  if ("overlap" %in% names(metrics))
    parts[[length(parts) + 1L]] <- mk(metrics$overlap, "Hypervolume overlap")
  for (j in seq_along(cov_cols))
    parts[[length(parts) + 1L]] <- mk(metrics[[cov_cols[j]]], cov_lab[j])
  if ("prop_true_abs" %in% names(metrics))
    parts[[length(parts) + 1L]] <- mk(metrics$prop_true_abs, "Proportion true absences")
  long <- do.call(rbind, parts)

  long$panel <- factor(long$panel, levels = panel_levels)
  lv <- .experiment_sampler_levels(long$sampler, sampler_labels)
  long$sampler <- factor(long$sampler, levels = lv$keys, labels = lv$labels)
  long$species_label <- factor(long$species_label,
                               levels = unique(metrics$species_label))
  long
}

# --- house theme for these figures (portable: no Lato/showtext dependency) ----
experiment_theme <- function(base_size = 9, base_family = "") {
  theme_bw(base_size = base_size, base_family = base_family) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      strip.background = element_rect(fill = "grey92", colour = NA),
      strip.text = element_text(size = base_size - 1),
      strip.text.y.left = element_text(angle = 90),
      legend.position = "bottom",
      legend.title = element_blank(),
      panel.spacing = unit(0.55, "lines"),
      plot.title = element_text(face = "bold", size = base_size + 1))
}

# --- the demanded figure ------------------------------------------------------
# Boxplots of every metric (rows) for every species (columns), one box per
# sampler, the box summarising the empirical distribution across the R Bernoulli
# realisations. Dimension-agnostic: coverage expands to one row per PC axis.
#   x               results list / metrics df (see collect_experiment_metrics)
#   jitter          overlay realisation points (TRUE for small R)
#   returns a ggplot (size it with experiment_panel_count(); see save helper)
plot_experiment_boxplots <- function(x,
                                     sampler_palette = EXPERIMENT_SAMPLER_COLORS,
                                     sampler_labels  = EXPERIMENT_SAMPLER_LABELS,
                                     base_size = 9, base_family = "",
                                     jitter = TRUE, title = NULL) {
  m    <- collect_experiment_metrics(x)
  long <- experiment_metrics_long(m, sampler_labels = sampler_labels)

  # Palette must key off the SAME samplers that survive the ok-filter inside
  # experiment_metrics_long(), so colours match the plotted factor levels.
  ok_keys <- if ("status" %in% names(m)) m$sampler[m$status == "ok"] else m$sampler
  lv   <- .experiment_sampler_levels(ok_keys, sampler_labels)
  cols <- ifelse(lv$keys %in% names(sampler_palette),
                 unname(sampler_palette[lv$keys]), "grey50")
  names(cols) <- lv$labels

  p <- ggplot(long, aes(.data$sampler, .data$value, fill = .data$sampler)) +
    geom_boxplot(width = 0.62, alpha = 0.85, linewidth = 0.3,
                 outlier.size = 0.5, outlier.alpha = 0.5)
  if (jitter)
    p <- p + geom_jitter(width = 0.12, height = 0, size = 0.5, alpha = 0.25,
                         show.legend = FALSE)
  p +
    facet_grid(panel ~ species_label, scales = "free_y", switch = "y") +
    scale_fill_manual(values = cols, name = NULL) +
    labs(x = NULL, y = NULL, title = title) +
    experiment_theme(base_size, base_family) +
    theme(strip.placement = "outside")
}

# Number of facet rows the combined figure will have (overlap + k coverage axes
# + true-absences) — use to size the output canvas.
experiment_panel_count <- function(x) {
  m <- collect_experiment_metrics(x)
  ("overlap" %in% names(m)) + length(.experiment_coverage_cols(m)) +
    ("prop_true_abs" %in% names(m))
}

# --- one metric in isolation (flexible per-panel layout) ----------------------
# metric: "overlap" | "prop_true_abs" | "coverage" (all axes, faceted by axis).
plot_experiment_metric <- function(x, metric = c("overlap", "prop_true_abs", "coverage"),
                                    sampler_palette = EXPERIMENT_SAMPLER_COLORS,
                                    sampler_labels  = EXPERIMENT_SAMPLER_LABELS,
                                    base_size = 9, base_family = "",
                                    jitter = TRUE, title = NULL) {
  metric <- match.arg(metric)
  m    <- collect_experiment_metrics(x)
  long <- experiment_metrics_long(m, sampler_labels = sampler_labels)
  sel  <- switch(metric,
    overlap       = long[long$panel == "Hypervolume overlap", , drop = FALSE],
    prop_true_abs = long[long$panel == "Proportion true absences", , drop = FALSE],
    coverage      = long[grepl("^Range coverage", as.character(long$panel)), , drop = FALSE])
  sel$panel <- droplevels(sel$panel)

  ok_keys <- if ("status" %in% names(m)) m$sampler[m$status == "ok"] else m$sampler
  lv   <- .experiment_sampler_levels(ok_keys, sampler_labels)
  cols <- ifelse(lv$keys %in% names(sampler_palette),
                 unname(sampler_palette[lv$keys]), "grey50")
  names(cols) <- lv$labels

  p <- ggplot(sel, aes(.data$sampler, .data$value, fill = .data$sampler)) +
    geom_boxplot(width = 0.62, alpha = 0.85, linewidth = 0.3,
                 outlier.size = 0.5, outlier.alpha = 0.5)
  if (jitter)
    p <- p + geom_jitter(width = 0.12, height = 0, size = 0.5, alpha = 0.25,
                         show.legend = FALSE)
  facets <- if (metric == "coverage") facet_grid(panel ~ species_label, scales = "free_y",
                                                  switch = "y")
            else facet_wrap(~ species_label, nrow = 1)
  p + facets +
    scale_fill_manual(values = cols, name = NULL) +
    labs(x = NULL,
         y = switch(metric, overlap = "Hypervolume overlap",
                    prop_true_abs = "Proportion true absences",
                    coverage = "Range coverage"),
         title = title) +
    experiment_theme(base_size, base_family) +
    theme(strip.placement = "outside")
}

# --- saver: vector PDF (cairo_pdf) by default, 600-dpi PNG for .png -----------
save_experiment_figure <- function(plot, file, width = 7.0, height = NULL,
                                    dpi = 600, base_height_per_panel = 1.7,
                                    panel_count = NULL) {
  if (is.null(height)) {
    n <- if (!is.null(panel_count)) panel_count else 4L
    height <- max(4.0, base_height_per_panel * n)
  }
  if (grepl("\\.pdf$", file, ignore.case = TRUE)) {
    ggsave(file, plot, device = grDevices::cairo_pdf,
           width = width, height = height, units = "in")
  } else {
    ggsave(file, plot, device = "png",
           width = width, height = height, units = "in", dpi = dpi)
  }
  invisible(file)
}

# =============================================================================
# Report-style per-species figures (ported from 6_compare_5d.R) — SAME layout and
# colours as the GaussNiche report. Used for the appendix PC-matrix figures and
# the final geographic species x sampler grid. Dimension-agnostic.
# =============================================================================

# Report sampler colours (steelblue / purple as in the report; buffer-out added).
EXPERIMENT_REPORT_COLORS <- c(random = "steelblue", mcmc = "purple",
                              buffer = "#D55E00", uniform = "darkorange")

# plot_pc_matrix(): k x k scatterplot matrix of the niche in PC space, for one
# species. Diagonal = per-PC marginal density (background grey, presences red,
# each sampler's PA density); lower triangle = samplers[1] pseudo-absences, upper
# triangle = samplers[2]; every cell = background KDE (grey) + presences (gold) +
# niche suitability contours (firebrick, 0.1/0.5/0.95). `res` is one
# virtualSpecies()/virtualSpecies_nd() result. Returns a patchwork ggplot.
plot_pc_matrix <- function(res, samplers = c("random", "mcmc"),
                           sampler_palette = EXPERIMENT_REPORT_COLORS,
                           title = NULL, subtitle = NULL, max_pts = 1200L, ng = 70L) {
  for (pk in c("mvtnorm", "MASS", "patchwork"))
    if (!requireNamespace(pk, quietly = TRUE)) stop("plot_pc_matrix requires '", pk, "'.")
  pcs <- grep("^PC[0-9]+$", names(res$background), value = TRUE)
  pcs <- pcs[order(as.integer(sub("PC", "", pcs)))]
  k   <- length(pcs); bg <- res$background
  s1  <- samplers[1]; s2 <- samplers[min(2L, length(samplers))]
  col_lower <- unname(sampler_palette[s1]); col_upper <- unname(sampler_palette[s2])
  sub_ <- function(d) if (!is.null(d) && nrow(d) > max_pts) d[sample(nrow(d), max_pts), ] else d
  pres     <- sub_(bg[bg$pa == 1L, , drop = FALSE])
  pa_lower <- sub_(res$samplers[[s1]]$pseudo_ref)
  pa_upper <- if (length(samplers) >= 2L) sub_(res$samplers[[s2]]$pseudo_ref) else NULL
  has_pa   <- !is.null(pa_lower) || !is.null(pa_upper)
  mu   <- res$niche$mu; Sigma <- res$niche$Sigma
  peak <- mvtnorm::dmvnorm(matrix(mu, nrow = 1L), mean = mu, sigma = Sigma)
  rng  <- setNames(lapply(pcs, function(p) range(bg[[p]])), pcs)

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
      scale_fill_grey(start = .92, end = .45, guide = "none") +
      geom_point(data = pres, aes(.data[[xc]], .data[[yc]]), colour = "#FFC400",
                 size = .4, alpha = .6)
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
  cells <- vector("list", k * k)
  for (i in seq_len(k)) for (j in seq_len(k))
    cells[[(i - 1L) * k + j]] <-
      if (i == j) diagonal(i)
      else if (i > j) scatter(i, j, pa_lower, col_lower)
      else if (has_pa) scatter(i, j, pa_upper, col_upper)
      else patchwork::plot_spacer()
  if (is.null(subtitle))
    subtitle <- sprintf("lower = %s (%s)  |  upper = %s (%s)  |  presences gold  |  niche red 0.1/0.5/0.95",
                        s1, col_lower, s2, col_upper)
  patchwork::wrap_plots(cells, nrow = k, ncol = k) +
    patchwork::plot_annotation(title = title, subtitle = subtitle,
                               theme = theme(plot.title = element_text(face = "bold")))
}

# plot_geo_sampler_grid(): base-R species x sampler grid of the suitability map
# (YlOrRd) with each sampler's pseudo-absences overplotted (coloured per sampler).
# Draws to the CURRENT graphics device (wrap in cairo_pdf/png in the caller).
# Ported from 6_compare_5d.R's "pseudo-absences on the suitability map" page.
plot_geo_sampler_grid <- function(res_list, samplers = c("random", "buffer", "mcmc"),
                                  sampler_palette = EXPERIMENT_REPORT_COLORS,
                                  species_labels = NULL, crs = "EPSG:4326",
                                  main = "Pseudo-absences on the suitability map (species x sampler)") {
  if (!requireNamespace("terra", quietly = TRUE)) stop("plot_geo_sampler_grid requires 'terra'.")
  sp_names <- names(res_list); n_sp <- length(sp_names); n_s <- length(samplers)
  lab <- function(nm) if (!is.null(species_labels) && nm %in% names(species_labels))
    unname(species_labels[[nm]]) else nm
  suit_r <- function(r) terra::rast(r$background[, c("x", "y", "suit")], type = "xyz", crs = crs)
  op <- graphics::par(mfrow = c(n_sp, n_s), mar = c(2, 2, 3, 1), oma = c(0, 0, 2.5, 0))
  on.exit(graphics::par(op), add = TRUE)
  for (nm in sp_names) {
    sr <- suit_r(res_list[[nm]])
    for (s in samplers) {
      terra::plot(sr, col = grDevices::hcl.colors(100, "YlOrRd", rev = TRUE),
                  main = sprintf("%s · %s", lab(nm), s))
      pts <- res_list[[nm]]$samplers[[s]]$pseudo_ref[, c("x", "y"), drop = FALSE]
      if (!is.null(pts) && nrow(pts))
        graphics::points(pts$x, pts$y, col = unname(sampler_palette[s]), cex = .35, pch = 16)
    }
  }
  graphics::mtext(main, outer = TRUE, line = .5, font = 2)
  invisible(NULL)
}
