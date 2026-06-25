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
