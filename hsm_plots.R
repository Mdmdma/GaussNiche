# hsm_plots.R --- downstream HSM comparison figures + significance tests
# =============================================================================
# Consumes the tidy HSM stats produced by hsm_eval.R / run_5d_hsm.R (one row per
# species x sampler x realisation x predictor_set x algorithm) and renders the
# USE-paper downstream comparison (Da Re et al. 2023, Analysis/2_ViolinPlots.R +
# 3_DunnTest.R), ported to GaussNiche's 4 controlled species and 3 samplers:
#
#   plot_hsm_violin()  per-species violin/box of each predictive metric by
#                      sampler, faceted by algorithm x metric (USE's ModelType x
#                      metric layout).
#   hsm_dunn_test()    Dunn pairwise test (uniform+ vs random, uniform+ vs
#                      buffer-out), holm-adjusted, per (species, predictor_set,
#                      algorithm, metric). WITHIN-species contrast (the 50 points
#                      are realisations of ONE species, not 50 species) -- no
#                      single omnibus p.
#   summarize_hsm()    per-cell median table for the appendix.
#   save_hsm_figure()  vector PDF / PNG saver.
#
# Reuses the sampler colour/label maps from experiment_plots.R (source it first).
# Modules return objects; saving is the caller's job (repo convention).
# Tolerant of ragged data: buffer-out yields zero pseudo-absences for tightly
# clustered presences, so some (species, realisation) cells are simply absent.
# =============================================================================

suppressPackageStartupMessages(library(ggplot2))

`%||%` <- function(a, b) if (is.null(a)) b else a

# Pretty labels / display order for the metrics we plot. Truth metrics are the
# headline (predicted vs the KNOWN suitability); the rest are internal validation.
HSM_METRIC_LABELS <- c(
  cor_truth = "COR (truth)", rmse_truth = "RMSE (truth)",
  auc = "AUC", tss = "TSS", boyce = "CBI",
  sensitivity = "Sensitivity", specificity = "Specificity",
  cor = "COR", kappa = "Kappa", deviance = "Deviance", rmse = "RMSE", r2 = "R2")

HSM_ALGO_LABELS <- c(glm = "GLM", gam = "GAM", rf = "RF", brt = "BRT", maxent = "Maxent")

# Default metric panel for the violin figures (USE's set + the truth pair).
HSM_DEFAULT_METRICS <- c("cor_truth", "rmse_truth", "auc", "tss", "boyce")

# Fallbacks if experiment_plots.R was not sourced (keep colours/labels in sync).
if (!exists("EXPERIMENT_SAMPLER_LABELS"))
  EXPERIMENT_SAMPLER_LABELS <- c(random = "RND", buffer = "buffer-out",
                                 uniform = "USE", mcmc = "uniform+", nn = "NN")
if (!exists("EXPERIMENT_SAMPLER_COLORS"))
  EXPERIMENT_SAMPLER_COLORS <- c(random = "#E69F00", buffer = "#D55E00",
                                 uniform = "#56B4E9", mcmc = "#009E73", nn = "#CC79A7")

# --- collect / reshape --------------------------------------------------------

#' Combine per-species virtualSpecies_nd() results' $hsm slots into one df, or
#' pass through a data.frame. `species_labels` optionally maps keys -> labels.
collect_hsm_metrics <- function(x, species_labels = NULL) {
  hsm <- if (is.data.frame(x)) x
         else do.call(rbind, lapply(x, function(el)
           if (!is.null(el$hsm)) el$hsm else NULL))
  if (is.null(hsm) || !nrow(hsm)) stop("no HSM rows found in `x`")
  if (is.null(hsm$species_label))
    hsm$species_label <- if (!is.null(species_labels))
      ifelse(hsm$species %in% names(species_labels),
             unname(species_labels[hsm$species]), hsm$species)
    else as.character(hsm$species)
  hsm
}

#' Long form for plotting: one row per (…, metric, value), ok rows + finite only.
hsm_metrics_long <- function(hsm, metrics = HSM_DEFAULT_METRICS,
                             ok_only = TRUE, species_labels = NULL) {
  hsm <- collect_hsm_metrics(hsm, species_labels)
  if (ok_only && "status" %in% names(hsm))
    hsm <- hsm[hsm$status == "ok", , drop = FALSE]
  metrics <- intersect(metrics, names(hsm))
  parts <- lapply(metrics, function(m) {
    d <- hsm[is.finite(hsm[[m]]), , drop = FALSE]
    if (!nrow(d)) return(NULL)
    data.frame(species = d$species, species_label = d$species_label,
               sampler = d$sampler, predictor_set = d$predictor_set,
               algorithm = d$algorithm, realization = d$realization,
               metric = m, value = d[[m]], stringsAsFactors = FALSE)
  })
  long <- do.call(rbind, parts)
  if (is.null(long) || !nrow(long)) stop("no finite metric values to plot")
  # factor orders: sampler (paper order), algorithm, metric
  s_keys <- intersect(names(EXPERIMENT_SAMPLER_LABELS), unique(long$sampler))
  long$sampler <- factor(long$sampler, levels = s_keys,
                         labels = unname(EXPERIMENT_SAMPLER_LABELS[s_keys]))
  a_keys <- intersect(names(HSM_ALGO_LABELS), unique(long$algorithm))
  long$algorithm <- factor(long$algorithm, levels = a_keys,
                           labels = unname(HSM_ALGO_LABELS[a_keys]))
  m_keys <- intersect(names(HSM_METRIC_LABELS), metrics)
  long$metric <- factor(long$metric, levels = m_keys,
                        labels = unname(HSM_METRIC_LABELS[m_keys]))
  long
}

# --- house theme --------------------------------------------------------------
hsm_theme <- function(base_size = 9) {
  theme_bw(base_size = base_size) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          strip.background = element_rect(fill = "grey92", colour = NA),
          legend.position = "bottom", legend.title = element_blank(),
          axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          plot.title = element_text(face = "bold", size = base_size + 1))
}

#' Per-(predictor_set) violin figure: metric (cols) x algorithm (rows), one
#' violin per sampler. `species` NULL -> facet by species too (columns become
#' metric x species); otherwise a single species.
plot_hsm_violin <- function(hsm, predictor_set, species = NULL,
                            metrics = HSM_DEFAULT_METRICS,
                            sampler_palette = EXPERIMENT_SAMPLER_COLORS,
                            species_labels = NULL, base_size = 9,
                            title = NULL) {
  long <- hsm_metrics_long(hsm, metrics, species_labels = species_labels)
  long <- long[long$predictor_set == predictor_set, , drop = FALSE]
  if (!is.null(species)) long <- long[long$species == species, , drop = FALSE]
  if (!nrow(long)) stop("no rows for predictor_set='", predictor_set, "'")
  s_keys <- intersect(names(sampler_palette), names(EXPERIMENT_SAMPLER_LABELS))
  cols <- setNames(unname(sampler_palette[s_keys]),
                   unname(EXPERIMENT_SAMPLER_LABELS[s_keys]))
  facets <- if (is.null(species)) facet_grid(algorithm ~ species_label + metric, scales = "free_y")
            else facet_grid(algorithm ~ metric, scales = "free_y")
  ggplot(long, aes(.data$sampler, .data$value, fill = .data$sampler)) +
    geom_violin(scale = "width", alpha = 0.85, linewidth = 0.25) +
    stat_summary(fun = median, geom = "point", size = 1.1, colour = "grey15") +
    facets +
    scale_fill_manual(values = cols, drop = FALSE) +
    labs(x = NULL, y = NULL,
         title = title %||% sprintf("Downstream HSM accuracy (%s predictors)", predictor_set)) +
    hsm_theme(base_size)
}

# --- Dunn pairwise test (uniform+ vs each other sampler) ----------------------

#' Holm-adjusted Dunn test per (species, predictor_set, algorithm, metric),
#' returning the comparisons that involve uniform+ (the mcmc sampler). The
#' "Auniform_plus" relabel exploits dunn.test's alphabetical group ordering so
#' the reported sign of Z is "uniform+ vs other". Within-species contrast.
hsm_dunn_test <- function(hsm, metrics = HSM_DEFAULT_METRICS, ref = "mcmc",
                          min_per_group = 5L) {
  if (!requireNamespace("dunn.test", quietly = TRUE))
    stop("hsm_dunn_test requires the 'dunn.test' package.")
  hsm <- collect_hsm_metrics(hsm)
  if ("status" %in% names(hsm)) hsm <- hsm[hsm$status == "ok", , drop = FALSE]
  metrics <- intersect(metrics, names(hsm))
  relabel <- function(s) ifelse(s == ref, "Auniform_plus", s)

  grid <- unique(hsm[, c("species", "predictor_set", "algorithm")])
  out <- list()
  for (i in seq_len(nrow(grid))) {
    g <- grid[i, ]
    sub0 <- hsm[hsm$species == g$species & hsm$predictor_set == g$predictor_set &
                hsm$algorithm == g$algorithm, , drop = FALSE]
    for (m in metrics) {
      d <- sub0[is.finite(sub0[[m]]), , drop = FALSE]
      grp <- relabel(d$sampler)
      tab <- table(grp)
      keep <- names(tab)[tab >= min_per_group]
      d <- d[grp %in% keep, , drop = FALSE]; grp <- grp[grp %in% keep]
      if (length(unique(grp)) < 2L) next
      res <- tryCatch(
        dunn.test::dunn.test(x = d[[m]], g = grp, method = "holm",
                             altp = TRUE, table = FALSE, list = FALSE, kw = FALSE),
        error = function(e) NULL)
      if (is.null(res)) next
      comp <- res$comparisons[grepl("Auniform_plus", res$comparisons)]
      idx  <- which(grepl("Auniform_plus", res$comparisons))
      if (!length(idx)) next
      out[[length(out) + 1L]] <- data.frame(
        species = g$species, predictor_set = g$predictor_set,
        algorithm = g$algorithm, metric = m,
        comparison = gsub("Auniform_plus", "uniform+", comp),
        Z = round(res$Z[idx], 4),
        p_adj = signif(res$altP.adjusted[idx], 4),
        significant = res$altP.adjusted[idx] < 0.05,
        stringsAsFactors = FALSE)
    }
  }
  if (!length(out)) stop("no Dunn comparisons could be computed")
  do.call(rbind, out)
}

#' Quick significance-summary plot (share of comparisons p<0.05), faceted by
#' metric, split by which sampler uniform+ is compared against.
plot_hsm_dunn_summary <- function(dunn_df, base_size = 10) {
  dunn_df$other <- sub("^uniform\\+ - ", "", sub(" - uniform\\+$", "", dunn_df$comparison))
  dunn_df$sig <- ifelse(dunn_df$significant, "p<0.05", "ns")
  ggplot(dunn_df, aes(.data$other, fill = .data$sig)) +
    geom_bar(position = "fill", width = 0.6) +
    facet_wrap(~ metric, nrow = 1) +
    scale_fill_manual(values = c("p<0.05" = "#ffce72", "ns" = "#98c1d9")) +
    labs(x = "uniform+ compared against", y = "Relative proportion", fill = NULL) +
    theme_light(base_size = base_size) +
    theme(legend.position = "bottom", panel.grid = element_blank())
}

# --- median summary table (for the appendix) ----------------------------------
summarize_hsm <- function(hsm, metrics = HSM_DEFAULT_METRICS) {
  hsm <- collect_hsm_metrics(hsm)
  if ("status" %in% names(hsm)) hsm <- hsm[hsm$status == "ok", , drop = FALSE]
  metrics <- intersect(metrics, names(hsm))
  agg <- aggregate(hsm[metrics],
                   by = hsm[c("species", "predictor_set", "algorithm", "sampler")],
                   FUN = function(v) round(stats::median(v, na.rm = TRUE), 4))
  agg[order(agg$species, agg$predictor_set, agg$algorithm, agg$sampler), ]
}

# =============================================================================
# CROSS-SPECIES AGGREGATE VIEWS (algorithm dimension collapsed)
# =============================================================================
# Four ordered shades + four point shapes for the four species, used by the
# single-view overlay plot.
HSM_SPECIES_SHADES <- c("#cfe1f2", "#7bb0d6", "#3182bd", "#08519c")
HSM_SPECIES_SHAPES <- c(21L, 22L, 23L, 24L)
.hsm_sampler_relabel <- c(random = "RND", buffer = "buffer-out", mcmc = "uniform+")

# Species-label levels in species-key order (sp1..sp4) rather than alphabetical.
.hsm_species_order <- function(df) {
  if (all(c("species", "species_label") %in% names(df))) {
    u <- unique(df[c("species", "species_label")])
    as.character(u$species_label[order(u$species)])
  } else sort(unique(as.character(df$species_label)))
}

# Average each metric over the HSM algorithms within a realisation -> long form
# with the `algorithm` dimension collapsed (one value per species x sampler x
# predictor_set x realisation x metric). This is the "averaged over the HSMs" step.
hsm_algo_average <- function(long) {
  aggregate(value ~ species + species_label + sampler + predictor_set + realization + metric,
            data = long, FUN = function(v) mean(v, na.rm = TRUE))
}

#' Single-view species comparison: all four species in one plot. Samplers are
#' differentiated by HUE (fill colour), species by SHADE (alpha) plus a central
#' symbol; samplers on x, faceted by metric. Central symbol = `center` (median)
#' of the HSM-averaged value over realisations; bars = IQR.
plot_hsm_species_overlay <- function(hsm, predictor_set = "pc5",
                                     metrics = HSM_DEFAULT_METRICS, center = median,
                                     sampler_palette = EXPERIMENT_SAMPLER_COLORS,
                                     species_labels = NULL, base_size = 10, title = NULL) {
  long <- hsm_metrics_long(hsm, metrics, species_labels = species_labels)
  long <- long[long$predictor_set == predictor_set, , drop = FALSE]
  if (!nrow(long)) stop("no rows for predictor_set='", predictor_set, "'")
  ra <- hsm_algo_average(long)
  sp_lvls <- .hsm_species_order(ra)
  ra$species_label <- factor(ra$species_label, levels = sp_lvls)
  agg <- aggregate(value ~ species_label + sampler + metric, data = ra,
                   FUN = function(v) c(mid = center(v),
                                       lo = stats::quantile(v, .25, names = FALSE),
                                       hi = stats::quantile(v, .75, names = FALSE)))
  agg <- data.frame(agg[c("species_label", "sampler", "metric")],
                    mid = agg$value[, "mid"], lo = agg$value[, "lo"], hi = agg$value[, "hi"])
  nsp    <- length(sp_lvls)
  alphas <- setNames(seq(0.45, 1, length.out = nsp), sp_lvls)
  s_keys <- intersect(names(sampler_palette), names(EXPERIMENT_SAMPLER_LABELS))
  cols   <- setNames(unname(sampler_palette[s_keys]), unname(EXPERIMENT_SAMPLER_LABELS[s_keys]))
  dodge  <- position_dodge(width = 0.7)
  ggplot(agg, aes(.data$sampler, .data$mid, group = .data$species_label)) +
    geom_errorbar(aes(ymin = .data$lo, ymax = .data$hi, colour = .data$sampler,
                      alpha = .data$species_label),
                  width = 0.25, position = dodge, linewidth = 0.5, show.legend = FALSE) +
    geom_point(aes(fill = .data$sampler, alpha = .data$species_label,
                   shape = .data$species_label),
               size = 2.9, colour = "grey25", position = dodge) +
    facet_wrap(~ metric, scales = "free_y", nrow = 1) +
    scale_fill_manual(values = cols, name = "Sampler") +
    scale_colour_manual(values = cols, guide = "none") +
    scale_alpha_manual(values = alphas, name = "Species") +
    scale_shape_manual(values = HSM_SPECIES_SHAPES[seq_len(nsp)], name = "Species") +
    guides(fill = guide_legend(order = 1, override.aes = list(shape = 21, alpha = 1, size = 3))) +
    labs(x = NULL, y = NULL,
         title = title %||% sprintf("Species comparison — %s predictors (hue = sampler, shade = species; median ± IQR, HSMs averaged)", predictor_set)) +
    theme_bw(base_size = base_size) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
          strip.background = element_rect(fill = "grey92", colour = NA),
          axis.text.x = element_text(angle = 30, hjust = 1),
          legend.position = "right",
          plot.title = element_text(face = "bold", size = base_size + 1))
}

#' Four species in a row (columns), HSMs averaged: boxplots of the HSM-averaged
#' value over realisations. Samplers differentiated by HUE (fill), species by
#' SHADE (alpha) as well as the facet column; faceted metric (rows) x species (cols).
plot_hsm_species_row <- function(hsm, predictor_set = "pc5",
                                 metrics = HSM_DEFAULT_METRICS,
                                 sampler_palette = EXPERIMENT_SAMPLER_COLORS,
                                 species_labels = NULL, base_size = 9, title = NULL) {
  long <- hsm_metrics_long(hsm, metrics, species_labels = species_labels)
  long <- long[long$predictor_set == predictor_set, , drop = FALSE]
  if (!nrow(long)) stop("no rows for predictor_set='", predictor_set, "'")
  ra <- hsm_algo_average(long)
  sp_lvls <- .hsm_species_order(ra)
  ra$species_label <- factor(ra$species_label, levels = sp_lvls)
  nsp    <- length(sp_lvls)
  alphas <- setNames(seq(0.5, 1, length.out = nsp), sp_lvls)
  s_keys <- intersect(names(sampler_palette), names(EXPERIMENT_SAMPLER_LABELS))
  cols   <- setNames(unname(sampler_palette[s_keys]), unname(EXPERIMENT_SAMPLER_LABELS[s_keys]))
  ggplot(ra, aes(.data$sampler, .data$value, fill = .data$sampler, alpha = .data$species_label)) +
    geom_boxplot(width = 0.62, linewidth = 0.3, outlier.size = 0.4, outlier.alpha = 0.4) +
    facet_grid(metric ~ species_label, scales = "free_y", switch = "y") +
    scale_fill_manual(values = cols, name = "Sampler") +
    scale_alpha_manual(values = alphas, guide = "none") +
    labs(x = NULL, y = NULL,
         title = title %||% sprintf("Species comparison, HSMs averaged — %s predictors (hue = sampler, shade = species)", predictor_set)) +
    hsm_theme(base_size) + theme(strip.placement = "outside")
}

#' Aggregate-metric table: one row per `by` group, one column per metric.
#' algo_average=TRUE first means over the HSM algorithms within a realisation,
#' then applies `fun` over realisations (so each algorithm is weighted equally).
hsm_aggregate_table <- function(hsm,
                                metrics = c("auc", "tss", "boyce", "sensitivity",
                                            "specificity", "cor_truth", "rmse_truth", "r2"),
                                by = c("predictor_set", "species", "sampler"),
                                algo_average = TRUE, fun = median, digits = 3,
                                relabel_sampler = TRUE) {
  hsm <- collect_hsm_metrics(hsm)
  if ("status" %in% names(hsm)) hsm <- hsm[hsm$status == "ok", , drop = FALSE]
  metrics <- intersect(metrics, names(hsm))
  if (algo_average) {
    grp <- c(setdiff(by, "algorithm"), "realization")
    hsm <- aggregate(hsm[metrics], by = hsm[grp], FUN = function(v) mean(v, na.rm = TRUE))
  }
  agg <- aggregate(hsm[metrics], by = hsm[by], FUN = function(v) round(fun(v, na.rm = TRUE), digits))
  agg <- agg[do.call(order, agg[by]), , drop = FALSE]
  if (relabel_sampler && "sampler" %in% names(agg))
    agg$sampler <- ifelse(agg$sampler %in% names(.hsm_sampler_relabel),
                          .hsm_sampler_relabel[agg$sampler], agg$sampler)
  rownames(agg) <- NULL
  agg
}

#' Difference table: every sampler minus the reference (default RND), per metric,
#' per `by` group -- makes the "who beats the benchmark, and by how much" explicit.
hsm_delta_table <- function(hsm, metrics = HSM_DEFAULT_METRICS, ref = "random",
                            by = c("predictor_set", "species"),
                            algo_average = TRUE, digits = 3) {
  tab <- hsm_aggregate_table(hsm, metrics, by = c(by, "sampler"),
                             algo_average = algo_average, fun = median, digits = 8,
                             relabel_sampler = FALSE)
  keys <- unique(tab[by]); out <- list()
  for (i in seq_len(nrow(keys))) {
    sub <- merge(keys[i, , drop = FALSE], tab, by = by)
    rr  <- sub[sub$sampler == ref, , drop = FALSE]
    if (!nrow(rr)) next
    for (s in setdiff(unique(sub$sampler), ref)) {
      sr <- sub[sub$sampler == s, , drop = FALSE]
      row <- keys[i, , drop = FALSE]
      slab <- ifelse(s %in% names(.hsm_sampler_relabel), .hsm_sampler_relabel[s], s)
      rlab <- ifelse(ref %in% names(.hsm_sampler_relabel), .hsm_sampler_relabel[ref], ref)
      row$comparison <- paste0(slab, " - ", rlab)
      for (m in metrics) row[[m]] <- round(sr[[m]] - rr[[m]], digits)
      out[[length(out) + 1L]] <- row
    }
  }
  do.call(rbind, out)
}

# =============================================================================
# TUNING-GRID HEATMAP (env_cutoff x species_cutoff)
# =============================================================================
#' Heatmap of a metric's median over the (environmental.cutof.percentile x
#' species.cutoff.threshold) grid, pooled over species + algorithms + realisations.
#' If `ref` (a baseline numeric, e.g. the fixed RND median) is supplied, the DELTA
#' value - ref is shown on a diverging scale (blue = uniform+ better for
#' higher_better metrics); otherwise the raw value on a viridis scale. `baseline`
#' = c(random=, buffer=) annotates the fixed references in the subtitle.
plot_tune_heatmap <- function(tune_df, metric, ref = NULL, baseline = NULL,
                              higher_better = TRUE, title = NULL) {
  d <- collect_hsm_metrics(tune_df)
  if ("status" %in% names(d)) d <- d[d$status == "ok", , drop = FALSE]
  d <- d[is.finite(d[[metric]]), , drop = FALSE]
  if (!all(c("env_cutoff", "species_cutoff") %in% names(d)))
    stop("tune_df needs env_cutoff + species_cutoff columns")
  agg <- aggregate(d[[metric]], by = list(env = d$env_cutoff, sp = d$species_cutoff),
                   FUN = function(v) stats::median(v, na.rm = TRUE))
  names(agg)[3] <- "val"
  agg$env <- factor(agg$env, levels = sort(unique(agg$env)))
  agg$sp  <- factor(agg$sp,  levels = sort(unique(agg$sp), decreasing = TRUE))
  if (!is.null(ref)) {
    agg$plotval <- agg$val - ref
    lim <- max(abs(agg$plotval), na.rm = TRUE); lim <- if (lim == 0) 1e-6 else lim
    lo <- if (higher_better) "#b2182b" else "#2166ac"
    hi <- if (higher_better) "#2166ac" else "#b2182b"
    fillscale <- scale_fill_gradient2(low = lo, mid = "grey95", high = hi,
                                      midpoint = 0, limits = c(-lim, lim), name = "Δ vs RND")
    lab <- sprintf("%+.3f", agg$plotval)
  } else {
    agg$plotval <- agg$val
    fillscale <- scale_fill_viridis_c(option = "viridis",
                                      direction = if (higher_better) 1 else -1, name = metric)
    lab <- sprintf("%.3f", agg$val)
  }
  sub <- if (!is.null(baseline))
    sprintf("fixed baselines: RND = %.3f,  buffer-out = %.3f",
            baseline["random"], baseline["buffer"]) else NULL
  ggplot(agg, aes(.data$env, .data$sp, fill = .data$plotval)) +
    geom_tile(colour = "white", linewidth = 0.6) +
    geom_text(aes(label = lab), size = 3) +
    fillscale +
    labs(x = "environmental.cutof.percentile", y = "species.cutoff.threshold",
         title = title %||% metric, subtitle = sub) +
    theme_minimal(base_size = 10) + theme(panel.grid = element_blank())
}

# =============================================================================
# FULL PER-SPECIES REPORT (all metrics; where uniform+ beats the naive samplers)
# =============================================================================
# Direction of "better" for every metric (TRUE = higher is better). Ordered
# truth-recovery first, then discrimination, then the rest.
HSM_METRIC_HIGHER <- c(cor_truth = TRUE, rmse_truth = FALSE, auc = TRUE, boyce = TRUE,
                       tss = TRUE, sensitivity = TRUE, specificity = TRUE, kappa = TRUE,
                       cor = TRUE, r2 = TRUE, rmse = FALSE, deviance = FALSE)
HSM_ALL_METRICS <- names(HSM_METRIC_HIGHER)

#' All-metrics violin ("flute") grid for ONE species x predictor_set: one panel
#' per metric, three sampler violins each (distribution over realisations x HSM
#' algorithms), median point overlaid. Free y-scales since metrics differ in range.
plot_hsm_species_metrics <- function(hsm, species, predictor_set,
                                     metrics = HSM_ALL_METRICS,
                                     sampler_palette = EXPERIMENT_SAMPLER_COLORS,
                                     species_labels = NULL, base_size = 9, ncol = 3, title = NULL) {
  long <- hsm_metrics_long(hsm, metrics, species_labels = species_labels)
  long <- long[long$species == species & long$predictor_set == predictor_set, , drop = FALSE]
  if (!nrow(long)) stop("no rows for ", species, " / ", predictor_set)
  s_keys <- intersect(names(sampler_palette), names(EXPERIMENT_SAMPLER_LABELS))
  cols   <- setNames(unname(sampler_palette[s_keys]), unname(EXPERIMENT_SAMPLER_LABELS[s_keys]))
  splab  <- as.character(long$species_label[1])
  ggplot(long, aes(.data$sampler, .data$value, fill = .data$sampler)) +
    geom_violin(scale = "width", alpha = 0.85, linewidth = 0.2) +
    stat_summary(fun = median, geom = "point", size = 1, colour = "grey15") +
    facet_wrap(~ metric, scales = "free_y", ncol = ncol) +
    scale_fill_manual(values = cols, name = NULL) +
    labs(x = NULL, y = NULL,
         title = title %||% sprintf("%s  —  %s predictors", splab, predictor_set),
         subtitle = "violins: distribution over 50 realisations x 5 HSM algorithms; point = median. Internal metrics = vs sampled PA; *_truth = vs known suitability") +
    theme_bw(base_size = base_size) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
          strip.background = element_rect(fill = "grey92", colour = NA),
          axis.text.x = element_text(angle = 25, hjust = 1),
          legend.position = "none", plot.title = element_text(face = "bold"))
}

#' Win matrix: per (species, predictor_set, metric), the focal sampler (uniform+)
#' median vs the naive samplers (RND + buffer-out), classified beats both / one /
#' none, with the signed margin vs the BETTER naive sampler. `tol` ignores ties.
hsm_win_matrix <- function(hsm, metrics = HSM_ALL_METRICS,
                           ref_samplers = c("random", "buffer"), focal = "mcmc", tol = 0.002) {
  hsm <- collect_hsm_metrics(hsm)
  if ("status" %in% names(hsm)) hsm <- hsm[hsm$status == "ok", , drop = FALSE]
  metrics <- intersect(metrics, names(hsm))
  med <- aggregate(hsm[metrics],
                   by = hsm[c("species", "species_label", "predictor_set", "sampler")],
                   FUN = function(v) median(v, na.rm = TRUE))
  keys <- unique(med[c("species", "species_label", "predictor_set")]); out <- list()
  for (i in seq_len(nrow(keys))) {
    sub <- merge(keys[i, , drop = FALSE], med)
    g1 <- function(s, m) { x <- sub[sub$sampler == s, m]; if (length(x)) x[1] else NA_real_ }
    for (m in metrics) {
      higher <- isTRUE(HSM_METRIC_HIGHER[[m]])
      u  <- g1(focal, m); rs <- vapply(ref_samplers, g1, numeric(1), m = m)
      better <- function(a, b) if (higher) (a > b + tol) else (a < b - tol)
      nbeat  <- sum(vapply(rs, function(rv) isTRUE(better(u, rv)), logical(1)))
      naive_best <- if (higher) max(rs, na.rm = TRUE) else min(rs, na.rm = TRUE)
      out[[length(out) + 1L]] <- data.frame(keys[i, , drop = FALSE], metric = m,
        uniform_plus = round(u, 3), RND = round(rs[["random"]], 3), buffer = round(rs[["buffer"]], 3),
        delta_vs_best_naive = round(if (higher) u - naive_best else naive_best - u, 3),
        category = c("beats none", "beats one", "beats both")[nbeat + 1L],
        stringsAsFactors = FALSE)
    }
  }
  do.call(rbind, out)
}

#' Tile plot of the win matrix for one predictor_set: metric (rows) x species
#' (cols), green/amber/red = uniform+ beats both / one / none naive samplers,
#' annotated with the signed margin vs the better naive sampler.
plot_hsm_win_matrix <- function(hsm, predictor_set, metrics = HSM_ALL_METRICS, title = NULL) {
  wm <- hsm_win_matrix(hsm, metrics)
  wm <- wm[wm$predictor_set == predictor_set, , drop = FALSE]
  if (!nrow(wm)) stop("no win-matrix rows for predictor_set='", predictor_set, "'")
  wm$metric <- factor(wm$metric, levels = rev(metrics),
                      labels = rev(unname(ifelse(metrics %in% names(HSM_METRIC_LABELS),
                                                 HSM_METRIC_LABELS[metrics], metrics))))
  sp_lvls <- .hsm_species_order(wm)
  wm$species_label <- factor(wm$species_label, levels = sp_lvls)
  wm$category <- factor(wm$category, levels = c("beats both", "beats one", "beats none"))
  catcols <- c("beats both" = "#1a9850", "beats one" = "#fee08b", "beats none" = "#d73027")
  ggplot(wm, aes(.data$species_label, .data$metric, fill = .data$category)) +
    geom_tile(colour = "white", linewidth = 0.6) +
    geom_text(aes(label = sprintf("%+.3f", .data$delta_vs_best_naive)), size = 2.7) +
    scale_fill_manual(values = catcols, name = "uniform+ vs naive\n(RND & buffer-out)", drop = FALSE) +
    labs(x = NULL, y = NULL,
         title = title %||% sprintf("Where uniform+ beats the naive samplers — %s predictors", predictor_set),
         subtitle = "cell = margin of uniform+ vs the BETTER naive sampler (>0 better); median over algos x realisations") +
    theme_minimal(base_size = 10) +
    theme(panel.grid = element_blank(), axis.text.x = element_text(face = "bold"),
          plot.title = element_text(face = "bold"))
}

# --- saver (cairo PDF / 600-dpi PNG) ------------------------------------------
save_hsm_figure <- function(plot, file, width = 9, height = 7, dpi = 600) {
  if (grepl("\\.pdf$", file, ignore.case = TRUE))
    ggsave(file, plot, device = grDevices::cairo_pdf, width = width, height = height, units = "in")
  else
    ggsave(file, plot, device = "png", width = width, height = height, units = "in", dpi = dpi)
  invisible(file)
}
