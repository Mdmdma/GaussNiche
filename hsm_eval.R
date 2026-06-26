# hsm_eval.R --- downstream habitat-suitability-model (HSM/SDM) evaluation
# =============================================================================
# The "downstream" half of the USE-paper experiment (Da Re et al. 2023,
# Analysis/1_Modelling_framework.R), ported to GaussNiche. For one Bernoulli
# realisation of one species under one pseudo-absence sampler, fit the SAME five
# model families USE used and score their predictive accuracy:
#
#   GLM     stats::glm     (binomial, LINEAR terms -- exactly USE's `Observed ~ .`)
#   GAM     mgcv::gam      (binomial, per-predictor thin-plate smooths)
#   RF      ranger::ranger (regression on 0/1 -> OOB r.squared, like USE's rf$rsq)
#   BRT     gbm::gbm       (distribution = "bernoulli")
#   Maxent  maxnet::maxnet (pure-R Maxent -- no rJava / maxent.jar / Java)
#
# Two metric families per fit (cf. USE):
#   * internal  -- on a held-out test split, presences vs SAMPLED pseudo-absences
#                  (discrimination against the sampler): AUC, COR, TSS, Kappa,
#                  Sensitivity, Specificity (at the max-TSS threshold), Deviance,
#                  RMSE, Boyce (CBI), and a per-algorithm R^2.
#   * truth     -- USE's rmse.distribution / cor.distribution: predicted vs the
#                  KNOWN Gaussian suitability `suit` over the whole background
#                  (no test split, no circularity -- the cleanest signal, only
#                  possible because GaussNiche simulates the truth).
#
# CRITICAL (faithfulness): the GLM is LINEAR, like USE. The niche is exactly
# Gaussian, so logit(suit) is exactly quadratic in the PCs; a quadratic GLM would
# be near-oracle and blind to the sampler, washing out the contrast. Keep it linear.
#
# Design for use inside virtualSpecies_nd()'s furrr workers:
#   * make_sdm_hook() returns a CLOSURE that fits + scores everything for one
#     (sampler, realisation) and returns a tidy data.frame of stat rows. It is
#     passed to virtualSpecies_nd(sdm_hook = ...) and called single-pass on the
#     PA set already drawn for the diagnostics (so the expensive pa_mcmc draw is
#     never re-paid).
#   * Every backend fit/predict is wrapped tryCatch -> NA so an unstable
#     rare-species fit cannot abort the whole (sampler, realisation) task.
#   * All package calls are fully-qualified (pkg::fun) so furrr workers resolve
#     them via the inherited .libPaths(); the hook also prepends `libpaths`.
#
# Source this file (assumes wd = repo root). Pairs with hsm_plots.R (figures +
# Dunn test) and run_5d_hsm.R (the driver).
# =============================================================================

# --- algorithm registry -------------------------------------------------------
HSM_ALGORITHMS <- c("glm", "gam", "rf", "brt", "maxent")

# Stat columns produced per (predictor_set, algorithm) row, in order.
HSM_METRIC_COLS <- c("auc", "cor", "tss", "kappa", "sensitivity", "specificity",
                     "deviance", "rmse", "boyce", "r2", "rmse_truth", "cor_truth")

# =============================================================================
# METRIC PRIMITIVES (base-R; no heavy deps)
# =============================================================================

# Rank-based ROC AUC (Mann-Whitney U); NA if a class is empty.
.hsm_auc <- function(obs, pred) {
  ok <- is.finite(pred) & !is.na(obs)
  obs <- obs[ok]; pred <- pred[ok]
  pos <- pred[obs == 1]; neg <- pred[obs == 0]
  if (!length(pos) || !length(neg)) return(NA_real_)
  r <- rank(c(pos, neg))
  (sum(r[seq_along(pos)]) - length(pos) * (length(pos) + 1) / 2) /
    (length(pos) * length(neg))
}

.hsm_rmse <- function(obs, pred) {
  ok <- is.finite(pred) & is.finite(obs)
  if (!any(ok)) return(NA_real_)
  sqrt(mean((obs[ok] - pred[ok])^2))
}

# Binomial deviance of probabilistic predictions (clamped away from 0/1).
.hsm_deviance <- function(obs, pred) {
  ok <- is.finite(pred) & !is.na(obs)
  obs <- obs[ok]; p <- pmin(pmax(pred[ok], 1e-6), 1 - 1e-6)
  if (!length(obs)) return(NA_real_)
  -2 * sum(obs * log(p) + (1 - obs) * log(1 - p))
}

# Threshold-based metrics at the threshold maximising TSS (= sens + spec - 1),
# mirroring sdm's threshold_based "max(se+sp)". Returns sens/spec/tss/kappa.
.hsm_threshold_metrics <- function(obs, pred, n_thr = 200L) {
  ok <- is.finite(pred) & !is.na(obs)
  obs <- obs[ok]; pred <- pred[ok]
  out <- c(sensitivity = NA_real_, specificity = NA_real_,
           tss = NA_real_, kappa = NA_real_)
  np <- sum(obs == 1); na <- sum(obs == 0)
  if (np == 0 || na == 0) return(out)
  thr <- sort(unique(c(0, pred, 1)))
  if (length(thr) > n_thr) thr <- stats::quantile(thr, seq(0, 1, length.out = n_thr), names = FALSE)
  best_tss <- -Inf
  for (t in thr) {
    pp <- pred >= t
    tp <- sum(pp & obs == 1); fn <- np - tp
    fp <- sum(pp & obs == 0); tn <- na - fp
    se <- tp / (tp + fn); sp <- tn / (tn + fp)
    tss <- se + sp - 1
    if (is.finite(tss) && tss > best_tss) {
      best_tss <- tss
      n <- np + na
      po <- (tp + tn) / n
      pe <- ((tp + fp) * (tp + fn) + (fn + tn) * (fp + tn)) / n^2
      kap <- if (pe == 1) NA_real_ else (po - pe) / (1 - pe)
      out <- c(sensitivity = se, specificity = sp, tss = tss, kappa = kap)
    }
  }
  out
}

#' Continuous Boyce Index (CBI), Hirzel et al. 2006 -- compact reimplementation
#' matching ecospat::ecospat.boyce's continuous (nclass = 0) form, so we avoid
#' the heavy ecospat dependency tree.
#'
#' @param pred_eval predicted suitability at all evaluation points.
#' @param obs       0/1 indicator aligned with pred_eval (1 = presence).
#' @param window    moving-window width as a fraction of the prediction range.
#' @param n_bins    number of moving windows.
#' @return Spearman correlation between window midpoint and predicted/expected
#'         frequency ratio; NA if degenerate.
boyce_index <- function(pred_eval, obs, window = 0.1, n_bins = 101L) {
  ok <- is.finite(pred_eval) & !is.na(obs)
  pred_eval <- pred_eval[ok]; obs <- obs[ok]
  pres <- pred_eval[obs == 1]
  if (length(pres) < 3L) return(NA_real_)
  rng <- range(pred_eval)
  if (diff(rng) <= 0) return(NA_real_)
  w  <- window * diff(rng)
  ctr <- seq(rng[1] + w / 2, rng[2] - w / 2, length.out = n_bins)
  pe <- vapply(ctr, function(c) {
    lo <- c - w / 2; hi <- c + w / 2
    f_all  <- mean(pred_eval >= lo & pred_eval <= hi)
    if (f_all == 0) return(NA_real_)
    mean(pres >= lo & pres <= hi) / f_all
  }, numeric(1))
  keep <- is.finite(pe)
  if (sum(keep) < 3L) return(NA_real_)
  suppressWarnings(stats::cor(ctr[keep], pe[keep], method = "spearman"))
}

# =============================================================================
# PER-ALGORITHM FIT + PREDICT  (each returns NULL/NA on failure)
# =============================================================================

# Pick a GAM basis dimension k so the total basis stays well below n (mgcv errors
# when total basis dim >= n); keeps GAM stable for the small rare-species fits.
.hsm_gam_k <- function(n, n_pred) {
  budget <- floor((0.5 * n) / max(1L, n_pred))
  as.integer(max(3L, min(10L, budget)))
}

# Fit one model family. `train` has the predictor columns + a 0/1 `pa` column.
# Warnings (glm fitted-probabilities-0/1, gbm OOB, etc.) are non-fatal -> muffled;
# only hard errors collapse to NULL (-> NA metrics, never an aborted task).
.hsm_fit_one <- function(algo, train, predictors, seed) {
  suppressWarnings(tryCatch({
    if (algo == "glm") {
      stats::glm(stats::reformulate(predictors, response = "pa"),
                 data = train, family = stats::binomial())
    } else if (algo == "gam") {
      k <- .hsm_gam_k(nrow(train), length(predictors))
      terms <- vapply(predictors, function(v) sprintf("s(%s, k = %d)", v, k), character(1))
      mgcv::gam(stats::reformulate(terms, response = "pa"),
                data = train, family = stats::binomial(), method = "REML")
    } else if (algo == "rf") {
      # Regression forest on the 0/1 response -> OOB r.squared (USE's rf$rsq).
      ranger::ranger(stats::reformulate(predictors, response = "pa"),
                     data = train, num.trees = 500L, num.threads = 1L, seed = seed)
    } else if (algo == "brt") {
      set.seed(seed)
      gbm::gbm(stats::reformulate(predictors, response = "pa"),
               data = train, distribution = "bernoulli",
               n.trees = 1500L, interaction.depth = 2L, shrinkage = 0.01,
               bag.fraction = 0.75, n.minobsinnode = 5L,
               verbose = FALSE, n.cores = 1L)
    } else if (algo == "maxent") {
      set.seed(seed)
      maxnet::maxnet(p = train$pa, data = train[, predictors, drop = FALSE])
    } else stop("unknown algorithm: ", algo)
  }, error = function(e) NULL))
}

# Best BRT iteration is a property of the fitted model, not of newdata, so compute
# it ONCE (at fit time) and reuse for every predict instead of re-scanning all
# trees on each call -- bit-identical, ~halves the gbm.perf work.
.hsm_brt_best_nt <- function(model) {
  if (is.null(model)) return(NULL)
  nt <- tryCatch(suppressWarnings(gbm::gbm.perf(model, method = "OOB", plot.it = FALSE)),
                 error = function(e) NA)
  if (is.na(nt) || nt < 1) model$n.trees else nt
}

# Predict suitability in [0,1] for one fitted model; NA vector on failure.
# `nt` = precomputed BRT best-iteration (from .hsm_brt_best_nt); ignored otherwise.
.hsm_predict_one <- function(algo, model, newdata, predictors, nt = NULL) {
  if (is.null(model)) return(rep(NA_real_, nrow(newdata)))
  p <- tryCatch({
    if (algo == "glm" || algo == "gam") {
      as.numeric(stats::predict(model, newdata = newdata, type = "response"))
    } else if (algo == "rf") {
      as.numeric(stats::predict(model, data = newdata, num.threads = 1L)$predictions)
    } else if (algo == "brt") {
      if (is.null(nt)) nt <- .hsm_brt_best_nt(model)   # fallback if not hoisted
      as.numeric(gbm::predict.gbm(model, newdata = newdata, n.trees = nt, type = "response"))
    } else if (algo == "maxent") {
      as.numeric(stats::predict(model, newdata = newdata, type = "cloglog", clamp = TRUE))
    } else rep(NA_real_, nrow(newdata))
  }, error = function(e) rep(NA_real_, nrow(newdata)))
  pmin(pmax(p, 0), 1)   # clamp to [0,1] (regression RF can overshoot slightly)
}

# Per-algorithm R^2, branched exactly as USE: Dsquared for glm/gam, ranger OOB
# r.squared for rf, cor(obs,pred)^2 for brt/maxent.
.hsm_r2 <- function(algo, model, obs_test, pred_test) {
  if (is.null(model)) return(NA_real_)
  tryCatch({
    if (algo == "glm" || algo == "gam") {
      (model$null.deviance - model$deviance) / model$null.deviance  # D^2
    } else if (algo == "rf") {
      r <- model$r.squared
      if (is.null(r)) NA_real_ else as.numeric(r)
    } else {
      cc <- suppressWarnings(stats::cor(obs_test, pred_test))
      if (is.na(cc)) NA_real_ else cc^2
    }
  }, error = function(e) NA_real_)
}

# =============================================================================
# ONE (predictor_set x algorithm) GRID FOR ONE TRAINING SET
# =============================================================================

# Fit + score every algorithm for one predictor set; returns a tidy data.frame
# with one row per algorithm (HSM_METRIC_COLS + status + counts).
.hsm_eval_predictor_set <- function(train, test, background, predictors,
                                    algorithms, truth, seed) {
  do.call(rbind, lapply(algorithms, function(algo) {
    s <- seed + as.integer(factor(algo, levels = HSM_ALGORITHMS)) * 101L
    model    <- .hsm_fit_one(algo, train, predictors, s)
    nt       <- if (algo == "brt") .hsm_brt_best_nt(model) else NULL  # once, reused below
    pred_te  <- .hsm_predict_one(algo, model, test, predictors, nt = nt)
    pred_bg  <- .hsm_predict_one(algo, model, background, predictors, nt = nt)
    thr      <- .hsm_threshold_metrics(test$pa, pred_te)
    row <- data.frame(
      algorithm   = algo,
      status      = if (is.null(model)) "fit_failed" else "ok",
      auc         = .hsm_auc(test$pa, pred_te),
      cor         = suppressWarnings(stats::cor(test$pa, pred_te)),
      tss         = unname(thr["tss"]),
      kappa       = unname(thr["kappa"]),
      sensitivity = unname(thr["sensitivity"]),
      specificity = unname(thr["specificity"]),
      deviance    = .hsm_deviance(test$pa, pred_te),
      rmse        = .hsm_rmse(test$pa, pred_te),
      boyce       = boyce_index(pred_te, test$pa),
      r2          = .hsm_r2(algo, model, test$pa, pred_te),
      rmse_truth  = .hsm_rmse(truth, pred_bg),
      cor_truth   = suppressWarnings(stats::cor(truth, pred_bg)),
      stringsAsFactors = FALSE)
    rm(model); row
  }))
}

# =============================================================================
# THE HOOK BUILDER  (passed to virtualSpecies_nd(sdm_hook = ...))
# =============================================================================

#' Build the SDM hook for one species.
#'
#' Returns a closure `function(pres_r, pseudo_r, background, pc_cols,
#' realization, sampler)` that virtualSpecies_nd() calls once per (sampler,
#' realisation), reusing the already-drawn PA set. It fits + scores all five
#' model families for every predictor set on a single train/test split and
#' returns a tidy data.frame of stat rows (NULL never returned for a usable PA
#' set; tiny-N sets yield status = "skip_few_pres" rows).
#'
#' @param species       Species key, stamped onto every row.
#' @param predictor_sets Named list of character vectors, e.g.
#'                       list(pc5 = c("PC1",...,"PC5"), raw12 = <12 layer names>).
#'                       Each must be columns present in pres_r/pseudo_r/background.
#' @param algorithms    Subset of HSM_ALGORITHMS to fit.
#' @param test_frac     Held-out test fraction (default 0.3, as USE's 30%).
#' @param min_pres      Minimum presences in TRAIN to attempt fitting; below this
#'                      the realisation is recorded as skip_few_pres.
#' @param seed          Base RNG seed; the realisation index is layered on.
#' @param libpaths      Library path(s) to prepend inside furrr workers so the
#'                      SDM packages resolve (defaults to the current .libPaths()).
#' @param source_file   Absolute path to this file (hsm_eval.R). Belt-and-suspenders
#'                      against furrr globals-shipping gaps: each worker re-sources
#'                      it once if the helpers are absent. hsm_eval.R has no
#'                      top-level library() calls, so sourcing is side-effect-free.
make_sdm_hook <- function(species, predictor_sets,
                          algorithms = HSM_ALGORITHMS,
                          test_frac = 0.3, min_pres = 12L, seed = 777L,
                          libpaths = .libPaths(), source_file = NULL) {
  force(species); force(predictor_sets); force(algorithms)
  force(test_frac); force(min_pres); force(seed); force(libpaths); force(source_file)
  all_pred <- unique(unlist(predictor_sets, use.names = FALSE))

  function(pres_r, pseudo_r, background, pc_cols, realization, sampler) {
    .libPaths(unique(c(libpaths, .libPaths())))   # workers: ensure SDM libs visible
    if (!is.null(source_file) && file.exists(source_file) &&
        !exists(".hsm_eval_predictor_set", mode = "function"))
      source(source_file)                          # workers: ensure helpers present

    set_names  <- names(predictor_sets)
    # regular (predictor_set x algorithm) skeleton, filled with NA + a status.
    skeleton <- function(status, np_tr = NA_integer_, np_te = NA_integer_,
                         n_tr = NA_integer_, n_te = NA_integer_) {
      grid <- expand.grid(predictor_set = set_names, algorithm = algorithms,
                          stringsAsFactors = FALSE)
      base <- data.frame(species = species, sampler = sampler,
                         realization = realization, grid, status = status,
                         n_pres_train = np_tr, n_pres_test = np_te,
                         n_train = n_tr, n_test = n_te, stringsAsFactors = FALSE)
      for (m in HSM_METRIC_COLS) base[[m]] <- NA_real_
      base[, c("species", "sampler", "realization", "predictor_set", "algorithm",
               "status", "n_pres_train", "n_pres_test", "n_train", "n_test",
               HSM_METRIC_COLS)]
    }

    # Assemble presence/pseudo-absence rows carrying every needed predictor col.
    need <- c(all_pred, "suit")
    if (!all(need %in% names(pres_r)) || !all(need %in% names(pseudo_r)))
      return(skeleton("missing_cols"))
    pres <- pres_r[, all_pred, drop = FALSE];   pres$pa <- 1L
    abs  <- pseudo_r[, all_pred, drop = FALSE]; abs$pa  <- 0L

    # Stratified train/test split (seeded per realisation; RNG restored by caller).
    set.seed(seed + realization)
    sp_idx <- function(n) { k <- max(1L, round((1 - test_frac) * n)); sample.int(n, k) }
    ip <- sp_idx(nrow(pres)); ia <- sp_idx(nrow(abs))
    train <- rbind(pres[ip, , drop = FALSE], abs[ia, , drop = FALSE])
    test  <- rbind(pres[-ip, , drop = FALSE], abs[-ia, , drop = FALSE])
    np_tr <- length(ip); np_te <- nrow(pres) - np_tr

    if (np_tr < min_pres || np_te < 3L || nrow(test) < 5L)
      return(skeleton("skip_few_pres", np_tr, np_te, nrow(train), nrow(test)))

    truth <- background$suit
    rows <- do.call(rbind, lapply(set_names, function(ps) {
      preds <- predictor_sets[[ps]]
      df <- .hsm_eval_predictor_set(train, test, background, preds,
                                    algorithms, truth, seed + realization)
      cbind(species = species, sampler = sampler, realization = realization,
            predictor_set = ps,
            n_pres_train = np_tr, n_pres_test = np_te,
            n_train = nrow(train), n_test = nrow(test), df,
            stringsAsFactors = FALSE)
    }))
    # reorder to the skeleton's column order
    rows[, c("species", "sampler", "realization", "predictor_set", "algorithm",
             "status", "n_pres_train", "n_pres_test", "n_train", "n_test",
             HSM_METRIC_COLS)]
  }
}
