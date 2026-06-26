---
name: publication-figures
description: Regenerate the appendix MCMC-diagnostic figures of the markov-chain-sampler paper ("Pseudo-absences generation through a Markov Chain sampler") — autocorrelation, Gelman-Rubin, trace/posterior — as publication-ready vector PDFs where the manuscript expects them. Use when those appendix figures look low-DPI / off-house-style / decoupled from the current code; the authoritative producer is the USE.MCMC-side publication-figures skill (the figures come from its vignette). NOTE: the former Results "distribution of PC values" figure (generate-combined-plot) was REMOVED from the paper; the current Results figures are the GaussNiche §2.3 sampler-comparison boxplots produced by run_2d_experiment.R / run_5d_experiment.R, NOT this skill. NOT for in-vignette exploratory plots or the GaussNiche per-species report PDFs.
---

# Publication figures for the markov-chain-sampler paper

This skill regenerates **only four figures** that the manuscript embeds, makes
them publication-grade, and writes them under the **exact filenames** the paper's
`\includegraphics` paths require. It is deliberately scoped: it does not touch the
many exploratory plots in the vignette or the GaussNiche report PDFs.

> **These figures are produced by the USE.MCMC vignette, not by GaussNiche.** The
> author's chosen approach is to improve them *in place* in the vignette chunks,
> and the authoritative producer-side guide lives at
> `../USE.MCMC/.claude/skills/publication-figures/SKILL.md` — **prefer it.** Do
> not rebuild the chains or figures from scratch in GaussNiche (the vignette
> already builds `coda.chain.list` and plots the diagnostics); edit the chunks.
> Final-paper figures are **vector PDFs**, produced by the vignette and copied via
> `make figures`. (`trace_plots` is now a real vignette chunk producing
> `trace_plots-1.pdf`, no longer a static hand-made PNG.) Any remaining 2-D/3-D or
> step-count caption notes are tracked in that skill's "Known caption mismatches"
> section.

The manuscript lives in a **sibling repo**:
`../markov-chain-sampler-paper/` (relative to the GaussNiche repo root, which is
your working directory). The graphics directory the paper reads from is:

```
../markov-chain-sampler-paper/graphics/
```

> Run everything with the working directory at the GaussNiche repo root (see the
> root `.claude/CLAUDE.md`). Sibling paths are relative to that root. Do not bake
> in any host- or user-specific absolute path.

> **Vignette references use knitr *chunk labels*** (e.g. `generate-combined-plot`),
> not line numbers. The MCMC vignette is regenerated/edited often, so line numbers
> drift; chunk labels are unique and stable. Find a chunk by searching its label in
> `../USE.MCMC/vignettes/insights-on-MCMC-pseudo-absence-sampling-vignette.Rmd`.

## The filename coupling — DO NOT RENAME (this breaks the paper)

The paper's `\includegraphics` calls hard-code these basenames in
`../markov-chain-sampler-paper/graphics/`. Renaming silently breaks a figure.
Several names are inherited from knitr chunk names where **the chunk name *is*
the filename**, misspelling included.

| In-scope figure | Required filename in `graphics/` | Paper ref |
| --- | --- | --- |
| Appendix — autocorrelation | `autocorreltation-plot-1.pdf` *(sic — keep the typo)* | `text/appendix/appendix.tex`, `fig:autocorrelation_analysis` |
| Appendix — Gelman-Rubin | `gelman-plot-1.pdf` | `text/appendix/appendix.tex`, `fig:gelman-rubin-appendix` |
| Appendix — trace + posterior | `trace_plots-1.pdf` | `text/appendix/appendix.tex`, `fig:traceplot` |

All three appendix figures are **vector PDFs** produced by the USE.MCMC vignette
and copied into `graphics/` by the paper's `make figures` target — all three are
listed in the paper `Makefile`, including `trace_plots-1.pdf`, which is now a real
vignette chunk (no longer a static hand-made PNG).

> **Results figure moved out of this skill's scope.** The former
> `generate-combined-plot-1.png` (`fig:posterior_comparison`, "distribution of PC
> values") has been **removed from the paper**. The current Results figures are the
> GaussNiche §2.3 sampler-comparison boxplots — `sampler-comparison-boxplots.pdf`
> (2-D, `fig:sampler_comparison`) and `sampler-comparison-boxplots-5d.pdf`
> (5-D, `fig:sampler_comparison_5d`) — produced by `run_2d_experiment.R` /
> `run_5d_experiment.R` and copied into `graphics/` directly, NOT by this skill.
> The per-figure recipe sections further below that still reference `.png`
> filenames or `generate-combined-plot` are **historical**; for the appendix
> diagnostics the USE.MCMC-side skill is authoritative.

If you decide to switch any figure to a vector PDF (recommended, see Shared
styling), you MUST also update the matching `\includegraphics{...}` extension in
`text/results.tex` / `text/appendix/appendix.tex` **and** the paper's `Makefile`.
That is a manuscript edit in the sibling repo — get the user's go-ahead first.

## Before you touch the figures — read these blockers

1. **No GaussNiche script retains MCMC chains.** Both `pa_mcmc` wrappers
   (`virtualSpecies_fn.R:241-255`, `virtualSpecies_nd_fn.R:88-99,184-190`) call
   `USE.MCMC::paSamplingMcmc(..., plot_proc = FALSE)`, which post-processes and
   `rbind`s chains into one pooled point cloud, discarding per-chain identity and
   the raw trace. A repo-wide grep for `coda` / `gelman` / `effectiveSize` /
   `traceplot` / `autocorr` / `as.mcmc` / `mcmc.list` returns **zero** matches in
   any GaussNiche driver or fn. **Therefore none of the three appendix
   diagnostics can be produced from a GaussNiche result object today.** You must
   build a `coda::mcmc.list` yourself with the lower-level
   `USE.MCMC::mcmcSampling()` (see "Chain builder").

2. **`precomputeMcmcEnvironment()` exists, but it caches the *environment*, not
   the *chains*.** As of USE.MCMC 0.0.4 it is exported and `paSamplingMcmc()` takes
   a real `precomputed.env` argument (`paSamplingMcmc.R:35`): it precomputes the
   PCA + environmental GMM + proposal covariance + distance threshold so many
   species sharing one environment can skip that setup. The GaussNiche `pa_mcmc`
   wrappers' `precomputed.env =` pass-through (`virtualSpecies_fn.R:252`,
   `virtualSpecies_nd_fn.R:98`) is therefore **correct** against the current
   package — it was only broken on older 0.0.2 checkouts, not now. But that cache
   is orthogonal to what you need here: `paSamplingMcmc()` still `rbind`s its
   chains into one pooled cloud (`paSamplingMcmc.R:191`) and dedups (`:213`),
   discarding per-chain identity and the raw trace (Blocker 1). So precompute does
   **not** hand you a `coda.chain.list` — for the appendix diagnostics you still
   build chains with `mcmcSampling()` (see "Chain builder"). (Other precompute
   hooks: the `precomputed.pca` arg and the internal `precompute_gmm_params()`.)

3. **Caption vs source dimensionality/length mismatches.** The current vignette
   convergence section runs **3D (PC1,PC2,PC3), 4 chains, 5000 steps**. But the
   paper captions say: autocorrelation = **50000 steps, 2D**; Gelman-Rubin = **3D**
   (matches); trace = **4 chains, 50000 steps, 2D**. So autocorrelation and trace
   need a regenerated **2D, 50000-step** ensemble to match the captions — OR the
   captions need correcting. **Do not silently pick one.** Resolve with the user
   which is the intended final setting, then generate accordingly.

4. **`trace_plots.png` has no scripted source at all.** It is hand-curated,
   outside both the Makefile and the vignette; there is no `traceplot`/
   `coda::traceplot` call anywhere in USE.MCMC. Any reproducible version is
   authored from scratch and will differ stylistically from the committed static
   PNG. Flag this to the user before overwriting it.

5. **Heavy runs go through `sbatch/`.** Building a 50000-step × 4-chain ensemble is
   not a login-node task. Use the existing `sbatch/` submission pattern (per the
   root CLAUDE.md and the `euler-rstudio-server` skill); on a plain workstation a
   direct `Rscript` is fine but expect minutes. `mcmcSampling()` auto-dispatches a
   C++ loop when the built-in density/proposal factories are used.

6. **Font registration.** The shared theme targets `base_family = "Lato"` (the
   journal body font). Lato must be installed and registered with the active
   graphics device (e.g. `systemfonts`/`showtext` for `cairo_pdf`) or text falls
   back silently. Confirm Lato is available on the machine; otherwise fall back to
   Helvetica/Arial and note the substitution to the user.

---

## Shared styling — ONE source of truth

Put this in a small helper at the GaussNiche repo root, e.g.
`publication_theme.R`, and `source()` it from every figure recipe below. It is
grounded in the WileyNJDv5 two-column (LATO2COL) journal spec.

```r
## publication_theme.R --- shared look for the markov-chain-sampler paper figures
suppressPackageStartupMessages({
  library(ggplot2)
})

## --- Journal geometry (inches), from WileyNJDv5.cls text-width probe ---
PUB_WIDTH_COL  <- 3.43   # single column (~87 mm; two-column LATO2COL, derived)
PUB_WIDTH_FULL <- 7.00   # full text width / figure* (177.8 mm, WileyNJDv5.cls:1076)
PUB_HEIGHT_MAX <- 9.00   # text block height; never exceed
PUB_DPI        <- 600    # default raster DPI if PDF is not used (>=300 min)

## --- Base font: Lato if registered, else a safe fallback (see Blocker 6) ---
pub_family <- tryCatch({
  fams <- if (requireNamespace("systemfonts", quietly = TRUE))
    systemfonts::system_fonts()$family else character(0)
  if ("Lato" %in% fams) "Lato" else "Helvetica"   # default fallback; note it
}, error = function(e) "Helvetica")               # default fallback

## --- Colourblind-safe + greyscale-robust categorical palette (Okabe-Ito) ---
## ONE named source of truth for the point-type / sampler series. Reconciles the
## scattered GaussNiche literals (Background grey50, Presences firebrick,
## Pseudo-absences darkorange) and the vignette's hardcoded green/orange/red/black.
PUB_COLORS <- c(
  "Virtual presence" = "#000000",  # black
  "USE-uniform"      = "#E69F00",  # orange
  "Random-geo"       = "#56B4E9",  # sky blue
  "MCMC-2d"          = "#009E73",  # bluish green
  "MCMC-3d"          = "#0072B2",  # blue
  "Background"       = "#999999",  # grey  (KDE / available E-space)
  "Presences"        = "#D55E00",  # vermillion
  "Pseudo-absences"  = "#CC79A7"   # reddish purple
)
## Redundant linetype encoding so series survive greyscale.
PUB_LINETYPES <- c(
  "Virtual presence" = "solid", "USE-uniform" = "longdash",
  "Random-geo" = "dotted", "MCMC-2d" = "dashed", "MCMC-3d" = "dotdash",
  "Background" = "solid", "Presences" = "dashed", "Pseudo-absences" = "dotted"
)

## --- The theme. base_size 8 for single-column/appendix; 9 for the full-width
## --- Results figure. Author at FINAL width so ggsave does not shrink fonts. ---
pub_theme <- function(base_size = 8, base_family = pub_family) {
  theme_minimal(base_size = base_size, base_family = base_family) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.25),  # >= ~0.25 pt
      axis.line        = element_line(linewidth = 0.3),
      legend.position  = "bottom",
      legend.title     = element_blank(),
      legend.key.size  = unit(0.8, "lines"),
      plot.title       = element_blank(),    # titles live in LaTeX \caption{}
      plot.subtitle    = element_blank(),
      strip.text       = element_text(face = "italic")  # lowercase (a)(b)(c)
    )
}

pub_color_scale    <- function(...) scale_colour_manual(values = PUB_COLORS, ...)
pub_fill_scale     <- function(...) scale_fill_manual(values = PUB_COLORS, ...)
pub_linetype_scale <- function(...) scale_linetype_manual(values = PUB_LINETYPES, ...)

## --- Standard saver. PREFER vector PDF (XeLaTeX embeds it; DPI irrelevant). ---
## Pass full=TRUE for the figure* Results plot, FALSE for single-column appendix.
pub_save <- function(plot, file, full = FALSE, height = NULL,
                     device = c("pdf", "png")) {
  device <- match.arg(device)
  width  <- if (full) PUB_WIDTH_FULL else PUB_WIDTH_COL
  if (is.null(height)) height <- if (full) 8.5 else 2.6  # defaults; set per fig
  if (device == "pdf") {
    ggplot2::ggsave(file, plot, device = grDevices::cairo_pdf,
                    width = width, height = height, units = "in")
  } else {
    ggplot2::ggsave(file, plot, device = "png",
                    width = width, height = height, units = "in", dpi = PUB_DPI)
  }
}
```

Styling rules every figure obeys:

- **Title/caption belong in LaTeX `\caption{}`**, never baked into the image.
  In particular, remove the cowplot `draw_label("Comparison to other PA sampling
  methods")` title from the Results chunk.
- **Author at the final printed width**; do not draw a big canvas and let `ggsave`
  shrink it (that shrinks the fonts below the 6 pt floor).
- **Data lines `linewidth` ~0.6–0.8**; axis/grid lines `>= 0.3`.
- **One collected legend** (`plot_layout(guides = "collect")` or
  `cowplot::get_legend`), not a legend on a single panel.
- **Vector PDF is the biggest quality win** over the current ~96-dpi PNGs. If
  raster is mandated, 600-dpi PNG/TIFF-LZW; never JPEG for plots.
- **Colour space RGB** (online/PDF-first); the typesetter converts to CMYK if
  required.

---

## Figure 1 — Results: "Distribution of the PC values using different sampling methods"

- **Filename (REQUIRED):** `generate-combined-plot-1.png` in
  `../markov-chain-sampler-paper/graphics/`.
- **Paper ref / caption:** `text/results.tex:2-7`, `\label{fig:posterior_comparison}`,
  `figure*` full-width (`width=1\linewidth`). Prose at `text/results.tex:8`.
- **What it must show:** a 2-column composite. **Left column = MCMC sampled in
  2D, right column = MCMC sampled in 3D.** Each column is three stacked per-PC-axis
  (PC1/PC2/PC3) kernel-density panels comparing the distribution of PC values from
  each sampling method — **MCMC, USE-grid `uniform`, random-geographical** —
  against the **virtual-presence** distribution. (Panel semantics come from the
  vignette, not the paper prose — see Open question 5; confirm with the user.)

### Producer (authoritative)

There is **no exact GaussNiche producer**. The canonical source is the USE.MCMC
vignette
`../USE.MCMC/vignettes/insights-on-MCMC-pseudo-absence-sampling-vignette.Rmd`:

- chunk `compare-to-uniform-models-2d` builds `plot.2d.comparisson`,
- chunk `compare-to-uniform-models-3d` builds `plot.3d.comparisson`,
- chunk `generate-combined-plot` `cowplot::plot_grid()`s them and
  adds the title → `figure-html/generate-combined-plot-1.png`, copied by
  `make figures`.

Data inputs (built earlier in the vignette): `env.with.pc.sf` (rastPCA scores),
`virtual.presence.points.pc` (virtual species presences),
`mapped.sampled.points.selected` (MCMC PAs, 2D),
`sampled.points.mcmc.higher.dim` (MCMC PAs, 3D),
`USE.MCMC::paSampling(grid.res = 10)` extracted to PC space (the USE-grid
`uniform`), and a uniform geographic subsample of `env.with.pc.sf` (random-geo).
Each must carry PC1/PC2/PC3 columns.

> The closest GaussNiche analogue is `unified_bias()` (`6_compare_5d.R:91-119`,
> per-PC `geom_density` overlay, `wrap_plots(ncol = 3)`) and
> `all_species_response()` (`6_compare_5d.R:141-162`). **They do not reproduce
> this figure** — they compare GaussNiche's own `random`/`mcmc` samplers against
> background + presences, omit the USE-grid `uniform` sampler, and have no
> 2D-vs-3D MCMC contrast. Reproducing Fig 1 inside GaussNiche would require adding
> (a) a USE-grid `uniform` run, (b) a paired MCMC-2d vs MCMC-3d run on the same
> virtual species, and (c) a generalised multi-point-type density-overlay — none
> exist today. **Default recommendation: keep this figure vignette-driven** (it
> preserves the chunk-name = filename coupling).

### Publication gap

`generate-combined-plot-1.png` is the worst figure in the paper (~96 dpi,
672×864 px, rendered at full width 7.0 in → soft/pixelated). Other gaps: raw
default `theme()`, asymmetric legend (only under PC3), non-colourblind-safe
hardcoded green/orange/red/black (the chunk's `scale_color_manual`), no Lato/8–9 pt control, a
baked-in `draw_label` title duplicating the LaTeX caption, bare PC1/PC2/PC3 axis
labels.

### Recipe (vignette-driven, recommended)

In the vignette chunks, with `source("publication_theme.R")` available:

1. Replace the hardcoded `scale_color_manual(green/orange/red/black)` with the
   shared scales: `pub_color_scale()` + `pub_linetype_scale()` keyed on
   `c("Virtual presence","USE-uniform","Random-geo","MCMC-2d","MCMC-3d")` so the
   series survive greyscale and red-green CB.
2. `theme_set(pub_theme(base_size = 9))` (full-width figure → base_size 9).
3. **Delete the cowplot `draw_label` title** in the `generate-combined-plot` chunk
   so the title lives only in `\caption{}`.
4. Build a single shared legend via `cowplot::get_legend` + `plot_grid` rel
   layout (not just on PC3).
5. Export at full width. Either set the chunk options to
   `dev = "png", dpi = 600, fig.width = 7.0, fig.height = 8.5` so
   `make figures` copies a 600-dpi PNG **named `generate-combined-plot-1.png`**;
   OR (preferred) save vector:
   ```r
   ggsave("generate-combined-plot-1.pdf", combined,
          device = cairo_pdf, width = 7.0, height = 8.5, units = "in")
   ```
   and switch the `\includegraphics` extension in `text/results.tex` + the paper
   Makefile to `.pdf` (manuscript edit — confirm with user).

`height = 8.5` in (≤ `PUB_HEIGHT_MAX`) is a **default** sized to a 3-row ×
2-column stack; adjust to the native aspect ratio so LaTeX does no rescaling.

> If the user insists on a GaussNiche-native reproduction: generalise
> `unified_bias()` to take an arbitrary named list of point-type data.frames
> (presence, USE-uniform, random-geo, MCMC-2d, MCMC-3d), emit a 2-col (2D|3D) ×
> 3-row (PC1/PC2/PC3) `wrap_plots`, apply `pub_theme(9)` + the shared scales, and
> `pub_save(p, ".../graphics/generate-combined-plot-1.png", full = TRUE,
> height = 8.5, device = "png")`. This additionally requires adding the USE-grid
> `uniform` run and the 2D-vs-3D MCMC pair to a driver. Verify any new USE.MCMC
> call signatures against the package before relying on them.

---

## Chain builder — prerequisite for all three appendix figures

All three appendix diagnostics consume one object: `coda.chain.list`, a
`coda::mcmc.list` of N independent chains. **No GaussNiche entry point yields a
chain** (Blocker 1), so build it with `USE.MCMC::mcmcSampling()` directly,
exactly as the vignette `highe-dim-convergence` chunk does
(chunk `highe-dim-convergence`). `mcmcSampling()` returns a `data.frame` with one column per
dimension plus a `density` column; one call = one chain.

`mcmcSampling()` signature (verify against USE.MCMC `R/mcmcSampling.R:18-22`):

```r
mcmcSampling(dataset, dimensions, densityFunction, proposalFunction,
             n.sample.points, burnIn = 1000, verbose = TRUE,
             covariance.correction = 1, max.burnin.cycles = 50,
             engine = c("auto", "R", "cpp"))
```

Standalone builder (the seed is layered per the GaussNiche convention; ship this
as a small script and run it via `sbatch/` for the 50000-step ensemble):

```r
## build_chains.R --- reproduce coda.chain.list for the appendix diagnostics.
## Mirrors USE.MCMC vignette chunk `highe-dim-convergence`.
suppressPackageStartupMessages({
  library(USE.MCMC); library(coda); library(parallel)
  library(terra); library(sf); library(mclust)
})
set.seed(42)  # vignette global seed

## --- env.rast MUST be the ORIGINAL rasters (root CLAUDE.md contract) ---
env.data.raster <- USE.MCMC::Worldclim_tmp |> terra::rast(type = "xyz") |> round(2)
rpc <- rastPCA(env.data.raster, stand = TRUE)
env.data.raster.with.pc <- c(env.data.raster, rpc$PCs)
env.data.sf  <- env.data.raster |> as.data.frame(xy = TRUE) |>
  sf::st_as_sf(coords = c("x", "y"))
env.with.pc.sf <- rpc$PCs |> as.data.frame(xy = TRUE) |> na.omit() |>
  sf::st_as_sf(coords = c("x", "y")) |> sf::st_join(env.data.sf)

vp <- getVirtualSpeciesPresencePoints(env.data = env.data.raster.with.pc,
                                      n.samples = 300)
virtual.presence.points.pc <- terra::extract(env.data.raster.with.pc,
                                             vp$sample.points, bind = TRUE) |>
  sf::st_as_sf()

## --- DIMENSIONS / LENGTH: set per Blocker 3 (resolve with the user). ---
## To MATCH the autocorrelation+trace captions: dimensions <- c("PC1","PC2"); n.steps <- 50000
## Vignette/Gelman-as-shipped:                 dimensions <- c("PC1","PC2","PC3"); n.steps <- 5000
dimensions <- c("PC1", "PC2", "PC3")   # CHANGE per the agreed caption
n.steps    <- 5000                     # CHANGE per the agreed caption (e.g. 50000)
n.chains   <- 4

env.sub <- env.with.pc.sf[
  stats::runif(min(nrow(env.with.pc.sf), 2000), 1, nrow(env.with.pc.sf)), ] |>
  sf::st_drop_geometry()
env.model <- mclust::densityMclust(env.sub[dimensions])
env.dens  <- mclust::predict.densityMclust(env.model,
                sf::st_drop_geometry(env.with.pc.sf[dimensions]))
env.threshold <- stats::quantile(env.dens, 0.04)

species.model <- mclust::densityMclust(
  sf::st_drop_geometry(virtual.presence.points.pc[dimensions]))
species.cutoff.threshold <- stats::quantile(species.model$density, 0.9)

cov.mat <- 0.075 * stats::cov(sf::st_drop_geometry(env.with.pc.sf)[dimensions])
proposalFunction <- addHighDimGaussian(cov.mat = cov.mat,
                                       dim = length(dimensions))
densityFunction  <- mclustDensityFunction(
  env.model = env.model, species.model = species.model, dim = dimensions,
  threshold = env.threshold, species.cutoff.threshold = species.cutoff.threshold)

chain.list <- mclapply(seq_len(n.chains), function(i) {
  utils::capture.output({
    s <- mcmcSampling(dataset = env.sub, dimensions = dimensions,
                      n.sample.points = n.steps,
                      densityFunction = densityFunction,
                      proposalFunction = proposalFunction,
                      burnIn = 1000, covariance.correction = 50, verbose = TRUE)
  })
  coda::as.mcmc(s[dimensions])
}, mc.cores = 1)                 # bump mc.cores under an sbatch allocation
coda.chain.list <- coda::mcmc.list(chain.list)
saveRDS(coda.chain.list, "coda_chain_list.rds")
```

Notes: `covariance.correction = 50` and `cov` scaling `0.075` are the vignette
values; `burnIn = 1000` uses Robbins-Monro adaptation toward target acceptance
0.234. `engine = "auto"` selects the C++ loop because the built-in factories
carry the `rcpp_spec` attribute. Names `addHighDimGaussian`,
`mclustDensityFunction`, `getVirtualSpeciesPresencePoints` are confirmed USE.MCMC
exports; if you adapt arguments, **verify against USE.MCMC**.

---

## Figure 2 — Appendix: autocorrelation

- **Filename (REQUIRED):** `autocorreltation-plot-1.png` *(keep the misspelling)*
  in `../markov-chain-sampler-paper/graphics/`.
- **Paper ref / caption:** `text/appendix/appendix.tex:5-10`,
  `\label{fig:autocorrelation_analysis}`, `width=0.9\linewidth`. Caption: 50000
  steps, **two** dimensions, ACF flat after ~30 lags (justifies thinning = 30).
- **What it must show:** ACF of **one** MCMC chain, one panel per dimension,
  flattening by ~lag 30.

### Producer

Vignette chunk `autocorreltation-plot`:
```r
autocorr.plot(coda.chain.list[1], lag.max = 50)   # FIRST chain only; coda
```
→ `figure-html/autocorreltation-plot-1.png` → `graphics/` via `make figures`.

> Caption mismatch (Blocker 3): the vignette runs **3D at 5000 steps**, the
> caption says **2D at 50000**. Build `coda.chain.list` with
> `dimensions = c("PC1","PC2")`, `n.steps = 50000` to match the caption — or have
> the user amend the caption — before rendering.

### Publication gap

672×480 px, ~96 dpi embedded. Base-R `coda` graphics: no ggplot theme, no Lato,
default black lines, off-house axis labels.

### Recipe (preferred: re-implement ACF in ggplot so the shared theme applies)

```r
source("publication_theme.R"); library(coda)
ccl <- readRDS("coda_chain_list.rds")
ac  <- coda::autocorr(ccl[1], lags = 0:50)         # first chain
## ac is [lag, var, var]; take the diagonal (per-dimension ACF), melt to long
dims <- coda::varnames(ccl)
acdf <- do.call(rbind, lapply(seq_along(dims), function(j)
  data.frame(lag = 0:50, acf = ac[, j, j], dimension = dims[j])))
n   <- coda::niter(ccl[1])
bnd <- 1.96 / sqrt(n)
p_acf <- ggplot(acdf, aes(lag, acf)) +
  geom_hline(yintercept = c(-bnd, bnd), linetype = "dashed", linewidth = 0.3) +
  geom_segment(aes(xend = lag, yend = 0), linewidth = 0.6) +
  facet_wrap(~ dimension, ncol = 1) +
  labs(x = "Lag", y = "Autocorrelation") +
  pub_theme(base_size = 8) + theme(legend.position = "none")
pub_save(p_acf, "../markov-chain-sampler-paper/graphics/autocorreltation-plot-1.png",
         full = FALSE, height = 2.6, device = "png")  # keep PNG + the typo name
```

Alternative if you must keep base `coda` graphics, wrap them at 600 dpi:
```r
png("../markov-chain-sampler-paper/graphics/autocorreltation-plot-1.png",
    width = 3.43, height = 2.6, units = "in", res = 600, family = pub_family)
coda::autocorr.plot(ccl[1], lag.max = 50); dev.off()
```
The appendix uses `width=0.9\linewidth`; authoring at single-column 3.43 in and
changing the LaTeX to `width=\linewidth` is the cleanest (manuscript edit —
confirm). `height = 2.6` in is a **default**.

---

## Figure 3 — Appendix: Gelman-Rubin

- **Filename (REQUIRED):** `gelman-plot-1.png` in
  `../markov-chain-sampler-paper/graphics/`.
- **Paper ref / caption:** `text/appendix/appendix.tex:13-18`,
  `\label{fig:gelman-rubin-appendix}`, `width=0.9\linewidth`. Caption: "chain run
  in **three** dimensions" (this one **matches** the vignette). Prose
  (`appendix.tex:12`): good convergence.
- **What it must show:** the Gelman-Rubin shrink-factor (PSRF) evolution for the
  multi-chain ensemble (4 chains), one panel per dimension.

### Producer

Vignette chunk `gelman-plot` — the whole chunk:
```r
coda::gelman.plot(coda.chain.list)
```
→ `figure-html/gelman-plot-1.png` → `graphics/`. Needs `>= 2` chains; build the
4-chain ensemble in **3D** (caption-consistent) with the chain builder.

### Publication gap

672×480 px, ~96 dpi. Base-R `coda` graphics, default colours (median + 97.5%
quantile), no Lato/8 pt, not CB-tuned.

### Recipe (preferred: ggplot re-implementation)

```r
source("publication_theme.R"); library(coda)
ccl <- readRDS("coda_chain_list.rds")   # 4 chains, 3D
gp  <- coda::gelman.plot(ccl, autoburnin = FALSE)  # returns $shrink + $last.iter
dims <- dimnames(gp$shrink)[[2]]
grdf <- do.call(rbind, lapply(seq_along(dims), function(j)
  data.frame(iteration = gp$last.iter,
             median    = gp$shrink[, j, 1],
             q975      = gp$shrink[, j, 2],
             dimension = dims[j])))
p_gr <- ggplot(grdf, aes(iteration)) +
  geom_line(aes(y = median),  linewidth = 0.6) +
  geom_line(aes(y = q975), linetype = "dashed", linewidth = 0.6) +
  geom_hline(yintercept = 1.1, linetype = "dotted", linewidth = 0.3) +
  facet_wrap(~ dimension, ncol = 1, scales = "free_y") +
  labs(x = "Last iteration in window", y = "Shrink factor") +
  pub_theme(base_size = 8) + theme(legend.position = "none")
pub_save(p_gr, "../markov-chain-sampler-paper/graphics/gelman-plot-1.png",
         full = FALSE, height = 2.6, device = "png")
```
`gelman.plot()` produces the base plot as a side effect AND returns the
shrink-factor data used above; **verify the returned structure against your coda
version**. Base-R fallback (600 dpi) is the same `png(...); coda::gelman.plot(ccl);
dev.off()` wrapper as Figure 2. `height = 2.6` in is a **default**.

---

## Figure 4 — Appendix: trace + posterior

- **Filename (REQUIRED):** `trace_plots.png` in
  `../markov-chain-sampler-paper/graphics/`.
- **Paper ref / caption:** `text/appendix/appendix.tex:19-24`,
  `\label{fig:traceplot}`, `width=0.9\linewidth`. Caption: "trace and posteriors
  of **4 chains** run for **50000 steps** in **two** dimensions", good mixing.
- **What it must show:** trace plots of the 4 chains alongside their marginal
  posterior densities.

### Producer — there is NONE (Blocker 4)

`trace_plots.png` is **static / hand-curated**: not in the Makefile, not produced
by any chunk, and there is no `traceplot`/`coda::traceplot` call anywhere in
USE.MCMC. This is the largest reproducibility gap of the four. A regenerated
version is authored from scratch over a freshly built `coda.chain.list` and will
differ stylistically from the committed PNG — **tell the user before overwriting**.

The caption says **2D, 50000 steps, 4 chains**, but the vignette ensemble is 3D
at 5000 (Blocker 3). Build a dedicated **2D, 50000-step, 4-chain** ensemble with
the chain builder (`dimensions = c("PC1","PC2")`, `n.steps = 50000`).

### Recipe (preferred: ggplot trace | posterior, composed)

```r
source("publication_theme.R")
library(coda); library(patchwork)
ccl <- readRDS("coda_chain_list.rds")    # 2D, 4 chains, 50000 steps
dims <- coda::varnames(ccl)
## long form: one row per (chain, iteration, dimension)
long <- do.call(rbind, lapply(seq_along(ccl), function(ci) {
  m <- as.matrix(ccl[[ci]])
  do.call(rbind, lapply(dims, function(d)
    data.frame(chain = factor(ci), iteration = seq_len(nrow(m)),
               value = m[, d], dimension = d)))
}))
chain_cols <- setNames(unname(PUB_COLORS[c("USE-uniform","Random-geo",
                                           "MCMC-2d","MCMC-3d")]), levels(long$chain))
p_trace <- ggplot(long, aes(iteration, value, colour = chain)) +
  geom_line(linewidth = 0.3, alpha = 0.8) +
  scale_colour_manual(values = chain_cols) +
  facet_wrap(~ dimension, ncol = 1, scales = "free_y") +
  labs(x = "Iteration", y = "Value") + pub_theme(base_size = 8)
p_post <- ggplot(long, aes(value, colour = chain)) +
  geom_density(linewidth = 0.6) +
  scale_colour_manual(values = chain_cols) +
  facet_wrap(~ dimension, ncol = 1, scales = "free") +
  labs(x = "Value", y = "Density") + pub_theme(base_size = 8)
fig <- (p_trace | p_post) + plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
pub_save(fig, "../markov-chain-sampler-paper/graphics/trace_plots.png",
         full = FALSE, height = 3.0, device = "png")  # keep the .png name
```

Base-R alternative (the canonical coda call, consistent with the package's coda
dependency) — note it has no scripted precedent and is purely a regeneration:
```r
png("../markov-chain-sampler-paper/graphics/trace_plots.png",
    width = 3.43, height = 3.0, units = "in", res = 600, family = pub_family)
coda::traceplot(ccl); coda::densplot(ccl); dev.off()
```
For the "posterior" half you may instead reuse the vignette `mcmc-posterior`
style (chunk `mcmc-posterior`): per-PC `geom_density` of sampled
pseudo-absences vs environment vs presence — descriptive posterior-marginal
checks, not coda diagnostics. `height = 3.0` in is a **default**.

If you make this a real scripted figure (recommended — it leaves the "static"
category), add it to the paper's Makefile so `make figures` regenerates it, and
confirm the final 4 chains / 50000 steps / 2D settings against the caption.

---

## End-to-end checklist

1. Resolve Blocker 3 with the user (autocorr/trace = 2D 50000? or amend
   captions?). Confirm Lato availability (Blocker 6).
2. Confirm whether Fig 1 stays vignette-driven (default) or is reproduced in
   GaussNiche.
3. `source("publication_theme.R")`. Build `coda.chain.list` once per required
   (dimensions, n.steps) setting via the chain builder — run heavy ensembles
   through `sbatch/`.
4. Render each figure to its **exact** required filename in
   `../markov-chain-sampler-paper/graphics/`. Do NOT rename; keep the
   `autocorreltation` typo and the `trace_plots.png` name.
5. If you switch any figure to `.pdf`, edit the matching `\includegraphics`
   extension in `text/results.tex` / `text/appendix/appendix.tex` AND the paper
   Makefile (sibling-repo manuscript edits — get user sign-off).
6. Rebuild the manuscript (`latexmk -xelatex main.tex` in the sibling repo) to
   confirm the figures resolve and are not rescaled.
```
