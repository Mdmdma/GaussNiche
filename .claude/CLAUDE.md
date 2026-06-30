# GaussNiche — agent guide

This file orients an LLM agent working on GaussNiche. It is **portable**: it
intentionally avoids any cluster-, host-, or user-specific instructions so it
can be shared, vendored, or read by anyone cloning the repo. Anything
machine-specific (Euler HPC, RStudio Server, R-package install pitfalls)
lives in `.claude/skills/` and is loaded by the agent on demand.

If you are on a non-trivial host (cluster login node, RStudio Server,
containerised R) consult the relevant skill before running anything:

- `.claude/skills/euler-rstudio-server/SKILL.md` — running inside an Euler
  JupyterHub-launched RStudio Server session (compute node + apptainer
  container). Read this first if `uname -n` returns anything other than a
  plain workstation hostname.
- `.claude/skills/euler-r-spack-setup/SKILL.md` — installing/repairing R
  packages with native deps (sf, terra, hypervolume, USE.MCMC, devtools)
  inside the Euler rocker/rstudio container.

If neither applies (plain Linux/macOS workstation), the rest of this file
is all you need.

Three **task** skills are host-independent (load on demand, not gated by host):

- `.claude/skills/publication-figures/SKILL.md` — authoritative guide to ALL
  figures in the markov-chain-sampler paper and its downstream-HSM report: Part A
  the appendix MCMC diagnostics (autocorrelation / Gelman-Rubin / trace, from the
  USE.MCMC vignette) and Part B the GaussNiche-produced figures (§2.3
  sampler-comparison boxplots, PC-matrix / geo-grid, hsm_full_report, the
  cross-species aggregates, the Dunn summary, the uniform+ ablation heatmaps),
  with the figure→producer→destination map + the copy-by-hand-into-`graphics/`
  convention.
- `.claude/skills/hsm-ablation/SKILL.md` — run / recompute / extend / interpret
  the uniform+ (pa_mcmc) settings ablation (environmental.cutof.percentile ×
  species.cutoff.threshold sweep vs the fixed RND/buffer baseline). Read it when
  the environment changes and the ablation must be recomputed, or to tune uniform+.
- `.claude/skills/adversarial-paper-review/SKILL.md` — heavy, non-routine
  multi-persona ADVERSARIAL review of the `../markov-chain-sampler-paper`
  manuscript. Fans out reader-base personas (ecologist / mathematician / student /
  interested public / computer-scientist) + specialists (reproducibility /
  implementation / abstract-vs-open-questions / Reviewer-2 skeptic) that review the
  RENDERED submission PDF first (render `main.pdf` → page PNGs — the login node has
  no rasterizer), adversarially refute every finding, then synthesize an editor
  report. Runs via the Workflow tool (`adversarial_review.workflow.js`). Read it to
  critique the paper from its target audiences before submission.

---

## What this project is

GaussNiche is an R framework for **simulating virtual species in
multivariate environmental space** using **bivariate Gaussian** suitability
functions fitted on PCA axes (PC1, PC2 of a stack of bioclimatic rasters),
plus diagnostics for the quality of pseudo-absence sampling strategies.

The full pipeline:

1. Stack environmental rasters → PCA → keep (PC1, PC2) as the analysis
   E-space. PC scores per non-NA cell form the **background**.
2. Define a niche = bivariate Gaussian with optimum **μ** and
   covariance **Σ = [[σ₁², ρ σ₁ σ₂], [ρ σ₁ σ₂, σ₂²]]** in PC space.
   Suitability is the Gaussian density normalised to peak = 1.0 at μ.
3. **Bernoulli draw** per background cell with probability = suitability →
   binary presence/absence. Repeated for `n_realizations` independent draws
   to quantify stochastic variation.
4. **Pseudo-absence sampling** via a swappable sampler interface
   (`pa_random`, `pa_uniform`, `pa_mcmc`). All samplers share the same
   call signature so they can be added/swapped without touching the main
   function.
5. **Diagnostics per realisation**: hypervolume intersection between
   presence and pseudo-absence point sets (class overlap), pseudo-absence
   PC range / background PC range (sampling bias), proportion of
   pseudo-absences that hit true-absence cells.
6. **Back-projection** of the suitability surface and the reference
   realisation to geographic rasters.

Why hypervolumes use a single, **pre-computed background bandwidth**:
estimate_bandwidth() depends on N and on the spread of the input — running
it per-realisation makes hypervolumes incomparable across species and
samplers. The single bandwidth is computed once from the full background
and fixed everywhere.

## Repository map

```
GaussNiche/
├── 1_developing_framework.R       linear, step-by-step development of the
│                                  pipeline (PCA → niche → Bernoulli → PA →
│                                  diagnostics → back-projection). Reference
│                                  for the methodology before abstraction.
├── virtualSpecies_fn.R            modular wrappers + samplers. Source this
│                                  from any analysis script. Public exports:
│                                    pa_random()       sampler (uniform from bg)
│                                    pa_uniform()      sampler (USE::paSampling)
│                                    pa_mcmc()         sampler (USE.MCMC::paSamplingMcmc)
│                                    compute_bandwidth() one-time bg bandwidth
│                                    sampler_effective_settings() resolves the
│                                                      effective tunables (call
│                                                      override vs sampler
│                                                      default) per sampler
│                                    summarize_parameters() one-page replication
│                                                      record (niche / sampler /
│                                                      run controls / env tables)
│                                                      built from species
│                                                      results' $parameters slot
│                                    virtualSpecies()  main pipeline; returns
│                                                      list(niche, parameters,
│                                                      suit_rast, pa_rast,
│                                                      response_curves,
│                                                      samplers, plots, bw,
│                                                      background)
├── 2_testing_wrapper_function.R   reproducible example exercising
│                                  virtualSpecies() with all three samplers
│                                  on USE.MCMC::Worldclim_tmp.
├── 3_ virtualSpecies_cases.R      additional case-study runs (note the
│                                  space in the filename — keep it; quote
│                                  it when sourcing).
├── 4_multi_species_comparison.R   2-D multi-species comparison driver (CLI
│                                  cutoff sweep); base for run_2d_experiment.R.
├── run_2d_experiment.R            paper §2.3 experiment: the 4 Table-1 species
│                                  in 2-D (RND/buffer-out/USE/uniform+) over R=50
│                                  Bernoulli realisations; saves results + the
│                                  demanded boxplots to <scratch>/.../results2d.
│                                  Run via sbatch/submit_2d_experiment.sh.
├── run_5d_experiment.R            5-D arm of the §2.3 experiment (env5d
│                                  'lean_natural', 5 PCs; samplers RND/buffer-out/
│                                  uniform+ — USE omitted, it is 2-D-only). Via
│                                  sbatch/submit_5d_experiment.sh → results5d_experiment.
├── run_5d_hsm.R                   DOWNSTREAM-HSM arm (appendix): same 4 species /
│                                  3 samplers / R=50 / IDENTICAL seeds as
│                                  run_5d_experiment.R, but fits the 5 USE model
│                                  families (GLM/GAM/RF/BRT/Maxent, Java-free) to
│                                  each realisation on TWO predictor sets (5 PCs +
│                                  12 raw layers) via the sdm_hook, scoring
│                                  internal (held-out) AND truth-based accuracy.
│                                  Writes hsm_metrics/dunn/violins to
│                                  results5d_hsm (no 90 MB result rds). Via
│                                  sbatch/submit_5d_hsm.sh (deps:
│                                  sbatch/submit_install_hsm_deps.sh first).
├── tune_5d_hsm.R                  uniform+ settings ABLATION sweep (appendix): a
│                                  grid of pa_mcmc knobs (environmental.cutof.percentile
│                                  x species.cutoff.threshold) scoring downstream HSM
│                                  metrics vs the fixed RND/buffer baseline. Speedups
│                                  A1-A4 baked in (uniform+ only / compute_hypervolume=
│                                  FALSE / pc5 only / 25 reals); sets the furrr plan
│                                  ONCE (manage_plan=FALSE) across all cell×species
│                                  calls (committed grid: 40 cells, env→0.25 × 5
│                                  species_cutoff incl =1 uniform → 160 calls). Emits
│                                  per-cell median table + value/Δ-vs-RND heatmaps to
│                                  results5d_tune. Via sbatch/submit_tune_5d_hsm.sh.
├── bench_hsm.R                    single-threaded profiler decomposing one
│                                  realisation's cost (sampler / hypervolume / SDM
│                                  hook / background-predict) — measured that the 5-D
│                                  hypervolume_gaussian is ~70% of a task. Via
│                                  sbatch/submit_bench_hsm.sh.
├── profile_samplers.R             single-threaded profiler decomposing the cost of
│                                  uniform (USE paSampling) vs uniform+ (paSamplingMcmc)
│                                  into MODEL-FITTING vs SAMPLING/CHAIN, in 2-D and 5-D.
│                                  Measured: the MCMC chain is ~0.007 s/10k steps (2-D),
│                                  the GMM EM fits dominate (~99% of a call); USE's
│                                  optimRes grid search ~380 s. Persists results/profile/
│                                  profile_samplers.{csv,rds}. Via
│                                  sbatch/submit_profile_samplers.sh.
├── experiment_plots.R             dimension-agnostic metric boxplots (overlap /
│                                  per-axis range coverage / prop true-absence).
│                                  Consumes virtualSpecies() AND
│                                  virtualSpecies_nd() results; discovers
│                                  samplers + rel_cov_* axes at runtime.
├── hsm_eval.R                     downstream HSM library: fit_hsms / per-algo
│                                  predict + metrics / boyce_index (CBI reimpl) /
│                                  make_sdm_hook() (the closure run_5d_hsm.R passes
│                                  to virtualSpecies_nd(sdm_hook=)). GLM is LINEAR
│                                  (USE-faithful; a quadratic GLM is oracle on a
│                                  Gaussian niche → sampler-blind). Maxent = maxnet
│                                  (no Java). Every backend tryCatch→NA.
├── hsm_plots.R                    downstream comparison figures (USE 2_ViolinPlots
│                                  + 3_DunnTest port): plot_hsm_violin (metric ×
│                                  algorithm by sampler), hsm_dunn_test (uniform+ vs
│                                  others, WITHIN-species, holm), summarize_hsm,
│                                  plot_hsm_species_overlay/_row (cross-species agg:
│                                  hue=sampler/shade=species), plot_tune_heatmap
│                                  (knob-sweep heatmaps), plot_hsm_species_metrics +
│                                  hsm_win_matrix/plot_hsm_win_matrix (all-metrics
│                                  violins + where-uniform+-beats-naive matrix).
│                                  Reuses experiment_plots.R sampler colours.
├── hsm_report.R                   builds a scrollable multi-page PDF from a 3-sampler
│                                  HSM CSV (results5d_hsm): overview medians, the win
│                                  matrix per predictor set, and per-species all-metrics
│                                  violin grids + hsm_win_matrix.csv.
├── pa_buffer.R                    geographic "buffer-out" PA sampler (excl. cells
│                                  within 50 km of presences; EPSG:3035 distance).
│                                  Dimension-agnostic spatial baseline; add to any
│                                  pa_samplers list (used in both experiments).
├── imports.R                     preflight package-load check (terra, USE.MCMC,
│                                  hypervolume, sf, future/furrr, …); run inside the
│                                  apptainer env before SLURM dispatch (used by
│                                  sbatch/submit_smoke_test.sh).
├── results/{2d,5d}/              persisted experiment figures + metrics CSV +
│                                  summary .rds (re-style figures without re-running;
│                                  survives the scratch purge).
├── vignettes/2d-experiment.Rmd    reproducible vignette for the §2.3 results
│                                  (loads the SLURM run, or a small inline demo).
├── README.md                      methodology summary (matches §"What this
│                                  project is" above).
├── LICENSE                        GPL-3 (see file).
├── .claude/
│   ├── CLAUDE.md                  this file.
│   └── skills/
│       ├── euler-rstudio-server/  runtime guidance for Euler RStudio.
│       ├── euler-r-spack-setup/   R-install recipe for Euler rocker/rstudio.
│       ├── publication-figures/   authoritative figure→producer map for the paper
│       │                          + downstream-HSM report (Part A vignette MCMC
│       │                          diagnostics, Part B GaussNiche figures).
│       ├── hsm-ablation/          run/recompute/interpret the uniform+ settings
│       │                          ablation sweep.
│       └── adversarial-paper-review/  heavy multi-persona adversarial review of the
│                                  markov-chain-sampler manuscript (PDF-first;
│                                  Workflow-driven; adversarial_review.workflow.js).
```

## Paper figures produced here

`run_2d_experiment.R` and `run_5d_experiment.R` (via `experiment_plots.R`) produce
the markov-chain-sampler paper's **Results** figures
(`sampler-comparison-boxplots.pdf` 2-D, `sampler-comparison-boxplots-5d.pdf` 5-D)
and the **appendix** per-species PC scatterplot matrices (`pcmatrix-5d-sp1..4.pdf`)
and geographic pseudo-absence grid (`geo-grid-5d.pdf`). They are written under
`<scratch>/GaussNiche/results{2d,5d_experiment}/` and **copied by hand into
`../markov-chain-sampler-paper/graphics/`** — NOT via the paper's `make figures`
(which pulls only the methods + MCMC-diagnostic figures from the USE.MCMC vignette).
Re-render without recomputing via `sbatch --export=ALL,MODE=figures
sbatch/submit_5d_experiment.sh`. The paper's `CLAUDE.md` holds the full
figure→producer map for both repos.

`run_5d_hsm.R` additionally produces the **appendix** downstream-HSM figures
(`hsm_violin_{pc5,raw12}_*.pdf`, `hsm_dunn_summary_5d_*.pdf`) + the
`hsm_metrics_5d_full.csv` / `hsm_dunn_5d_full.csv` tables under
`<scratch>/GaussNiche/results5d_hsm/`, hand-copied into the paper's `graphics/`
for `text/appendix/appendix.tex` (same copy-by-hand convention as the other
GaussNiche figures).

## Higher-dimensional environment (the `5d-niche` extension)

These files keep the 2-D pipeline above intact (as the reference) and add a
**k-dimensional** sibling so the *same* analysis (4 species: generalist/specialist
× common/rare) runs in a 5-D environment where ≥5 PCs are needed for >80% variance.
The richer E-space is built by adding climate-orthogonal datasets (soil, terrain,
vegetation) to the WorldClim baseline, all via the `geodata` package on the
WorldClim 10 arc-min Central/W-Europe grid. Files added (portable R; the
`sbatch/submit_*.sh` are the Euler/apptainer runners — see the skills):

- `build_env_stack.R`    download + harmonise a broad candidate pool (clim/soil/
                         ter/veg/anth) onto the baseline grid; cache to scratch.
                         INCREMENTAL: each block is cached to block_<name>.tif keyed
                         by its layer set, so a SLIGHT env change re-harmonises only
                         the changed block; terrain is derived on a buffered Europe
                         crop (not the globe). Both bit-identical to a clean build
                         (validated: 41 s vs the old multi-hour build). Overrides:
                         FORCE_REBUILD=1, ENV_OUT_DIR=<sandbox>.
- `analyze_env_stack.R`  per-variable stats, correlation/VIF, per-block
                         orthogonality-to-climate, full PCA with cumulative-80 /
                         Kaiser / broken-stick "useful PC" criteria, forward block
                         inclusion, curated candidate configs.
- `build_final_stack.R`  subset the LOCKED layer set (`lean_natural`: 5 climate +
                         soil pH/clay/bdod/SOC + terrain TRI/TPI + tree cover →
                         5 PCs ≈ 82.5%), run PCA, save `background_5d.rds`.
- `virtualSpecies_nd_fn.R`  k-D engine. Exports `pa_random`, `pa_mcmc` (forwards
                         `dimensions = c("PC1",…,"PCk")` to the dimension-agnostic
                         `USE.MCMC::paSamplingMcmc`), `pa_nn` (k-D via the now
                         dimension-general `USE.MCMC::paSamplingNn`: uniform-box
                         proposals + NN-remap with a √(d/2) support-threshold
                         correction; acceptance falls with d, so it is
                         compute-bound in high d — raise `n.candidates`/`n.tr`),
                         `compute_bandwidth_nd()`, `virtualSpecies_nd()`. `pa_nn`
                         is provided but NOT in the default sampler list — add it
                         via `pa_samplers = list(random = pa_random, mcmc =
                         pa_mcmc, nn = pa_nn)`. Only the KDE backend
                         `USE.MCMC::paSampling` remains 2-D-only.
                         `virtualSpecies_nd()` also accepts an optional
                         `sdm_hook=` closure (default NULL → byte-identical to
                         before): called once per (sampler, realisation) on the
                         already-drawn PA set to fit downstream HSMs, its tidy
                         rows collected into the result's `$hsm` slot. The hook's
                         RNG is isolated (save/restore `.Random.seed`) so the
                         hypervolume diagnostics are unperturbed. See hsm_eval.R.
                         Perf flags (all default to the original behaviour ->
                         run_5d_* unchanged): `compute_hypervolume=FALSE` skips the
                         2 hypervolume_gaussian() calls (~70% of a 5-D task;
                         overlap->NA); `make_plots=FALSE` skips the per-call
                         ggplot/KDE/projection objects; `compute_reference=FALSE`
                         skips the one-off reference-PA draw per sampler;
                         `manage_plan=FALSE` lets the caller own the future plan.
                         pa_mcmc also takes `precomputed.env` as a FILE PATH
                         (inlined per-worker cache) so a sweep ships a string, not
                         the multi-MB bundle, to each worker. See the "Downstream
                         HSM ... running it fast" section below for the full why.
- `5_highdim_species.R`  the 4 species in 5-D; mu/sigma as multiples of each PC's
                         background sd (generalist σ=1.0·sd / specialist σ=0.5·sd;
                         common μ = KDE mode of the 5-D background / rare μ = fixed
                         peripheral point c(1.55,1.13,1.4,0.8,1.2)). `smoke` mode for
                         validation. (run_5d_experiment.R is the §2.3 5-D arm built
                         on this engine — see the main repository map.)
- `6_compare_5d.R`       cross-species/sampler report (combined metrics + PDF).
- `bench_cache.R`        benchmark + correctness check for the
                         `USE.MCMC::precomputeMcmcEnvironment()` cache that
                         `pa_mcmc` forwards via `precomputed.env=` (cached vs
                         uncached timing + bit-identity). Run with
                         `sbatch/submit_bench_cache.sh`.

Key k-D differences vs the 2-D module: `mu`/`sigma` are length-k vectors,
`Σ = diag(sigma)·Cor·diag(sigma)`; suitability/Bernoulli/`hypervolume_gaussian`/
bandwidth operate on the k PC columns; niche & PA plots are drawn on the PC1×PC2
*projection* (suitability sliced at μ); coverage is reported per axis
(`rel_cov_PC1…PCk`). All heavy steps run on SLURM via the `sbatch/submit_*.sh` scripts
(apptainer + the rocker SIF), never the login node.

## Downstream HSM, the uniform+ ablation, and running it fast (what & why)

Read this to understand *what* the recent work is and *why* it is shaped the way
it is — it is the fast on-ramp for anyone picking this up.

**The research arc (what we're doing).** Beyond the sampling-quality diagnostics,
the goal is to reproduce the USE paper's *downstream* test in 5-D: fit the same 5
HSM families (GLM/GAM/RF/BRT/Maxent, Java-free) to each Bernoulli realisation of
the 4 species under each sampler, and ask whether **uniform+ (`pa_mcmc`) yields
better species-distribution models** than the naive samplers (RND, buffer-out).
Settled finding: **uniform+ beats the naive samplers on discrimination
(AUC/TSS/Sensitivity/Kappa) but loses on truth-recovery (cor/RMSE-to-known-
suitability, CBI)** — random pseudo-absences reconstruct the true surface best.
The **ablation** (`tune_5d_hsm.R`) sweeps the two `pa_mcmc` cutoffs to test whether
tuning closes that gap: it does not — `environmental.cutof.percentile` is the lever
(`species.cutoff.threshold`, incl. the `=1` pure-uniform endpoint, is flat), and
raising it only lifts truth-recovery to buffer-out's level, never to RND, before the
sampler degenerates (env-cutoff ≳ 0.25, a borderline/stochastic onset). So the
deficit is **structural to environmental-uniform sampling**, framed as an
objective-dependent trade-off. (NB: the `=1`/unify/CRAN USE.MCMC rebuild did **not**
change the sampler — verified byte-identical, `a57e34e`=HEAD; the once-noted "env=0.2
degenerate" was a transient development-build artifact, not the committed mechanism.) Full detail + how to recompute on an env change:
the **hsm-ablation** skill. The figures it produces: the **publication-figures** skill.

**Where the time goes (cost model, measured by `bench_hsm.R`).** Per 5-D
realisation, single-threaded: `hypervolume_gaussian` ×2 ≈ 15.8 s (the dominant
cost — 5-D Monte-Carlo); `pa_mcmc` ≈ 2 s (~1.74 s is a fixed `densityMclust` GMM
fit + NN remap, the chain itself ~free); the SDM hook ≈ 4.5 s (fit + predict over
the ~16.5k-cell background). For a downstream-HSM **sweep** the hypervolume is pure
overhead, so the engine has flags to strip everything the sweep does not read.

**Engine perf flags** (`virtualSpecies_nd`; every one defaults to the original
behaviour, so `run_5d_experiment.R` / `run_5d_hsm.R` are byte-for-byte unchanged):
- `compute_hypervolume = FALSE` — skip the 2 hypervolume calls (~70 % of a task).
- `make_plots = FALSE` — skip the per-call ggplot / KDE / projection objects.
- `compute_reference = FALSE` — skip the one-off reference-PA draw per sampler.
- `manage_plan = FALSE` — the caller sets the `future` plan once (no per-call respawn).
- `pa_mcmc(precomputed.env = <path>)` — ship the MCMC bundle to workers by FILE
  PATH (a string, cached per worker), not the multi-MB object.
`tune_5d_hsm.R` turns them all on; each is bit-identical when off.

**Single-node parallelism (ONE node, ONE job, ≤64 cores — why this shape).** The
ablation runs as a single `sbatch` on one node: `tune_5d_hsm.R` fans the
independent `(cell × species)` jobs across the node's cores in ONE `future` pool
(each job's realisations serial inside its worker). The default 4×4 grid is
16×4 = 64 jobs — a clean fit for 64 cores (one wave). **Deliberately no SLURM job
arrays**: a single node / single job is the easiest thing for someone else to
reproduce (`sbatch sbatch/submit_tune_5d_hsm.sh`). Results depend only on the
layered seeds inside `eval_realization`, so this is bit-identical to a serial run.

**Env-build is incremental** (`build_env_stack.R`). Each environmental block
(soil/terrain/veg/anth) is cached to its own `block_<name>.tif` keyed by its layer
set, so editing one block's variables — a *slight* env change — re-harmonises only
that block. Terrain is derived on a buffered European crop, not the globe. Both are
bit-identical to a clean build. `FORCE_REBUILD=1` ignores the caches; `ENV_OUT_DIR=
<dir>` builds to a sandbox (used to bit-identity-validate a change without touching
the live `env5d`).

**Discipline — do this whenever you change anything here:**
- **Bit-identity:** an "optimization" must NOT change the numbers. Validate by
  re-running on the unchanged input and diffing against the committed reference:
  the full per-realisation rows live in the committed summaries
  `results/5d_tune/summary_tune_grid.rds` (`$tune`, the 20k-row grid = 40 cells:
  8 env × 5 species_cutoff incl the `=1` pure-uniform endpoint) and
  `results/5d_hsm/summary_5d_hsm_full.rds` (`$hsm`, the 6k-row table), plus the
  `env5d` rasters. The raw `*_metrics_5d_*.csv` dumps are **gitignored** (multi-MB,
  redundant with those rds) — diff against `readRDS(...)$tune`/`$hsm`, or
  regenerate a CSV with `write.csv(readRDS(f)$tune, …, row.names=FALSE)`.
  `bench_cache.R` is the template for a with/without bit-identity check.
- **Never edit a script while a SLURM job is running it** — `Rscript` parses the
  file incrementally, so the running job dies with "unexpected end of input".
- **Regenerate, don't recompute:** every figure/table comes from a saved CSV
  (`hsm_report.R`, `hsm_aggregate_report.R`, `plot_tune_heatmap`), so a
  post-processing fix never re-runs the sweep.

## Sampler interface

Every pseudo-absence sampler is a function with this signature:

```r
sampler(background, N_pa, pres = NULL, seed = 123, ...)
# returns a data.frame with the same columns as `background`
```

`background` is the full PC-space data.frame (x, y, PC1, PC2, suit).
`pres` is the subset where the Bernoulli draw was 1 — required by samplers
that exclude presence-like environments.
`...` is forwarded by `virtualSpecies()`, which is how sampler-specific
options reach individual samplers (e.g. `grid.res`, `thres` for `pa_uniform`;
`chain.length`, `burnIn`, `engine` for `pa_mcmc`). Samplers should swallow
unknown args via `...` so they remain interchangeable.

To add a new sampler:

1. Write it to the interface above.
2. Add it to the named list passed as `pa_samplers = list(name = fn, ...)`
   when calling `virtualSpecies()`.
3. Anything sampler-specific goes through `...`. Don't widen
   `virtualSpecies()`'s formal arguments.

## Critical contract — `env.rast` must be the ORIGINAL rasters

Both `pa_uniform` and `pa_mcmc` accept `env.rast` (forwarded as
`pa_env_rast` from `virtualSpecies()`). It **must** be the original
environmental SpatRaster (e.g. `envData`), **not** `rpc$PCs`. The USE /
USE.MCMC samplers run `rastPCA()` internally; re-PCA-ing the already
orthogonal PC scores produces a rotated E-space where the
presence-exclusion filter no longer corresponds to the analysis PC axes —
making `pa_uniform` and `pa_random` numerically indistinguishable. Both
samplers re-match results back to the analysis PC scores by **geographic
coordinates** (rounded to 4 d.p. ≈ 11 m).

## Conventions

- **Reproducibility**: seeds are layered — `seed_base + r` for Bernoulli
  realisation `r`, `seed_pseudo_base + r` for the PA draw of realisation
  `r`. Don't add `set.seed()` calls inside helpers; rely on the layered
  seeds so changes to one sampler don't perturb another.
- **`max_pres` cap**: presences are capped (default 500) before
  `hypervolume_gaussian()` and the sampler call to keep wall-time
  bounded. The full Bernoulli draw is still used for prevalence and the
  true-absence match — never replace the full draw with the subsample.
- **No `setwd()` in tracked code**. `1_developing_framework.R` has a
  legacy `setwd()` line; leave it as a working-directory hint but don't
  copy that pattern elsewhere. Sourced scripts should assume the working
  directory IS the repo root.
- **Plot generation**: every diagnostic returns its ggplot in
  `result$samplers[[name]]$plots` rather than printing as a side-effect.
  Side-effects are only printed by the example scripts (1_*, 2_*, 3_*).
- **Pseudo-absence equality**: `N_pa` is derived from the **subsampled**
  presence count so the pres:PA ratio (`bgk_prev`) stays consistent
  across realisations.
- **Downstream HSM arm**: the 50 Bernoulli realisations ARE the replication
  (one held-out 70/30 split per realisation, no USE-style internal `n=5`
  subsampling). HSM stats live in their OWN tidy CSV (`results5d_hsm/`), never
  inside the 90 MB per-species result rds. Maxent = `maxnet` (pure-R, never
  Java). The GLM is **linear** (USE's `Observed ~ .`): a quadratic GLM is
  near-oracle on the exactly-Gaussian niche and blind to the sampler. Dunn
  comparisons are **within-species** (the realisations are one species on one
  fixed background, not 50 independent species) — no single omnibus p.

## Running the example

```r
source("virtualSpecies_fn.R")
source("2_testing_wrapper_function.R")  # full reproducible run
```

Heavy bits and their wall-time:

- **`hypervolume_gaussian()`** is called 2× per realisation per sampler.
  ~0.5–2 s each at the default `max_pres = 500`.
- **`pa_mcmc`** with `chain.length = 20000` takes ~10–20 s per call.
- **`USE.MCMC::optimRes()`** (§3b of `2_testing_wrapper_function.R`) is
  the slowest one-off step — minutes. Pre-compute and cache its return
  value if you'll re-run.

For fast iteration set `n_realizations = 10` and `chain.length = 5000`.

## Self-update protocol

This file should reflect the **current** state of the code, not a snapshot.
Update it when you (the agent) observe any of:

1. **A new R file is added** at the repo root (`*.R`). Add a line to the
   "Repository map" describing it.
2. **A public function is added, renamed, or removed** in
   `virtualSpecies_fn.R`. Update the "Repository map" entry for that
   file and, if it changes the sampler call signature, the "Sampler
   interface" section.
3. **A new sampler** is added. Add it to the list of samplers in the
   Repository map and confirm the "Sampler interface" still applies.
4. **A new convention** emerges from the user's feedback (e.g. "always
   use seed X", "never call Y in helpers"). Add it under "Conventions".
5. **A new skill** is added under `.claude/skills/`. Add a one-line
   pointer at the top of this file.
6. **An assumption you find written here turns out to be wrong** (e.g.
   `max_pres` default changed, an interface drifted). Fix it immediately
   — outdated guidance is worse than no guidance.

Update procedure: edit this file in place, in the smallest patch that
captures the change. Don't rewrite sections that didn't change. Don't
add timestamps or changelogs — `git log` is the source of truth for
"when did this change".

What **not** to put in this file:

- Anything machine-/cluster-/host-specific → goes in
  `.claude/skills/<topic>/SKILL.md`.
- One-off conversation context, plans, or to-dos.
- Re-statements of what is already obvious from `README.md`.
- Author/user identity, API keys, paths under `/cluster/`, anything in
  `$HOME`.

The litmus test: can a stranger who cloned this repo onto a Mac read
this file and start contributing? If yes, it belongs here. If they'd
need to be on Euler to make sense of it, it belongs in a skill.

## Things outside the agent's purview (don't change without asking)

- `LICENSE`. GPL-3 — leave it alone.
- `README.md` methodology section. The user maintains the prose
  description; this file (CLAUDE.md) is the agent-facing twin. Keep them
  consistent in *content*, but the README is authoritative for *wording*.
- The legacy `setwd()` line in `1_developing_framework.R`. The user
  treats that script as a frozen reference; don't refactor it.
