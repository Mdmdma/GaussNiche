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
│                                    virtualSpecies()  main pipeline; returns
│                                                      list(niche, suit_rast,
│                                                      pa_rast, response_curves,
│                                                      samplers, plots, bw,
│                                                      background)
├── 2_testing_wrapper_function.R   reproducible example exercising
│                                  virtualSpecies() with all three samplers
│                                  on USE.MCMC::Worldclim_tmp.
├── 3_ virtualSpecies_cases.R      additional case-study runs (note the
│                                  space in the filename — keep it; quote
│                                  it when sourcing).
├── README.md                      methodology summary (matches §"What this
│                                  project is" above).
├── LICENSE                        GPL-3 (see file).
├── .claude/
│   ├── CLAUDE.md                  this file.
│   └── skills/
│       ├── euler-rstudio-server/  runtime guidance for Euler RStudio.
│       └── euler-r-spack-setup/   R-install recipe for Euler rocker/rstudio.
```

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
