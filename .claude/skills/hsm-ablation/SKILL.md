---
name: hsm-ablation
description: Run, recompute, extend, and interpret the uniform+ (pa_mcmc) settings ABLATION — a grid sweep of the two MCMC tunables (environmental.cutof.percentile × species.cutoff.threshold) measuring downstream-HSM accuracy against the fixed RND / buffer-out baselines, to test whether tuning can close uniform+'s truth-recovery gap. For the markov-chain-sampler paper appendix. Use whenever the user changes the 5-D environment and wants the ablation recomputed, asks to extend/sweep the pa_mcmc cutoffs, tune uniform+, or interpret the ablation result. Driver = tune_5d_hsm.R; baseline (all 3 samplers) = run_5d_hsm.R; outputs in results/5d_tune. Pairs with the publication-figures skill (the ablation heatmaps) and the euler-* skills (run on SLURM, never the login node).
---

# uniform+ settings ablation (the §2.3 downstream-HSM tuning sweep)

## What it is and why
The downstream-HSM result (baseline `run_5d_hsm.R`) showed uniform+ (the MCMC
sampler, `pa_mcmc`) **beats RND/buffer-out on discrimination** (AUC, TSS,
Sensitivity, Kappa) but **loses on truth-recovery** (`cor_truth`, `rmse_truth`,
CBI) — random pseudo-absences reconstruct the true suitability surface better.
This ablation asks: *can we tune uniform+'s two cutoffs to close that
truth-recovery gap without giving back the discrimination edge?* It sweeps a grid
of the cutoffs, fits the same 5 HSM families per cell, and compares the
downstream metrics against the **fixed** RND/buffer-out baseline (those samplers
are knob-invariant, so they are NOT re-run — reused from `results/5d_hsm`).

## The two knobs (and which is the lever)
Both are `pa_mcmc` / `USE.MCMC::paSamplingMcmc` args, forwarded through
`virtualSpecies_nd()`'s `...`:

- **`environmental.cutof.percentile`** — the percentile of the environmental GMM
  density **below which the chain is forbidden** (`paSamplingMcmc.R:144`). Raising
  it trims the low-density / rare-environment tail that uniform+ over-samples →
  pseudo-absences become more availability-like → better truth-recovery. **This is
  THE lever.** It also shrinks the samplable space (see degeneracy below).
- **`species.cutoff.threshold`** — the presence-GMM percentile setting the
  presence-exclusion strength (`paSamplingMcmc.R`, the `stats::quantile(species.densities, …)`
  line). **Direction matters: HIGHER = WEAKER exclusion** (the target is
  `1 − sp_density / quantile(sp_densities, p)`, so a high quantile `p` excludes only
  the densest presences; a low `p` excludes almost everywhere). **A non-lever** here:
  flat effect on every metric across the `{0.9,0.75,0.6,0.5}` grid — **and flat at
  the `1.0` endpoint too** (see below).
  **Endpoint — `species.cutoff.threshold = 1`** skips the presence GMM entirely and
  samples the environment **uniformly** (presence-model OFF; `USE.MCMC` commit
  `6903ab1`, C++ `Inf`-cutoff sentinel). This is the *continuous limit* of the above,
  NOT a reversal: at `p = 1` only the single max-density point would be excluded
  (measure-zero ≈ no exclusion), so pure uniform is the natural endpoint. **The
  committed grid now INCLUDES `1.0`** as the fifth `species_cutoff` column: it
  confirms the non-lever finding extends to the no-exclusion limit (its cell medians
  sit inside the `{0.9..0.5}` band, e.g. cor_truth 0.618 at env 0.001, 0.695 at env
  0.15). So disabling presence exclusion entirely does not change downstream-HSM
  accuracy.

The swept knobs only set thresholds **after** the per-realisation `densityMclust`
species-GMM fit, which is why the GMM fit is identical across cells (a known
redundancy — see the performance note at the bottom).

## How to run
Heavy → SLURM only (per the root CLAUDE.md + euler skills). Driver `tune_5d_hsm.R`,
submit `sbatch/submit_tune_5d_hsm.sh`.

```bash
# install HSM deps first if the env is fresh:  sbatch sbatch/submit_install_hsm_deps.sh
# PARALLEL smoke (validates engine flags + grid plumbing; REQUIRED before full):
sbatch --export=ALL,MODE=smoke --cpus-per-task=2 --time=00:30:00 sbatch/submit_tune_5d_hsm.sh
# full grid (script DEFAULT is 4x4: env {0.001,0.005,0.01,0.05} x species {0.9,0.75,0.6,0.5}):
sbatch sbatch/submit_tune_5d_hsm.sh
# the COMMITTED canonical grid = 8 env x 5 species_cutoff (incl the 1.0 uniform
# endpoint) = 40 cells. A comma in a value breaks --export=<list>, so set the vars
# in the submitting env and use --export=ALL:
MODE=grid ENV_CUTOFFS="0.001,0.005,0.01,0.05,0.1,0.15,0.2,0.25" \
  SPECIES_CUTOFFS="1.0,0.9,0.75,0.6,0.5" \
  sbatch --export=ALL sbatch/submit_tune_5d_hsm.sh
```
Env overrides: `ENV_CUTOFFS`, `SPECIES_CUTOFFS`, `N_REALIZATIONS` (default 25),
`MAX_PRES` (300), `MIN_PRES` (12). Inputs: `env5d/{final_stack_lean_natural.tif,
background_5d.rds}` + the baseline `results5d_hsm/hsm_metrics_5d_full.csv`
(RND/buffer reference). The committed grid: env ∈ {0.001,0.005,0.01,0.05,0.1,0.15,0.2,0.25},
species ∈ {1.0,0.9,0.75,0.6,0.5}. env=0.25 is ~24% random-fallback-contaminated and
env=0.3 degenerates further — see "the result" / pitfalls.

## Single-node parallelism + speedups (why the grid runs in ~minutes)
The sweep runs as ONE `sbatch` on ONE node (no job arrays — easiest to reproduce).
`tune_5d_hsm.R` fans the independent `(cell × species)` jobs across the node's
cores in a SINGLE `future` pool, each job's realisations serial INSIDE its worker
(`parallel = FALSE` on the `virtualSpecies_nd` call). The default 4×4 grid is
16 cells × 4 species = 64 jobs → a clean fit for a 64-core node (one wave,
`submit_tune_5d_hsm.sh` requests `--cpus-per-task=64`); the committed 40-cell grid
is 160 jobs → 3 waves, ~11 min wall. A `SpatRaster` cannot
cross a worker boundary, so `envData` is `terra::wrap()`ed once and `unwrap()`ed
inside each job. A worker also runs under `future`'s L'Ecuyer RNG, but the
Bernoulli `pa_matrix` draw uses `set.seed()` with no kind — so `virtualSpecies_nd`
pins `RNGkind("Mersenne-Twister")` at entry, making a worker bit-identical to the
main process (without this, the draws silently diverge under parallelism).

**The recent `USE.MCMC` rebuild did NOT change the sampler — proven byte-identical.**
The `=1`/unify/CRAN commits look scary but none touch the `species_cutoff<1`
mechanism: the `=1` branch is additive and inert for `<1`, the unify commit is
interface-only (arg rename + guards), and `predict.densityMclust`→`stats::predict`
is a verified no-op (`identical=TRUE`). Confirmed empirically: rebuilding `a57e34e`
(0.0.4, pre-`=1`/pre-unify) and diffing vs HEAD (0.0.5) gives **byte-identical
pseudo-absences** (same points, same env GMM densities/threshold) at env 0.05 and
0.2. `env5d` + packages are unchanged too. *To re-verify after any future
`USE.MCMC` change:* `git -C ~/USE.MCMC worktree add <tmp> <old-sha>`; in apptainer
`R CMD INSTALL --no-multiarch --library=<tmplib> <tmp>` (deps from the main lib);
then run `paSamplingMcmc(precomputed.env=<bundle>, seed.number=s, …)` with fixed
seeds under each lib (`.libPaths(<tmplib>)` vs main) and diff the sampled
coordinates — identical fingerprints ⇒ mechanism unchanged.
**Why the OLD pre-`b716c0c` grid still doesn't reproduce bit-for-bit:** it was
generated *mid-development* (the jobs ran the night of Jun 26, before the build/RNG
state was finalized + committed; the `=1` work landed hours later). So it carries a
different random *stream*, not a different mechanism — every per-realisation value
differs, yet per-cell medians match the clean grid within 0.02 (cor_truth 0.021,
auc 0.006, tss 0.022, rmse_truth 0.014; mixed sign). The **current** 40-cell grid
is clean and fully reproducible: re-confirm a GaussNiche optimization by a
*matched-key* diff (key on `species × sampler × realization × predictor_set ×
algorithm × env_cutoff × species_cutoff`, max abs diff over shared rows) at a fixed
build — when env 0.2/0.25 were added the env≤0.15 sub-grid re-ran bit-for-bit
(max |diff| 0).

The sweep also strips the engine to just the downstream-HSM numbers via flags:
- **A1** `pa_samplers = list(mcmc = pa_mcmc)` — uniform+ only (RND/buffer reused).
- **A2** `compute_hypervolume = FALSE` — skip the 2 `hypervolume_gaussian()`
  calls (the dominant ~70%-of-task 5-D cost; overlap → NA).
- **A3** `predictor_sets = list(pc5 = pc_cols)` — pc5 only (raw12 same story).
- **A4** `n_realizations = 25` (ranks cells; re-confirm a winner at 50).
- **`make_plots = FALSE`** — no per-call ggplot/KDE/projection objects.
- **`compute_reference = FALSE`** — no one-off reference-PA draw per sampler.
- **`precomputed.env = bundle_path`** — the MCMC bundle is `saveRDS`-ed once and
  passed BY FILE PATH; `pa_mcmc` loads it per-worker-cached (inlined, ships a
  string not the multi-MB object). The hooks are built once per species and reused
  across cells (identical train/test split per realisation → fair comparison).
All flags live in `virtualSpecies_nd_fn.R`; each defaults to the original behaviour
(`make_plots`/`compute_hypervolume`/`compute_reference`/`manage_plan` = TRUE,
object `precomputed.env`), so `run_5d_hsm.R` / `run_5d_experiment.R` are untouched.

## The result (settled — do not re-litigate)
- **`environmental.cutof.percentile` is the lever**: raising it 0.001 → 0.25
  monotonically lifts `cor_truth` ~0.62 → ~0.72 and lowers `rmse_truth`,
  confirming the rare-environment-trimming mechanism.
- **`species_cutoff` is a non-lever — including the `1.0` pure-uniform endpoint.**
  At every env level the five `species_cutoff` columns `{1.0,0.9,0.75,0.6,0.5}` are
  flat (spread ≤0.01 in cor_truth), so turning the presence model OFF entirely
  (`= 1.0`) lands inside the band: disabling presence exclusion does not change
  downstream-HSM accuracy.
- **The gap never closes.** cor_truth reaches ~0.70 (env 0.15–0.20) and ~0.72 at
  env 0.25 — i.e. it only catches up to *buffer-out* (0.725), still **~−0.04 below
  RND (0.759)**, with diminishing returns; and the 0.25 climb is partly artefact
  (fallback, below).
- **Degeneracy onset is ~0.25–0.3, borderline/stochastic — and ALWAYS was in the
  committed mechanism.** The clean committed sampler (`a57e34e`=HEAD) samples env=0.2
  robustly (probe 10/10 at sp=1; grid ~4% fallback) — verified by rebuilding `a57e34e`
  (300 pts at env=0.2, byte-identical to HEAD). The earlier "env=0.2 degenerate" note
  was a **development-build artifact** (the old grid was generated mid-development),
  NOT a property of the committed code — disregard it. Soft degeneracy onsets at
  env≈0.25 (~24% of grid units fall back) and env≈0.3 (~20% at sp=1, probe). The
  committed grid runs to env=0.25; env≥0.3 is not worth committing.
- **Discrimination erodes** as env rises (AUC 0.97 → 0.94, TSS 0.88 → 0.80 by
  env 0.25) — not a free lunch.
- **Random-fallback contamination grows with env.** Units whose 5 `species_cutoff`
  rows are identical (incl. uniform `=1` — a knob-blind random draw, which a genuine
  chain never produces): 0.001 1% / 0.005 0 / 0.01 1% / 0.05 0 / 0.1 0 / 0.15 4% /
  0.2 4% / **0.25 24%** (~1% baseline = sparse-presence rare-species reals,
  env-independent). **The env=0.25 row is materially contaminated** — its cor_truth
  (~0.72) is inflated toward RND (random recovers the true surface better), so treat
  it as an UPPER BOUND; env≤0.2 is essentially clean. The genuine-MCMC ceiling is at
  or below ~0.72 — which only strengthens "the gap is structural".

**Conclusion: uniform+'s truth-recovery deficit is structural to
environmental-uniform sampling, not a settings artifact — tuning narrows it ~50%
then the sampler degenerates before closing it.** The honest appendix framing is
**objective-dependent**: uniform+/buffer-out for discrimination, random for
suitability-surface recovery.

## Outputs (durable in the repo; survive the scratch purge)
`results/5d_tune/` — canonical = the `_grid` generation (env 0.001→0.25):
- `summary_tune_grid.rds` — **the committed source of truth.** `$tune` is the full
  per-(cell, species, realisation, algorithm) table (20k rows = 40 cells × 4 species
  × 25 reals × 5 algos, tagged `env_cutoff`/`species_cutoff`); `$cell` the 40 per-cell
  medians; `$ref` the baseline.
- `tune_cell_medians_grid.csv` — per-cell median table (+ `n` survivorship column).
- `tune_heat_<metric>_{value,delta}_grid.pdf` — value + Δ-vs-RND heatmaps (the
  ablation figures; see the publication-figures skill).
- `TUNE_ablation_config.txt` — the run config.
- `tune_hsm_metrics_5d_grid.csv` — the raw CSV is **gitignored** (3.4 MB, redundant
  with `summary_tune_grid.rds$tune`). A sweep rewrites it to the scratch dir; or
  regenerate it from the rds:
  `write.csv(readRDS("results/5d_tune/summary_tune_grid.rds")$tune, "…grid.csv", row.names=FALSE)`.
`results/5d_hsm/` is the same shape: figures + `summary_5d_hsm_full.rds` (`$hsm` =
the full 6k-row table) are committed; `hsm_metrics_5d_full.csv` is gitignored.
Scratch run dir: `<scratch>/GaussNiche/results5d_tune/`. Keep only the latest
generation (drop superseded `_full` / `_full24` and any `*_smoke*`).

## Recompute on an environment change (the user's main use case)
1. Rebuild the 5-D environment (`build_final_stack.R` etc. → new
   `env5d/{final_stack_lean_natural.tif, background_5d.rds}`).
2. **Re-run the baseline** `run_5d_hsm.R` (full) → fresh `results5d_hsm/
   hsm_metrics_5d_full.csv` (the RND/buffer reference the ablation reads).
3. Run the ablation: `sbatch sbatch/submit_tune_5d_hsm.sh` (smoke first).
4. Persist the `_grid` outputs into `results/5d_tune/` and update the config note.
The species niches in `tune_5d_hsm.R` are constructed from the background sd / KDE
mode exactly as `run_5d_hsm.R` (so they match); confirm prevalences are sane.

## Reusable post-processing (regenerate tables/figures WITHOUT recompute)
The committed source of truth is `summary_tune_grid.rds$tune` (the raw CSV is a
gitignored derivative — regenerate it from the rds if a tool wants a file path). To
re-style or recover after a post-processing crash, regenerate from it (no fitting):
- `hsm_aggregate_report.R <csv> <outdir>` — aggregate figures + tables.
- `plot_tune_heatmap(tune, metric, ref=, baseline=)` (hsm_plots.R) — one heatmap.
- combine generations by `rbind`-ing the per-run CSVs, then re-aggregate.

## Pitfalls / contamination checks
- **Random fallback** — at aggressive `env_cutoff`, `pa_mcmc` may return 0 PAs and
  fall back to `pa_random`; those realisations are then *random*, inflating
  `cor_truth` for the wrong reason. **In the single-node (outer-`furrr`) design the
  fallback `warning()` fires inside a worker and is swallowed at the worker boundary
  — it never reaches the `.err`, so grepping the log misses it** (it worked only in
  the old inner-parallel loop). Detect from the DATA: a fallback draw is knob-blind,
  so for a fixed `(species, realisation, algorithm, env_cutoff)` its `cor_truth` is
  **identical across ALL 5 `species_cutoff` — including the uniform `=1` column, which
  a genuine chain never matches** (uniform skips the GMM that `<1` fits). Count those
  units per env for the contamination rate (measured on the 40-cell grid: 4% at
  env 0.2, **24% at env 0.25**; ~1% baseline at low env = sparse rare-species reals).
  Exclude them for a clean number, or treat a high-fallback cell (env≥0.25) as an
  upper bound.
- **Survivorship vs fallback** — two faces of the same boundary degeneracy. *Old*
  builds **skipped** a `pa_mcmc`-returns-<5-PAs realisation (no rows → cell `n` < 100,
  survivorship bias). The *current* build instead **recovers** it via the
  `pa_mcmc`→`pa_random` fallback (cell `n` = 100, but that row is random). Check the
  `n` column AND the fallback signature (above); high-env cells lean mildly optimistic
  either way.
- **NEVER edit a script while a SLURM job is running it.** `Rscript` parses the
  file incrementally, so an edit shifts the not-yet-parsed tail and the job dies
  with "unexpected end of input" after running for a while. The sweep data is
  written before the post-processing, so on such a crash just regenerate the
  tables/figures from the saved raw CSV (above) — do not re-run the sweep.

## Performance note (the one remaining big redundancy)
The grid loops cells OUTER, so for each (species, realisation) the ~1.7 s
`densityMclust` species-GMM is refit identically for all cells (~23× redundant on
a 24-cell grid). Eliminating it needs a `precomputed.species.model` API in
USE.MCMC (mirroring `precomputed.env`) so the GMM is fit once per
(species, realisation) and reused across cells — a USE.MCMC edit, not GaussNiche.
Evaluate before a large re-sweep.
