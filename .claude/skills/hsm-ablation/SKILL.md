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
  presence-exclusion strength (`paSamplingMcmc.R:161`). **A non-lever** here:
  flat effect on every metric across the grid.

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
# full grid (default 4x4):
sbatch sbatch/submit_tune_5d_hsm.sh
# custom / extended grid:
sbatch --export=ALL,ENV_CUTOFFS=0.1,0.15,SPECIES_CUTOFFS=0.9,0.75,0.6,0.5 sbatch/submit_tune_5d_hsm.sh
```
Env overrides: `ENV_CUTOFFS`, `SPECIES_CUTOFFS`, `N_REALIZATIONS` (default 25),
`MAX_PRES` (300), `MIN_PRES` (12). Inputs: `env5d/{final_stack_lean_natural.tif,
background_5d.rds}` + the baseline `results5d_hsm/hsm_metrics_5d_full.csv`
(RND/buffer reference). The default grid: env ∈ {0.001,0.005,0.01,0.05}, species
∈ {0.9,0.75,0.6,0.5}; the usable range was extended to env {0.1,0.15} (0.2 is
degenerate — see below).

## Single-node parallelism + speedups (why a 24-cell grid runs in ~minutes)
The sweep runs as ONE `sbatch` on ONE node (no job arrays — easiest to reproduce).
`tune_5d_hsm.R` fans the independent `(cell × species)` jobs across the node's
cores in a SINGLE `future` pool, each job's realisations serial INSIDE its worker
(`parallel = FALSE` on the `virtualSpecies_nd` call). The default 4×4 grid is
16 cells × 4 species = 64 jobs → a clean fit for a 64-core node (one wave,
`submit_tune_5d_hsm.sh` requests `--cpus-per-task=64`). A `SpatRaster` cannot
cross a worker boundary, so `envData` is `terra::wrap()`ed once and `unwrap()`ed
inside each job. A worker also runs under `future`'s L'Ecuyer RNG, but the
Bernoulli `pa_matrix` draw uses `set.seed()` with no kind — so `virtualSpecies_nd`
pins `RNGkind("Mersenne-Twister")` at entry, making a worker bit-identical to the
main process (without this, the draws silently diverge under parallelism).

**Bit-identity is on VALUES, not row counts (validated).** After an env or
parallelism change, re-confirm by a *matched-key* diff (key on `species × sampler ×
realization × predictor_set × algorithm × env_cutoff × species_cutoff`, max abs diff
over shared rows) — NOT by row count. Every genuine-MCMC row matches the committed
grid bit-for-bit (max abs diff `0`, checked at env 0.001–0.05 and env 0.1). Row
counts differ slightly because the engine's `pa_mcmc`→`pa_random` fallback (shipped
with `precomputed.env`) recovers a few degenerate-boundary units the older committed
grid omits — ≤0.6 % of rows at the env-ceiling, shifting medians by <0.006. Identical
values + a few extra boundary fallbacks = the restructure is clean.

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
- **`environmental.cutof.percentile` is the lever**: raising it 0.001 → 0.15
  monotonically lifts `cor_truth` ~0.63 → ~0.70 and lowers `rmse_truth`,
  confirming the rare-environment-trimming mechanism. `species_cutoff` is flat.
- **The gap never closes.** Best usable cell (env=0.15) `cor_truth` ≈ 0.696, still
  **−0.063 below RND (0.759)** and below buffer-out (0.725), with **diminishing
  returns**.
- **env = 0.2 is DEGENERATE** — the top-20%-density exclusion shrinks the
  samplable space so much that `pa_mcmc` can't form ≥5 pseudo-absences and **every
  realisation is skipped** (0 usable rows). This is NOT a random fallback
  (fallback count = 0); ~0.1–0.15 is the hard ceiling.
- **Discrimination erodes slightly** as env rises (AUC 0.97 → 0.96, TSS 0.87 →
  0.84) — not a free lunch.
- **Survivorship bias grows** near the ceiling: surviving (species × realisation)
  units per cell drop 100 → ~96 (env 0.1) → ~91 (0.15) → 0 (0.2), so high-env
  `cor_truth` is mildly optimistic.

**Conclusion: uniform+'s truth-recovery deficit is structural to
environmental-uniform sampling, not a settings artifact — tuning narrows it ~50%
then the sampler degenerates before closing it.** The honest appendix framing is
**objective-dependent**: uniform+/buffer-out for discrimination, random for
suitability-surface recovery.

## Outputs (durable in the repo; survive the scratch purge)
`results/5d_tune/` — canonical = the `_grid` generation (env 0.001→0.15):
- `summary_tune_grid.rds` — **the committed source of truth.** `$tune` is the full
  per-(cell, species, realisation, algorithm) table (11.7k rows, tagged
  `env_cutoff`/`species_cutoff`); `$cell` the per-cell medians; `$ref` the baseline.
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
  the old inner-parallel loop). Detect from the DATA instead: a fallback draw ignores
  *both* MCMC knobs, so its rows are **identical across `species_cutoff`** for a fixed
  `(species, realisation, env_cutoff)` and carry a `cor_truth` far above the cell's
  MCMC median (≈0.85 vs ≈0.68 at env=0.1). Flag any high-`env_cutoff`
  `(species, realisation)` whose rows repeat verbatim across `species_cutoff` as a
  random fallback and exclude it for a clean number (at env=0.1 ≈12/2000 units).
- **Survivorship** — a different degeneracy: `pa_mcmc` returns <5 PAs → the
  realisation is skipped (no rows, NOT a fallback). Check the `n` column in
  `tune_cell_medians_grid.csv`; a cell well below 100 is survivorship-biased.
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
