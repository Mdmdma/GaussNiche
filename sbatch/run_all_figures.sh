#!/usr/bin/env bash
# run_all_figures.sh --- one command to rebuild ALL Part-B paper figures from scratch.
#
# LOGIN-NODE LAUNCHER (not a batch job): submits the 4 compute jobs with NO
# inter-dependency so they run FULLY IN PARALLEL, then a final PLOT job that renders
# every figure from the saved data (it depends on all 4). Run on the login node:
#   bash sbatch/run_all_figures.sh
#
# Structure (why):
#   env build (optional) ─▶ 4 compute jobs, all parallel, all compute+save only:
#       run_2d_experiment   (env-independent; bundled Worldclim_tmp)
#       run_5d_experiment   (env5d)
#       run_5d_hsm          (env5d)
#       tune_5d_hsm         (env5d)  ← no longer depends on run_5d_hsm: the ablation
#                                      Δ-vs-RND baseline merge moved to make_figures.R.
#   ─▶ submit_make_figures.sh  (afterok: all 4)  renders 2-D/5-D boxplots, PC-matrix,
#       geo-grid, HSM report + aggregates + violins + Dunn, ablation value+Δ heatmaps.
#
# Options (env vars):
#   MODE=full|smoke   forwarded to every job (default full).
#   SKIP_ENV=1        env5d already on scratch -> skip the two env-build jobs.
#   RUN_INSTALL=1     prepend the one-off HSM dep install as a prereq of the compute jobs.
#   DRY_RUN=1         print the sbatch commands without submitting.
set -euo pipefail
cd "$(dirname "$0")/.."
MODE="${MODE:-full}"
EXPORTS="ALL,MODE=$MODE"

sub() {   # sub "<colon-joined dep jobids or empty>" <sbatch script...> -> jobid
  local dep="$1"; shift
  local cmd=(sbatch --parsable --export="$EXPORTS")
  [[ -n "$dep" ]] && cmd+=(--dependency="afterok:$dep")
  cmd+=("$@")
  if [[ "${DRY_RUN:-0}" == "1" ]]; then echo "DRYRUN: ${cmd[*]}" >&2; echo "0"; else "${cmd[@]}"; fi
}
colon_join() { local out=""; local x; for x in "$@"; do [[ -n "$x" ]] && out="${out:+$out:}$x"; done; echo "$out"; }

echo "[orchestrate] MODE=$MODE SKIP_ENV=${SKIP_ENV:-0} RUN_INSTALL=${RUN_INSTALL:-0} DRY_RUN=${DRY_RUN:-0}" >&2

# --- env build (unless cached on scratch) -----------------------------------
ENVDEP=""
if [[ "${SKIP_ENV:-0}" != "1" ]]; then
  jid_es=$(sub "" sbatch/submit_build_env_stack.sh)
  jid_fs=$(sub "$jid_es" sbatch/submit_final_stack.sh)
  ENVDEP="$jid_fs"
  echo "  env:    build_env_stack=$jid_es -> build_final_stack=$jid_fs" >&2
else
  echo "  env:    SKIP_ENV=1 -> using existing env5d on scratch" >&2
fi

# --- optional one-off HSM dep install ---------------------------------------
INSTDEP=""
if [[ "${RUN_INSTALL:-0}" == "1" ]]; then
  INSTDEP=$(sub "" sbatch/submit_install_hsm_deps.sh)
  echo "  deps:   install_hsm_deps=$INSTDEP" >&2
fi
CDEP="$(colon_join "$ENVDEP" "$INSTDEP")"    # compute-job prereqs (env + install only)

# --- the 4 compute jobs: NO inter-dependency -> fully parallel ---------------
jid_2d=$(sub    ""      sbatch/submit_2d_experiment.sh)    # env-independent
jid_5dexp=$(sub "$CDEP" sbatch/submit_5d_experiment.sh)
jid_hsm=$(sub   "$CDEP" sbatch/submit_5d_hsm.sh)
jid_tune=$(sub  "$CDEP" sbatch/submit_tune_5d_hsm.sh)      # parallel with hsm now

# --- plot job: render everything once all 4 compute jobs succeed -------------
PLOTDEP="$(colon_join "$jid_2d" "$jid_5dexp" "$jid_hsm" "$jid_tune")"
jid_fig=$(sub "$PLOTDEP" sbatch/submit_make_figures.sh)

cat >&2 <<EOF
[orchestrate] submitted (the 4 compute jobs run in parallel; deps shown are env/install prereqs only):
  2d_experiment    = $jid_2d      (dep: none)
  5d_experiment    = $jid_5dexp   (dep: ${CDEP:-none})
  5d_hsm           = $jid_hsm     (dep: ${CDEP:-none})
  tune_5d_hsm      = $jid_tune    (dep: ${CDEP:-none})   <- now parallel with 5d_hsm
  make_figures     = $jid_fig     (dep: all 4 compute jobs)
Watch:   squeue -u \$USER
EOF
