#!/usr/bin/env bash
# GaussNiche harness — compute-context detector (sourceable + standalone).
#
# Detects whether the current process is running inside a SLURM compute
# allocation (an interactive `srun --pty` / `salloc` session, or a batch job)
# vs. on an Euler login node vs. a generic host, and registers the resources
# the allocation holds. Consumed by:
#   * .claude/hooks/session_start_compute.sh  (injects the context into a session)
#   * .claude/bin/gnrun                        (routes jobs: run-here vs sbatch)
#
# SOURCE it to get the GN_* variables, or run it standalone for a one-line report:
#   source .claude/hooks/compute_env.sh && gn_detect_compute && gn_report
#   bash   .claude/hooks/compute_env.sh          # prints the report
#
# Defensive everywhere: on a Mac/workstation with no SLURM it reports
# "workstation" and exits 0. It must NEVER fail a session start, so no `set -e`.

gn_detect_compute() {
  GN_NODE="$(hostname -s 2>/dev/null || echo unknown)"
  GN_JOBID="${SLURM_JOB_ID:-${SLURM_JOBID:-}}"

  # --- classify the host -------------------------------------------------
  # compute : inside a SLURM allocation on a non-login node (use it fully)
  # login   : an Euler login node (eu-login-*) — submit heavy work
  # host    : anything else (generic Linux/macOS workstation)
  #
  # The `! ^eu-login` guard matters: `salloc` can leave you in a shell that
  # still lives ON the login node with SLURM_JOB_ID set — that is NOT a compute
  # node, so we must not treat it as one.
  if [[ -n "$GN_JOBID" && ! "$GN_NODE" =~ ^eu-login ]]; then
    GN_CLASS="compute"
  elif [[ "$GN_NODE" =~ ^eu-login ]]; then
    GN_CLASS="login"
  else
    GN_CLASS="host"
  fi

  # --- CPUs allocated ----------------------------------------------------
  # SLURM_CPUS_PER_TASK is what our --cpus-per-task launches request; fall back
  # to the node-total, then to nproc for the non-SLURM case.
  GN_CPUS="${SLURM_CPUS_PER_TASK:-${SLURM_CPUS_ON_NODE:-}}"
  if [[ -z "$GN_CPUS" || ! "$GN_CPUS" =~ ^[0-9]+$ ]]; then
    GN_CPUS="$(nproc 2>/dev/null || echo 1)"
  fi

  # --- memory allocated (whole GB) --------------------------------------
  # SLURM reports memory in MB. --mem sets SLURM_MEM_PER_NODE; --mem-per-cpu
  # sets SLURM_MEM_PER_CPU (multiply by the core count for the node total).
  local mem_mb=""
  if [[ "${SLURM_MEM_PER_NODE:-}" =~ ^[0-9]+$ ]]; then
    mem_mb="$SLURM_MEM_PER_NODE"
  elif [[ "${SLURM_MEM_PER_CPU:-}" =~ ^[0-9]+$ ]]; then
    mem_mb=$(( SLURM_MEM_PER_CPU * GN_CPUS ))
  fi
  if [[ -n "$mem_mb" ]]; then
    GN_MEM_GB=$(( mem_mb / 1024 ))
  else
    GN_MEM_GB=""
  fi

  # --- GPUs --------------------------------------------------------------
  GN_GPUS="${SLURM_GPUS_ON_NODE:-${SLURM_JOB_GPUS:-${SLURM_GPUS:-0}}}"
  # SLURM_JOB_GPUS can be a comma-separated id list (e.g. "0,1"); count entries.
  if [[ "$GN_GPUS" == *,* ]]; then
    GN_GPUS="$(printf '%s' "$GN_GPUS" | tr ',' '\n' | grep -c . || true)"
  fi
  [[ "$GN_GPUS" =~ ^[0-9]+$ ]] || GN_GPUS=0

  # --- wall-time remaining ----------------------------------------------
  # `|| GN_TIME_LEFT=""` matters: under a caller's `set -euo pipefail` (gnrun),
  # a non-zero squeue or a SIGPIPE from `head` closing the pipe (141 with
  # pipefail) would otherwise abort the caller mid-detection.
  GN_TIME_LEFT=""
  if [[ -n "$GN_JOBID" ]] && command -v squeue >/dev/null 2>&1; then
    GN_TIME_LEFT="$(squeue -h -j "$GN_JOBID" -o '%L' 2>/dev/null | head -n1 | tr -d '[:space:]')" || GN_TIME_LEFT=""
  fi

  GN_ON_COMPUTE=0
  [[ "$GN_CLASS" == "compute" ]] && GN_ON_COMPUTE=1
  export GN_CLASS GN_ON_COMPUTE GN_NODE GN_JOBID GN_CPUS GN_MEM_GB GN_GPUS GN_TIME_LEFT
}

gn_report() {
  case "$GN_CLASS" in
    compute)
      echo "compute node  ${GN_NODE}  job ${GN_JOBID}"
      echo "  cpus=${GN_CPUS}  mem=${GN_MEM_GB:-?}G  gpus=${GN_GPUS}  time_left=${GN_TIME_LEFT:-?}"
      ;;
    login)
      echo "login node    ${GN_NODE}  (no allocation — submit heavy work with sbatch/srun)"
      ;;
    *)
      echo "host          ${GN_NODE}  (no SLURM — run as usual)"
      ;;
  esac
}

# Standalone invocation prints the report; sourcing defines the functions only.
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
  gn_detect_compute
  gn_report
fi
