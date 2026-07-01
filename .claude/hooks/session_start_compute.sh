#!/usr/bin/env bash
# SessionStart hook — register the live SLURM compute context into every Claude
# Code session started in this repo, and tell the model how to route jobs
# (run-here vs. sbatch). Wired via .claude/settings.json.
#
# Emits the SessionStart hook payload on stdout:
#   {"hookSpecificOutput":{"hookEventName":"SessionStart","additionalContext":"..."}}
# On a generic (non-Euler) host it prints nothing — a clean no-op, so this is
# safe to commit and safe for anyone who clones the repo onto a workstation.
#
# It must NEVER block a session start, so: no `set -e`, always exit 0, and the
# JSON is built with a pure-bash escaper (no python3/jq dependency).
set -uo pipefail

HOOK_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=/dev/null
source "$HOOK_DIR/compute_env.sh" 2>/dev/null || exit 0
gn_detect_compute 2>/dev/null || exit 0

# Build the context block for this host class. Unquoted heredoc so ${GN_*}
# expand; a literal "$USER" is written as \$USER (no shell expansion of it).
CTX=""
if [[ "$GN_ON_COMPUTE" == "1" ]]; then
  read -r -d '' CTX <<EOF
COMPUTE ALLOCATION DETECTED (injected by .claude/hooks/session_start_compute.sh)

This Claude Code session is running INSIDE a SLURM compute allocation, not on a
login node:
  node        ${GN_NODE}
  slurm job   ${GN_JOBID}
  cpus        ${GN_CPUS}
  memory      ${GN_MEM_GB:-?} GB
  gpus        ${GN_GPUS}
  time left   ${GN_TIME_LEFT:-unknown}

How to run compute this session:
- The machine-wide "never compute on a login node, always sbatch" rule is
  SATISFIED: you are already on allocated compute. Run jobs DIRECTLY here. Do
  NOT sbatch from inside this allocation (it double-allocates and queues the
  work behind your own job).
- Preferred: run \`.claude/bin/gnrun <script.R> [mode] [workers]\` — it executes
  the script directly here via apptainer using up to ${GN_CPUS} workers, e.g.
  \`.claude/bin/gnrun run_5d_hsm.R smoke\`.
- Stay inside the allocation: cap parallel workers at ${GN_CPUS}; keep total
  memory under ${GN_MEM_GB:-?} GB or SLURM kills the session; keep BLAS/OMP
  threads pinned to 1 (APPTAINERENV_OMP_NUM_THREADS=1 ...) so furrr/PSOCK
  workers do not oversubscribe the cores.
- This allocation ends in ${GN_TIME_LEFT:-unknown}; anything longer must go to a
  separate sbatch job. See the euler-compute-session skill for the manual
  apptainer form and how this interactive session was launched.
EOF
elif [[ "$GN_CLASS" == "login" ]]; then
  read -r -d '' CTX <<EOF
LOGIN NODE (injected by .claude/hooks/session_start_compute.sh)

This session is on an Euler login node (${GN_NODE}) — no SLURM allocation. Run
NO heavy compute here. Submit it with \`.claude/bin/gnrun <script.R> [mode]\`
(auto-sbatch) or the sbatch/submit_*.sh scripts. For immediate interactive
compute, launch a compute session — see the euler-compute-session skill.
EOF
fi

# Nothing to inject on a generic host → clean no-op.
[[ -n "$CTX" ]] || exit 0

# --- pure-bash JSON string escaping (order matters: backslash first) --------
esc=$CTX
esc=${esc//\\/\\\\}     # \  -> \\
esc=${esc//\"/\\\"}     # "  -> \"
esc=${esc//$'\t'/\\t}   # tab
esc=${esc//$'\r'/}      # strip CR
esc=${esc//$'\n'/\\n}   # newline -> \n

printf '{"hookSpecificOutput":{"hookEventName":"SessionStart","additionalContext":"%s"}}\n' "$esc"
exit 0
