---
name: euler-compute-session
description: Launch a long interactive SLURM compute session on ETH Euler and run the `claude` CLI INSIDE it, so a Claude Code session has immediate compute (no per-job sbatch round-trip). Covers the exact `srun --pty` launch recipe (account es_schin, resources, 6h), the ETH-proxy setup that lets `claude` reach the API from a compute node, and how the repo's SessionStart hook (.claude/hooks/session_start_compute.sh) auto-registers the allocation while `.claude/bin/gnrun` routes each R job (run-here on compute vs. sbatch on login). Use when the user wants an interactive compute session / dev session with compute in hand, mentions srun/salloc/interactive node, or when a session banner says "COMPUTE ALLOCATION DETECTED". Pairs with [[euler-rstudio-server]] (the JupyterHub/apptainer runtime) and [[euler-r-spack-setup]] (installing R native-dep packages).
---

# Running Claude Code on an Euler compute node (immediate-compute dev session)

The idea: instead of editing on the login node and `sbatch`-ing every run (queue
wait per iteration), grab **one interactive SLURM allocation** and run `claude`
*inside* it. Then Claude Code has cores/RAM in hand and runs jobs directly —
fast iteration — while the machine-wide "never compute on the login node" rule is
still honoured (you're on allocated compute, not the login node).

The repo's harness makes this automatic:

- **`.claude/hooks/session_start_compute.sh`** (SessionStart hook, wired in
  `.claude/settings.json`) detects on startup whether the session is inside a
  SLURM allocation and **injects the live resources** (node, job id, cpus, mem,
  gpus, time-left) + the routing rule into context. On a login node / plain
  workstation it's a clean no-op.
- **`.claude/bin/gnrun <script.R> [mode] [workers]`** routes a job by where it
  runs: **on a compute allocation → `apptainer exec … Rscript` directly** with
  the allocation's cores; **on a login node → `sbatch sbatch/submit_gnrun.sh`**.
- Both share the detector **`.claude/hooks/compute_env.sh`** (`gn_detect_compute`
  → `GN_CLASS` / `GN_ON_COMPUTE` / `GN_CPUS` / `GN_MEM_GB` / `GN_GPUS` /
  `GN_TIME_LEFT`). Run `bash .claude/hooks/compute_env.sh` for a one-line report.

## Launch recipe (copy-paste, from the login node)

Resources default to **2× the standard heavy run** (`run_5d_experiment.R` /
`run_5d_hsm.R` use 16 cpu × 6 GB) → **32 cpu × 6 GB = 192 GB**, for 6 h. Scale
`-c` / `--mem-per-cpu` / `--time` as needed.

```bash
# 0) tmux ON THE LOGIN NODE first — an interactive srun job dies if your SSH
#    connection drops. tmux keeps the allocation (and the claude session) alive.
tmux new -s claude          # reattach later with:  tmux attach -t claude

# 1) 6h interactive shell ON a compute node: 32 physical cores, 6 GB/core, es_schin.
#    Use srun --pty (NOT salloc: on this cluster salloc strands you on eu-login,
#    because use_interactive_step is not set). Do NOT name a partition — Euler
#    routes by --time (6h -> normal.24h). AMD EPYC nodes are HT-off, so 32 = 32
#    physical cores.
srun --account=es_schin \
     --ntasks=1 --cpus-per-task=32 --mem-per-cpu=6G \
     --time=06:00:00 \
     --pty bash

# --- you are now on eu-a*/eu-c*/eu-g* --- sanity check:
hostname; echo "job=$SLURM_JOB_ID cpus=$SLURM_CPUS_PER_TASK"; nproc     # expect eu-*, 32, 32

# 2) Outbound HTTPS for the claude CLI via the ETH proxy. eth_proxy.sh self-inits
#    Lmod (the srun shell has no `module` function), sets http(s)_proxy (lower +
#    UPPER) = http://proxy.service.consul:3128, and forwards the APPTAINERENV_
#    twins so in-container R downloads are proxied too.
cd ~/GaussNiche
source sbatch/eth_proxy.sh
curl -sS -o /dev/null -w 'api.anthropic.com -> %{http_code}\n' https://api.anthropic.com   # 401/200 = reachable

# 3) claude is a self-contained ELF at ~/.local/bin/claude (shared $HOME, no node
#    needed). --pty gives it a real TTY, so the TUI works. If PATH is missing it
#    (non-login shell), use the full path or export it:
export PATH="$HOME/.local/bin:$PATH"
claude
```

Inside that `claude`, the SessionStart hook fires and you'll see a **COMPUTE
ALLOCATION DETECTED** block — then just run jobs directly, e.g.

```
.claude/bin/gnrun run_5d_hsm.R smoke          # runs here, now, 32 workers
.claude/bin/gnrun tune_5d_hsm.R full          # the ablation sweep, in-allocation
```

`gnrun` on a compute node is equivalent to the manual form the `sbatch/submit_*.sh`
use (sources `eth_proxy.sh` then `config_r_studio` for spack native libs, pins
BLAS/OMP to 1):

```bash
apptainer exec --bind "/cluster/home/$USER,/cluster/software,/cluster/scratch/$USER" \
  /cluster/scratch/$USER/rocker_rstudio_4.5.sif Rscript <script.R> <mode> <workers>
```

## Why these choices (verified on this cluster)

- **`srun --pty bash`, not `salloc`.** `LaunchParameters` lacks
  `use_interactive_step`, so `salloc` only grabs the allocation and leaves you on
  `eu-login-*` (a second `srun --pty bash` hop would be needed). `srun --pty bash`
  lands you ON the compute node with `SLURM_*` set — which is what the hook and
  `gnrun` key off.
- **Proxy = `proxy.service.consul:3128`.** `module load eth_proxy` resolves to
  this on compute nodes (`proxy.ethz.ch:3128` is 403'd from batch nodes — the
  `config_r_studio` comment saying otherwise is stale). `claude` reads the
  lowercase `https_proxy`; `eth_proxy.sh` sets both cases and its `no_proxy` does
  not cover `api.anthropic.com`, so the API traffic is correctly proxied.
- **Workers get all 32 cores.** With `--ntasks=1 --cpus-per-task=32`, the `claude`
  process, `apptainer`, and every forked furrr/PSOCK `Rscript` worker are children
  of the one step and inherit its cgroup cpuset = all 32 cores. `SLURM_CPUS_PER_TASK`
  is what `gnrun` / the R drivers use as the worker count. Call `apptainer exec`
  **directly** — do NOT wrap it in a nested `srun` (fights for step resources). If
  you ever see everything pinned to one core, add `--cpu-bind=none` to the `srun`.

## Caveats

- **Disconnect kills an interactive job.** Always launch inside `tmux`/`screen`
  on the login node (step 0). Detach (`Ctrl-b d`) before dropping SSH.
- **Don't `sbatch` from inside the allocation.** You're already on compute —
  nesting a batch job double-allocates and queues behind your own job. `gnrun`
  already refuses to (it runs directly when `GN_ON_COMPUTE=1`).
- **Stay within the allocation.** Cap workers at the allocated cpus and keep total
  memory under the allocated GB (OOM kills the whole session, `claude` included).
  Keep `APPTAINERENV_OMP_NUM_THREADS=1` (gnrun/sbatch set this) so workers don't
  oversubscribe.
- **The clock is real.** A 6 h allocation ends in 6 h and takes the `claude`
  session with it. For a run longer than the time left, `sbatch` a separate job
  (e.g. `submit_tune_5d_hsm.sh`) instead. `git commit` your work before the wall.
- **`--mem-per-cpu=6G × 32 = 192 GB`** exceeds the per-core RAM ratio of the
  249 GB/128-core nodes, so SLURM places you on a node with ≥192 GB free — may add
  a little queue time. Lower `--mem-per-cpu` if you want to schedule faster.
- For anything bigger than one node or longer than a session, go back to the
  `sbatch/submit_*.sh` scripts — this skill is for *interactive iteration*, not
  for the final production sweeps.
```
