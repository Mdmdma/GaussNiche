---
name: euler-rstudio-server
description: Operate inside an RStudio Server session on the ETH Euler HPC cluster — sessions are launched by JupyterHub onto a SLURM-allocated compute node and run R inside a rocker/rstudio apptainer container. Use when the user mentions Euler, RStudio Server, JupyterHub, the rocker container, compute-node networking (ETH proxy), or when the working hostname looks like `eu-a*-NNN` / `eu-g*-NNN` (compute node) rather than `eu-login-*`. Pairs with the [[euler-r-spack-setup]] skill, which handles installing R packages that need native libs.
---

# Operating inside RStudio Server on ETH Euler

This skill is the **runtime/environmental** companion to [[euler-r-spack-setup]]
(which is the install-time recipe). It tells the agent what the RStudio
Server session actually is, what restrictions apply, and how to use it
correctly without falling back on either "login node" rules (too strict) or
"random Linux box" rules (too permissive).

## Am I in one?

Check the hostname:

```bash
uname -n
```

- `eu-login-NN` — **login node**. No heavy compute. Submit jobs with
  `sbatch` / `srun`. Pure-R unit tests (seconds, no real data) are OK.
- `eu-a*-NNN`, `eu-g*-NNN`, `eu-c*-NNN`, anything else — **compute node**.
  Either: (a) you SSH'd into an `srun` allocation, or (b) you are inside
  the terminal pane of a JupyterHub-launched RStudio Server session.
  Either way: this node is **your** SLURM allocation, with the resources
  you asked for at launch. You can use those resources fully.

If you don't know, **assume login node** and submit. The cost of a
needless `sbatch` is zero; clobbering a login node is real (see the
machine-wide CLAUDE.md at `~/.claude/CLAUDE.md`).

## What an Euler RStudio Server session actually is

1. **JupyterHub** at `https://jupyter.euler.hpc.ethz.ch/` is the entry
   point. The user picks "RStudio" as the software, requests CPUs/mem/time
   (and optionally GPUs), and clicks start.
2. JupyterHub issues an **`sbatch`** under the hood. SLURM places the job
   on a compute node — that is where the rest happens.
3. On the compute node, JupyterHub:
   - sources `~/.config/euler/jupyterhub/config_r_studio` (a bash script
     the user owns) to set up modules, env vars, and proxy settings;
   - runs `rserver.sh`, which launches an **apptainer** container
     (image: `docker://rocker/rstudio:${RSTUDIO_TAG}`, e.g. `4.5`);
   - inside that container, `rserver` (the RStudio Server daemon) starts
     and serves on a private port, proxied back through JupyterHub to
     the user's browser.
4. Random auto-generated credentials are at `~/.rstudio/.password`.

Consequences the agent must know:

- The R process you talk to runs **inside an apptainer container**.
  Host paths under `/cluster/home/$USER`, `/cluster/scratch/$USER`,
  `/cluster/project/*` are bind-mounted in, but `/usr/local/...`,
  `/lib/...`, `/etc/...` come from the container's debian/ubuntu base, not
  from Euler's host OS.
- The container has **its own R library** at
  `/usr/local/lib/R/site-library`, which is a **64 MB writable overlay**
  that fills fast. **Never install packages there.** See [[euler-r-spack-setup]].
- Host modules (`module load proj/9.2.1`) only take effect for processes
  spawned from `config_r_studio` and only if `APPTAINERENV_LD_LIBRARY_PATH`
  / `APPTAINERENV_LD_PRELOAD` are exported (apptainer wipes them
  otherwise). The user's existing `config_r_studio` already does this.
- The terminal pane in RStudio runs **outside** the apptainer container,
  on the bare compute node — so a `module load r/4.5.1` there gives you a
  different R than the one served by RStudio (host R vs. container R).
  Don't confuse the two when debugging "library not found" errors.

## Resource policy on a compute node

You hold a SLURM allocation. Within it, **use what you asked for**:

- Multi-core builds, large compiles, real model runs, full pipelines on
  real data — all fine.
- Don't *exceed* the allocation (memory OOM kills the session; CPU
  oversubscription just slows it down).
- If the user wants to do something bigger than the live allocation,
  recommend they re-launch the JupyterHub session with more resources, or
  submit a separate `sbatch` job (the latter goes to a different node and
  can be much larger).
- Scratch — `/cluster/scratch/$USER` — for all bulk data, model outputs,
  rendered plots, intermediate `.rds` files. Auto-purged ~14 days; not
  backed up. `$HOME` has a tight quota — fill it and tooling breaks.

## Networking from a compute node

Compute nodes have **no direct outbound internet**. Two ways to escape:

1. **`module load eth_proxy`** — sets `http_proxy` / `https_proxy` /
   `HTTP_PROXY` / `HTTPS_PROXY` to `http://proxy.ethz.ch:3128`. This is
   what the user's `config_r_studio` already does, so any process that
   inherits the JupyterHub environment (including the in-container R) is
   already proxied. Verify with:

   ```bash
   echo "$http_proxy $https_proxy"
   ```

2. **Manually in R** if the env var didn't propagate:

   ```r
   Sys.setenv(http_proxy  = "http://proxy.ethz.ch:3128")
   Sys.setenv(https_proxy = "http://proxy.ethz.ch:3128")
   ```

   (Older docs use `proxy.service.consul:3128`; both resolve to the same
   service. Prefer `proxy.ethz.ch:3128` — it's the documented public name.)

Anything that talks to the internet is affected: `install.packages()`,
`geodata::worldclim_country()` and friends, `git clone` from external
hosts (GitHub is reachable), `WebFetch`/`WebSearch` in this very agent.

Be considerate: the ETH proxy has limited bandwidth and rate limits.
Don't fan-out hundreds of concurrent downloads. If a download throws
`Could not resolve host` or a 403, the proxy is the first suspect.

## Sharing files with the host (the RStudio terminal vs. R)

The RStudio R session sees:

- `~/...` — the host `$HOME`, mounted in.
- `/cluster/scratch/...` — the host scratch, mounted in.
- `/cluster/project/...` — group/project storage (if the user has any).

So `setwd("~/GaussNiche")` from RStudio R behaves the same as
`cd ~/GaussNiche` from the host terminal. **But do not assume the same
binaries are available** — the container's `Sys.which("gcc")` returns
the container's gcc, not the spack one. If a build needs spack tooling,
see the cc_wrap/ wrappers in [[euler-r-spack-setup]].

## Common failure modes (in roughly the order you'll hit them)

### `library(terra)` fails with `libproj.so.25: cannot open shared object`

The apptainer container wiped `LD_LIBRARY_PATH` from the
`module load proj` you did outside it. Confirm by running, **inside an
RStudio terminal pane** (i.e. inside the container):

```bash
echo "$LD_LIBRARY_PATH" | tr ':' '\n' | grep -E 'proj|gdal' || echo "NO SPACK PATHS"
```

If `NO SPACK PATHS`: the `APPTAINERENV_LD_LIBRARY_PATH` line is missing
from `config_r_studio`. See [[euler-r-spack-setup]] §B0.

### `install.packages()` hangs or `Could not resolve host cloud.r-project.org`

Proxy is not set inside the R session. Either `module load eth_proxy`
didn't run (check `config_r_studio`), or apptainer dropped the env vars
on container entry. Fix: `Sys.setenv(http_proxy=..., https_proxy=...)`
inside R, or re-launch JupyterHub after editing `config_r_studio`.

### Package "installed successfully" but `library(...)` fails with `readRDS: unknown type 0`

The 64 MB container overlay at `/usr/local/lib/R/site-library` filled
during install and truncated the package's lazyload DB. Reinstall with
`lib = "~/R/rocker-rstudio/4.5"`. See [[euler-r-spack-setup]] §A.

### Edits to `config_r_studio` "don't take effect"

That file is sourced **at JupyterHub session start**. The user must end
the current session and start a fresh one from `jupyter.euler.hpc.ethz.ch`.
Restarting R within RStudio is not sufficient.

### `git push` works but `install.packages` doesn't (or vice versa)

`git`/`ssh` to GitHub goes over port 22, which is allowed without the
proxy. `install.packages()` goes over HTTPS, which needs the proxy.
They have independent failure modes — don't infer one from the other.

## Submitting an external SLURM job from inside an RStudio session

Even with an interactive session, sometimes you want to fire off a
heavy job to a different node (more cores, GPU, longer runtime). The
RStudio terminal pane has the full SLURM CLI:

```bash
sbatch --time=04:00:00 --mem-per-cpu=4G --cpus-per-task=8 my_script.sh
squeue -u $USER
sacct -j <jobid>
```

This is the standard escape valve for "I need more than this session
gives me." Don't try to scale up by running heavier inside the current
session — its limits are fixed at launch.

## Sources

- `~/.config/euler/jupyterhub/config_r_studio` (user's actual setup)
- ETH SciComp wiki: JupyterHub, Accessing the clusters, Using the batch system
- ETH HPC docs: JupyterHub services, firewall/proxy
- Cross-references: [[euler-r-spack-setup]] for package compile/install
  recipes; machine-wide `~/.claude/CLAUDE.md` for login-node rules.
