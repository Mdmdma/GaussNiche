---
name: euler-r-spack-setup
description: Install or repair R packages on the ETH Euler RStudio Server stack (rocker/rstudio:4.5 in apptainer with spack-provided GDAL/GEOS/PROJ etc.) — including USE.MCMC + GaussNiche, sf, terra, hypervolume, devtools. Use when the user reports R install failures from RStudio on Euler ("library(terra) not found", "configure: error: GDALAllRegister", "GNU MP not found", "hb-ft.h", "krb5int_strlcpy", "no package called 'X'"), wants to add R packages requiring native libs from spack, or needs the install_block.R / GaussNiche stack working.
---

# Euler RStudio + spack R-package setup

This skill encodes the working recipe (and pitfalls) for installing source-built R packages inside the **rocker/rstudio:4.5 apptainer container** that ETH Euler launches via JupyterHub. The native libraries (GDAL, GEOS, PROJ, krb5, fontconfig, harfbuzz, ...) come from the cluster's **spack stack 2024-06**, not from apt.

## When this is relevant

- The user is on Euler and `~/.config/euler/jupyterhub/config_r_studio` exists.
- They want to install R packages that need to compile C/C++ code linking against GDAL / sf / terra / hypervolume / devtools / fs / fontconfig / harfbuzz / etc.
- Or `library(...)` fails in RStudio with `there is no package called 'X'` after a previous install attempt.

If those conditions don't hold (e.g. user is on a generic Linux box, or only needs pure-R packages), most of this is overkill — fall back to plain `install.packages()`.

## Persistent state from the previous install

Two categories, by lifetime:

**Runtime — under `$HOME`, persistent forever:**
- `~/R/rocker-rstudio/4.5/` — compiled R packages (terra, sf, USE.MCMC, hypervolume, …).
- `~/R/lib_shim/` — single-file shim dir with `libstdc++.so.6` and `libkrb5support.so.0` symlinks (see gotchas B-libstdc + B). Loaded by `~/.config/euler/jupyterhub/config_r_studio` on every JupyterHub session start. Lose this dir → `library(Rcpp)` and PSOCK workers break.
- `~/.config/euler/jupyterhub/config_r_studio` — the bash script JupyterHub sources before launching rserver; sets all the module loads + the apptainer LD_LIBRARY_PATH / LD_PRELOAD / R_LD_LIBRARY_PATH passthroughs.
- `~/.Rprofile` — prepends `~/R/rocker-rstudio/4.5` to `.libPaths()`.

**Compile-time scratch — `/cluster/scratch/merler/USE.MCMC_setup/`, auto-purged ~14 days idle:**
- `spack_env.sh` — full `PATH`/`LD_LIBRARY_PATH`/`PKG_CONFIG_PATH`/`CPATH`/`LIBRARY_PATH` setup
- `cc_wrap/` — `gcc`, `g++` wrappers + `pkg-config` symlink
- `config_r_studio.before_2026-05-12` / `config_r_studio.after_2026-05-12` — old/new versions
- `Rprofile` — copy of `~/.Rprofile`
- `install_logs/` — all `block_*.log` from the original session
- `libstdcxx_shim/` — **legacy** location of the shim dir, superseded by `~/R/lib_shim/`. Safe to delete after confirming `~/R/lib_shim/` is on `config_r_studio`'s `LIBSTDCXX_SHIM` line.

Compile-time scratch is only needed when **adding or rebuilding** R packages with native deps. For runtime (loading already-installed packages), only the `$HOME` items matter. If scratch gets purged, recreate `spack_env.sh` + `cc_wrap/` only if you need to install something new.

## Where R packages live

- Library: `~/R/rocker-rstudio/4.5/` (this is `R_LIBS_USER` inside the rocker container).
- Always install with `lib = "~/R/rocker-rstudio/4.5"` and prepend it to `.libPaths()` because `R_LIBS` is set in the container env (overrides default `R_LIBS_USER` placement).
- **Never install into `/usr/local/lib/R/site-library`** inside the container — that's a writable overlay capped at **64 MB** and fills up after ~3 packages, silently truncating files (e.g. terra's lazyload DB → `readRDS: unknown type 0`).

## Standard install sequence

### 1. Source the env (compile-time)

```bash
source /cluster/scratch/merler/USE.MCMC_setup/spack_env.sh   # or /tmp/spack_env.sh
export PATH=/cluster/scratch/merler/USE.MCMC_setup/cc_wrap:$PATH
```

The env script exports paths for GDAL/GEOS/PROJ/sqlite/udunits/openssl/krb5/curl/json-c/hdf5/libjpeg/libpng/libtiff/libgeotiff/netcdf-c/openjpeg/zstd/xz/libxml2/libiconv/zlib-ng/abseil-cpp/cmake **plus** the libs needed for the wider package set: gmp/fontconfig/freetype/harfbuzz/fribidi/libwebp/glib/bzip2/graphite2/pcre2/util-linux-uuid.

The wrappers inject `-L $KRB5/lib -Wl,-rpath-link=$KRB5/lib` into every gcc/g++ invocation; the `pkg-config` symlink points at spack `pkgconf` (the rocker container has no pkg-config).

It also exports `USE_BUNDLED_LIBUV=1` so `fs` builds its bundled libuv (no system libuv-dev).

### 2. Run R install with explicit lib

```r
.libPaths(c("~/R/rocker-rstudio/4.5", .libPaths()))
install.packages(<pkgs>,
  lib    = "~/R/rocker-rstudio/4.5",
  repos  = "https://cloud.r-project.org",
  Ncpus  = parallel::detectCores())
```

Run R via `nohup R --no-save --no-restore -e '...' > LOG 2>&1 &` and tail the log. Heavy compiles (terra, sf, s2 with vendored abseil, stringi with ICU) take 5-10 min each.

### 3. Make it persistent for future RStudio sessions

Already done in this user's setup (see `/cluster/scratch/merler/USE.MCMC_setup/`):

- **`~/.Rprofile`** — prepends `~/R/rocker-rstudio/4.5` to `.libPaths()`:
  ```r
  local({
    user_lib <- "~/R/rocker-rstudio/4.5"
    if (dir.exists(user_lib) && !(normalizePath(user_lib) %in% normalizePath(.libPaths()))) {
      .libPaths(c(user_lib, .libPaths()))
    }
  })
  ```
- **`~/.config/euler/jupyterhub/config_r_studio`** — adds these `module load` lines after the existing terra/sf set:
  ```
  module load gmp/6.2.1 \
              fontconfig/2.13.1-b4vx4ls freetype/2.11.1-fy5fkou \
              harfbuzz/7.3.0-v2ihsoq fribidi/1.0.12-f3zew33 \
              libwebp/1.2.3
  module load krb5/1.20.1-7iahlck
  export LD_PRELOAD="$(dirname "$(readlink -f "$(command -v krb5-config)")")/../lib/libkrb5support.so.0${LD_PRELOAD:+:$LD_PRELOAD}"
  ```
  The `LD_PRELOAD` is essential — see krb5 gotcha below.

After editing config_r_studio, the user must restart the **JupyterHub** session (not just R) to pick it up.

## Gotchas to know about (in roughly the order you'll hit them)

### A. Container `/usr/local/lib/R/site-library` is a 64 MB overlay
- It's writable but tiny. R will install there by default (R_LIBS sets it as first writable path).
- When it fills, install partially writes files (e.g. `Lazyload.rds` truncated → `readRDS unknown type 0` when loading the package later).
- Symptom: `Error in readRDS(mapfile)` from a package that "installed successfully" earlier.
- Fix: pass `lib = "~/R/rocker-rstudio/4.5"` to `install.packages()`, and if anything is already in site-library, `mv` it to the user lib.

### B0. apptainer drops host `LD_LIBRARY_PATH` and `LD_PRELOAD`
- `module load proj/9.2.1` (and the rest) only sets these in the **outer** shell. The rocker-rstudio apptainer container that hosts RStudio Server wipes them at launch — inside, you see only the rocker defaults (`/usr/local/lib/R/lib`, `/.singularity.d/libs`, `/usr/lib/x86_64-linux-gnu`, ...).
- Symptom: `library(terra)` → `unable to load shared object .../terra.so: libproj.so.25: cannot open shared object file`. Same pattern for libgdal, libgeos, libharfbuzz, anything spack-provided.
- Fix: at the **bottom** of `config_r_studio` (after every `module load` has run), export the apptainer pass-through twins:
  ```bash
  export APPTAINERENV_LD_LIBRARY_PATH="$LD_LIBRARY_PATH"
  export SINGULARITYENV_LD_LIBRARY_PATH="$LD_LIBRARY_PATH"
  export APPTAINERENV_LD_PRELOAD="$LD_PRELOAD"
  export SINGULARITYENV_LD_PRELOAD="$LD_PRELOAD"
  ```
  apptainer translates `APPTAINERENV_X` to `X` inside the container. Both names exported because apptainer accepts either spelling.

### B-libstdc. spack gcc-12.2.0 libstdc++ shadows the container's gcc-13 libstdc++
- Packages built inside the rocker/rstudio:4.5 container link against **gcc 13.3.0**'s libstdc++ (`libstdc++.so.6.0.33`, exports up to `GLIBCXX_3.4.33`). The whole Rcpp ecosystem — Rcpp, terra, sf, hypervolume, USE.MCMC — has `GLIBCXX_3.4.32` in its NEEDED symbols.
- `module load stack/2024-06 r/4.5.1` (and most other modules) prepend **spack gcc-12.2.0**'s runtime path to `LD_LIBRARY_PATH`. That path holds `libstdc++.so.6.0.30`, max `GLIBCXX_3.4.30`. Spack stack 2024-06 ships no newer gcc.
- `rserver` itself is a C++ binary that links against `libstdc++.so.6`. At session start the dynamic linker resolves it from `LD_LIBRARY_PATH` → picks the spack copy → the forked `rsession` (the R interpreter) inherits the wrong libstdc++ in its address space.
- Symptom (in RStudio R console, fresh session): `library(Rcpp)` / `library(terra)` / `library(sf)` →
  ```
  unable to load shared object '.../Rcpp/libs/Rcpp.so':
    /cluster/.../gcc-12.2.0-bj2twcn.../lib64/libstdc++.so.6: version `GLIBCXX_3.4.32' not found
  ```
- Confirm: `grep libstdc /proc/$(pgrep -u $USER rsession | head -1)/maps` shows `gcc-12.2.0-bj2twcn.../libstdc++.so.6.0.30`.
- **Fix: a single-file shim directory on `APPTAINERENV_LD_LIBRARY_PATH`.** `LD_PRELOAD` looks tempting but does **not** work for `rsession`: it has `NoNewPrivs: 1` set (check `grep NoNewPrivs /proc/<rsession_pid>/status`), which puts the dynamic linker into secure-execution mode and silently drops any `LD_PRELOAD` entry with slashes. `LD_LIBRARY_PATH` is kept in secure mode. Recipe:
  ```bash
  # one-time, on the host (outside the container) — keep under $HOME so
  # scratch auto-purge doesn't blow it away
  mkdir -p ~/R/lib_shim
  ln -sfn /usr/lib/x86_64-linux-gnu/libstdc++.so.6.0.33 ~/R/lib_shim/libstdc++.so.6
  ```
  The symlink target `/usr/lib/x86_64-linux-gnu/libstdc++.so.6.0.33` exists only *inside* the rocker container; the symlink stores the target as a string and ld.so resolves it at runtime, so this works fine.
  ```bash
  # in config_r_studio, AFTER all module loads
  LIBSTDCXX_SHIM=~/R/lib_shim
  export APPTAINERENV_LD_LIBRARY_PATH="${LIBSTDCXX_SHIM}:${LD_LIBRARY_PATH}"
  export SINGULARITYENV_LD_LIBRARY_PATH="${LIBSTDCXX_SHIM}:${LD_LIBRARY_PATH}"
  ```
  This handles the main rsession process. **It does NOT handle PSOCK parallel workers** spawned by `parallel::makeCluster()` — those need the extra fix in gotcha B-psock below.
  Why a single-file shim and not just prepending `/usr/lib/x86_64-linux-gnu` directly: that dir also holds older `libcurl.so.4`, `libssl.so.3`, `libxml2.so.2`, etc. that would shadow the spack copies the R packages were compiled against. The shim isolates libstdc++ alone.
  Editing `config_r_studio` requires a fresh **JupyterHub** session — restarting R inside RStudio is insufficient because `rserver`'s libstdc++ is locked in at process start.

### B. krb5 versioned symbol mismatch
- The system `/usr/lib/x86_64-linux-gnu/libkrb5support.so.0` is missing the `krb5int_strlcpy@krb5support_0_MIT` symbol that the spack `libcom_err.so.3` (pulled in transitively by libcurl → libgdal) needs.
- At configure-time: `/usr/bin/ld: ... libcom_err.so.3: undefined reference to 'krb5int_strlcpy@krb5support_0_MIT'`
- At runtime: same message but `undefined symbol`, often surfaced through `parallel::makeCluster` PSOCK workers when they run `library(sf)` (e.g. `USE.MCMC::optimRes()` spins up 5 workers and they all fail). The main rsession may have already imported sf fine — the failure is worker-specific.
- Fix at compile-time: the gcc wrappers add `-Wl,-rpath-link=$KRB5/lib`.
- Fix at runtime — **same shim dir as gotcha B-libstdc**, because in `rsession` `LD_PRELOAD` is silently stripped (NoNewPrivs=1 → secure mode → entries with slashes dropped). Add a symlink for `libkrb5support.so.0` so it wins via the shim:
  ```bash
  KRB5_LIB=/cluster/software/stacks/2024-06/spack/opt/spack/linux-ubuntu22.04-x86_64_v3/gcc-12.2.0/krb5-1.20.1-7iahlcki6fuy6bcolamzqtojxptqlian/lib
  ln -sfn "$KRB5_LIB/libkrb5support.so.0" ~/R/lib_shim/libkrb5support.so.0
  ```
  The shim dir is already first on `APPTAINERENV_LD_LIBRARY_PATH` after gotcha B-libstdc was set up. The legacy `LD_PRELOAD=$KRB5/lib/libkrb5support.so.0` line in `config_r_studio` is a no-op inside rsession but harmless — keep it for the rserver parent and for host-side terminal sessions where LD_PRELOAD does work.
  Adding the symlink takes effect **without restarting JupyterHub** for the main rsession process, but **not** for already-running PSOCK workers (they cache the resolution). For the optimRes / parallel workflow, also apply gotcha B-psock and restart once.

### B-psock. PSOCK parallel workers spawn `Rscript`, which re-prepends `/usr/lib/x86_64-linux-gnu` to LD_LIBRARY_PATH
- Symptom: `library(sf)` works fine in the main R session, but `USE.MCMC::optimRes()` (or anything that calls `parallel::makeCluster(cr, type="PSOCK")`) crashes with:
  ```
  Error in checkForRemoteErrors(val) :
    5 nodes produced errors; first error: unable to load shared object .../sf.so:
    .../krb5-.../lib/libcom_err.so.3: undefined symbol: krb5int_strlcpy, version krb5support_0_MIT
  ```
- Diagnosis: PSOCK workers are spawned by `parallel:::.slaveRSOCK` via `system("Rscript ...")`. `Rscript`'s startup wrapper sources `${R_HOME}/etc/ldpaths`, which **prepends** `R_HOME/lib:/usr/local/lib:/usr/lib/x86_64-linux-gnu:JAVA_HOME/lib/server` to `LD_LIBRARY_PATH`. That shoves our `~/R/lib_shim` from position 2 down to position ~9 — behind `/usr/lib/x86_64-linux-gnu`, which has the broken system `libkrb5support.so.0.1` (no `krb5int_strlcpy@krb5support_0_MIT`). The system copy then "satisfies" libcom_err.so.3's NEEDED before the dynamic linker ever sees the shim.
- The main rsession R process doesn't hit this because `rserver` embeds libR directly without going through R's startup wrapper.
- Fix: pre-export `R_LD_LIBRARY_PATH` with the shim in front. `etc/ldpaths` uses `: ${R_LD_LIBRARY_PATH=...}`, which leaves an already-set value alone:
  ```bash
  # in config_r_studio, alongside the LD_LIBRARY_PATH passthroughs
  export APPTAINERENV_R_LD_LIBRARY_PATH="${LIBSTDCXX_SHIM}:/usr/local/lib/R/lib:/usr/local/lib:/usr/lib/x86_64-linux-gnu"
  export SINGULARITYENV_R_LD_LIBRARY_PATH="${LIBSTDCXX_SHIM}:/usr/local/lib/R/lib:/usr/local/lib:/usr/lib/x86_64-linux-gnu"
  ```
  After a fresh JupyterHub session, PSOCK workers inherit `R_LD_LIBRARY_PATH` with the shim first, so the final worker `LD_LIBRARY_PATH` is `~/R/lib_shim:/usr/local/lib/R/lib:/usr/local/lib:/usr/lib/x86_64-linux-gnu:<JAVA>:<rest>` — shim wins for libstdc++ and libkrb5support.
- Verify: in RStudio R console, after JupyterHub restart, run
  ```r
  cl <- parallel::makeCluster(2L, type="PSOCK")
  parallel::clusterEvalQ(cl, {
    library(sf)
    grep("libkrb5support|libstdc", readLines("/proc/self/maps"), value=TRUE)[c(1,5)]
  })
  parallel::stopCluster(cl)
  ```
  Expect paths under `~/R/lib_shim/...` and `/cluster/.../krb5-.../lib/...` — not `/usr/lib/x86_64-linux-gnu/...`.

### C. s2 build needs `pkg-config`
- The rocker container has **no** pkg-config. s2's configure script (and several others — ragg, textshaping, systemfonts) fail with `pkg-config: not found`.
- Fix: symlink spack `pkgconf` as `pkg-config` and put it on PATH (the `cc_wrap/` dir does this).

### D. `fs` needs libuv or USE_BUNDLED_LIBUV=1
- Container has no libuv-dev; spack has none either.
- Fix: `export USE_BUNDLED_LIBUV=1` (already in `spack_env.sh`). `fs` then builds the bundled libuv source.

### E. `rcdd` (hypervolume → geometry → rcdd) needs GNU MP
- `configure: error: GNU MP not found`. Add `gmp/6.2.1` from spack — already in `spack_env.sh`.

### F. `systemfonts` / `ragg` / `textshaping` need a chain of pkg-config-discoverable libs
- `fontconfig/2.13.1`, `freetype/2.11.1`, `harfbuzz/7.3.0`, `fribidi/1.0.12`, `libwebp/1.2.3` (the 1.2.3 install is the one with `libwebpmux.pc`; 1.2.4 doesn't have it).
- `fontconfig.pc` requires `uuid.pc` (provided by `util-linux-uuid/2.38.1`).
- All wired up in `spack_env.sh`.

### G. `devtools::install_local()` needs `remotes`
- After devtools 2.5.0 the function dispatches to remotes but doesn't install it transitively. If you skipped block (c)'s remotes install, manually `install.packages("remotes")` first.

### H. terra corruption after overlay-fill
- If terra "installed" while the overlay was full, its `Lazyload.rds` is truncated and `library(terra)` fails with `readRDS unknown type 0`. Reinstall terra cleanly into the user lib. Anything depending on terra (raster, geodata, tidyterra) will also need reinstalling.

## Verification commands

```bash
source /cluster/scratch/merler/USE.MCMC_setup/spack_env.sh
export PATH=/cluster/scratch/merler/USE.MCMC_setup/cc_wrap:$PATH
R --no-save -e '
  .libPaths(c("~/R/rocker-rstudio/4.5", .libPaths()))
  for (p in c("terra","sf","s2","USE.MCMC","tidyterra","geodata","hypervolume","devtools","ragg")) {
    cat(sprintf("%-12s ", p))
    tryCatch({ suppressMessages(library(p, character.only = TRUE)); cat("OK\n") },
             error = function(e) cat("FAIL:", conditionMessage(e), "\n"))
  }
'
```

For the parallel/MCMC code path that surfaces the krb5 issue, source `~/GaussNiche/2_testing_wrapper_function.R` (uses `library(USE.MCMC)`, not `library(USE)`).

## Editing the install_block.R script

`~/install_block.R` is the original recipe but doesn't include the runtime/build-time gotchas above. If asked to "re-run install_block.R", source `spack_env.sh` and put `cc_wrap/` on PATH FIRST, and override its `install.packages()` calls to set `lib = "~/R/rocker-rstudio/4.5"`. Otherwise the overlay fills and installs corrupt silently.

## Out of scope

- This skill is for the in-container RStudio R 4.5 stack. **Not** for the host-side `module load r/4.5.1` (terminal/sbatch). Those use a totally different lib path and don't go through the apptainer overlay.
- This skill is **not** for installing system libraries — it just wires up existing spack-installed ones.
- Don't `apt install` inside the container — overlay full + non-persistent.
