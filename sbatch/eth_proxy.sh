# Enable outbound internet from a SLURM compute node via the ETH proxy, and make it
# reachable from inside apptainer. SOURCE this in an sbatch script BEFORE any
# `apptainer exec` that needs the network (install.packages, geodata downloads, R CMD
# check --as-cran URL checks, git clone from external hosts).
#
# Why this and not a hardcoded proxy: `module load eth_proxy` is the ETH-maintained
# mechanism and currently resolves to http://proxy.service.consul:3128 (the public
# name proxy.ethz.ch:3128 is refused with a 403 from batch compute nodes). A
# non-interactive sbatch shell has no `module` function until the Lmod init is sourced.
# apptainer only forwards host env vars that are prefixed APPTAINERENV_, so we export
# the proxy under those names too. Verified: host curl + in-container R both reach CRAN.

# 1) make `module` available in this non-interactive shell, then load the proxy module
if [ -n "${MODULESHOME:-}" ] && [ -f "$MODULESHOME/init/bash" ]; then
  # shellcheck disable=SC1090
  . "$MODULESHOME/init/bash"
fi
module load eth_proxy 2>/dev/null || true

# 2) fall back to the known proxy if the module was unavailable
: "${http_proxy:=http://proxy.service.consul:3128}"
: "${https_proxy:=$http_proxy}"
export http_proxy https_proxy
export HTTP_PROXY="$http_proxy" HTTPS_PROXY="$https_proxy"

# 3) forward the proxy INTO apptainer (only APPTAINERENV_-prefixed vars are passed through)
export APPTAINERENV_http_proxy="$http_proxy"   APPTAINERENV_https_proxy="$https_proxy"
export APPTAINERENV_HTTP_PROXY="$http_proxy"   APPTAINERENV_HTTPS_PROXY="$https_proxy"
export APPTAINERENV_no_proxy="${no_proxy:-localhost,127.0.0.1}"
export APPTAINERENV_NO_PROXY="$APPTAINERENV_no_proxy"

echo "[eth_proxy] http_proxy=$http_proxy (forwarded into apptainer)"
