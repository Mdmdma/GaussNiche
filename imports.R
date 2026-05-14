# Smallest possible test: load every package 4_multi_species_comparison.R
# uses, in the same order. No computation, no side effects beyond library().
# Run as the verification step before dispatching anything to SLURM.

pkgs <- c("terra", "USE.MCMC", "mvtnorm", "MASS", "ggplot2", "viridis",
          "hypervolume", "patchwork", "sf", "tictoc", "future", "furrr")

failed <- character()
for (p in pkgs) {
  cat(sprintf("%-12s ", p))
  ok <- tryCatch({
    suppressMessages(library(p, character.only = TRUE))
    cat("OK", as.character(packageVersion(p)), "\n")
    TRUE
  }, error = function(e) {
    cat("FAIL:", conditionMessage(e), "\n")
    FALSE
  })
  if (!ok) failed <- c(failed, p)
}

cat("\nlibPaths:\n"); print(.libPaths())
if (length(failed)) {
  stop(sprintf("FAILED packages: %s", paste(failed, collapse = ", ")))
} else {
  cat("\nALL OK\n")
}
