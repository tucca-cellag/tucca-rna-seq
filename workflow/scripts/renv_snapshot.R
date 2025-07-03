# tucca-rna-seq/workflow/scripts/renv_snapshot.R
#
# This script is executed by the `renv_snapshot` rule in `renv.smk`.
#
# Purpose:
# To create a `renv.lock` file that captures all R package dependencies
# required by the tucca-rna-seq pipeline. `renv::snapshot()` automatically scans
# the R scripts (`.R`, `.Rmd`) in the project directory to find required
# packages and records their versions in `renv.lock`.
#
# How it works:
# This script runs inside the pipeline's dedicated R Conda environment.
# It uses `renv::snapshot(type = "all")` to ensure that it captures not only
# explicit `library()` calls but also packages used in other ways, providing
# a comprehensive lockfile. The `force = TRUE` argument is used to bypass the
# interactive prompt, which is essential for non-interactive execution within
# the Snakemake workflow.
#
# The resulting `renv.lock` file serves as the blueprint for recreating the
# interactive R environment, ensuring it matches the pipeline's environment.

# --- 1. Setup and Logging ---
log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")
date()

# --- 2. Load renv ---
library(renv)

# --- 3. Snapshot Dependencies ---
message("Starting renv snapshot...")

# The `type = "all"` option ensures we capture all packages used, not just
# those explicitly attached with `library()`. `force = TRUE` is critical for
# non-interactive use in the Snakemake pipeline.
renv::snapshot(
  project = getwd(),
  type = "all",
  prompt = FALSE
)

message("renv snapshot completed successfully.")
date()
