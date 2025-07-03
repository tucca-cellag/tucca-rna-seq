# tucca-rna-seq/workflow/scripts/renv_restore.R
#
# This script is executed by the `renv_restore` rule in `renv.smk`.
#
# Purpose:
# To install all R packages specified in the `renv.lock` file into the
# project's local `renv` library. This step pre-populates the interactive
# R environment, so the user does not have to wait for `renv::restore()` to
# run when they first open RStudio.
#
# How it works:
# The script calls `renv::restore()`, which reads `renv.lock` and installs
# the exact versions of all required packages. The `prompt = FALSE` argument
# is essential for non-interactive execution within the Snakemake workflow,
# as it prevents the script from asking for user confirmation.
#
# By automating this step, the user is provided with a ready-to-use,
# fully provisioned R environment as soon as the pipeline finishes.

# --- 1. Setup and Logging ---
log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")
date()

# --- 2. Load renv ---
library(renv)

# --- 3. Restore Project Library ---
message("Restoring renv library from lockfile...")

# This command installs packages from the lockfile into the project library.
# `prompt = FALSE` is critical for non-interactive use. `clean = TRUE` ensures
# the library is an exact match to the lockfile.
# We use `snakemake@workflow$basedir` to ensure the project path is always
# correct, regardless of the script's execution context.
renv::restore(
  project = base::getwd(),
  lockfile = base::file.path(base::getwd(), "renv.lock")
  prompt = FALSE,
  clean = TRUE
)

message("renv library restored successfully.")
date()
