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
log <- base::file(snakemake@log[[1]], open = "wt")
base::sink(log)
base::sink(log, type = "message")
base::date()

# --- 2. Set Build Environment Variables ---
# This is a robust method to ensure all R packages that compile C/C++ code
# can find the libraries and headers provided by the conda environment. We
# explicitly set the environment variables that the compiler and linker use,
# which is necessary when the job scheduler (e.g. SLURM) sanitizes them.
conda_prefix <- base::Sys.getenv("CONDA_PREFIX")
if (nchar(conda_prefix) > 0) {
  base::Sys.setenv(
    PKG_CONFIG_PATH = base::file.path(conda_prefix, "lib", "pkgconfig"),
    CPPFLAGS = paste0("-I", base::file.path(conda_prefix, "include")),
    LDFLAGS = paste0("-L", base::file.path(conda_prefix, "lib"))
  )
}

# --- 3. Load renv and Restore ---
base::message("Restoring renv library from lockfile...")
base::library(renv)
renv::restore(
  project = base::getwd(),
  prompt = FALSE,
  clean = TRUE
)

base::message("renv library restored successfully.")
base::date()
