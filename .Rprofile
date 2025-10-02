# This file configures the R session when the project is opened.
#
# The `if (interactive())` check is crucial for our hybrid Conda/renv setup.
#
# - For INTERACTIVE sessions (like RStudio), it activates `renv`, giving the
#   user access to the project's pre-built package library for analysis.
#
# - For NON-INTERACTIVE sessions (like scripts run by the Snakemake pipeline),
#   it does nothing. This allows the script to use the isolated Conda
#   environment provided by Snakemake, which is essential for pipeline
#   reproducibility.
#
# This logic prevents the user's interactive environment from interfering with
# the pipeline's execution environment.
if (interactive()) {
  if (file.exists("renv/activate.R")) {
    source("renv/activate.R")
  }
}
