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

# --- 2. Diagnostic Information ---
base::message("--- STARTING BUILD DIAGNOSTICS ---")

base::message("Querying environment variables...")
env_vars <- c(
  "PATH", "LD_LIBRARY_PATH", "PKG_CONFIG_PATH", "CFLAGS", "LDFLAGS",
  "CPPFLAGS", "CONDA_PREFIX"
)
for (v in env_vars) {
  base::message(paste0(v, ": ", base::Sys.getenv(v)))
}

base::message("\nRunning build tool diagnostics...")
base::message("Output of: pkg-config --cflags --libs libxml-2.0")
try(base::system("pkg-config --cflags --libs libxml-2.0"))

base::message("\nOutput of: xml2-config --cflags --libs")
try(base::system("xml2-config --cflags --libs"))

conda_prefix <- base::Sys.getenv("CONDA_PREFIX")
if (nchar(conda_prefix) > 0) {
  base::message("\nChecking for libxml2 in $CONDA_PREFIX/lib...")
  try(base::system(paste("ls -lh", base::shQuote(base::file.path(conda_prefix, "lib")), "| grep xml")))

  base::message("\nChecking for xml2-config in $CONDA_PREFIX/bin...")
  try(base::system(paste("ls -lh", base::shQuote(base::file.path(conda_prefix, "bin")), "| grep xml")))
}

base::message("--- END OF BUILD DIAGNOSTICS ---")

# --- 3. Load renv ---
base::library(renv)

# --- 4. Restore Project Library ---
base::message("Restoring renv library from lockfile...")

# Set configure arguments for packages that have trouble finding libraries
# in a conda environment. This explicitly tells the build system where to look.
lib_path <- base::file.path(base::Sys.getenv("CONDA_PREFIX"), "lib")
include_path <- base::file.path(base::Sys.getenv("CONDA_PREFIX"), "include")

base::options(
  configure.args = base::c(
    XML = base::paste0(
      "--with-xml-config=", base::file.path(base::Sys.getenv("CONDA_PREFIX"), "bin", "xml2-config"),
      " --with-libxml-include=", include_path,
      " --with-libxml-lib=", lib_path
    ),
    curl = base::paste0(
      "--with-curl-config=", base::file.path(base::Sys.getenv("CONDA_PREFIX"), "bin", "curl-config")
    ),
    openssl = base::paste0(
      "--with-ssl-include=", include_path,
      " --with-ssl-lib=", lib_path
    )
  )
)

# This command installs packages from the lockfile into the project library.
# `prompt = FALSE` is critical for non-interactive use. `clean = TRUE` ensures
# the library is an exact match to the lockfile.
# We use `snakemake@workflow$basedir` to ensure the project path is always
# correct, regardless of the script's execution context.
renv::restore(
  project = base::getwd(),
  lockfile = base::file.path(base::getwd(), "renv.lock"),
  prompt = FALSE,
  clean = TRUE
)

base::message("renv library restored successfully.")
base::date()
