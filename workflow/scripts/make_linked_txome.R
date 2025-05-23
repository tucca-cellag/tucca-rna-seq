log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")
date()

# Using the following approach to set a few environment variables because
# Singularity/Apptainer containers on GitHub runners do not offer the necessary
# write permissions for setTximetaBFC() to create
# /home/runner/.config/R/tximeta/bfcloc.json

# 1. Define and set XDG_CACHE_HOME to a writable path
# This is where rappdirs::user_cache_dir() (used by tximeta/BiocFileCache)
# will try to create user-specific config/cache subdirectories like 'tximeta'.
# The problematic 'bfcloc.json' should end up in a subdirectory here if tximeta
# still tries to write it.

# This path is relative to the Snakemake working directory
r_user_cache_base_dir <- file.path(getwd(), ".r_user_cache_for_tximeta")
if (!dir.exists(r_user_cache_base_dir)) {
  dir.create(r_user_cache_base_dir, recursive = TRUE, showWarnings = FALSE)
}
# Set R_USER_CACHE_DIR, which rappdirs should respect.
Sys.setenv(R_USER_CACHE_DIR = r_user_cache_base_dir)
# Also set XDG_CACHE_HOME for broader compatibility,
# R_USER_CACHE_DIR might take precedence for rappdirs in R
Sys.setenv(XDG_CACHE_HOME = r_user_cache_base_dir)

# 2. Define your main tximeta data cache path
# This is where the actual large cache files (e.g., SQLite DBs) will go.
# This was 'resources/tximeta' in your Snakemake rule.
# Get the desired cache path from Snakemake parameters
tximeta_data_cache_path <- snakemake@params[["tximeta_cache"]]

# Ensure the target cache directory exists.
# This path is relative to the Snakemake working directory
# (e.g., "resources/tximeta") and should be writable within containers
if (!dir.exists(tximeta_data_cache_path)) {
  dir.create(tximeta_data_cache_path, recursive = TRUE, showWarnings = FALSE)
}

# 3. Set environment variables for tximeta and BiocFileCache data
Sys.setenv(TXIMETA_HUB_CACHE = tximeta_data_cache_path)
Sys.setenv(BIOCFILECACHE_CACHE = tximeta_data_cache_path)

suppressPackageStartupMessages({
  library(devtools)
})

# --- START DEBUGGING BLOCK ---
# Check what rappdirs resolves to AFTER setting env vars and BEFORE loading tximeta.
# This requires the 'rappdirs' package to be available in the Conda environment.
# If 'rappdirs' is not found, these messages will indicate that, but the script will continue.
message("--- Pre-tximeta Load Debug ---")
message(paste("Timestamp:", Sys.time()))
message(paste("R_USER_CACHE_DIR set to:", Sys.getenv("R_USER_CACHE_DIR")))
message(paste("XDG_CACHE_HOME set to:", Sys.getenv("XDG_CACHE_HOME")))
message(paste("TXIMETA_HUB_CACHE set to:", Sys.getenv("TXIMETA_HUB_CACHE")))
message(paste("BIOCFILECACHE_CACHE set to:", Sys.getenv("BIOCFILECACHE_CACHE")))
if (requireNamespace("rappdirs", quietly = TRUE)) {
  message("rappdirs package is available.")
  resolved_tximeta_config_dir <- tryCatch(
    {
      rappdirs::user_cache_dir("tximeta")
    },
    error = function(e) {
      paste("Error calling rappdirs::user_cache_dir('tximeta'):", e$message)
    }
  )
  message(paste("rappdirs::user_cache_dir('tximeta') resolves to:", resolved_tximeta_config_dir))

  resolved_bfc_config_dir <- tryCatch(
    {
      rappdirs::user_cache_dir("BiocFileCache")
    },
    error = function(e) {
      paste("Error calling rappdirs::user_cache_dir('BiocFileCache'):", e$message)
    }
  )
  message(paste("rappdirs::user_cache_dir('BiocFileCache') resolves to:", resolved_bfc_config_dir))
} else {
  message("rappdirs package NOT available for pre-check. Ensure it's in the Conda env if issues persist.")
}
message("--- End Pre-tximeta Load Debug ---")
# --- END DEBUGGING BLOCK ---

devtools::session_info()

message("Attempting to load tximeta...")
suppressPackageStartupMessages({
  library(tximeta)
})
message("tximeta loaded successfully.")

# Explicitly set the tximeta cache using its function.
# The bfcloc.json will point to tximeta_data_cache_path.
message(paste("Attempting to set tximeta BFC to:", tximeta_data_cache_path))
tximeta::setTximetaBFC(dir = tximeta_data_cache_path, quiet = FALSE)
message("setTximetaBFC called.")

organism_split <- strsplit(snakemake@params[["organism"]], "_")[[1]]
organism_reformat <- paste(paste(organism_split[1], organism_split[2]))

if (snakemake@params[["source"]] %in% c("Ensembl", "GENCODE")) {
  # Enforce creation of a TxDb object for Ensembl and GENCODE when
  # is called makeLinkedTxome
  source <- paste0("Local", snakemake@params[["source"]])
  # TODO: If support for makeLinkedTxome(source = c("Ensembl", "GENCODE")) is
  # added instead of forcing "LocalEnsembl" and "LocalGENCODE" the logic in
  # workflow/scripts/tximeta.R will need to be refactored to work with EnsDb
  # objects
  # See: https://github.com/thelovelab/tximeta/blob/devel/R/tximeta.R
  # Specifically getTxDb() call and definition
} else {
  source <- snakemake@params[["source"]]
}

message("Calling makeLinkedTxome...")
tximeta::makeLinkedTxome(
  # index_dir is a list of files, select the first file's dirname
  indexDir = dirname(snakemake@input[["index_dir"]])[1],
  source = source,
  organism = organism_reformat,
  release = snakemake@params[["release"]],
  genome = snakemake@params[["genome"]],
  fasta = snakemake@input[["fasta"]],
  gtf = snakemake@input[["gtf"]],
  write = TRUE,
  jsonFile = snakemake@output[["jsonFile"]]
)
message("makeLinkedTxome finished.")
