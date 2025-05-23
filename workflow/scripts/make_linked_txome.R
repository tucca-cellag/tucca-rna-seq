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
workflow_xdg_cache_dir <- file.path(getwd(), ".workflow_xdg_cache")
if (!dir.exists(workflow_xdg_cache_dir)) {
  dir.create(workflow_xdg_cache_dir, recursive = TRUE, showWarnings = FALSE)
}
Sys.setenv(XDG_CACHE_HOME = workflow_xdg_cache_dir)

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
devtools::session_info()

suppressPackageStartupMessages({
  library(tximeta)
})

# Explicitly set the tximeta cache using its function.
# The bfcloc.json will point to tximeta_data_cache_path.
tximeta::setTximetaBFC(dir = tximeta_data_cache_path, quiet = TRUE)

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
