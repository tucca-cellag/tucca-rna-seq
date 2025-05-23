log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")
date()
suppressPackageStartupMessages({
  library(devtools)
})
devtools::session_info()

# Configure tximeta to use cache and config directories within the project
# workspace. This ensures that tximeta's files (like bfcloc.json and the main
# data cache) are stored in predictable, project-local paths. This complements
# system-level configurations or container environment setups (e.g., GitHub
# Actions bind mounts that make system-level config paths like
# /home/runner/.config/R writable). See issue #11 for more context.

# 1. Define and set R_USER_CACHE_DIR/XDG_CACHE_HOME to a writable path in the
# workspace. This directs where packages like tximeta (via rappdirs) might
# place user-specific configuration files (e.g., bfcloc.json) if not managed
# by setTximetaBFC's primary mechanism.
r_user_cache_base_dir <- file.path(getwd(), ".r_user_cache_for_tximeta")
if (!dir.exists(r_user_cache_base_dir)) {
  dir.create(r_user_cache_base_dir, recursive = TRUE, showWarnings = FALSE)
}
Sys.setenv(R_USER_CACHE_DIR = r_user_cache_base_dir)
Sys.setenv(XDG_CACHE_HOME = r_user_cache_base_dir)
message(paste(
  "Set R_USER_CACHE_DIR and XDG_CACHE_HOME to:", r_user_cache_base_dir
))

# 2. Define the main tximeta data cache path (where large SQLite DBs, etc.,
# will go). This path is typically provided via Snakemake params
# (e.g., "resources/tximeta").
tximeta_data_cache_path <- snakemake@params[["tximeta_cache"]]
if (!dir.exists(tximeta_data_cache_path)) {
  dir.create(tximeta_data_cache_path, recursive = TRUE, showWarnings = FALSE)
}
message(paste(
  "Main tximeta data cache path defined as:", tximeta_data_cache_path
))

# 3. Set environment variables for tximeta and BiocFileCache to use the defined
# data cache path.
Sys.setenv(TXIMETA_HUB_CACHE = tximeta_data_cache_path)
Sys.setenv(BIOCFILECACHE_CACHE = tximeta_data_cache_path)
message(paste(
  "Set TXIMETA_HUB_CACHE and BIOCFILECACHE_CACHE to:", tximeta_data_cache_path
))

# Load tximeta
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

# Prepare parameters for makeLinkedTxome
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

message(paste(
  "Parameters for makeLinkedTxome:",
  "\n  Index Dir (first):", dirname(snakemake@input[["index_dir"]])[1],
  "\n  Source:", source,
  "\n  Organism:", organism_reformat,
  "\n  Release:", snakemake@params[["release"]],
  "\n  Genome:", snakemake@params[["genome"]],
  "\n  Fasta:", snakemake@input[["fasta"]],
  "\n  GTF:", snakemake@input[["gtf"]],
  "\n  JSON Output:", snakemake@output[["jsonFile"]]
))

# Call makeLinkedTxome
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
message(paste("Script completed at:", Sys.time()))
