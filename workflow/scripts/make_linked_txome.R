log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")
date()
suppressPackageStartupMessages({
  library(devtools)
})
devtools::session_info()

suppressPackageStartupMessages({
  library(tximeta)
})

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

setTximetaBFC(snakemake@params[["tximeta_cache"]])

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
