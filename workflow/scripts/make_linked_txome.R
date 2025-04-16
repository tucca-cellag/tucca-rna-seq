log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

suppressPackageStartupMessages({
  library(tximeta)
  library(devtools)
})

organism_split <- strsplit(snakemake@params[["organism"]], "_")[[1]]
organism_no_underscore <- paste(paste(organism_split[1], organism_split[2]))

makeLinkedTxome(
  indexDir = dirname(snakemake@input[["indexDir"]]),
  source = snakemake@params[["source"]],
  organism = organism_no_underscore,
  release = snakemake@params[["release"]],
  genome = snakemake@params[["genome"]],
  fasta = snakemake@input[["fasta"]],
  gtf = snakemake@input[["gtf"]],
  write = TRUE,
  jsonFile = snakemake@output
)

date()
session_info()
