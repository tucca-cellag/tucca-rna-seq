log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

suppressPackageStartupMessages({
  library(tximeta)
  library(devtools)
})

ss <- strsplit(organism, "_")[[1]]
organism <- paste(paste(ss[1], ss[2]))

makeLinkedTxome(
  indexDir = dirname(indexDir),
  source = snakemake@params[["source"]],
  organism = snakemake@params[["organism"]],
  release = snakemake@params[["release"]],
  genome = snakemake@params[["genome"]],
  fasta = snakemake@input[["fasta"]],
  gtf = snakemake@input[["gtf"]],
  write = TRUE,
  jsonFile = snakemake@output
)

date()
session_info()
