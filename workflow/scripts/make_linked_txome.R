log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

suppressPackageStartupMessages({
  library(tximeta)
  library(devtools)
})

org_split <- strsplit(organism, "_")[[1]]
organism <- paste(paste(org_split[1], org_split[2]))

makeLinkedTxome(
  indexDir = dirname(snakemake@input[["indexDir"]]),
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
