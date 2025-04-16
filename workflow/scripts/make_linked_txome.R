log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

suppressPackageStartupMessages({
  library(tximeta)
  library(devtools)
})

organism_split <- strsplit(snakemake@params[["organism"]], "_")[[1]]
organism_reformat <- paste(paste(organism_split[1], organism_split[2]))

# Check the value of snakemake@params[["source"]]
if (snakemake@params[["source"]] %in% c("Ensembl", "GENCODE")) {
  source <- paste0("Local", snakemake@params[["source"]])
} else {
  source <- snakemake@params[["source"]]
}

makeLinkedTxome(
  indexDir = dirname(snakemake@input[["indexDir"]]),
  source = source,
  organism = organism_reformat,
  release = snakemake@params[["release"]],
  genome = snakemake@params[["genome"]],
  fasta = snakemake@input[["fasta"]],
  gtf = snakemake@input[["gtf"]],
  write = TRUE,
  jsonFile = snakemake@output[[1]]
)

date()
session_info()
