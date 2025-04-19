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
  library(tidyverse)
})

# Load linkedTxome
# TODO: If support for makeLinkedTxome(source = c("Ensembl", "GENCODE")) is
# added instead of forcing "LocalEnsembl" and "LocalGENCODE" the logic in
# workflow/scripts/tximeta.R will need to be refactored to work with EnsDb
# objects
# See: https://github.com/thelovelab/tximeta/blob/devel/R/tximeta.R
# Specifically getTxDb() definition where either EnsDb() or loadDb() is called
tximeta::loadLinkedTxome(snakemake@input[["linkedTxome"]])

# Create coldata
files <- snakemake@input[["quant"]]
files
names <- basename(dirname(files))
names
coldata <- data.frame(files, names)
coldata

extra <- snakemake@params[["extra"]]

# Create summarized experiment using tximeta
se <- tximeta(coldata)

## Summarize to gene level
gse <- summarizeToGene(
  se,
  assignRanges = "abundant",
  countsFromAbundance = "lengthScaledTPM"
)

# Retrieve TxDb or EnsDb database
db <- retrieveDb(se)

class(db)

saveRDS(se, file = snakemake@output[["se"]])
saveRDS(gse, file = snakemake@output[["gse"]])
