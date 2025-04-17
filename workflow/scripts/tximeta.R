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

# Load linkedTxome
tximeta::loadLinkedTxome(snakemake@input[["linkedTxome"]])

# Create coldata
files <- snakemake@input[["files"]]
files
names <- basename(sub(basename(files), "", files))
names
coldata <- data.frame(files, names)
coldata

# Create summarized experiment using tximeta
se <- tximeta::tximeta(coldata)

## Summarize to gene level
gse <- summarizeToGene(se)

saveRDS(se, file = snakemake@output[["se"]])
saveRDS(gse, file = snakemake@output[["gse"]])
