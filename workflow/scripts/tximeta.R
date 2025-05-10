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
# this script will need to be refactored to work with EnsDb objects
# See: https://github.com/thelovelab/tximeta/blob/devel/R/tximeta.R
# specifically getTxDb() call and definition
tximeta::loadLinkedTxome(snakemake@input[["linkedTxome"]])

# load required metadata dependencies
sample_names <- snakemake@params[["sample_names"]]
samples <- read.delim(snakemake@config$samples)
units <- read.delim(snakemake@config$units)

# create coldata
coldata <- left_join(samples, units, by = "sample_name") %>%
  mutate(names = paste(sample_name, unit_name, sep = "_")) %>%
  mutate(files = paste0("results/salmon/", names, "/quant.sf")) %>%
  select(names, sample_name, unit_name, everything(), -sra, -fq1, -fq2)
rownames(coldata) <- coldata$names
coldata


if (snakemake@config[["ref_assembly"]][["source"]] == "RefSeq") {
  skipSeqinfo <- TRUE
} else {
  skipSeqinfo <- FALSE
}

extra <- snakemake@params[["tximeta_extra"]]

# Create summarized experiment using tximeta
se <- tximeta(coldata,
  skipSeqinfo = skipSeqinfo,
  # TODO: skipSeqinfo = FALSE triggers gtf2RefSeq() which errors out due to
  # an incorrectly formatted assembly_report.txt. Might be impossible to fix
  # without using a different download service for the RefSeq data or tximeta
  # is outdated in how it deals with RefSeq inputs wrt skipSeqinfo
  extra
)

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

save.image(file = snakemake@output[["image"]])
