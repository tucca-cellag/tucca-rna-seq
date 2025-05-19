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

factor_specs_list <- snakemake@params[["factors"]]

# Factorize columns and set reference levels
if (!is.null(factor_specs_list) && length(factor_specs_list) > 0) {
  for (i in seq_along(factor_specs_list)) {
    # Get the factor and ref level
    factor <- factor_specs_list[[i]]
    factor_name <- factor$name
    ref_level <- factor$reference_level

    # TODO: These if-else blocks can probably be dealt w/ in the config schema
    if (factor_name %in% colnames(coldata)) {
      coldata[[factor_name]] <- as.factor(coldata[[factor_name]])

      if (ref_level %in% levels(coldata[[factor_name]])) {
        coldata[[factor_name]] <- relevel(
          coldata[[factor_name]],
          ref = ref_level
        )
        message(paste0(
          "Column '", factor_name, "' converted to factor with",
          "reference level '", ref_level, "'."
        ))
      } else {
        warning(paste0(
          "Reference level '", ref_level, "' not found in column '",
          factor_name, "'. Using default reference level for this factor."
        ))
      }
    } else {
      warning(paste0(
        "Column '", factor_name, "' specified for factorization (in",
        "specification #", i, ") not found in coldata."
      ))
    }
  }
}

print("Final coldata structure after factorization:")
str(coldata)

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
