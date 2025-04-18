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
files <- snakemake@input[["quant"]]
files
names <- basename(dirname(files))
names
coldata <- data.frame(files, names)
coldata

extra <- snakemake@params[["extra"]]

# Create summarized experiment using tximeta
# Evaluate the following
se <- base::eval(
  # ... parsed expression
  base::parse(
    # ... of tximeta and its arguments
    text = base::paste0(
      "tximeta::tximeta(", coldata, extra, ");"
    )
  )
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
