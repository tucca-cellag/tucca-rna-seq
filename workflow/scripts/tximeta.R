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

snakemake@input[["files"]]
files <- list.files(
  "results/salmon",
  full = TRUE,
  recursive = TRUE,
  pattern = "^+quant.sf$"
)
print("")
files

saveRDS(files, file = snakemake@output[[1]])
