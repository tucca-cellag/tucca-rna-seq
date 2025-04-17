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

files <- snakemake@input[["files"]]
files
names <- basename(sub(basename(files), "", files))
names
coldata <- data.frame(files, names)
coldata

saveRDS(files, file = snakemake@output[[1]])
