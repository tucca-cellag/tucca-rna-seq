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

print(snakemake@inputs[["quants_path"]])

files <- file.path(list.files("results/salmon", full = TRUE), "quant.sf")
all(file.exists(files))
list.files("results/salmon")
names <- list.files("results/salmon")
names <- sub(".quant", "", names)
names
coldata <- data.frame(files, names)
coldata

saveRDS(coldata, file = snakemake@output[[1]])
