log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = c("output", "message"))
date()
suppressPackageStartupMessages({
  library(devtools)
})
devtools::session_info()

suppressPackageStartupMessages({
  library(DESeq2)
})

# Get Snakemake inputs and params
dds <- readRDS(snakemake@input[["dds"]])
transformation_method <- tolower(snakemake@params[["method"]])
extra <- snakemake@params[["extra"]]

dst_cmd <- base::paste0(
  "DESeq2::", transformation_method, "(dds", extra, ")"
)
base::message("Transformation command:")
base::message(dst_cmd)

# Perform the transformation
tryCatch(
  {
    dst <- base::eval(base::parse(text = dst_cmd))
  },
  error = function(msg) {
    stop(
      "Unsupported transformation_method: '", transformation_method,
      "'. Choose 'vst' or 'rlog'."
    )
  }
)

# Save the DESeqTransform object
saveRDS(dst, file = snakemake@output[["dst"]])

save.image(file = snakemake@output[["image"]])
