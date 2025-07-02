# -----------------------------------------------------------------------------
# tucca-rna-seq/workflow/scripts/deseqtransform.R
#
# Purpose:
# This script applies a variance-stabilizing transformation to a DESeqDataSet
# object. The transformation method (e.g., rlog, vst) is specified via
# Snakemake parameters. This is a common step before downstream visualization
# analyses like PCA or sample clustering, which assume homoscedastic data
# (variance is constant across the range of mean values).
#
# Inputs:
#   - dds: An RDS file containing the raw DESeqDataSet object.
#
# Outputs:
#   - dst: An RDS file containing the transformed DESeqTransform object.
#   - image: An .RData file saving the entire R session image for debugging.
#
# Parameters:
#   - method: The transformation method to use (e.g., "rlog", "vst").
#   - extra: A string of extra arguments to be passed to the transformation
#            function.
#
# -----------------------------------------------------------------------------
log <- base::file(snakemake@log[[1]], open = "wt")
base::sink(log)
base::sink(log, type = "message")
base::date()

base::library(package = "devtools", character.only = TRUE)
devtools::session_info()

base::library(package = "DESeq2", character.only = TRUE)

# Get Snakemake inputs and params
dds <- base::readRDS(snakemake@input[["dds"]])
transformation_method <- base::tolower(snakemake@params[["method"]])
extra_params_str <- base::as.character(snakemake@params[["extra"]])

cmd_str <- "dds"
if (extra_params_str != "") {
  cmd_str <- base::paste(cmd_str, extra_params_str, sep = ",")
}

dst_cmd <- base::paste0(
  "DESeq2::", transformation_method, "(", cmd_str, ")"
)
base::message("Transformation command:")
base::message(dst_cmd)

# Perform the transformation
base::tryCatch(
  {
    dst <- base::eval(base::parse(text = dst_cmd))
  },
  error = function(e) {
    detailed_error_message <- base::paste0(
      "Error executing DESeq2 transformation ('", transformation_method,
      "'). \n", "Command attempted: '", dst_cmd, "'\n",
      "Original error message: ", e$message
    )
    base::stop(detailed_error_message, call. = FALSE)
  }
)

# Save the DESeqTransform object
base::saveRDS(dst, file = snakemake@output[["dst"]])

base::save.image(file = snakemake@output[["image"]])
