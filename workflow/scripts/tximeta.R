log <- base::file(snakemake@log[[1]], open = "wt")
base::sink(log)
base::sink(log, type = "message")
base::date()
base::library(devtools)
devtools::session_info()

base::library(tximeta)
base::library(tidyverse)

# Load linkedTxome
# TODO: If support for makeLinkedTxome(source = c("Ensembl", "GENCODE")) is
# added instead of forcing "LocalEnsembl" and "LocalGENCODE" the logic in
# this script will need to be refactored to work with EnsDb objects
# See: https://github.com/thelovelab/tximeta/blob/devel/R/tximeta.R
# specifically getTxDb() call and definition
tximeta::loadLinkedTxome(snakemake@input[["linkedTxome"]])

# load required metadata dependencies
sample_names <- snakemake@params[["sample_names"]]
samples <- utils::read.delim(snakemake@config$samples)
units <- utils::read.delim(snakemake@config$units)

# create coldata
coldata <- dplyr::left_join(samples, units, by = "sample_name") %>%
  dplyr::mutate(names = base::paste(sample_name, unit_name, sep = "_")) %>%
  dplyr::mutate(files = base::paste0("results/salmon/", names, "/quant.sf")) %>%
  dplyr::select(
    names, sample_name, unit_name, tidyselect::everything(), -sra, -fq1, -fq2
  )
base::rownames(coldata) <- coldata$names

factor_specs_list <- snakemake@params[["factors"]]

# Factorize columns and set reference levels
if (!base::is.null(factor_specs_list) && base::length(factor_specs_list) > 0) {
  for (i in base::seq_along(factor_specs_list)) {
    # Get the factor and ref level
    factor_spec <- factor_specs_list[[i]]
    factor_name <- factor_spec$name
    ref_level <- factor_spec$reference_level

    # TODO: These if-else blocks can probably be dealt w/ in the config schema
    if (factor_name %in% base::colnames(coldata)) {
      coldata[[factor_name]] <- base::as.factor(coldata[[factor_name]])

      if (ref_level %in% base::levels(coldata[[factor_name]])) {
        coldata[[factor_name]] <- stats::relevel(
          coldata[[factor_name]],
          ref = ref_level
        )
        base::message(base::paste0(
          "Column '", factor_name, "' converted to factor with ",
          "reference level '", ref_level, "'."
        ))
      } else {
        base::warning(base::paste0(
          "Reference level '", ref_level, "' not found in column '",
          factor_name, "'. Using default reference level for this factor."
        ))
      }
    } else {
      base::warning(base::paste0(
        "Column '", factor_name, "' specified for factorization (in",
        "specification #", i, ") not found in coldata."
      ))
    }
  }
}

base::print("Final coldata structure after factorization:")
utils::str(coldata)

if (snakemake@config[["ref_assembly"]][["source"]] == "RefSeq") {
  skip_seq_info <- TRUE
} else {
  skip_seq_info <- FALSE
}

extra_params_str <- base::as.character(snakemake@params[["extra"]])

# TODO: skipSeqinfo = FALSE triggers gtf2RefSeq() which errors out due to
# an incorrectly formatted assembly_report.txt. Might be impossible to fix
# without using a different download service for the RefSeq data or tximeta
# is outdated in how it deals with RefSeq inputs wrt skipSeqinfo
cmd_str <- "coldata = coldata, skipSeqinfo = skip_seq_info"
if (extra_params_str != "") {
  cmd_str <- base::paste(cmd_str, extra_params_str, sep = ",")
}

tximeta_cmd <- base::paste0(
  "tximeta::tximeta(", cmd_str, ")"
)
base::message("tximeta command:")
base::message(tximeta_cmd)

# Create se object via tximeta
base::tryCatch(
  {
    se <- base::eval(base::parse(text = tximeta_cmd))
  },
  error = function(e) {
    detailed_error_message <- base::paste0(
      "Error executing tximeta::tximeta() \n",
      "Command attempted: '", tximeta_cmd, "'\n",
      "Original error message: ", e$message
    )
    base::stop(detailed_error_message, call. = FALSE)
  }
)

## Summarize to gene level
gse <- tximeta::summarizeToGene(
  se,
  assignRanges = "abundant",
  countsFromAbundance = "lengthScaledTPM"
)

# Retrieve TxDb or EnsDb database
db <- tximeta::retrieveDb(se)

base::class(db)

base::saveRDS(se, file = snakemake@output[["se"]])
base::saveRDS(gse, file = snakemake@output[["gse"]])

base::save.image(file = snakemake@output[["image"]])
