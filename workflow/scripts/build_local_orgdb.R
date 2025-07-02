# workflow/scripts/build_local_orgdb.R
#
# Purpose:
# This script uses AnnotationForge to create a local OrgDb package from a GTF
# file when the pre-built package is not available on Bioconductor.
#
# Inputs:
#   - gtf: Path to the genome annotation GTF file.
#
# Outputs:
#   - out_dir: A directory where the source package will be built.
#
# Parameters:
#   - From Snakemake params: species, genus, version, author, extra

log <- base::file(snakemake@log[[1]], open = "wt")
base::sink(log)
base::sink(log, type = "message")
base::date()

base::library(package = "devtools", character.only = TRUE)
devtools::session_info()

base::library(package = "AnnotationForge", character.only = TRUE)


# --- Parameters ---
gtf_file <- snakemake@input[["gtf"]]
tax_id_file <- snakemake@input[["tax_id_file"]]
output_dir <- snakemake@output[["out_dir"]]
species <- snakemake@params[["species"]]
genus <- snakemake@params[["genus"]]
version <- snakemake@params[["version"]]
author <- snakemake@params[["author"]]
extra_params_str <- snakemake@params[["extra"]]

# Read the taxonomy ID from the input file
tax_id <- base::readLines(tax_id_file, warn = FALSE)[1]
base::message("Using Taxonomy ID: ", tax_id)


base::message("Starting AnnotationForge process...")
# ... (rest of the messages are the same) ...

# Build command string dynamically
cmd_list <- base::list(
  version = base::sprintf("'%s'", version),
  author = base::sprintf("'%s'", author),
  maintainer = base::sprintf("'%s'", author),
  outputDir = base::sprintf("'%s'", output_dir),
  tax_id = base::sprintf("'%s'", tax_id), # Use the fetched tax_id
  genus = base::sprintf("'%s'", genus),
  species = base::sprintf("'%s'", species),
  NCBIFiles = base::sprintf("'%s'", gtf_file)
)

# ... (The rest of the script for building and evaluating the command is unchanged) ...

# Convert list to a named string of arguments
default_args <- base::paste(
  base::names(cmd_list), cmd_list,
  sep = " = ", collapse = ", "
)

final_args <- default_args
if (extra_params_str != "") {
  final_args <- base::paste(final_args, extra_params_str, sep = ", ")
}

cmd <- base::paste0("AnnotationForge::makeOrgPackageFromNCBI(", final_args, ")")
base::message("Command to be executed:")
base::message(cmd)

# Execute the command
base::eval(base::parse(text = cmd))

base::message("Local OrgDb package built successfully in: ", output_dir)
base::date()
