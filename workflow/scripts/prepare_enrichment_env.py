#!/usr/bin/env python3
"""
Generates a post-deployment shell script for the enrichment Conda environment.

This script creates a shell script that, when run, will install the correct
organism-specific Bioconductor database package into the already-created
base environment. All print statements are redirected to the log file
specified by the Snakemake rule.
"""
import sys
from snakemake.script import snakemake

# --- Log setup ---
# This is the correct way to capture all script output to the log file
# specified in the Snakefile's `log` directive.
try:
    log_file = open(snakemake.log[0], "w", encoding="utf-8")
    sys.stdout = sys.stderr = log_file
except IOError:
    print(f"Warning: Could not open log file {snakemake.log[0]} for writing.")
    # The script will continue and print to the console if the log can't be
    # opened.

# --- Get parameters from Snakemake ---
# The desired package name (e.g., "org.Hs.eg.db") from the rule's params
org_db_pkg = snakemake.params.org_db_pkg
output_script_path = snakemake.output.script

# The full conda package name, e.g., "bioconductor-org.hs.eg.db"
conda_pkg_name = f"bioconductor-{org_db_pkg.lower()}"

print("--- Preparing post-deployment script ---")
print(f"Log file: {snakemake.log[0]}")
print(f"Output script path: {output_script_path}")
print(f"Dynamic package to install: {conda_pkg_name}")

# --- Main logic ---
# Write the shell script that will perform the installation
with open(output_script_path, "w", encoding="utf-8") as f:
    f.write("#!/usr/bin/env bash\n")
    f.write("set -e\n\n")
    f.write(f"# Dynamically installing: {conda_pkg_name}\n")
    f.write("mamba install -y --quiet -c bioconda -c " +
            f"conda-forge {conda_pkg_name}\n")

print("\nSuccessfully generated post-deployment script.")
print("---")

# Close the log file if it was opened
if 'log_file' in locals() and not log_file.closed:
    log_file.close()
