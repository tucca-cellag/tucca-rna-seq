#!/usr/bin/env python3
"""
Generates a post-deployment shell script for the enrichment Conda environment.

This script creates a shell script that, when run, will install the correct
organism-specific Bioconductor database package into the already-created
base environment. All print statements are redirected to the log file
specified by the Snakemake rule.
"""
import contextlib
from snakemake.script import snakemake

# The 'with' statement now manages the log file handle, satisfying the linter.
try:
    with open(snakemake.log[0], "w", encoding="utf-8") as log_file:
        # The redirect_stdout/stderr context managers handle the redirection.
        with contextlib.redirect_stdout(log_file), contextlib.redirect_stderr(log_file):
            # --- All code within this block will be logged ---

            # --- Get parameters from Snakemake ---
            org_db_pkg = snakemake.params.org_db_pkg
            output_script_path = snakemake.output.script

            # --- Main logic ---
            # Ensure the package name is lowercase for consistency.
            conda_pkg_name = f"bioconductor-{org_db_pkg.lower()}"

            print("--- Preparing post-deployment script ---")
            print(f"Log file: {snakemake.log[0]}")
            print(f"Output script path: {output_script_path}")
            print(f"Dynamic package to install: {conda_pkg_name}")

            # Write the shell script that will perform the installation.
            with open(output_script_path, "w", encoding="utf-8") as f:
                f.write("#!/usr/bin/env bash\n")
                f.write("set -e\n\n")
                f.write(f"# Dynamically installing: {conda_pkg_name}\n")
                f.write("mamba install -y --quiet -c bioconda -c "
                        f"conda-forge {conda_pkg_name}\n")

            print("\nSuccessfully generated post-deployment script.")
            print("---")

except IOError as e:
    print(f"Warning: Could not open log file {snakemake.log[0]} " +
          f"for writing. Error: {e}")
