#!/usr/bin/env python3
"""
Prepares a Conda environment file for enrichment analysis.

This script takes a template Conda environment YAML, adds a specified
organism-specific Bioconductor database package to its dependencies,
and writes the result to a new YAML file.

File: workflow/scripts/prepare_enrichment_env.py
"""

import datetime
import sys

import yaml
from snakemake.script import snakemake

# --- Log setup ---
# The log file path is provided by Snakemake
log_file_path = snakemake.log[0]
with open(log_file_path, "w", encoding="utf-8") as log_f:
    # Redirect stdout and stderr to the log file
    sys.stdout = sys.stderr = log_f

    print(f"--- Log generated on {datetime.date.today()} ---")

    # --- Get parameters from Snakemake ---
    template_path = snakemake.input.template
    org_db_pkg = snakemake.params.org_db_pkg

    print(f"Input template: {template_path}")
    print(f"Organism DB package to add: {org_db_pkg}")

    # --- Main logic ---
    conda_pkg_name = f"bioconductor-{org_db_pkg.lower()}"

    with open(template_path, "r", encoding="utf-8") as f:
        env_config = yaml.safe_load(f)

    if conda_pkg_name not in env_config["dependencies"]:
        print(f"Adding '{conda_pkg_name}' to dependencies.")
        env_config["dependencies"].append(conda_pkg_name)
    else:
        print(f"'{conda_pkg_name}' is already in the dependency list.")

    with open(template_path, "w", encoding="utf-8") as f:
        yaml.dump(env_config, f, sort_keys=False)

    print("Successfully wrote updated environment file.")
