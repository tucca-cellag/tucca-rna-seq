#!/bin/bash
# .test/slurm-hpc-testing/run_snakemake_tests.sh
#
# This script mimics the jobs in .github/workflows/test.yml by invoking
# Snakemake with the targets and configuration the workflow uses.
#
# Usage:
#   ./run_snakemake_tests.sh {lint|env-init|local-reads|sra-reads}
#
# IMPORTANT: Replace "YOUR_NCBI_API_KEY" with your actual key.
#
# This script assumes that you are in the repository root where the Snakefile
# (workflow/Snakefile) and the profiles directory exist, and that both the
# singularity and snakemake modules are loaded.
#
# Example:
#   ./run_snakemake_tests.sh sc-genome

# Hard-coded API key—replace this placeholder with your actual NCBI API key.
API_KEY="YOUR_NCBI_API_KEY"

# Check that singularity is in the PATH.
if ! command -v singularity &>/dev/null; then
  echo "Error: Singularity is not available. Please load its module."
  exit 1
fi

# Check that snakemake is in the PATH.
if ! command -v snakemake &>/dev/null; then
  echo "Error: Snakemake is not available. Please load its module."
  exit 1
fi

# Check that snakemake is in the PATH.
if ! command -v conda &>/dev/null; then
  echo "Error: Conda is not available. Please load its module (miniforge)."
  exit 1
fi

if [ -z "${1}" ]; then
  echo "Usage: $0 {lint|env-init|local-reads|sra-reads}"
  exit 1
fi

# The first positional parameter is the task.
TASK=$1

# Global settings for the Snakemake call.
PROFILE="profiles/slurm-dev"

case $TASK in
lint)
  echo "Running Snakemake lint..."
  # The --lint flag will check for problems in your Snakefile.
  snakemake --lint --profile ${PROFILE}
  ;;
local-reads)
  echo "Running Snakemake workflow on local reads ..."
  snakemake all --profile ${PROFILE} \
    --configfile .test/singularity/local_reads/config/config.yaml \
    --config api_keys="{\"ncbi\": \"${API_KEY}\"}"
  ;;
sra-reads)
  echo "Running Snakemake workflow on SRA reads ..."
  snakemake all --profile ${PROFILE} \
    --configfile .test/singularity/sra_reads/config/config.yaml \
    --config api_keys="{\"ncbi\": \"${API_KEY}\"}"
  ;;
*)
  echo "Invalid task provided: $TASK"
  echo "Usage: $0 {lint|env-init|local-reads|sra-reads}"
  exit 1
  ;;
esac
