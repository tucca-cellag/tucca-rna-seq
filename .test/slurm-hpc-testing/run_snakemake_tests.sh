#!/bin/bash
set -e
# .test/slurm-hpc-testing/run_snakemake_tests.sh
#
# This script mimics the jobs in .github/workflows/test.yml by invoking
# Snakemake with the targets and configuration the workflow uses.
#
# Usage:
#   sh ./run_snakemake_tests.sh {lint|target-rule|local-reads-refseq|local-reads-ensembl|sra-reads}
#
# IMPORTANT: Replace "YOUR_NCBI_API_KEY" with your actual key.
#
# This script assumes that you are in the repository root where the Snakefile
# (workflow/Snakefile) and the profiles directory exist, and that both the
# singularity and snakemake modules are loaded.
#
# Example:
#   sh ./run_snakemake_tests.sh lint

# Hard-coded API key—replace this placeholder with your actual NCBI API key.
API_KEY="YOUR_NCBI_API_KEY"

# Check that the first positional parameter is set.
if [ -z "${1}" ]; then
  echo "Usage: $0 {lint|target-rule|local-reads-refseq|local-reads-ensembl|sra-reads}"
  exit 1
fi

echo "Purging all loaded modules for a clean environment ..."
module purge

# Check that singularity is in the PATH.
if ! command -v singularity &>/dev/null; then
  echo -e "Singularity's module is not loaded.\nNow loading singularity/3.8.4 ..."
  module load singularity/3.8.4
fi

# Check that snakemake is in the PATH.
if ! command -v snakemake &>/dev/null; then
  echo -e "Singularity's module is not loaded.\nNow loading snakemake/8.27.1 ..."
  module load snakemake/8.27.1
fi

# Check that snakemake is in the PATH.
if ! command -v conda &>/dev/null; then
  echo -e "No Conda module is loaded.\nNow loading miniforge/24.11.2-py312 ..."
  module load miniforge/24.11.2-py312
fi

# The first positional parameter is the task.
TASK=$1
TARGET_RULE=$2
SOURCE=$3
EXTRA_ARG=$4

# Global settings for the Snakemake call.
PROFILE="profiles/slurm-dev"

case $TASK in
lint)
  echo "Running Snakemake lint..."
  # The --lint flag will check for problems in your Snakefile.
  snakemake --lint --profile ${PROFILE}
  ;;
target-rule)
  echo "Dry-run on local reads using a ${SOURCE} assembly..."
  snakemake ${TARGET_RULE} -np --workflow-profile ${PROFILE} \
    --configfile .test/singularity/local_reads_${SOURCE}/config/config.yaml \
    --config api_keys="{\"ncbi\": \"${API_KEY}\"}" ${EXTRA_ARG}
  echo "The dry-run was successful!!"
  echo "Running Snakemake workflow on local reads using a ${SOURCE} assembly..."
  snakemake ${TARGET_RULE} --workflow-profile ${PROFILE} \
    --configfile .test/singularity/local_reads_${SOURCE}/config/config.yaml \
    --config api_keys="{\"ncbi\": \"${API_KEY}\"}" ${EXTRA_ARG}
  ;;
local-reads-refseq)
  echo "Dry-run on local reads using a RefSeq assembly..."
  snakemake all -np --workflow-profile ${PROFILE} \
    --configfile .test/singularity/local_reads_refseq/config/config.yaml \
    --config api_keys="{\"ncbi\": \"${API_KEY}\"}"
  echo "The dry-run was successful!!"
  echo "Running Snakemake workflow on local reads using a RefSeq assembly..."
  snakemake all --workflow-profile ${PROFILE} \
    --configfile .test/singularity/local_reads_refseq/config/config.yaml \
    --config api_keys="{\"ncbi\": \"${API_KEY}\"}"
  ;;
local-reads-ensembl)
  echo "Dry-run on local reads using an Ensembl assembly ..."
  snakemake all -np --workflow-profile ${PROFILE} \
    --configfile .test/singularity/local_reads_ensembl/config/config.yaml \
    --config api_keys="{\"ncbi\": \"${API_KEY}\"}"
  echo "The dry-run was successful!!"
  echo "Running Snakemake workflow on local reads using an Ensembl assembly ..."
  snakemake all --workflow-profile ${PROFILE} \
    --configfile .test/singularity/local_reads_ensembl/config/config.yaml \
    --config api_keys="{\"ncbi\": \"${API_KEY}\"}"
  ;;
sra-reads)
  echo "Dry-run on SRA reads ..."
  snakemake all -np --workflow-profile ${PROFILE} \
    --configfile .test/singularity/sra_reads/config/config.yaml \
    --config api_keys="{\"ncbi\": \"${API_KEY}\"}"
  echo "The dry-run was successful!!"
  echo "Running Snakemake workflow on SRA reads ..."
  snakemake all --workflow-profile ${PROFILE} \
    --configfile .test/singularity/sra_reads/config/config.yaml \
    --config api_keys="{\"ncbi\": \"${API_KEY}\"}"
  ;;
*)
  echo "Invalid task provided: $TASK"
  echo "Usage: $0 {lint|target-rule|local-reads-refseq|local-reads-ensembl|sra-reads}"
  exit 1
  ;;
esac
