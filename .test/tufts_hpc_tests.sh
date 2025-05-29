#!/bin/bash
set -euo pipefail
# .test/tufts_hpc_tests.sh
#
# This script mimics the jobs in .github/workflows/test.yml by invoking
# Snakemake with the targets and configuration the workflow uses.
#
# IMPORTANT: Replace "YOUR_NCBI_API_KEY" with your actual key.
#
# This script assumes that you are in the repository root where the Snakefile
# (workflow/Snakefile) and the profiles directory exist, and that both the
# singularity and snakemake modules are loaded.

# Hard-coded API key—replace this placeholder with your actual NCBI API key.
API_KEY="YOUR_NCBI_API_KEY"
DEFAULT_CONDA_PREFIX="/cluster/tufts/kaplanlab/bbromb01/tucca-rna-seq-dev"

# Function to print usage information
print_usage() {
  echo "Usage: $0 <task> [options]"
  echo ""
  echo "This script runs various Snakemake workflow tests."
  echo ""
  echo "Tasks:"
  echo "  lint                             Lint the Snakemake workflow."
  echo "  conda-create-envs-only           Create Conda environments only (no other options)."
  echo "  test                             Run a specific target rule with local reads."
  echo "    Options for test:"
  echo "      -t <target_name>             Optional Snakemake target (defaults to 'all')."
  echo "      -r                           Use 'refseq' as source (defaults to 'ensembl')."
  echo "      -c                           Use complex configuration (defaults to basic)."
  echo "  sra-reads                        Run workflow with SRA reads."
  echo "    Options for sra-reads:"
  echo "      -t <target_name>             Optional Snakemake target (defaults to 'all')."
  echo "                                   Config is from: .test/data/sra_reads/config/config.yaml"
  echo ""
  echo "Global Options (can be used with any task, or alone for help):"
  echo "  -h                               Show this help message and exit."
  echo "  -p [<path>]                      Set the --conda-prefix for Snakemake (for tasks: conda-create-envs-only, test, sra-reads)."
  echo "                                     - If -p is omitted, no --conda-prefix is used."
  echo "                                     - If -p is used without <path>, defaults to: ${DEFAULT_CONDA_PREFIX}"
  echo "                                     - If -p <path> is used, the specified <path> is used."
  echo ""
  echo "Examples:"
  echo "  sh $0 test -t my_custom_rule -r -x \"--cores 1\"  # Uses refseq source"
  echo "  sh $0 test -t my_custom_rule      # Uses ensembl source (default)"
  echo "  sh $0 test -p                   # Uses --conda-prefix ${DEFAULT_CONDA_PREFIX}"
  echo "  sh $0 test -p /my/custom/path   # Uses --conda-prefix /my/custom/path"
  echo "  sh $0 test -t my_custom_rule -r -x \"--cores 1\" -p # Uses refseq, default prefix"
  echo "  sh $0 sra-reads -t specific_sra_target -x \"--debug\""
  echo "  sh $0 -h"
}

# --- Initial Checks and Task Definition ---
if [ -z "${1}" ] || [ "${1}" == "-h" ]; then
  print_usage
  exit 0
fi

TASK=$1
shift # Remove task from arguments, so getopts processes the rest for options

# --- Option Parsing ---
# Initialize variables for options
TARGET_OPT="all"
EXTRA_SNAKEMAKE_ARGS_OPT=""
CONFIG_SUBDIR="config_basic"
USE_REFSEQ_SOURCE=false
P_OPTION_USED=false
CONDA_PREFIX_ARG="" # Will store the actual prefix path if -p is used

# The leading colon in getopts string enables silent error handling.
# p:: means -p can take an optional argument.
while getopts ":t:x:crhp::" opt; do
  case ${opt} in
  t) TARGET_OPT="${OPTARG}" ;;
  x) EXTRA_SNAKEMAKE_ARGS_OPT="${OPTARG}" ;;
  c) CONFIG_SUBDIR="config_complex" ;;
  r) USE_REFSEQ_SOURCE=true ;;
  p)
    P_OPTION_USED=true
    if [ -z "${OPTARG}" ]; then # -p was used without a value
      CONDA_PREFIX_ARG="${DEFAULT_CONDA_PREFIX}"
    else # -p was used with a value
      CONDA_PREFIX_ARG="${OPTARG}"
    fi
    ;;
  h)
    print_usage
    exit 0
    ;;
  \?)
    echo "Invalid option: -${OPTARG}" >&2
    print_usage
    exit 1
    ;;
  :) # Handles missing arguments for options that require them (e.g., t, x)
    echo "Option -${OPTARG} requires an argument." >&2
    print_usage
    exit 1
    ;;
  esac
done
shift $((OPTIND - 1)) # Remove parsed options and their arguments from $@

# Check for any remaining non-option arguments (unexpected)
if [ "$#" -gt 0 ]; then
  echo "Error: Unexpected arguments found after options: $*" >&2
  print_usage
  exit 1
fi

# --- Module Loading ---
echo "Purging all loaded modules for a clean environment ..."
module purge

# Check that singularity is in the PATH.
if ! command -v singularity &>/dev/null; then
  echo -e "Singularity's module is not loaded.\nNow loading singularity/3.8.4 ..."
  module load singularity/3.8.4
fi

# Check that snakemake is in the PATH.
if ! command -v snakemake &>/dev/null; then
  echo -e "Snakemake's module is not loaded.\nNow loading snakemake/8.27.1 ..."
  module load snakemake/8.27.1
fi

# Check that conda is in the PATH.
if ! command -v conda &>/dev/null; then
  echo -e "No Conda module is loaded.\nNow loading miniforge/24.11.2-py312 ..."
  module load miniforge/24.11.2-py312
fi

# Global settings for the Snakemake call.
PROFILE="profiles/slurm-dev"

# --- Task Execution ---
# Prepare dynamic arguments for snakemake calls that use conda prefix
SNAKEMAKE_CONDA_ARGS=()
CONDA_PREFIX_MESSAGE="" # For logging
if [ "${P_OPTION_USED}" = true ]; then
  SNAKEMAKE_CONDA_ARGS+=("--conda-prefix" "${CONDA_PREFIX_ARG}")
  CONDA_PREFIX_MESSAGE=" with conda prefix: ${CONDA_PREFIX_ARG}"
else
  CONDA_PREFIX_MESSAGE=" (no --conda-prefix specified)"
fi

case $TASK in
lint)
  echo "Linting the Snakemake workflow..."
  # Lint task does not use --conda-prefix
  snakemake --lint --verbose --workflow-profile ${PROFILE} ${EXTRA_SNAKEMAKE_ARGS_OPT}
  ;;
conda-create-envs-only)
  echo "Running snakemake --conda-create-envs-only${CONDA_PREFIX_MESSAGE}"
  snakemake --conda-create-envs-only --workflow-profile ${PROFILE} "${SNAKEMAKE_CONDA_ARGS[@]}" ${EXTRA_SNAKEMAKE_ARGS_OPT}
  ;;
test)
  SMK_TARGET=${TARGET_OPT:-all}
  SOURCE=$([ "${USE_REFSEQ_SOURCE}" = true ] && echo "refseq" || echo "ensembl")

  echo "Dry-run for target '${SMK_TARGET}' using a '${SOURCE}' assembly with '${CONFIG_SUBDIR}' config${CONDA_PREFIX_MESSAGE}..."
  if snakemake ${SMK_TARGET} -np --workflow-profile ${PROFILE} \
    --configfile ".test/local_reads/${SOURCE}/${CONFIG_SUBDIR}/config.yaml" \
    --config api_keys="{\"ncbi\": \"${API_KEY}\"}" "${SNAKEMAKE_CONDA_ARGS[@]}" ${EXTRA_SNAKEMAKE_ARGS_OPT}; then
    echo "The dry-run for target '${SMK_TARGET}' on '${SOURCE}' assembly with '${CONFIG_SUBDIR}' was successful."
    echo "Running Snakemake workflow for target '${SMK_TARGET}' using a '${SOURCE}' assembly with '${CONFIG_SUBDIR}' config${CONDA_PREFIX_MESSAGE}..."
    snakemake ${SMK_TARGET} --workflow-profile ${PROFILE} \
      --configfile ".test/local_reads/${SOURCE}/${CONFIG_SUBDIR}/config.yaml" \
      --config api_keys="{\"ncbi\": \"${API_KEY}\"}" "${SNAKEMAKE_CONDA_ARGS[@]}" ${EXTRA_SNAKEMAKE_ARGS_OPT}
  else
    DRY_RUN_EXIT_CODE=$?
    echo "The dry-run for target '${SMK_TARGET}' on '${SOURCE}' assembly with '${CONFIG_SUBDIR}' failed with exit code ${DRY_RUN_EXIT_CODE}. Skipping actual run." >&2
    exit ${DRY_RUN_EXIT_CODE}
  fi
  ;;
sra-reads)
  SMK_TARGET=${TARGET_OPT:-all}

  echo "Dry-run on SRA reads (target: '${SMK_TARGET}')${CONDA_PREFIX_MESSAGE}..."
  if snakemake ${SMK_TARGET} -np --workflow-profile ${PROFILE} \
    --configfile ".test/sra_reads/config/config.yaml" \
    --config api_keys="{\"ncbi\": \"${API_KEY}\"}" "${SNAKEMAKE_CONDA_ARGS[@]}" ${EXTRA_SNAKEMAKE_ARGS_OPT}; then
    echo "The dry-run on SRA reads (target: '${SMK_TARGET}') was successful!!"
    echo "Running Snakemake workflow on SRA reads (target: '${SMK_TARGET}')${CONDA_PREFIX_MESSAGE}..."
    snakemake ${SMK_TARGET} --workflow-profile ${PROFILE} \
      --configfile ".test/sra_reads/config/config.yaml" \
      --config api_keys="{\"ncbi\": \"${API_KEY}\"}" "${SNAKEMAKE_CONDA_ARGS[@]}" ${EXTRA_SNAKEMAKE_ARGS_OPT}
  else
    DRY_RUN_EXIT_CODE=$?
    echo "The dry-run on SRA reads (target: '${SMK_TARGET}') failed with exit code ${DRY_RUN_EXIT_CODE}. Skipping actual run." >&2
    exit ${DRY_RUN_EXIT_CODE}
  fi
  ;;
*)
  echo "Invalid task provided: $TASK" >&2
  print_usage
  exit 1
  ;;
esac

echo "Script finished for task: $TASK"
exit 0
