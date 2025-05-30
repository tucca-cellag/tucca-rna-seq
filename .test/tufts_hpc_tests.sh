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
DEFAULT_CONDA_PREFIX="/cluster/tufts/kaplanlab/bbromb01/tucca-rna-seq-dev/envs"
DEFAULT_SINGULARITY_PREFIX="/cluster/tufts/kaplanlab/bbromb01/tucca-rna-seq-dev/envs"
PROFILE_PROD_PATH="profiles/slurm"
PROFILE_DEV_PATH="profiles/slurm-dev"

# Function to print usage information
print_usage() {
  echo "Usage: bash $0 <task> [options]  OR  ./$0 <task> [options] (if executable)"
  echo ""
  echo "This script runs various Snakemake workflow tests."
  echo ""
  echo "Tasks:"
  echo "  lint                             Lint the Snakemake workflow."
  echo "  conda-create-envs-only           Create Conda environments only."
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
  echo "  -d                               Use development Slurm profile (${PROFILE_DEV_PATH})."
  echo "                                   Default is production profile (${PROFILE_PROD_PATH})."
  echo "  -p                               Use default prefixes for Snakemake (for tasks: conda-create-envs-only, test, sra-reads):"
  echo "                                     --conda-prefix ${DEFAULT_CONDA_PREFIX}"
  echo "                                     --singularity-prefix ${DEFAULT_SINGULARITY_PREFIX}"
  echo "                                   If custom paths are needed for these, use the -x option instead of/in addition to -p."
  echo "  -x \"<extra_args>\"            Optional extra arguments for Snakemake (quote if multiple)."
  echo "                                   Example: -x \"--cores 1 --conda-prefix /custom/conda --singularity-prefix /custom/singularity\""
  echo ""
  echo "Examples:"
  echo "  bash $0 test -t rule1             # No default prefixes used by -p"
  echo "  bash $0 test -p                   # Uses default conda and singularity prefixes"
  echo "  bash $0 test -d                   # Uses development Slurm profile"
  echo "  bash $0 test -x \"--conda-prefix /my/custom/conda\" # Uses custom conda prefix via -x"
  echo "  bash $0 test -p -x \"--conda-prefix /my/custom/conda\" # -p ignored for conda, default singularity used, custom conda from -x"
  echo "  bash $0 -h"
}

# --- Initial Checks ---
# Handle -h globally before task/option parsing
if [ -z "${1}" ] || [ "${1}" == "-h" ]; then
  print_usage
  exit 0
fi

# --- Task Definition ---
if [[ "${1}" =~ ^- ]]; then
  echo "Error: Task must be the first argument." >&2
  print_usage
  exit 1
fi
TASK=$1
shift # Remove task from arguments, so $@ contains only options now

# --- Option Parsing ---
# Initialize variables for options
TARGET_OPT="all"
EXTRA_SNAKEMAKE_ARGS_OPT=""
CONFIG_SUBDIR="config_basic"
USE_REFSEQ_SOURCE=false
USE_DEFAULT_P_PREFIXES=false # Changed variable name for clarity
USE_DEV_PROFILE=false

OPTIND=1 # Reset OPTIND
while getopts ":t:x:crdhp" opt; do
  case ${opt} in
  t) TARGET_OPT="${OPTARG}" ;;
  x) EXTRA_SNAKEMAKE_ARGS_OPT="${OPTARG}" ;;
  c) CONFIG_SUBDIR="config_complex" ;;
  r) USE_REFSEQ_SOURCE=true ;;
  d) USE_DEV_PROFILE=true ;;
  p) USE_DEFAULT_P_PREFIXES=true ;; # -p is a flag to use default for both
  h)
    print_usage
    exit 0
    ;;
  \?)
    echo "Invalid option: -${OPTARG}" >&2
    print_usage
    exit 1
    ;;
  :)
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
# Set PROFILE based on the -d flag
PROFILE="${PROFILE_PROD_PATH}" # Default to production profile
if [ "${USE_DEV_PROFILE}" = true ]; then
  PROFILE="${PROFILE_DEV_PATH}"
  echo "Using development Slurm profile: ${PROFILE}"
else
  echo "Using production Slurm profile: ${PROFILE}"
fi

# --- Prepare Snakemake dynamic arguments and informational message ---
SNAKEMAKE_DYNAMIC_ARGS=()
PATH_MESSAGE_SEGMENTS=()

if [ "${USE_DEFAULT_P_PREFIXES}" = true ]; then
  # Conda Prefix
  if [[ "${EXTRA_SNAKEMAKE_ARGS_OPT}" != *"--conda-prefix"* ]]; then
    SNAKEMAKE_DYNAMIC_ARGS+=("--conda-prefix" "${DEFAULT_CONDA_PREFIX}")
    PATH_MESSAGE_SEGMENTS+=("default conda prefix: ${DEFAULT_CONDA_PREFIX}")
  else
    PATH_MESSAGE_SEGMENTS+=("--conda-prefix from -x")
  fi
  # Singularity Prefix
  if [[ "${EXTRA_SNAKEMAKE_ARGS_OPT}" != *"--singularity-prefix"* ]]; then
    SNAKEMAKE_DYNAMIC_ARGS+=("--singularity-prefix" "${DEFAULT_SINGULARITY_PREFIX}")
    PATH_MESSAGE_SEGMENTS+=("default singularity prefix: ${DEFAULT_SINGULARITY_PREFIX}")
  else
    PATH_MESSAGE_SEGMENTS+=("--singularity-prefix from -x")
  fi
else
  PATH_MESSAGE_SEGMENTS+=("no default prefixes via -p")
  if [[ "${EXTRA_SNAKEMAKE_ARGS_OPT}" == *"--conda-prefix"* ]]; then
    PATH_MESSAGE_SEGMENTS+=("--conda-prefix from -x detected")
  fi
  if [[ "${EXTRA_SNAKEMAKE_ARGS_OPT}" == *"--singularity-prefix"* ]]; then
    PATH_MESSAGE_SEGMENTS+=("--singularity-prefix from -x detected")
  fi
fi

FINAL_PATH_MESSAGE=""
if [ ${#PATH_MESSAGE_SEGMENTS[@]} -gt 0 ]; then
  FINAL_PATH_MESSAGE=" ("
  for i in "${!PATH_MESSAGE_SEGMENTS[@]}"; do
    FINAL_PATH_MESSAGE+="${PATH_MESSAGE_SEGMENTS[$i]}"
    if [ $i -lt $((${#PATH_MESSAGE_SEGMENTS[@]} - 1)) ]; then
      FINAL_PATH_MESSAGE+="; "
    fi
  done
  FINAL_PATH_MESSAGE+=")"
fi

# --- Task Execution ---
case $TASK in
lint)
  echo "Linting the Snakemake workflow..."
  # Lint task does not use --conda-prefix or --singularity-prefix from -p
  snakemake --lint --verbose --workflow-profile ${PROFILE} ${EXTRA_SNAKEMAKE_ARGS_OPT}
  ;;
conda-create-envs-only)
  echo "Running snakemake --conda-create-envs-only${FINAL_PATH_MESSAGE}"
  snakemake --conda-create-envs-only --workflow-profile ${PROFILE} "${SNAKEMAKE_DYNAMIC_ARGS[@]}" ${EXTRA_SNAKEMAKE_ARGS_OPT}
  ;;
test)
  SMK_TARGET=${TARGET_OPT:-all}
  SOURCE=$([ "${USE_REFSEQ_SOURCE}" = true ] && echo "refseq" || echo "ensembl")

  echo "Dry-run for target '${SMK_TARGET}' using a '${SOURCE}' assembly with '${CONFIG_SUBDIR}' config${FINAL_PATH_MESSAGE}..."
  if snakemake ${SMK_TARGET} -np --workflow-profile ${PROFILE} \
    --configfile ".test/local_reads/${SOURCE}/${CONFIG_SUBDIR}/config.yaml" \
    --config api_keys="{\"ncbi\": \"${API_KEY}\"}" "${SNAKEMAKE_DYNAMIC_ARGS[@]}" ${EXTRA_SNAKEMAKE_ARGS_OPT}; then
    echo "The dry-run for target '${SMK_TARGET}' on '${SOURCE}' assembly with '${CONFIG_SUBDIR}' was successful."
    echo "Running Snakemake workflow for target '${SMK_TARGET}' using a '${SOURCE}' assembly with '${CONFIG_SUBDIR}' config${FINAL_PATH_MESSAGE}..."
    snakemake ${SMK_TARGET} --workflow-profile ${PROFILE} \
      --configfile ".test/local_reads/${SOURCE}/${CONFIG_SUBDIR}/config.yaml" \
      --config api_keys="{\"ncbi\": \"${API_KEY}\"}" "${SNAKEMAKE_DYNAMIC_ARGS[@]}" ${EXTRA_SNAKEMAKE_ARGS_OPT}
  else
    DRY_RUN_EXIT_CODE=$?
    echo "The dry-run for target '${SMK_TARGET}' on '${SOURCE}' assembly with '${CONFIG_SUBDIR}' failed with exit code ${DRY_RUN_EXIT_CODE}. Skipping actual run." >&2
    exit ${DRY_RUN_EXIT_CODE}
  fi
  ;;
sra-reads)
  SMK_TARGET=${TARGET_OPT:-all}

  echo "Dry-run on SRA reads (target: '${SMK_TARGET}')${FINAL_PATH_MESSAGE}..."
  if snakemake ${SMK_TARGET} -np --workflow-profile ${PROFILE} \
    --configfile ".test/sra_reads/config/config.yaml" \
    --config api_keys="{\"ncbi\": \"${API_KEY}\"}" "${SNAKEMAKE_DYNAMIC_ARGS[@]}" ${EXTRA_SNAKEMAKE_ARGS_OPT}; then
    echo "The dry-run on SRA reads (target: '${SMK_TARGET}') was successful!!"
    echo "Running Snakemake workflow on SRA reads (target: '${SMK_TARGET}')${FINAL_PATH_MESSAGE}..."
    snakemake ${SMK_TARGET} --workflow-profile ${PROFILE} \
      --configfile ".test/sra_reads/config/config.yaml" \
      --config api_keys="{\"ncbi\": \"${API_KEY}\"}" "${SNAKEMAKE_DYNAMIC_ARGS[@]}" ${EXTRA_SNAKEMAKE_ARGS_OPT}
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
