#!/bin/bash
# launch.sh: Environment initializer for the tucca-rna-seq-1.0.0 workflow.
#
# This script will:
#   • Load the required miniforge module (version 24.11.2)
#   • Source the conda initialization script so that the conda command is
#       available
#   • Create the workflow conda environment (tucca-rna-seq) if it doesn’t exist
#   • Activate the tucca-rna-seq environment
#
# Usage (run these commands in your shell):
#   $ source launch.sh
#   $ snakemake all --workflow-profile profiles/slurm --verbose
#

# Important: This script is meant to be sourced.
# If it’s not sourced, the environment changes (e.g., conda activation) will
#   not persist.

# 1. Ensure that the 'module' command is available.
if ! command -v module &>/dev/null; then
  echo "Error: 'module' command not found. Are you on a system that uses environment modules?"
  return 1 2>/dev/null || exit 1
fi

# 2. Load the required miniforge module with error checking.
echo "Loading miniforge/24.11.2-py312 module..."
if module load miniforge/24.11.2-py312; then
  echo "Successfully loaded miniforge/24.11.2-py312"
else
  echo "Error: Failed to load miniforge/24.11.2-py312 module."
  echo "Please check that the module is available and that you typed the correct version."
  return 1 2>/dev/null || exit 1
fi

# 3. Source the conda initialization script from the miniforge installation.
CONDA_PROFILE="/cluster/tufts/hpc/tools/miniforge3/24.11.2/etc/profile.d/conda.sh"
if [ ! -f "$CONDA_PROFILE" ]; then
  echo "Error: Conda initialization script not found at ${CONDA_PROFILE}"
  return 1 2>/dev/null || exit 1
fi

echo "Sourcing conda initialization from ${CONDA_PROFILE}..."
if ! . "$CONDA_PROFILE"; then
  echo "Error: Failed to source conda initialization from ${CONDA_PROFILE}"
  return 1 2>/dev/null || exit 1
fi

# 4. Check if the workflow conda environment (tucca-rna-seq) exists.
ENV_NAME="tucca-rna-seq-1.0.0"
if ! conda env list | awk '{print $1}' | grep -q "^${ENV_NAME}$"; then
  echo "Workflow conda environment '${ENV_NAME}' not found."
  echo "Configuring base conda environment before creating '${ENV_NAME}'..."

  # Activate base so that the following configuration applies there.
  echo "Activating base environment..."
  conda activate base

  echo "Adding bioconda channel to base conda configuration..."
  conda config --add channels bioconda

  echo "Adding conda-forge channel to base conda configuration..."
  conda config --add channels conda-forge

  echo "Setting channel_priority to strict in base conda configuration..."
  conda config --set channel_priority strict

  echo "Current conda channels configuration:"
  conda config --show channels

  echo "Creating environment '${ENV_NAME}' from install/tucca-rna-seq-1.0.0.yaml..."
  conda env create -f install/tucca-rna-seq-1.0.0.yaml -y
else
  echo "Workflow conda environment '${ENV_NAME}' already exists."

  echo "Activating base environment..."
  conda activate base
fi

# 5. Activate the workflow environment.
echo "Activating conda environment: ${ENV_NAME}"
if ! conda activate "${ENV_NAME}"; then
  echo "Error: Failed to activate conda environment '${ENV_NAME}'"
  return 1 2>/dev/null || exit 1
fi

# 6. Print confirmation.
CURRENT_ENV=$(conda info --envs | grep '*' | awk '{print $1}')
echo "Environment '${CURRENT_ENV}' is now active."
echo "tucca-rna-seq Version 1.0.0 is now loaded."
echo " "
echo "You may now run the tucca-rna-seq-1.0.0 workflow using snakemake."
echo " "
echo "If you have not yet configured the workflow for your project, please read the"
echo "   README.md file in the 'config' subdirectory for instructions on how to begin."
echo " "
echo "To RUN the STANDARD version of the workflow use the following command, run:"
echo "   snakemake all"
echo " "
echo "To run the workflow so that it integrates nicely with the TUFTS HPC, run:"
echo "   snakemake all --workflow-profile profiles/slurm"
echo "Note: This profile was last confirmed to work with the TUFTS HPC as of 99/99/9999"
echo " "
echo "For a DRY-RUN, optionally add the -n (dry-run) and -p (printshellcmds) tags:"
echo "   snakemake all -np"
echo "      OR"
echo "   snakemake all -n -p"
echo " "
echo "To use more command line options, add their tags after the target rule of the workflow ('all'):"
echo "   snakemake all <additional-command-line-options>"
echo "For further documentation of additional command line options visit the following:"
echo "    https://snakemake.readthedocs.io/en/stable/executing/cli.html#all-options"
echo " "
echo "When DEBUGGING, it is recommended to add the --verbose tag to your standard or dry-runs:"
echo "   snakemake all --workflow-profile profiles/slurm --verbose"
echo " "
echo "Graph Visualization Instructions:"
echo "  To generate a graphical representation of the workflow rules (rulegraph), run:"
echo "      snakemake --rulegraph | dot -Tpng > images/rulegraph.png"
echo "  This will create an image showing the relationships between the workflow rules."
echo " "
echo "  To generate the Directed Acyclic Graph (DAG) of the workflow, run:"
echo "      snakemake --dag | dot -Tpng > images/dag.png"
echo "  This image shows the dependency tree of the workflow and can be useful for debugging or understanding the structure."
echo " "
echo "If you need further documentation for this workflow please visit:"
echo "   https://tucca-cellag.github.io/tucca-rna-seq/introduction"
echo " "
echo "Happy Snakemaking!"
