#!/bin/bash

# Ensure the script is called with at least one argument
if [ "$#" -lt 1 ]; then
  echo "Usage: $0 <target_rule> [additional_params...]"
  exit 1
fi

TARGET_RULE=$1
shift       # Shift the command-line arguments to the left, so $2 becomes $1, $3 becomes $2, and so on
PARAMS="$@" # Capture all remaining arguments

# Update the repository
git pull

# Update the repository again to ensure that the latest changes are pulled
git pull

# Create symlink to bbromb01/raw_data directory if not already made
if [ ! -L data/raw_data ]; then
  ln -s /cluster/tufts/kaplanlab/bbromb01/raw_data data/raw_data
fi

# Load conda
module load miniconda/23.10
source activate base
source activate snakemake

# Run snakemake pipeline
snakemake $TARGET_RULE $PARAMS \
  --workflow-profile profiles/slurm
