#!/bin/bash

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
snakemake star_index \
  --workflow-profile profiles/slurm \
  -np
