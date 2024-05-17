#!/bin/bash

git pull
module load miniconda/23.10
source activate base
source activate snakemake
snakemake salmon_decoys \
  --cores 12 \
  --conda-frontend conda \
  --sdm conda \
  --executor slurm \
  --profile workflow/profiles/slurm/config.yaml \
  -np --rerun-incomplete
