#!/bin/bash

git pull
module load miniconda/23.10
source activate base
source activate snakemake
snakemake salmon_decoys \
  --workflow-profile profiles/slurm \
  -p --rerun-incomplete
