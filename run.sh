#!/bin/bash

git pull
module load miniconda/23.10
source activate base
source activate snakemake
snakemake star_index \
  --workflow-profile profiles/slurm \
  -p --rerun-incomplete
