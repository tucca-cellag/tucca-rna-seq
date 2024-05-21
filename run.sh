#!/bin/bash

git pull
module load miniconda/23.10
source activate base
source activate snakemake
snakemake unzip_genome \
  --workflow-profile profiles/slurm \
  -np --rerun-incomplete
