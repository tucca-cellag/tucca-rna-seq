#!/bin/bash

git pull
module load miniconda/23.10
source activate base
source activate snakemake
snakemake salmon_decoys \
  --cores 12 \
  --sdm conda \
  --conda-frontend conda \
  -np --rerun-incomplete
