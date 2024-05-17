#!/bin/bash

git pull
module load miniconda/23.10
source activate snakemake
snakemake salmon_decoys \
  --cores 12 \
  --sdm conda \
  -np --rerun-incomplete
