#!/bin/bash

module load miniconda/23.10
conda activate snakemake
snakemake salmon_decoys \
  --cores 12 \
  --sdm conda \
  -np --rerun-incomplete
