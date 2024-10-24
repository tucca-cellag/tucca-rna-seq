# Snakemake workflow: `tucca-rna-seq`

# `THIS REPO IS STILL UNDER CONSTRUCTION AND DOES NOT REPRESENT A COMPLETED PIPELINE`

[![Snakemake](https://img.shields.io/badge/snakemake-≥6.3.0-brightgreen.svg)](https://snakemake.github.io)
[![GitHub actions status](https://github.com/<owner>/<repo>/workflows/Tests/badge.svg?branch=main)](https://github.com/<owner>/<repo>/actions?query=branch%3Amain+workflow%3ATests)


TUCCA's Automated RNA-Seq Snakemake Pipeline using Salmon and DESeq2

# Installation for Tufts HPC Users

## 1. Get off of a login node

```bash
srun -p preempt -t 1-2:30:00 -n 2 --mem=32G --gres=gpu:1 --pty bash
```

## 2. Clone the repository

Go to the desired directory/folder on your file system on the Tufts HPC, then
clone/get the repository and move into the respective directory with:

```bash
git clone https://github.com/benjibromberg/tucca-rna-seq.git
```

## 3. Dependencies installation

For improved reproducibility and reusability of the workflow,
each individual step of the workflow runs either in its own [Conda][conda]
virtual environemnt. As a consequence, running this workflow has very few 
individual dependencies.

To install the necessary dependencies run:

```bash
module load miniforge/24.7.1-py312
conda env create -f tucca-rna-seq/install/base.yaml
conda env create -f tucca-rna-seq/install/snakemake.yaml
```

# Set-Up the `config` Directory

TODO: add this section

# Running the Pipeline

To run the pipeline:

```bash
cd tucca-rna-seq
module load miniforge/24.7.1-py312
sh run.sh <snakemake_target_rule> <extra_snakemake_parameters>
```

To perform a dry-run of the pipeline:

```bash
cd tucca-rna-seq
module load miniforge/24.7.1-py312
sh run.sh <snakemake_target_rule> -np <extra_snakemake_parameters>
```

When performing a dry-run using the -n (or --dry-run) flag, Snakemake will only
show the execution plan instead of actually performing the steps. The -p flag
instructs Snakemake to also print the resulting shell command for illustration.

# TODO: necessary Fixes 

* Hello world
* Hello world 2

[conda]: <https://docs.conda.io/projects/conda/en/latest/index.html>
[rule-graph]: images/rule_graph.svg
[sample-doc]: pipeline_documentation.md#read-sample-table
[snakemake]: <https://snakemake.readthedocs.io/en/stable/>
[dry-run]: <>
[slurm]: <https://slurm.schedmd.com/documentation.html>
