# Snakemake workflow: `tucca-rna-seq`

# `THIS REPO IS STILL UNDER CONSTRUCTION AND DOES NOT REPRESENT A COMPLETED PIPELINE`

[![Snakemake](https://img.shields.io/badge/snakemake-≥6.3.0-brightgreen.svg)](https://snakemake.github.io)
[![GitHub license](https://img.shields.io/github/license/benjibromberg/tucca-rna-seq?color=orange)](https://github.com/benjibromberg/tucca-rna-seq/blob/main/LICENSE)

The Tufts University Center for Cellular Agriculture's ([TUCCA][tucca])
RNA-Seq Snakemake Workflow for Cellular Agriculture Projects

* Currently only supporting reads2counts

# Installation

## Installation For Tufts HPC Users

### 1. Login to the Tufts HPC

* [Tufts HPC]
  * Main page for the Tufts HPC with helpful cheat sheets for the cluster at
    the bottom of the page
* [Tufts HPC's Welcome Page]
  * Contains links to many introductory resources for the Tufts HPC
* [Tufts HPC's New User Support]
  * Box folder with many helpful introductory PDFs. The most recent version of
    the `Intro to Tufts HPC` workshop is highly recommend.

### 2. Get off of a login node

Once logged into the Tufts HPC run the following to get off a login node and
onto a compute node so that you can start performing computations.

```bash
srun -p interactive --pty bash
```

This will allocate a compute node with bash shell, 1 CPU core (default),
4 hours (interactive partition default and maxmimum value), and 2GB of CPU
memory (default). If desired, modify with number of CPU cores `-n` , memory
`--mem=` , an alternative partition `-p` , and/or runtime `--time=` .

Depending on the current usage of the Tufts HPC, you may need to wait in a
queue to be allocated enough resources.

### 3. Clone the repository

Go to the desired directory/folder on your file system on the Tufts HPC, then
clone/get the repository and move into the respective directory with:

```bash
git clone https://github.com/benjibromberg/tucca-rna-seq.git
```

### 4. Dependencies installation

For improved reproducibility and reusability of the workflow,
each individual step of the workflow runs either in its own [Conda][conda]
virtual environment. As a consequence, running this workflow has very few
individual dependencies.

You will need to add the `bioconda` channel to your `base` conda environment
which loads everytime you load a conda manager like [miniforge][miniforge].
If this is problematic for some reason, you will need to load another conda
environment that includes the `conda-forge`, `bioconda`, and `default` channels
before loading the `snakemake` conda environment which will require manually
modifying the `run.sh` runner script.

To add the `bioconda` channel to your `base` conda environment run:

```bash
module load miniforge/24.7.1-py312
conda activate base
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict
conda config --show channels
```

The result of the final command should match the following:

```bash
channels:
  - conda-forge
  - bioconda
  - defaults
```

Then run the following to create the `snakemake` conda environment:

```bash
conda env create -f tucca-rna-seq/install/snakemake.yaml
```

# Set-Up the `config` Directory

TODO: add this section

# Running the Workflow

## Running the Workflow For Tufts HPC Users

### Load the Required Dependencies

```bash
cd tucca-rna-seq
module load miniforge/24.7.1-py312
conda activate base
conda activate snakemake
```

### To run the entire workflow

```bash
snakemake all <extra-snakemake-parameters> --workflow-profile profiles/slurm
```

### To perform a dry-run of the entire workflow

```bash
snakemake all -np <extra-snakemake-parameters> --workflow-profile profiles/slurm
```

When performing a [dry-run][dry-run] using the `-n` (or `--dry-run`) flag,
Snakemake will only show the execution plan instead of actually performing the
steps. The `-p` flag instructs Snakemake to also print the resulting shell
command for illustration.

### To run the workflow targeting a specific target rule

```bash
snakemake <snakemake-target-rule> <extra-snakemake-parameters> --workflow-profile profiles/slurm
```

### Workflow Visualization

To visualize the workflow, one can use the option `--dag`. This creates a
representation of the [DAG][dag-snake] in the graphviz dot language which has
to be postprocessed by the graphviz tool dot. This yields a
[directed acyclic graph][dag-wiki] (DAG) of jobs where the edges represent
dependencies. To visualize the DAG that would be executed running the workflow,
you can execute the following command optionally targeting a specific target
rule:

```bash
snakemake <optional-target-rule> --dag | dot -Tpdf > dag.pdf
```

## Help
If something is unclear or broken in the workflow, don't hesitate to
[file an issue in the `tucca-rna-seq` GitHub repository]. We will try to be
prompt with our assistance.

[conda]: <https://docs.conda.io/projects/conda/en/latest/index.html>
[sample-doc]: pipeline_documentation.md#read-sample-table
[snakemake]: <https://snakemake.readthedocs.io/en/stable/>
[dry-run]: <https://learn.flowdeploy.com/snakemake-dry-run>
[slurm]: <https://slurm.schedmd.com/documentation.html>
[Tufts HPC's Welcome Page]: <https://it.tufts.edu/high-performance-computing/hpc-welcome-page>
[Tufts HPC's New User Support]: <https://tufts.app.box.com/v/HPC-New-User/folder/46988375653>
[Tufts HPC]: <https://it.tufts.edu/high-performance-computing>
[miniforge]: <https://github.com/conda-forge/miniforge>
[dag-snake]: <https://snakemake.readthedocs.io/en/stable/executing/cli.html#visualization>
[dag-wiki]: <https://en.wikipedia.org/wiki/Directed_acyclic_graph>
[tucca]: <https://cellularagriculture.tufts.edu/>
[rna-seq-star-deseq2 Snakemake workflow]: https://github.com/snakemake-workflows/rna-seq-star-deseq2/tree/b3998c158a87cc9096f7cda8ae913adf2ac6da9d
[rna-seq-kallisto-sleuth Snakmake workflow]: https://github.com/snakemake-workflows/rna-seq-kallisto-sleuth/tree/main
[file an issue in the `tucca-rna-seq` GitHub repository]: https://github.com/benjibromberg/tucca-rna-seq/issues/new