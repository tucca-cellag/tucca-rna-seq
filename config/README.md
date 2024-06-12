# General settings

To configure this workflow, modify the following files to reflect your dataset and differential expression analysis model:

* `config/samples.tsv`: samples sheet with covariates and conditions
* `config/units.tsv`: (sequencing) units sheet with raw data paths
* `config/config.yaml`: general workflow configuration and differential expression model setup

## samples sheet

For each biological sample, add a line to the sample sheet in `config/samples.tsv`.
The column `sample` is required and gives the sample name.

Missing values can be specified by empty columns or by writing `NA`.

## units sheet

For each sample, add one or more sequencing unit lines (runs, lanes or replicates) to the unit sheet in `config/units.tsv`.
For each unit, provide either of the following:

* The path to two paired read FASTQ files in the columns `fq1`, `fq2`.

Missing values can be specified by empty columns or by writing `NA`.

## config.yaml

This file contains the general workflow configuration.
Configurable options should be explained in the comments above the respective entry or right here in this `config/README.md` section.
If something is unclear, don't hesitate to [file an issue in the `tucca-rna-seq` GitHub repository](https://github.com/benjibromberg/tucca-rna-seq/issues/new).
