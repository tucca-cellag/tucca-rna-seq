# Snakemake workflow: ngs-test-data

[![Snakemake](https://img.shields.io/badge/snakemake-≥4.2-brightgreen.svg)](https://snakemake.bitbucket.io)
[![Build Status](https://travis-ci.org/snakemake-workflows/ngs-test-data.svg?branch=master)](https://travis-ci.org/snakemake-workflows/ngs-test-data)

This workflow creates small test datasets for NGS data analyses. The generated data is available in the folders `ref` and `reads`, such that the repository can be directly used as a git submodule for continuous integration tests.

## Authors

* Johannes Köster (@johanneskoester), https://koesterlab.github.io

## Usage

### Step 1: Install workflow

If you simply want to use this workflow, download and extract the [latest release](https://github.com/snakemake-workflows/ngs-test-data/releases).
If you intend to modify and further develop this workflow, fork this reposity. Please consider providing any generally applicable modifications via a pull request.

In any case, if you use this workflow in a paper, don't forget to give credits to the authors by citing the URL of this repository and, if available, its DOI (see above).

### Step 2: Configure workflow

Configure the workflow according to your needs via editing the file `config.yaml`.

### Step 3: Execute workflow

Test your configuration by performing a dry-run via

    snakemake -n

Execute the workflow locally via

    snakemake --cores $N

using `$N` cores or run it in a cluster environment via

    snakemake --cluster qsub --jobs 100

or

    snakemake --drmaa --jobs 100

See the [Snakemake documentation](https://snakemake.readthedocs.io) for further details.

## Attribution

This workflow uses resources from the [EMBL-EBI](https://www.ebi.ac.uk) [Ensembl Biomart](https://www.ensembl.org/info/data/biomart/index.html) and the [Ensembl FTP downloads](https://www.ensembl.org/info/data/ftp/index.html).
Neither of these resources seems to have explicit licensing information mentioned, so the [general EMBL-EBI licensing](https://www.ebi.ac.uk/licencing) suggests that the [EMBL-EBI Terms of Use](https://www.ebi.ac.uk/about/terms-of-use/) apply, with no restrictions on re-use.

Additionally, data from the [`E-GEUV-1` ArrayExpress / BioStudy](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-GEUV-1) is used.
As no `License` field is included for this dataset, it also defaults to the [EMBL-EBI Terms of Use](https://www.ebi.ac.uk/about/terms-of-use/) according to [the BioStudies Help/FAQ](https://www.ebi.ac.uk/biostudies/help#).
