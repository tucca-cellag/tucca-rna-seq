# General settings

To configure this workflow, modify `config/config.yaml` according to your needs,
following the explanations provided in the file.

## `DESeq2` differential expression analysis

To run the differential expression analysis, you must tell DESeq2 which sample
annotations to use (annotations are columns in the samples.tsv file described
below). This is done in the `config.yaml` file with the entries under
`diffexp:`.

### Sample and Unit Sheets

The sample and unit sheet setup is specified via tab-separated tabular files
(`.tsv`). Each sample refers to an actual physical sample, and replicates (both
biological and technical) may be specified as separate samples. Each sample, may
correspond to one or more sequencing units (for example if you have several runs
or lanes per sample).

Missing values can be specified by empty columns or by writing `NA`.

## Sample Sheet (`samples.tsv`)

For each sample, add a line to the sample sheet in `config/samples.tsv`. For
each sample, you will always have to specify a `sample_name`.

In addition, all `variables_of_interest` and `batch_effects` specified in the
`config/config.yaml` under the `diffexp:` entry, will have to have corresponding
columns in the `config/samples.tsv`.

Finally, the sample sheet can contain any number of additional columns, so if
you are in doubt about whether you might at some point need some metadata you
already have at hand, just add it to the sample sheet. Your future self will
thank you.

### `samples.tsv` Example

| sample_name | treatment_1 | treatment_2 | sequencing_batch |
|-------------|-------------|-------------|------------------|
| A           | untreated   | untreated   | 1                |
| B           | untreated   | treated     | 1                |
| C           | treated     | untreated   | 1                |
| D           | treated     | untreated   | 2                |
| E           | treated     | treated     | 2                |

## Unit Sheet (`units.tsv`)

For each sample, add one or more sequencing unit lines (runs, lanes or
replicates) to the unit sheet in `config/units.tsv`.

### `.fastq` file source

For each unit, you will have to define a source for your `.fastq` files.
This can be done via the columns `fq1`, `fq2` and `sra`, with either of:

1. A single `.fastq` file for single-end reads (`fq1` column only; `fq2` and
  `sra` columns present, but empty). The entry can be any path on your system,
  but we suggest something like a `raw/` data directory within your analysis
  directory.
2. Two `.fastq` files for paired-end reads (columns `fq1` and `fq2`; column
  `sra` present, but empty). As for the `fq1` column, the `fq2` column can also
  point to anywhere on your system.
3. [NOT IMPLEMENTED YET]
  A sequence read archive (SRA) accession number (`sra` column only; `fq1` and
  `fq2` columns present, but empty). The workflow will automatically download
  the corresponding `.fastq` data (currently assumed to be paired-end). The
  accession numbers usually start with SRR or ERR and you can find accession
  numbers for studies of interest with the [SRA Run Selector]. If both local
  files and an SRA accession are specified for the same unit, the local files
  will be used.

### `units.tsv` Example

| sample_name | unit_name | fq1        | fq2        | sra |
|-------------|-----------|------------|------------|-----|
| A           | lane1     | A.1.fq.gz  | A.2.fq.gz  |     |
| A           | lane2     | A2.1.fq.gz | A2.2.fq.gz |     |
| B           | lane1     | B.1.fq.gz  | B.2.fq.gz  |     |
| C           | lane1     | C.1.fq.gz  | C.2.fq.gz  |     |
| D           | lane1     | D.1.fq.gz  | D.2.fq.gz  |     |
| E           | lane1     | E.1.fq.gz  | E.2.fq.gz  |     |

## config.yaml

This file contains the general workflow configuration. Configurable options
should be explained in the comments above the respective entry or right here in
this `config/README.md` section.

[SRA Run Selector]: https://trace.ncbi.nlm.nih.gov/Traces/study/