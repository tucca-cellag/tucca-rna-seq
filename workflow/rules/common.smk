import glob
import os
import pandas as pd

####### load config and sample sheets #######

# from snakemake.utils import validate

# validate(config, schema="../schemas/config.schema.yaml")

samples = (
    pd.read_csv(config["samples"], sep="\t", dtype={"sample_name": str})
    .set_index("sample_name", drop=False)
    .sort_index()
)

# validate(samples, schema="../schemas/samples.schema.yaml")

units = (
    pd.read_csv(config["units"], sep="\t", dtype={"sample_name": str, "unit_name": str})
    .set_index(["sample_name", "unit_name"], drop=False)
    .sort_index()
)

# validate(units, schema="../schemas/units.schema.yaml")


####### wildcard constraints #######


# Snakemake documentation on wildcard constraints (as of Snakemake 8.20.5)
# https://snakemake.readthedocs.io/en/stable/tutorial/additional_features.html#constraining-wildcards
wildcard_constraints:
    # Constrain the 'sample' wildcard to match any of the sample names listed
    # in the 'samples' DataFrame.
    # This ensures that the 'sample' wildcard can only take values from the
    # specified sample names.
    sample="|".join(samples["sample_name"]),
    # Constrain the 'unit' wildcard to match any of the unit names listed in the
    # 'units' DataFrame.
    # This ensures that the 'unit' wildcard can only take values from the
    # specified unit names.
    unit="|".join(units["unit_name"]),
    # Constrain the 'read' wildcard to match either "R1" or "R2".
    # This ensures that the 'read' wildcard can only take the values "R1" or
    # "R2", representing the read direction in paired-end sequencing.
    read="R1|R2",


####### helper functions #######


def get_fq_files(wildcards):
    """
    Retrieve the file path for a specified read direction (R1 or R2) for a
    given sample and unit.

    This function looks up the `units` DataFrame using the sample and unit
    specified in the `wildcards` object. It then returns the file path for the
    read direction specified in `wildcards.read`.

    Parameters:
    wildcards (object): An object containing the following attributes:
        - sample (str): The name of the sample.
        - unit (str): The name of the unit.
        - read (str): The read direction, either "R1" or "R2".

    Returns:
    str: The file path for the specified read direction.

    Raises:
    ValueError: If the read direction is not "R1" or "R2".

    Example:
    >>> wildcards = type('obj', (object,), {'sample': 'sample1', 'unit': 'unit1', 'read': 'R1'})
    >>> units = pd.DataFrame({
    ...     'sample': ['sample1', 'sample1'],
    ...     'unit': ['unit1', 'unit2'],
    ...     'fq1': ['sample1_unit1_R1.fastq.gz', 'sample1_unit2_R1.fastq.gz'],
    ...     'fq2': ['sample1_unit1_R2.fastq.gz', 'sample1_unit2_R2.fastq.gz']
    ... }).set_index(['sample', 'unit'])
    >>> get_fq_files(wildcards)
    'sample1_unit1_R1.fastq.gz'
    """
    u = units.loc[(wildcards.sample, wildcards.unit)]
    print(f"Getting fq files for {wildcards.sample} {wildcards.unit}")
    if wildcards.read == "R1":
        return u.fq1
    elif wildcards.read == "R2":
        return u.fq2
    else:
        raise ValueError("Invalid read direction: {}".format(wildcards.read))


def get_paired_reads(wildcards):
    """
    Given a sample name, return a list of dictionaries containing paired reads.

    This function retrieves the units associated with the specified sample from
    the `units` DataFrame. It then checks if the sample is paired-end using the
    `is_paired_end` function. If the sample is paired-end, it constructs a list
    of dictionaries, each containing the unit name and the paths to the paired
    read files (fq1 and fq2). If the sample is not paired-end, a ValueError is
    raised.

    Parameters:
    wildcards (object): An object containing the sample name as an attribute
                        (wildcards.sample).

    Returns:
    list: A list of dictionaries, each containing:
        - "unit" (str): The unit name.
        - "fq1" (str): The file path to the first read.
        - "fq2" (str): The file path to the second read.

    Raises:
    ValueError: If a single-end read is encountered for the specified sample.

    Example:
    >>> wildcards = type('obj', (object,), {'sample': 'sample1'})
    >>> units = pd.DataFrame({
    ...     'sample': ['sample1', 'sample1'],
    ...     'unit': ['unit1', 'unit2'],
    ...     'fq1': ['sample1_unit1_R1.fastq.gz', 'sample1_unit2_R1.fastq.gz'],
    ...     'fq2': ['sample1_unit1_R2.fastq.gz', 'sample1_unit2_R2.fastq.gz']
    ... }).set_index(['sample', 'unit'])
    >>> get_paired_reads(wildcards)
    [
        {'unit': 'unit1', 'fq1': 'sample1_unit1_R1.fastq.gz', 'fq2': 'sample1_unit1_R2.fastq.gz'},
        {'unit': 'unit2', 'fq1': 'sample1_unit2_R1.fastq.gz', 'fq2': 'sample1_unit2_R2.fastq.gz'}
    ]
    """
    sample_units = units.loc[wildcards.sample]

    paired_reads = []
    for unit_name, unit_info in sample_units.iterrows():
        if is_paired_end(wildcards.sample):
            fq1, fq2 = unit_info.fq1, unit_info.fq2
            print(
                f"Adding paired reads for {wildcards.sample}, unit {unit_name}: fq1={fq1}, fq2={fq2}"
            )
            paired_reads.append(fq1, fq2)
        else:
            raise ValueError(
                f"""
                Single-end read encountered for sample {wildcards.sample},
                unit {unit_name}. This pipeline does not currently surrort
                single-end reads. Paired-end reads are required for the
                pipeline to function.
                """
            )

    print(f"Completed getting paired reads for {wildcards.sample}: {paired_reads}")
    return paired_reads


def is_paired_end(sample):
    """
    Determine if a given sample is paired-end or single-end based on the units
    DataFrame.

    This function checks the `fq2` column in the `units` DataFrame for the
    specified sample. If all entries in the `fq2` column are non-null, the
    sample is considered paired-end. If all entries in the `fq2` column are
    null, the sample is considered single-end. If there is a mix of non-null
    and null entries, an assertion error is raised.

    Parameters:
    sample (str): The name of the sample to be checked.

    Returns:
    bool: True if the sample is paired-end, False if the sample is single-end.

    Raises:
    AssertionError: If the sample contains a mix of paired-end and single-end
    reads.

    Example:
    >>> units = pd.DataFrame({
    ...     'fq1': ['sample1_R1.fastq.gz', 'sample2_R1.fastq.gz'],
    ...     'fq2': ['sample1_R2.fastq.gz', None]
    ... }, index=['sample1', 'sample2'])
    >>> is_paired_end('sample1')
    True
    >>> is_paired_end('sample2')
    False
    """
    sample_units = units.loc[sample].dropna()
    paired = sample_units["fq2"].notna()
    all_paired = paired.all()
    all_single = (~paired).all()
    assert (
        all_single or all_paired
    ), "Mixed paired-end and single-end reads found for sample {}.".format(sample)
    return all_paired


def get_read_from_filename(filename, convention):
    """
    Determine the read type (R1 or R2) from a given filename based on the
    specified file naming convention.

    Parameters:
    filename (str): The name of the file to be checked.
    convention (str): The naming convention to be used. It can be either
                        "standard" or "numeric".

    Returns:
    str: "R1" if the filename corresponds to the first read, "R2" if it
            corresponds to the second read.

    Raises:
    ValueError: If the filename does not match any known read convention.

    Conventions:
    - "standard": The filename contains "_R1_" for the first read and "_R2_"
        for the second read. Common convention used by Genewiz.
    - "numeric": The filename ends with "_1.fq.gz" or "_1.fastq.gz" for the
        first read and "_2.fq.gz" or "_2.fastq.gz" for the second read. Common
        convention used by Novogene.

    Example:
    >>> get_read_from_filename("sample_R1_.fastq.gz", "standard")
    'R1'
    >>> get_read_from_filename("sample_1.fastq.gz", "numeric")
    'R1'
    """
    if convention == "standard":
        if "_R1_" in filename:
            return "R1"
        elif "_R2_" in filename:
            return "R2"
    elif convention == "numeric":
        if filename.endswith(("_1.fq.gz", "_1.fastq.gz")):
            return "R1"
        elif filename.endswith(("_2.fq.gz", "_2.fastq.gz")):
            return "R2"
    raise ValueError(
        "Filename does not match any known read convention: {}".format(filename)
    )


def get_final_output():
    """
    Generate a list of final output file paths for FastQC results based on the
    units DataFrame.

    This function iterates over each row in the `units` DataFrame and constructs
    the paths for the FastQC HTML and ZIP files for both read 1 and read 2,
    based on the naming convention specified in each row.

    Returns:
    list: A list of strings representing the file paths for the FastQC HTML and
            ZIP files.

    Example:
    >>> units = pd.DataFrame({
    ...     'sample_name': ['sample1', 'sample2'],
    ...     'unit_name': ['unit1', 'unit2'],
    ...     'fq1': ['sample1_unit1_R1_001.fastq.gz', 'sample2_unit2_R1_001.fastq.gz'],
    ...     'fq2': ['sample1_unit1_R2_001.fastq.gz', 'sample2_unit2_R2_001.fastq.gz'],
    ...     'convention': ['standard', 'standard']
    ... })
    >>> get_final_output()
    [
        'results/fastqc/sample1_unit1_R1.html',
        'results/fastqc/sample1_unit1_R2.html',
        'results/fastqc/sample1_unit1_R1_fastqc.zip',
        'results/fastqc/sample1_unit1_R2_fastqc.zip',
        'results/fastqc/sample2_unit2_R1.html',
        'results/fastqc/sample2_unit2_R2.html',
        'results/fastqc/sample2_unit2_R1_fastqc.zip',
        'results/fastqc/sample2_unit2_R2_fastqc.zip'
    ]
    """
    final_output = []

    # Ask for FastQC output files
    for index, row in units.iterrows():
        convention = row["convention"]
        read1_fq_html = "results/fastqc/{}_{}_{}.html".format(
            row.sample_name, row.unit_name, get_read_from_filename(row.fq1, convention)
        )
        read2_fq_html = "results/fastqc/{}_{}_{}.html".format(
            row.sample_name, row.unit_name, get_read_from_filename(row.fq2, convention)
        )
        read1_fq_zip = "results/fastqc/{}_{}_{}_fastqc.zip".format(
            row.sample_name, row.unit_name, get_read_from_filename(row.fq1, convention)
        )
        read2_fq_zip = "results/fastqc/{}_{}_{}_fastqc.zip".format(
            row.sample_name, row.unit_name, get_read_from_filename(row.fq2, convention)
        )
        bam = "results/star/{}_{}_Aligned.sortedByCoord.out.bam".format(
            row.sample_name, row.unit_name
        )

        print(f"Adding outputs for sample {row.sample_name}, unit {row.unit_name}:")
        print(f"  FastQC HTML (R1): {read1_fq_html}")
        print(f"  FastQC HTML (R2): {read2_fq_html}")
        print(f"  FastQC ZIP (R1): {read1_fq_zip}")
        print(f"  FastQC ZIP (R2): {read2_fq_zip}")
        print(f"  BAM: {bam}")

        final_output.extend(
            [read1_fq_html, read2_fq_html, read1_fq_zip, read2_fq_zip, bam]
        )

    print(f"Final output list: {final_output}")
    return final_output
