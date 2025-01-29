# common.smk

import glob
import os
import shutil
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
    # print(f"Getting fq files for {wildcards.sample} {wildcards.unit}")
    # Check if sample is an SRA read
    if pd.isna(u["fq1"]):
        accession = u["sra"]
        return "data/pe/{accession}_{read}.fastq".format(
            accession=accession, read=wildcards.read[1]
        )
    else:
        if wildcards.read == "R1":
            return u.fq1
        elif wildcards.read == "R2":
            return u.fq2
        else:
            raise ValueError("Invalid read direction: {}".format(wildcards.read))


def get_paired_reads(wildcards):
    """
    Retrieve the file paths for paired-end reads for a given sample and unit.
    """
    u = units.loc[(wildcards.sample, wildcards.unit)]

    # Check if sample is an SRA read
    if pd.isna(u["fq1"]) and pd.isna(u["fq2"]):
        # SRA-based sample; link to the download_sra rule
        accession = u["sra"]
        fq1 = f"data/pe/{accession}_1.fastq".format(accession=u.sra)
        fq2 = f"data/pe/{accession}_2.fastq"
        return [fq1, fq2]
    else:
        # Regular paired-end sample
        return [u.fq1, u.fq2]


def is_paired_end(sample):
    """
    TODO: This function is deprecated

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
    """ 
    sample_units = units.loc[sample].dropna()
    paired = sample_units["fq2"].notna()
    all_paired = paired.all()
    all_single = (~paired).all()
    assert (
        all_single or all_paired
    ), "Mixed paired-end and single-end reads found for sample {}.".format(sample)
    return all_paired """
    pass


def get_read_direction(filename, convention):
    """
    Determine the read direction (R1 or R2) from a given filename based on the
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
    >>> get_read_direction("sample_R1_.fastq.gz", "standard")
    'R1'
    >>> get_read_direction("sample_1.fastq.gz", "numeric")
    'R1'
    """
    # Check if sample is an SRA read
    # TODO: SRA reads break entire logic of this function. Refactor
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
    else:
        raise ValueError(
            "Filename does not match any known read convention: {}".format(filename)
        )


def cp_config_to_res_dir():
    """
    Copy the contents of the 'config' directory to the
    'results/last_run_config_snapshot' directory.

    This function creates the 'results/last_run_config_snapshot' directory if
    it does not already exist. It then iterates over all items in the 'config'
    directory. For each item, it checks if the item is a directory or a file.
    If it is a directory, it uses `shutil.copytree` to copy the entire directory
    tree to the destination. If it is a file, it uses `shutil.copy2` to copy
    the file to the destination.

    The `dirs_exist_ok=True` parameter in `shutil.copytree` ensures that
    existing directories are not overwritten, and the `exist_ok=True` parameter
    in `os.makedirs` ensures that no error is raised if the directory already
    exists.

    Example:
    >>> cp_config_to_res_dir()
    This will copy all files and directories from 'config' to
    'results/last_run_config_snapshot'.

    Raises:
    OSError: If an error occurs while creating the directory or copying files.
    """
    os.makedirs("results/last_run_config_snapshot", exist_ok=True)
    for item in os.listdir("config"):
        s = os.path.join("config", item)
        d = os.path.join("results/last_run_config_snapshot", item)
        if os.path.isdir(s):
            shutil.copytree(s, d, dirs_exist_ok=True)
        else:
            shutil.copy2(s, d)


def get_final_output():
    """
    Generate a list of final output file paths for FastQC results and BAM files
    based on the units DataFrame.

    This function iterates over each row in the `units` DataFrame and constructs
    the paths for the FastQC HTML and ZIP files for both read 1 and read 2,
    as well as the BAM file, based on the naming convention specified in
    each row.

    Returns:
    list: A list of strings representing the file paths for the FastQC HTML,
            ZIP files, and BAM files.

    TODO: finish this function def
    """
    # Copy config files to the results directory
    cp_config_to_res_dir()

    final_output = []

    # Iterate over each unit to collect expected outputs
    for index, row in units.iterrows():
        convention = row["convention"]

        # Define FastQC outputs
        read1_fq_html = f"results/fastqc/{row.sample_name}_{row.unit_name}_R1.html"
        read2_fq_html = f"results/fastqc/{row.sample_name}_{row.unit_name}_R2.html"
        read1_fq_zip = f"results/fastqc/{row.sample_name}_{row.unit_name}_R1_fastqc.zip"
        read2_fq_zip = f"results/fastqc/{row.sample_name}_{row.unit_name}_R2_fastqc.zip"

        # Define STAR outputs
        bam_file = f"results/star/{row.sample_name}_{row.unit_name}_Aligned.sortedByCoord.out.bam"
        log_final = f"results/star/{row.sample_name}_{row.unit_name}_Log.final.out"
        log_out = f"results/star/{row.sample_name}_{row.unit_name}_Log.out"
        log_progress = (
            f"results/star/{row.sample_name}_{row.unit_name}_Log.progress.out"
        )
        sj_out = f"results/star/{row.sample_name}_{row.unit_name}_SJ.out.tab"
        star_tmp = f"results/star/{row.sample_name}_{row.unit_name}__STARtmp"

        # Define Qualimap outputs
        qualimapReport = f"results/qualimap/{row.sample_name}_{row.unit_name}.qualimap/qualimapReport.html"
        qualimap_qc_results = f"results/qualimap/{row.sample_name}_{row.unit_name}.qualimap/rnaseq_qc_results.txt"

        # Define Salmon outputs
        salmon_quant = (
            f"results/salmon/{row.sample_name}_{row.unit_name}.salmon/quant.sf"
        )

        # Aggregate all unit-based outputs
        final_output.extend(
            [
                read1_fq_html,
                read2_fq_html,
                read1_fq_zip,
                read2_fq_zip,
                bam_file,
                log_final,
                log_out,
                log_progress,
                sj_out,
                star_tmp,
                qualimapReport,
                qualimap_qc_results,
                salmon_quant,
            ]
        )

    # Define non-unit-based outputs
    multiqc = f"results/multiqc/{config['params']['multiqc']['report_name']}.html"

    # Aggregate all non-unit-based outputs
    final_output.extend([multiqc])

    return final_output
