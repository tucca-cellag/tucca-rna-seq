# workflow/rules/common.smk

import os
from pathlib import Path
import pandas as pd
import shutil
from typing import Protocol, List
from snakemake.utils import validate
import re


# Define a protocol for wildcards so that we can expect
# them to have attributes: sample, unit, and read.
class Wildcard(Protocol):
    sample: str
    unit: str
    read: str


####### load config and sample sheets #######


validate(config, schema="../schemas/config.schema.yaml")

samples = pd.read_csv(config["samples"], sep="\t", dtype={"sample_name": str})
samples["sample_name"] = samples["sample_name"].str.strip()
samples = samples.set_index("sample_name", drop=False).sort_index()
validate(samples, schema="../schemas/samples.schema.yaml")

units = pd.read_csv(
    config["units"], sep="\t", dtype={"sample_name": str, "unit_name": str}
)
units["sample_name"] = units["sample_name"].str.strip()
units["unit_name"] = units["unit_name"].str.strip()
units["sample_unit"] = units["sample_name"] + "_" + units["unit_name"]
units = units.set_index(["sample_unit"], drop=False).sort_index()
validate(units, schema="../schemas/units.schema.yaml")

# Check that each (sample, unit) combination is unique.
if not units.index.is_unique:
    raise ValueError("Each (sample, unit) combination must be unique in units.tsv")

####### wildcard constraints #######


# Snakemake documentation on wildcard constraints
# https://snakemake.readthedocs.io/en/stable/tutorial/additional_features.html#constraining-wildcards
wildcard_constraints:
    # Constrain the 'sample' wildcard to match any of the sample names listed
    # in the 'samples' DataFrame.
    # This ensures that the 'sample' wildcard can only take values from the
    # specified sample names.
    #sample="|".join(re.escape(s) for s in samples.index.unique()),
    # Constrain the 'unit' wildcard to match any of the unit names listed in the
    # 'units' DataFrame.
    # This ensures that the 'unit' wildcard can only take values from the
    # specified unit names.
    # Escape each unit name to avoid regex meta-character issues.
    sample_unit="|".join(re.escape(u) for u in units.index.unique()),
    # Constrain the 'read' wildcard to match either "R1" or "R2".
    # This ensures that the 'read' wildcard can only take the values "R1" or
    # "R2", representing the read direction in paired-end sequencing.
    read="R1|R2",


####### helper functions #######


def get_fq_files(wildcards: Wildcard) -> str:
    """
    Retrieve the file path for a given read direction (R1 or R2) for a
    sample/unit.

    Parameters:
        wildcards (Wildcard): An object with attributes: sample, unit, read.

    Returns:
        str: The file path for the specified read.
    """
    u: pd.Series = get_unit_record(wildcards)
    if is_sra_read(u):
        return str(get_sra_filepath(str(u.sra), wildcards.read))
    else:
        if wildcards.read == "R1":
            return str(u.fq1)
        elif wildcards.read == "R2":
            return str(u.fq2)
        else:
            raise ValueError(f"Invalid read direction: {wildcards.read}")


def get_unit_record(wildcards: Wildcard) -> pd.Series:
    # Clean the wildcard values.
    sample_unit = wildcards.sample_unit.strip()
    try:
        record = units.loc[(sample_unit)]
    except KeyError:
        raise ValueError(
            f"Combination of sample_name and unit_name ({sample_unit}) not found in units.tsv"
        )

    # If the lookup returns a DataFrame instead of a Series then more than one
    # match was found.
    if isinstance(record, pd.DataFrame):
        raise ValueError(
            f"Multiple entries found for this combination of sample_name and unit_name ({sample_unit}). Ensure that the metadata has one unique entry per sample-unit pair."
        )
    return record


def is_sra_read(u: pd.Series) -> bool:
    """
    Determine if a unit record indicates an SRA-based read.

    Conditions:
        - u.sra is a non-empty string (after stripping whitespace)
        - u.fq1 and u.fq2 are NaN

    Parameters:
        u (pd.Series): A unit row with keys "sra", "fq1" and "fq2".

    Returns:
        bool: True if the record represents an SRA read, False otherwise.
    """
    return str(u.sra).strip() != "" and pd.isna(u.fq1) and pd.isna(u.fq2)


def get_sra_filepath(accession: str, read: str) -> Path:
    """
    Generate the file path for an SRA read given its accession and read
    direction.

    Parameters:
        accession (str): The SRA accession number.
        read (str): The read direction (e.g., "R1" or "R2").

    Returns:
        Path: The path to the SRA fastq file.
    """
    # Use read number (assumes read is like "R1" or "R2")
    read_num = read[1]
    return f"data/sra_reads/{accession}_{read_num}.fastq"


def get_paired_reads(wildcards: Wildcard) -> List[str]:
    """
    Retrieve file paths for paired-end reads for a given sample/unit.

    Parameters:
        wildcards (Wildcard): An object with attributes: sample, unit, read.

    Returns:
        List[str]: A list containing the paths for R1 and R2 reads.
    """
    u: pd.Series = get_unit_record(wildcards)
    if is_sra_read(u):
        fq1: str = str(get_sra_filepath(str(u.sra), "R1"))
        fq2: str = str(get_sra_filepath(str(u.sra), "R2"))
        return [fq1, fq2]
    else:
        return [str(u.fq1), str(u.fq2)]


#### get_final_output() helpers and definition ####


# Helper function for FastQC output paths.
def get_fastqc_paths(row: pd.Series) -> List[str]:
    sample_unit: str = row.sample_unit
    paths: List[str] = [
        f"results/fastqc/{sample_unit}_R1.html",
        f"results/fastqc/{sample_unit}_R2.html",
        f"results/fastqc/{sample_unit}_R1_fastqc.zip",
        f"results/fastqc/{sample_unit}_R2_fastqc.zip",
    ]
    return paths


# Helper function for STAR output paths.
def get_star_paths(row: pd.Series) -> List[str]:
    sample_unit: str = row.sample_unit
    paths: List[str] = [
        f"results/star/{sample_unit}_Aligned.sortedByCoord.out.bam",
        f"results/star/{sample_unit}_Log.final.out",
        f"results/star/{sample_unit}_Log.out",
        f"results/star/{sample_unit}_Log.progress.out",
        f"results/star/{sample_unit}_SJ.out.tab",
    ]
    return paths


# Helper function for Qualimap output paths.
def get_qualimap_paths(row: pd.Series) -> List[str]:
    sample_unit: str = row.sample_unit
    paths: List[str] = [
        f"results/qualimap/{sample_unit}.qualimap/qualimapReport.html",
        f"results/qualimap/{sample_unit}.qualimap/rnaseq_qc_results.txt",
    ]
    return paths


# Helper function for Salmon output paths.
def get_salmon_paths(row: pd.Series) -> List[str]:
    sample_unit: str = row.sample_unit
    return [f"results/salmon/{sample_unit}/quant.sf"]


# Main function that aggregates all expected outputs. Called by rule all.
def get_final_output() -> List[str]:
    """
    Aggregate a list of final output file paths (for FastQC, STAR, Qualimap,
    Salmon, etc.) and a MultiQC report.

    Returns:
        List[str]: A list of file paths representing the expected outputs.
    """
    final_output: List[str] = []

    # Iterate over each unit to collect outputs using helper functions.
    for _, row in units.iterrows():
        final_output.extend(get_fastqc_paths(row))
        final_output.extend(get_star_paths(row))
        final_output.extend(get_qualimap_paths(row))
        final_output.extend(get_salmon_paths(row))

    # Define non-unit-based output (e.g., MultiQC report)
    multiqc: str = f"results/multiqc/{config['params']['multiqc']['report_name']}.html"
    final_output.append(multiqc)

    return final_output
