################################################################################
# workflow/rules/common.smk - Shared configurations and helper functions for
# the workflow
#
# TODO: Table of Contents:
# ------------------
# 1. HEADER NAME (Line ~XX)
#    - description
#    - more description
################################################################################

################################################################################
#                                   IMPORTS                                    #
################################################################################

import os
from pathlib import Path
import pandas as pd
import shutil
from typing import Protocol, List
from snakemake.utils import validate
import re
import yaml
import hashlib

################################################################################
#                                GLOBAL VARIABLES                              #
################################################################################


# Define a protocol for wildcards so that we can expect
# them to have attributes: sample, unit, and read.
class Wildcard(Protocol):
    sample_unit: str
    read: str


################################################################################
#                      CONFIGURATION PARSING FUNCTIONS                         #
#         (e.g., for parsing samples.tsv, units.tsv, config.yaml complexities) #
################################################################################


# TODO: Determine why the following line triggers a linting error
validate(config, schema="../schemas/config.schema.yaml")

# Load in a pandas samples df from samples.tsv
samples = pd.read_csv(config["samples"], sep="\t", dtype={"sample_name": str})
samples["sample_name"] = samples["sample_name"].str.strip()
samples = samples.set_index("sample_name", drop=False).sort_index()
validate(samples, schema="../schemas/samples.schema.yaml")

# Load in a pandas units df from samples.tsv
units = pd.read_csv(
    config["units"], sep="\t", dtype={"sample_name": str, "unit_name": str}
)
units["sample_name"] = units["sample_name"].str.strip()
units["unit_name"] = units["unit_name"].str.strip()
# Create a combined sample_unit index for units
units["sample_unit"] = (
    units["sample_name"].astype(str).str.cat(units["unit_name"].astype(str), sep="_")
)
units = units.set_index(["sample_unit"], drop=False).sort_index()
validate(units, schema="../schemas/units.schema.yaml")

# Check that all sample_name values in units.tsv exist in samples.tsv
missing_samples = set(units["sample_name"]) - set(samples["sample_name"])
if missing_samples:
    raise ValueError(
        f"The following sample_name values in units.tsv are not found in "
        f"samples.tsv: {missing_samples}"
    )

# Check that all sample_name values in samples.tsv are used in units.tsv
unused_samples = set(samples["sample_name"]) - set(units["sample_name"])
if unused_samples:
    raise ValueError(
        f"The following sample_name values in samples.tsv are not used in "
        f"units.tsv: {unused_samples}"
    )

# Check that each (sample, unit) combination is unique.
if not units.index.is_unique:
    raise ValueError("Each (sample, unit) combination must be unique in units.tsv")

####### wildcard constraints #######


# Define global wildcard_constraints
# To learn more see:
# https://snakemake.readthedocs.io/en/stable/tutorial/additional_features.html#constraining-wildcards
wildcard_constraints:
    # Constrain the 'sample_unit' wildcard to match any of the sample_unit
    # names listed in the 'units' DataFrame.
    # This ensures that the 'sample_unit' wildcard can only take values
    # from the unique sample-unit combinations provided in `units.tsv`
    sample_unit="|".join(re.escape(u) for u in units.index.unique()),
    # Constrain the 'read' wildcard to "R1" or "R2"
    # The 'read' wildcard is only able to take these values as they represent
    # the two possible read directions in paired-end sequencing (forward and
    # reverse, respectively).
    read="R1|R2",


################################################################################
#                            UTILITY FUNCTIONS                                 #
#                (e.g., for file path manipulation, string formatting)         #
################################################################################


def get_fq_files(wildcards: Wildcard) -> str:
    """
    Used by rule fastqc to retrieve the file path for a given read direction
    (R1 or R2) for a sample/unit.

    Parameters:
        wildcards (Wildcard): An object with attributes: sample_unit, read.

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


def get_paired_reads(wildcards: Wildcard) -> List[str]:
    """
    Used by rules salmon_quant and star to retrieve file paths for paired-end
    reads for a given sample/unit.

    Parameters:
        wildcards (Wildcard): An object with attributes: sample_unit, read.

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


def get_unit_record(wildcards: Wildcard) -> pd.Series:
    """
    Retrieve a single unit record from the units DataFrame using the sample_unit
    wildcard.

    Used by get_fq_files and get_paired_reads to access unit-specific
    information.

    Parameters:
        wildcards (Wildcard): An object with the attribute: sample_unit.

    Returns:
        pd.Series: A single row from the units DataFrame corresponding to the
            sample_unit.

    Raises:
        ValueError: If the sample_unit is not found in units.tsv or if multiple
            matching entries are found.
    """
    # Clean the wildcard values.
    sample_unit = wildcards.sample_unit.strip()
    try:
        record = units.loc[(sample_unit)]
    except KeyError:
        raise ValueError(
            f"Combination of sample_name and unit_name ({sample_unit}) not "
            f"found in units.tsv"
        )

    # If the lookup returns a DataFrame instead of a Series then more than one
    # match was found.
    if isinstance(record, pd.DataFrame):
        raise ValueError(
            f"Multiple entries found for this combination of sample_name and "
            f"unit_name ({sample_unit}). Ensure that the metadata has one "
            f"unique entry per sample-unit pair."
        )
    return record


def is_sra_read(u: pd.Series) -> bool:
    """
    Determines if a unit record indicates an SRA-based read.

    Conditions:
        - u.sra is a non-empty string (after stripping whitespace)
        - u.fq1 and u.fq2 are NaN or empty strings

    Parameters:
        u (pd.Series): A unit row with keys "sra", "fq1" and "fq2".

    Returns:
        bool: True if the record represents an SRA read, False otherwise.
    """
    return (
        str(u.sra).strip() != ""
        and (pd.isna(u.fq1) or u.fq1 == "")
        and (pd.isna(u.fq2) or u.fq2 == "")
    )


def get_sra_filepath(accession: str, read: str) -> Path:
    """
    Generates the file path for an SRA read given its accession and read
    direction.

    Parameters:
        accession (str): The SRA accession number.
        read (str): The read direction (e.g., "R1" or "R2").

    Returns:
        Path: The path to the SRA fastq file.
    """
    # Use read number (assumes read is like "R1" or "R2")
    read_num = read[1]
    return Path(f"data/sra_reads/{accession}_{read_num}.fastq")


def get_md5_file(wildcards: Wildcard) -> str:
    """
    Finds the .md5 checksum file for a given FASTQ file.
    Assumes the .md5 file is co-located with the FASTQ file,
    with a .md5 extension appended.
    e.g., for /path/to/file.fq.gz, it looks for /path/to/file.fq.gz.md5
    """
    fq_path_str = get_fq_files(wildcards)
    md5_path = Path(f"{fq_path_str}.md5")
    if not md5_path.exists():
        raise FileNotFoundError(f"Required checksum file not found: {md5_path}")
    return str(md5_path)


def get_checksum_dependency(wildcards: Wildcard) -> List[str]:
    """
    Returns a list with the checksum validation flag file for a given read,
    but only if the read is from a local file (not SRA). Returns an empty list
    otherwise.
    """
    unit = get_unit_record(wildcards)
    if not is_sra_read(unit):
        return [f"results/checksums/{wildcards.sample_unit}_{wildcards.read}.valid"]
    return []


def get_paired_checksum_dependency(wildcards: Wildcard) -> List[str]:
    """
    Returns a list of checksum validation flag files for paired-end reads, but only
    if the reads are from local files (not SRA). Returns an empty list otherwise.
    """
    unit = get_unit_record(wildcards)
    if not is_sra_read(unit):
        return expand(
            "results/checksums/{sample_unit}_{read}.valid",
            sample_unit=wildcards.sample_unit,
            read=["R1", "R2"],
        )
    return []


def is_sra_subsampling_enabled() -> bool:
    """
    Determines if SRA subsampling is enabled in the configuration.

    Returns:
        bool: True if subsampling is enabled, False otherwise.
    """
    return get_sra_subsample_params()["enabled"]


def get_sra_subsample_params() -> dict:
    """
    Returns the subsampling parameters from the configuration.

    Returns:
        dict: Dictionary containing spot ID range parameters for SRA subsampling.
    """
    # Use .get() with default values to handle missing configuration
    sra_tools_config = config.get("params", {}).get("sra_tools", {})
    subsample_config = sra_tools_config.get("subsample", {})

    # Return default values if subsample config is missing
    return {
        "enabled": subsample_config.get("enabled", False),
        "min_spot_id": subsample_config.get("min_spot_id", 1),
        "max_spot_id": subsample_config.get("max_spot_id", 1000),
    }


def is_qualimap_enabled() -> bool:
    """
    Returns whether Qualimap is enabled in the configuration.

    Returns:
        bool: True if Qualimap is enabled, False otherwise.
    """
    qualimap_config = config.get("params", {}).get("qualimap_rnaseq", {})
    return qualimap_config.get("enabled", True)


def get_sra_download_rule_name() -> str:
    """
    Returns the name of the SRA download rule to use based on configuration.

    Returns:
        str: Either 'download_sra_pe_reads' or 'subsample_sra_pe_reads'
    """
    return (
        "subsample_sra_pe_reads"
        if is_sra_subsampling_enabled()
        else "download_sra_pe_reads"
    )


################################################################################
#                        GENOME DOWNLOAD HELPER FUNCTIONS                      #
#          Helper functions for single chromosome downloads and genome         #
#          configuration access                                                #
################################################################################


def get_chromosome_param():
    """
    Returns chromosome parameter if specified in config, otherwise None.

    Used by genome download rules to support single chromosome downloads
    for CI/CD testing and resource optimization.

    Returns:
        List[str] or None: List of chromosome names if specified, None otherwise.
    """
    return config["ref_assembly"].get("chromosome", None)


################################################################################
#                      DESEQ2 MULTI-ANALYSIS HELPERS                           #
#          Parses and provides access to DESeq2 analysis configurations        #
#          defined in config.yaml under 'diffexp.deseq2.analyses'              #
################################################################################

# Set global variables for differential expression analysis
DESEQ_ANALYSES_LIST = config["diffexp"]["deseq2"]["analyses"]


# Helper function to get the full configuration for a specific analysis by its name
def get_analysis_config_by_name(wildcards_analysis_name) -> dict:
    """
    Retrieve the configuration for a specific DESeq2 analysis by its name.

    Parameters:
        wildcards_analysis_name (str): The name of the DESeq2 analysis to
            retrieve.

    Returns:
        dict: The configuration for the specified analysis.

    Raises:
        ValueError: If the analysis name is not found in DESEQ_ANALYSES_LIST.

    Used by: DESeqDataSet_from_ranged_se_per_analysis
    """
    for analysis in DESEQ_ANALYSES_LIST:
        if analysis["name"] == wildcards_analysis_name:
            return analysis
    # This error should ideally be caught by schema validation or earlier config checks
    raise ValueError(
        f"Configuration for analysis '{wildcards_analysis_name}' not found in "
        f"DESEQ_ANALYSES_LIST."
    )


def get_dds_threads(wildcards: Wildcard) -> int:
    """Returns the thread count for DESeqDataSet rule."""
    return get_analysis_config_by_name(wildcards.analysis_name)["deseqdataset"]["threads"]


def get_dds_formula(wildcards: Wildcard) -> str:
    """Returns the formula for DESeqDataSet rule."""
    return get_analysis_config_by_name(wildcards.analysis_name)["deseqdataset"]["formula"]


def get_dds_min_counts(wildcards: Wildcard) -> int:
    """Returns the minimum counts for DESeqDataSet rule."""
    return get_analysis_config_by_name(wildcards.analysis_name)["deseqdataset"][
        "min_counts"
    ]


def get_dds_extra(wildcards: Wildcard) -> str:
    """Returns extra parameters for DESeqDataSet rule."""
    return get_analysis_config_by_name(wildcards.analysis_name)["deseqdataset"]["extra"]


# Set global variables for deseq2_wald_per_analysis
DESEQ_ANALYSES_NAMES = [analysis["name"] for analysis in DESEQ_ANALYSES_LIST]

# Create a list of all contrast jobs to run
# Used by deseq2_wald_per_analysis
CONTRAST_JOBS = []
for i, deseq_analysis in enumerate(DESEQ_ANALYSES_LIST):
    deseq_analysis_name = deseq_analysis["name"]
    # Each deseq_analysis MUST have a "contrasts" key (even if empty list)
    # as per the explicit definition requirement.
    for contrast_conf_item in deseq_analysis["contrasts"]:  # Direct access
        contrast_name_item = contrast_conf_item["name"]
        contrast_elements_item = contrast_conf_item["elements"]
        CONTRAST_JOBS.append(
            {
                "analysis_name": deseq_analysis_name,
                "contrast_name": contrast_name_item,
                "elements": contrast_elements_item,
                "config_index": i,  # Store index to the DESEQ_ANALYSES_LIST list
            }
        )


def get_analysis_config_by_index(index) -> dict:
    """
    Retrieve the full configuration for a specific analysis by its index in
    DESEQ_ANALYSES_LIST.

    Parameters:
        index (int): The index of the analysis in DESEQ_ANALYSES_LIST.

    Returns:
        dict: The configuration for the specified analysis.

    Used by: deseq2_wald_per_analysis
    """
    return DESEQ_ANALYSES_LIST[index]


def get_contrast_job_details(wildcards_analysis_name, wildcards_contrast_name) -> dict:
    """
    Retrieve contrast-specific job details, including the config_index, for a
    given analysis name and contrast name combination.

    Parameters:
        wildcards_analysis_name (str): The name of the analysis.
        wildcards_contrast_name (str): The name of the contrast within the
            analysis.

    Returns:
        dict: A dictionary containing contrast job details including
            analysis_name, contrast_name, elements, and config_index.

    Raises:
        ValueError: If the specified analysis and contrast combination is not
            found.

    Used by: deseq2_wald_per_analysis
    """
    for job in CONTRAST_JOBS:
        if (
            job["analysis_name"] == wildcards_analysis_name
            and job["contrast_name"] == wildcards_contrast_name
        ):
            return job
    raise ValueError(
        f"Details for contrast '{wildcards_contrast_name}' in analysis "
        f"'{wildcards_analysis_name}' not found in CONTRAST_JOBS."
    )


def get_wald_threads(wildcards: Wildcard) -> int:
    """Returns the thread count for DESeq2 Wald test rule."""
    config_idx = get_contrast_job_details(wildcards.analysis_name, wildcards.contrast_name)[
        "config_index"
    ]
    return get_analysis_config_by_index(config_idx)["wald"]["threads"]


def get_wald_deseq_extra(wildcards: Wildcard) -> str:
    """Returns DESeq extra parameters for DESeq2 Wald test rule."""
    config_idx = get_contrast_job_details(wildcards.analysis_name, wildcards.contrast_name)[
        "config_index"
    ]
    return get_analysis_config_by_index(config_idx)["wald"]["deseq_extra"]


def get_wald_shrink_extra(wildcards: Wildcard) -> str:
    """Returns shrink extra parameters for DESeq2 Wald test rule."""
    config_idx = get_contrast_job_details(wildcards.analysis_name, wildcards.contrast_name)[
        "config_index"
    ]
    return get_analysis_config_by_index(config_idx)["wald"]["shrink_extra"]


def get_wald_results_extra(wildcards: Wildcard) -> str:
    """Returns results extra parameters for DESeq2 Wald test rule."""
    config_idx = get_contrast_job_details(wildcards.analysis_name, wildcards.contrast_name)[
        "config_index"
    ]
    return get_analysis_config_by_index(config_idx)["wald"]["results_extra"]


def get_wald_contrast_elements(wildcards: Wildcard) -> List[str]:
    """Returns contrast elements for DESeq2 Wald test rule."""
    return get_contrast_job_details(wildcards.analysis_name, wildcards.contrast_name)[
        "elements"
    ]


################################################################################
#                   ENRICHMENT ENVIRONMENT HELPER FUNCTIONS                    #
################################################################################


# A mapping of species names (from config.yaml) to their official OrgDb package
# names. A user can add to this list if a new package becomes available.
ORGDB_MAPPING = {
    "Homo_sapiens": "org.Hs.eg.db",
    "Mus_musculus": "org.Mm.eg.db",
    "Rattus_norvegicus": "org.Rn.eg.db",
    "Saccharomyces_cerevisiae": "org.Sc.sgd.db",
    "Caenorhabditis_elegans": "org.Ce.eg.db",
    "Drosophila_melanogaster": "org.Dm.eg.db",
    "Danio_rerio": "org.Dr.eg.db",
    "Bos_taurus": "org.Bt.eg.db",
    "Sus_scrofa": "org.Ss.eg.db",
    "Canis_lupus_familiaris": "org.Cf.eg.db",
    "Gallus_gallus": "org.Gg.eg.db",
    "Arabidopsis_thaliana": "org.At.tair.db",
    "Macaca_mulatta": "org.Mmu.eg.db",
    "Xenopus_laevis": "org.Xl.eg.db",
    "Anopheles_gambiae": "org.Ag.eg.db",
    "Pan_troglodytes": "org.Pt.eg.db",
    "Escherichia_coli_K_12_strain_MG1655": "org.EcK12.eg.db",
    "Plasmodium_falciparum": "org.Pf.plasmo.db",
}


def get_orgdb_install_info(config):
    """
    Looks up the species in the ORGDB_MAPPING.
    Returns install info for Bioconda if found, otherwise for a local build.
    """
    species = config["ref_assembly"]["species"]
    if species in ORGDB_MAPPING:
        # Package exists, instruct to install from Bioconda
        pkg_name = ORGDB_MAPPING[species]
        return {
            "method": "bioconda",
            "source": pkg_name,
            "pkg_name": pkg_name,
            "local_build_needed": False,
        }
    else:
        # Package not found, instruct to build locally
        # The package name will be derived by AnnotationForge
        return {
            "method": "local",
            "source": "resources/enrichment/local_orgdb_build",
            "pkg_name": "built.locally",  # Placeholder, not used for install
            "local_build_needed": True,
        }


def get_kegg_organism_code(config):
    """
    Constructs the 3-letter KEGG organism code from the species name.
    Example: "Saccharomyces_cerevisiae" -> "sce"
    """
    try:
        species_name = config["ref_assembly"]["species"]
        parts = species_name.split("_")
        return (parts[0][0] + parts[1][:2]).lower()
    except (KeyError, IndexError) as e:
        raise ValueError(
            "Could not determine KEGG organism code. Ensure config['ref_assembly']['species'] "
            f"is in 'Genus_species' format. Error: {e}"
        )


def get_enrichment_outputs():
    """Generates a list of all enrichment analysis output files."""
    outputs = []
    # CONTRAST_JOBS is defined earlier in common.smk
    for job in CONTRAST_JOBS:
        analysis_name = job["analysis_name"]
        contrast_name = job["contrast_name"]
        outputs.append(
            f"resources/enrichment/{analysis_name}/{contrast_name}/gsea_results.RDS"
        )
        outputs.append(
            f"resources/enrichment/{analysis_name}/{contrast_name}/ora_results.RDS"
        )
        # Only include SPIA outputs if SPIA is enabled
        if config["enrichment"]["spia"]["enabled"]:
            outputs.append(
                f"resources/enrichment/{analysis_name}/{contrast_name}/spia_results.RDS"
            )
            outputs.append(
                f"resources/enrichment/{analysis_name}/{contrast_name}/spia_results_readable.RDS"
            )

    # Include SPIA data directory if SPIA is enabled
    if config["enrichment"]["spia"]["enabled"]:
        outputs.append("resources/enrichment/spia_data")

    return outputs


# Helper function to get dynamic inputs based on the dictionary lookup
def get_enrichment_deps(wildcards):
    info = get_orgdb_install_info(config)
    deps = {
        "dge_tsv": f"resources/deseq2/{wildcards.analysis}/{wildcards.contrast}/dge.tsv",
        "orgdb_install_flag": get_orgdb_install_flag(config),
    }
    if info["local_build_needed"]:
        deps["local_build"] = "resources/enrichment/local_orgdb_build"
        deps["tax_id"] = "resources/enrichment/tax_id.txt"
    return deps


# Helper function to get dynamic params
def get_enrichment_params(wildcards):
    """
    Constructs a dictionary of parameters for enrichment analysis rules.
    This function combines static configuration from config.yaml with
    dynamically generated values (like package names and install methods).
    """
    # Start with a copy of the enrichment section from config
    params = config["enrichment"].copy()

    # Get dynamic info for OrgDb installation
    info = get_orgdb_install_info(config)
    params["install_method"] = info["method"]
    params["install_source"] = info["source"]
    params["org_db_pkg"] = info["pkg_name"]

    # Get dynamic info for KEGG
    params["kegg_organism"] = get_kegg_organism_code(config)

    # Get species name for WikiPathways
    params["species"] = config["ref_assembly"]["species"]

    return params


def get_orgdb_install_flag(config):
    """Returns the path to the installation flag for the OrgDb package."""
    pkg_name = get_orgdb_install_info(config)["pkg_name"]
    return f"resources/enrichment/{pkg_name}.installed"


################################################################################
#                            TARGET RULE ALL HELPERS                           #
#          defines get_final_output() and it's helpers for use in              #
#        the workflow's target rule (rule all) in workflow/Snakefile           #
################################################################################


# Helper function for FastQC output paths.
def get_fastqc_paths(row: pd.Series) -> List[str]:
    sample_unit: str = row.sample_unit
    paths: List[str] = [
        f"results/fastqc/{sample_unit}_R1/{sample_unit}_R1.html",
        f"results/fastqc/{sample_unit}_R2/{sample_unit}_R2.html",
        f"results/fastqc/{sample_unit}_R1/{sample_unit}_R1_fastqc.zip",
        f"results/fastqc/{sample_unit}_R2/{sample_unit}_R2_fastqc.zip",
    ]
    return paths


# Helper function for Qualimap output paths.
def get_qualimap_paths(row: pd.Series) -> List[str]:
    """
    Returns Qualimap output paths only if Qualimap is enabled.
    Returns empty list if Qualimap is disabled to avoid requiring these files.
    """
    if not is_qualimap_enabled():
        return []

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
        final_output.extend(get_qualimap_paths(row))
        final_output.extend(get_salmon_paths(row))

    # Define non-unit-based output (e.g., DESeq2 results, MultiQC report)
    final_output.append("resources/deseq2/deseq2_analyses_complete.done")
    final_output.append("resources/enrichment/enrichment_analyses_complete.done")

    multiqc: str = f"results/multiqc/multiqc.html"
    final_output.append(multiqc)

    return final_output
