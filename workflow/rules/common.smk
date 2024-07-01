import glob
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


wildcard_constraints:
    sample="|".join(samples["sample_name"]),
    unit="|".join(units["unit_name"]),
    read="R1|R2",


####### helper functions #######


def is_paired_end(sample):
    sample_units = units.loc[sample].dropna()
    paired = sample_units["fq2"].notna()
    all_paired = paired.all()
    all_single = (~paired).all()
    assert (
        all_single or all_paired
    ), "Mixed paired-end and single-end reads found for sample {}.".format(sample)
    return all_paired


def get_fq_files(wildcards):
    u = units.loc[(wildcards.sample, wildcards.unit)]
    if wildcards.read == "R1":
        return u.fq1
    elif wildcards.read == "R2":
        return u.fq2
    else:
        raise ValueError("Invalid read direction: {}".format(wildcards.read))


def get_read_from_filename(filename, convention):
    if convention == "standard":
        if "_R1_" in filename:
            return "R1"
        elif "_R2_" in filename:
            return "R2"
    elif convention == "numeric":
        if filename.endswith("_1.fq.gz"):
            return "R1"
        elif filename.endswith("_2.fq.gz"):
            return "R2"
    raise ValueError(
        "Filename does not match any known read convention: {}".format(filename)
    )


def get_final_output():
    final_output = []
    for index, row in units.iterrows():
        convention = row["convention"]
        read1_html = "results/fastqc/{}_{}_{}.html".format(
            row.sample_name, row.unit_name, get_read_from_filename(row.fq1, convention)
        )
        read2_html = "results/fastqc/{}_{}_{}.html".format(
            row.sample_name, row.unit_name, get_read_from_filename(row.fq2, convention)
        )
        read1_zip = "results/fastqc/{}_{}_{}_fastqc.zip".format(
            row.sample_name, row.unit_name, get_read_from_filename(row.fq1, convention)
        )
        read2_zip = "results/fastqc/{}_{}_{}_fastqc.zip".format(
            row.sample_name, row.unit_name, get_read_from_filename(row.fq2, convention)
        )

        final_output.extend([read1_html, read2_html, read1_zip, read2_zip])
    return final_output


def get_paired_reads(wildcards):
    """
    Given a sample name, return a list of tuples containing paired reads.
    """
    sample_units = units.loc[wildcards.sample]

    paired_reads = []
    for unit_name, unit_info in sample_units.iterrows():
        if is_paired_end(wildcards.sample):
            paired_reads.append(
                {"unit": unit_name, "fq1": unit_info.fq1, "fq2": unit_info.fq2}
            )
        else:
            raise ValueError(
                f"Single-end read encountered for sample {wildcards.sample}, unit {unit_name}. Paired-end reads are required."
            )

    return paired_reads


""" 
def get_strandedness(units):
    if "strandedness" in units.columns:
        return units["strandedness"].tolist()
    else:
        strand_list = ["none"]
        return strand_list * units.shape[0]


def get_deseq2_threads(wildcards=None):
    # https://twitter.com/mikelove/status/918770188568363008
    few_coeffs = False if wildcards is None else len(get_contrast(wildcards)) < 10
    return 1 if len(samples) < 100 or few_coeffs else 6


def get_bioc_species_name():
    first_letter = config["ref"]["species"][0]
    subspecies = config["ref"]["species"].split("_")[1]
    return first_letter + subspecies


def get_contrast(wildcards):
    return config["diffexp"]["contrasts"][wildcards.contrast] """
