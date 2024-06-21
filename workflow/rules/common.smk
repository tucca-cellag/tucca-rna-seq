import glob
import pandas as pd

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


wildcard_constraints:
    sample="|".join(samples["sample_name"]),
    unit="|".join(units["unit_name"]),


def get_fq(wildcards):
    u = units.loc[(wildcards.sample, wildcards.unit)]
    if not is_paired_end(wildcards.sample):
        return {"fq1": f"{u.fq1}"}
    else:
        return {"fq1": f"{u.fq1}", "fq2": f"{u.fq2}"}


def is_paired_end(sample):
    sample_units = units.loc[sample]
    fq2_null = sample_units["fq2"].isnull()
    paired = ~fq2_null
    all_paired = paired.all()
    all_single = (~paired).all()
    assert (
        all_single or all_paired
    ), "invalid units for sample {}, must be all paired end or all single end".format(
        sample
    )
    return all_paired


def get_final_output():
    final_output = []
    final_output.extend(
        expand(
            ["results/fastqc/{sample}.html", "results/fastqc/{sample}_fastqc.zip"],
            sample=wildcards.sample,
        )
    )
    return final_output


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
