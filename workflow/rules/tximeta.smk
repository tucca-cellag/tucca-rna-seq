# workflow/rules/tximeta.smk


rule tximeta:
    input:
        quant=expand(
            "results/salmon/{sample_unit}/quant.sf",
            sample_unit=units.sample_unit.values.tolist(),
        ),
        lib=expand(
            "results/salmon/{sample_unit}/lib_format_counts.json",
            sample_unit=units.sample_unit.values.tolist(),
        ),
        aux_info=expand(
            "results/salmon/{sample_unit}/aux_info/",
            sample_unit=units.sample_unit.values.tolist(),
        ),
        cmd_info=expand(
            "results/salmon/{sample_unit}/cmd_info.json",
            sample_unit=units.sample_unit.values.tolist(),
        ),
        libparams=expand(
            "results/salmon/{sample_unit}/libParams/",
            sample_unit=units.sample_unit.values.tolist(),
        ),
        logs=expand(
            "results/salmon/{sample_unit}/logs/",
            sample_unit=units.sample_unit.values.tolist(),
        ),
        linkedTxome="results/salmon/transcriptome_index.json",
        fasta=branch(
            config["ref_assembly"]["source"] == "RefSeq",
            then="resources/datasets/ncbi_dataset/data/{genome_asc}/rna.fna".format(
                genome_asc=config["ref_assembly"]["accession"],
            ),
            otherwise="resources/ensembl/{species}.{genome_name}.cdna.fa".format(
                species=config["ref_assembly"]["species"],
                genome_name=config["ref_assembly"]["name"],
            ),
        ),
        gtf=branch(
            config["ref_assembly"]["source"] == "RefSeq",
            then="resources/datasets/ncbi_dataset/data/{genome_asc}/genomic.gtf".format(
                genome_asc=config["ref_assembly"]["accession"],
            ),
            otherwise="resources/ensembl/{species}.{genome_name}.gtf".format(
                species=config["ref_assembly"]["species"],
                genome_name=config["ref_assembly"]["name"],
            ),
        ),
    output:
        se="resources/tximeta/tximeta_se.RDS",
        gse="resources/tximeta/tximeta_gse.RDS",
        image="resources/tximeta/.RData",
    params:
        sample_names=samples.index.unique().tolist(),
        factors=config["diffexp"]["tximeta"]["factors"],
        extra=config["diffexp"]["tximeta"]["extra"],
    conda:
        "../envs/r_env.yaml"
    log:
        "logs/tximeta/tximeta.log",
    script:
        "../scripts/tximeta.R"
