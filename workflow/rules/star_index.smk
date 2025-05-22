# workflow/rules/star_index.smk


rule star_index:
    input:
        fasta=branch(
            config["ref_assembly"]["source"] == "RefSeq",
            then="resources/datasets/ncbi_dataset/data/{genome_asc}/{genome_asc}_{genome_name}_genomic.fna".format(
                genome_asc=config["ref_assembly"]["accession"],
                genome_name=config["ref_assembly"]["name"],
            ),
            otherwise="resources/ensembl/{species}.{genome_name}.dna.fa".format(
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
        directory(
            "resources/star/{genome_asc}_index".format(
                genome_asc=config["ref_assembly"]["accession"],
            )
        ),
    threads: 12
    params:
        sjdbOverhang=config["params"]["star_index"]["sjdbOverhang"],
        extra=config["params"]["star_index"]["extra"],
    log:
        "logs/star/star_index.log",
    wrapper:
        "v6.2.0/bio/star/index"
