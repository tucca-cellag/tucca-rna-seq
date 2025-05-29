# workflow/rules/make_linked_txome.smk


rule make_linked_txome:
    input:
        index_dir=multiext(
            "results/salmon/transcriptome_index/",
            "complete_ref_lens.bin",
            "ctable.bin",
            "ctg_offsets.bin",
            "duplicate_clusters.tsv",
            "info.json",
            "mphf.bin",
            "pos.bin",
            "pre_indexing.log",
            "rank.bin",
            "refAccumLengths.bin",
            "ref_indexing.log",
            "reflengths.bin",
            "refseq.bin",
            "seq.bin",
            "versionInfo.json",
        ),
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
        jsonFile="results/salmon/transcriptome_index.json",
    params:
        source=branch(
            config["ref_assembly"]["source"] == "RefSeq",
            then=config["ref_assembly"]["source"],
            otherwise="Local{source}".format(
                source=config["ref_assembly"]["source"],
            ),
        ),
        organism=config["ref_assembly"]["species"],
        release=config["ref_assembly"]["release"],
        genome=config["ref_assembly"]["name"],
        tximeta_cache="resources/tximeta",
    conda:
        "../envs/r_env.yaml"
    log:
        "logs/tximeta/make_linked_txome.log",
    script:
        "../scripts/make_linked_txome.R"
