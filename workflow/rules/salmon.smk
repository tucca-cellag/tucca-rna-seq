# Why do this?
# See the note on mapping using selective alignment here:
# https://salmon.readthedocs.io/en/latest/salmon.html

# TLDR
# "we recommend using selective alignment with a decoy-aware transcriptome, to
# mitigate potential spurious mapping of reads that actually arise from some
# unannotated genomic locus that is sequence-similar to an annotated
# transcriptome."

# For more information on this see the documentation for Salmon
# https://github.com/COMBINE-lab/salmon
# https://salmon.readthedocs.io/en/latest/salmon.html#preparing-transcriptome-indices-mapping-based-mode
# https://combine-lab.github.io/alevin-tutorial/2019/selective-alignment/
# Thread on "How does salmon deal with decoy?"
# https://www.biostars.org/p/456231/

# workflow/rules/salmon.smk


rule salmon_decoys:
    input:
        transcriptome=branch(
            using_refseq_assembly,
            then="resources/datasets/ncbi_dataset/data/{genome_asc}/rna.fna".format(,
            genome_asc=config["ref_assembly"]["accession"],
            ),
            otherwise="resources/ensembl/{species}.{genome_name}.cdna.fa".format(
                species=config["ref_assembly"]["species"],
                genome_name=config["ref_assembly"]["name"],
            ),
        ),
        genome=branch(
            using_refseq_assembly,
            then="resources/datasets/ncbi_dataset/data/{genome_asc}/{genome_asc}_{genome_name}_genomic.fna".format(
                genome_asc=config["ref_assembly"]["accession"],
                genome_name=config["ref_assembly"]["name"],
            ),
            otherwise="resources/ensembl/{species}.{genome_name}.dna.fa".format(
                species=config["ref_assembly"]["species"],
                genome_name=config["ref_assembly"]["name"],
            ),
        ),
    output:
        gentrome="resources/salmon/gentrome.fasta",
        decoys="resources/salmon/decoys.txt",
    threads: 2
    log:
        "logs/salmon/decoys.log",
    wrapper:
        "v5.10.0/bio/salmon/decoys"


rule salmon_index:
    input:
        sequences="resources/salmon/gentrome.fasta",
        decoys="resources/salmon/decoys.txt",
    output:
        multiext(
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
    params:
        extra=config["params"]["salmon_index"]["extra"],
    threads: 2
    log:
        "logs/salmon/transcriptome_index.log",
    wrapper:
        "v5.10.0/bio/salmon/index"


rule salmon_quant:
    input:
        # If you have multiple fastq files for a single sample
        # (e.g. technical replicates) use a list for r1 and r2.
        r1=lambda wildcards: get_paired_reads(wildcards)[0],
        r2=lambda wildcards: get_paired_reads(wildcards)[1],
        index=multiext(
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
    output:
        quant="results/salmon/{sample}_{unit}/quant.sf",
        lib="results/salmon/{sample}_{unit}/lib_format_counts.json",
        aux_info=directory("results/salmon/{sample}_{unit}/aux_info"),
        cmd_info="results/salmon/{sample}_{unit}/cmd_info.json",
        libparams=directory("results/salmon/{sample}_{unit}/libParams"),
        logs=directory("results/salmon/{sample}_{unit}/logs"),
    params:
        libtype=config["params"]["salmon_quant"]["libtype"],
        extra=config["params"]["salmon_quant"]["extra"],
    threads: 12
    log:
        "logs/salmon/salmon_quant_{sample}_{unit}.log",
    message:
        """
        Running Salmon Quant for:
            sample = {wildcards.sample},
            unit = {wildcards.unit}
        Running Salmon Quant with the inputs:
            {input.r1}
            {input.r2}
        """
    wrapper:
        "v5.10.0/bio/salmon/quant"
