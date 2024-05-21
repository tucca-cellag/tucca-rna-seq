import glob


rule salmon_decoys:
    input:
        transcriptome=expand(
            "results/datasets/ncbi_dataset/data/{genome}/rna.fna",
            genome=config["ncbi_genome_accession"],
        )[0],
        genome=glob.glob(
            (
                "results/datasets/ncbi_dataset/data/{genome}/{genome}_"
                + "*"
                + "_genomic.fna"
            ).format(genome=config["ncbi_genome_accession"])
        ),
    output:
        gentrome="results/salmon/gentrome.fasta.gz",
        decoys="results/salmon/decoys.txt",
    threads: 12
    conda:
        "../envs/salmon.yaml"
    log:
        "logs/salmon/decoys.log",
    params:
        kmer_len=31,
    shell:
        """
        (#!/bin/sh \
        grep "^>" {input.genome} | cut -d " " -f 1 > {output.decoys} \
        sed -i.bak -e 's/>//g' {output.decoys} \
        cat {input.transcriptome} {input.genome} > {output.gentrome} \
        salmon index -t {output.gentrome} \
        -i salmon_index_k{params.kmer_len} \
        -d {output.decoys} \
        -p {threads} \
        -k {params.kmer_len}) &> {log}
        """
