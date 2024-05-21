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


rule salmon_index:
    input:
        sequences=expand(
            "results/datasets/ncbi_dataset/data/{genome}/rna.fna",
            genome=config["ncbi_genome_accession"],
        )[0],
        gentrome="results/salmon/gentrome.fasta.gz",
        decoys="results/salmon/decoys.txt",
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
    threads: 12
    conda:
        "../envs/salmon.yaml"
    log:
        "logs/salmon/transcriptome_index.log",
    params:
        kmer_len=31,
    shell:
        """
        (# Index the concatenated transcriptome and genome
        salmon index -t {input.gentrome} \
        -i results/salmon/transcriptome_index \
        -d {input.decoys} \
        -p {threads} \
        -k {params.kmer_len}) &> {log}
        """
