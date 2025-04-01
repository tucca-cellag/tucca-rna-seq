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

# workflow/rules/salmon_decoys.smk

import glob


rule salmon_decoys:
    input:
        transcriptome=expand(
            "results/datasets/ncbi_dataset/data/{genome_asc}/rna.fna",
            genome_asc=config["genome"]["assembly_accession"],
        )[0],
        genome="results/datasets/ncbi_dataset/data/{genome_asc}/{genome_asc}_{genome_name}_genomic.fna".format(
            genome_asc=config["genome"]["assembly_accession"],
            genome_name=config["genome"]["assembly_name"],
        ),
    output:
        gentrome="results/salmon/gentrome.fasta.gz",
        decoys="results/salmon/decoys.txt",
    threads: 12
    conda:
        "../envs/salmon.yaml"
    log:
        "logs/salmon/decoys.log",
    shell:
        """
        (# Preparing decoy metadata (the full genome is used as decoy)
        grep \"^>\" {input.genome} | cut -d \" \" -f 1 > {output.decoys}
        sed -i.bak -e \'s/>//g\' {output.decoys}

        # Concatenate genome to end of transcriptome to make ref file for index
        cat {input.transcriptome} {input.genome} > {output.gentrome}) &> {log}
        """
