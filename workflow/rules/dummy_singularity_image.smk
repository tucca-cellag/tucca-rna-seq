# workflow/rules/dummy_singularity_image.smk


rule dummy_all_images:
    input:
        "results/singularity/dummy_fastqc_init.txt",
        "results/singularity/dummy_p7zip_init.txt",
        "results/singularity/dummy_ncbi_datasets_init.txt",
        "results/singularity/dummy_multiqc_init.txt",
        "results/singularity/dummy_qualimap_init.txt",
        "results/singularity/dummy_salmon_init.txt",
        "results/singularity/dummy_python_init.txt",
        "results/singularity/dummy_sra_tools_init.txt",
        "results/singularity/dummy_star_init.txt",
    output:
        touch("results/singularity/dummy_all_images_init.txt"),
    container:
        config["containers"]["ubuntu"]
    log:
        "logs/singularity/dummy_all_images.log",
    shell:
        "(echo 'All container images have been pulled.' && touch {output}) &> {log}"


rule dummy_fastqc:
    input:
        [],
    output:
        touch("results/singularity/dummy_fastqc_init.txt"),
    container:
        config["containers"]["fastqc"]
    log:
        "logs/singularity/fastqc.log",
    shell:
        "(echo 'Pulling container images...' && touch {output}) &> {log}"


rule dummy_p7zip:
    input:
        [],
    output:
        touch("results/singularity/dummy_p7zip_init.txt"),
    container:
        config["containers"]["p7zip"]
    log:
        "logs/singularity/p7zip.log",
    shell:
        "(echo 'Pulling container images...' && touch {output}) &> {log}"


rule dummy_ncbi_datasets:
    input:
        [],
    output:
        touch("results/singularity/dummy_ncbi_datasets_init.txt"),
    container:
        config["containers"]["ncbi_datasets"]
    log:
        "logs/singularity/ncbi_datasets.log",
    shell:
        "(echo 'Pulling container images...' && touch {output}) &> {log}"


rule dummy_multiqc:
    input:
        [],
    output:
        touch("results/singularity/dummy_multiqc_init.txt"),
    container:
        config["containers"]["multiqc"]
    log:
        "logs/singularity/multiqc.log",
    shell:
        "(echo 'Pulling container images...' && touch {output}) &> {log}"


rule dummy_qualimap:
    input:
        [],
    output:
        touch("results/singularity/dummy_qualimap_init.txt"),
    container:
        config["containers"]["qualimap"]
    log:
        "logs/singularity/qualimap.log",
    shell:
        "(echo 'Pulling container images...' && touch {output}) &> {log}"


rule dummy_salmon:
    input:
        [],
    output:
        touch("results/singularity/dummy_salmon_init.txt"),
    container:
        config["containers"]["salmon"]
    log:
        "logs/singularity/salmon.log",
    shell:
        "(echo 'Pulling container images...' && touch {output}) &> {log}"


rule dummy_python:
    input:
        [],
    output:
        touch("results/singularity/dummy_python_init.txt"),
    container:
        config["containers"]["python"]
    log:
        "logs/singularity/python.log",
    shell:
        "(echo 'Pulling container images...' && touch {output}) &> {log}"


rule dummy_sra_tools:
    input:
        [],
    output:
        touch("results/singularity/dummy_sra_tools_init.txt"),
    container:
        config["containers"]["sra_tools"]
    log:
        "logs/singularity/sra_tools.log",
    shell:
        "(echo 'Pulling container images...' && touch {output}) &> {log}"


rule dummy_star:
    input:
        [],
    output:
        touch("results/singularity/dummy_star_init.txt"),
    container:
        config["containers"]["star"]
    log:
        "logs/singularity/star.log",
    shell:
        "(echo 'Pulling container images...' && touch {output}) &> {log}"
