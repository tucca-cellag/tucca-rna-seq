# workflow/rules/dummy_singularity_image.smk


rule dummy_all_containers:
    input:
        expand(
            "results/singularity/dummy_{cont}_init.txt",
            cont=config["containers"].keys(),
        ),
    output:
        touch("results/singularity/dummy_all_images_init.txt"),
    container:
        config["containers"]["ubuntu"]
    log:
        "logs/singularity/dummy_all_containers.log",
    shell:
        "(echo 'All container images have been pulled.' && touch {output}) &> {log}"


rule dummy_container:
    output:
        "results/singularity/dummy_{cont}_init.txt",
    container:
        # Choose the container image from the config based on wildcards.cont
        lambda wc: config["containers"][wc.cont]
    log:
        "logs/singularity/{cont}.log",
    shell:
        """
        (echo "Pulling container: {wildcards.cont}"
        touch {output}) &> {log}
        """
