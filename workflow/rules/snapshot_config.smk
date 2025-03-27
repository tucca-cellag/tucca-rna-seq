# workflow/rules/snapshot_config.smk


rule snapshot_config:
    input:
        multiext(
            "config/",
            "config.yaml",
            "multiqc_config.yaml",
            "README.md",
            "samples.tsv",
            "units.tsv",
        ),
    output:
        "results/last_run_config_snapshot/snapshot_done.done",
    container:
        config["containers"]["python"]
    log:
        "logs/snapshot_config/snapshot_config.log",
    shell:
        "(python3 workflow/scripts/snapshot_config.py {output}) &> {log}"
