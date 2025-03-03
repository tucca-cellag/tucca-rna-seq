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
        "results/last_run_config_snapshot/snapshot_done.txt",
    container:
        "docker://python:latest"
    log:
        "logs/snapshot_config/snapshot_config.log",
    script:
        "workflow/scripts/snapshot_config.py"
