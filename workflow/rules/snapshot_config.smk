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
    run:
        cp_config_to_res_dir()
        # Write a marker file to signal completion of the snapshot
        with open(output[0], "w") as f:
            f.write("Configuration snapshot completed.")
