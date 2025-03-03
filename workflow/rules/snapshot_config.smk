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
    log:
        "logs/snapshot_config/snapshot_config.log",
    run:
        from contextlib import redirect_stdout

        with open(log[0], "w") as logfile:
            with redirect_stdout(logfile):
                cp_config_to_res_dir()
                print("Configuration snapshot completed.")
                # Also write a marker file for Snakemake’s output requirement
        with open(output[0], "w") as f:
            f.write("Configuration snapshot completed.")
