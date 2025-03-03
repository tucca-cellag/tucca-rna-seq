#!/usr/bin/env python3
"""
This script snapshots the configuration directory by copying its contents
into the results snapshot directory. It uses the cp_config_to_res_dir()
function from common.smk.
"""

from contextlib import redirect_stdout
import sys

# Import the cp_config_to_res_dir function from common.smk.
# Ensure that the workflow/rules directory is set as a module (i.e. has an __init__.py file)
from workflow.rules.common import cp_config_to_res_dir


def main():
    # Write log output if a log file is specified.
    if snakemake.log:
        # If snakemake.log is a list, use the first element.
        log_file = snakemake.log[0] if isinstance(
            snakemake.log, list) and snakemake.log else snakemake.log
        with open(log_file, "w") as logfile:
            with redirect_stdout(logfile):
                cp_config_to_res_dir()
                print("Configuration snapshot completed.")
    else:
        cp_config_to_res_dir()
        print("Configuration snapshot completed.")

    # Write the marker file to indicate the snapshot is complete.
    with open(snakemake.output[0], "w") as marker:
        marker.write("Configuration snapshot completed.")


if __name__ == "__main__":
    main()
