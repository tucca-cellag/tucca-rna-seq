#!/usr/bin/env python3
"""
This script snapshots the configuration directory by copying its contents
into the results snapshot directory. It uses the cp_config_to_res_dir()
function from common.smk.

File: workflow/scripts/snapshot_config.py
"""

import shutil
import sys
from pathlib import Path

CONFIG_DIR = Path("config")
RESULT_SNAPSHOT = Path("results/last_run_config_snapshot")


def cp_config_to_res_dir() -> None:
    """
    Copy configuration directory contents to a results snapshot directory.
    """
    RESULT_SNAPSHOT.mkdir(parents=True, exist_ok=True)
    for item in CONFIG_DIR.iterdir():
        destination = RESULT_SNAPSHOT / item.name
        if item.is_dir():
            shutil.copytree(item, destination, dirs_exist_ok=True)
        else:
            shutil.copy2(item, destination)


def main():
    cp_config_to_res_dir()
    print("Configuration snapshot completed.")
    # The marker file is provided as the first command-line argument.
    marker_file = sys.argv[1]
    with open(marker_file, "w", encoding="utf-8") as f:
        f.write("Configuration snapshot completed.")


if __name__ == "__main__":
    main()
