import yaml
import sys
import datetime
from snakemake.script import snakemake

# --- Log setup ---
# The log file path is provided by Snakemake
log_file_path = snakemake.log[0]
with open(log_file_path, "w") as log_f:
    # Redirect stdout and stderr to the log file
    sys.stdout = sys.stderr = log_f

    print(f"--- Log generated on {datetime.date.today()} ---")

    # --- Get parameters from Snakemake ---
    template_path = snakemake.input.template
    output_path = snakemake.output.env_file
    org_db_pkg = snakemake.params.org_db_pkg

    print(f"Input template: {template_path}")
    print(f"Output file: {output_path}")
    print(f"Organism DB package to add: {org_db_pkg}")

    # --- Main logic ---
    conda_pkg_name = f"bioconductor-{org_db_pkg.lower()}"

    with open(template_path) as f:
        env_config = yaml.safe_load(f)

    if conda_pkg_name not in env_config["dependencies"]:
        print(f"Adding '{conda_pkg_name}' to dependencies.")
        env_config["dependencies"].append(conda_pkg_name)
    else:
        print(f"'{conda_pkg_name}' is already in the dependency list.")

    with open(output_path, "w") as f:
        yaml.dump(env_config, f, sort_keys=False)

    print("Successfully wrote updated environment file.")
