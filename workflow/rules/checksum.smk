# workflow/rules/checksum.smk
#
# This rule validates the checksum of a FASTQ file against a provided .md5
# manifest file. It uses md5sum -c.
# The rule produces a flag file upon successful validation. If the checksum
# does not match, the rule will fail, halting the workflow.


rule validate_checksum:
    input:
        fq=get_fq_files,
        md5=get_md5_file,
    output:
        flag=touch("results/checksums/{sample_unit}_{read}.valid"),
    log:
        "logs/checksums/{sample_unit}_{read}.log",
    shell:
        # We cd into the directory of the fastq file because the .md5 file
        # typically contains relative paths. --status suppresses normal output
        # and returns an exit code, which is ideal for scripting.
        # All output (stdout and stderr) is redirected to the log file.
        "(cd $(dirname {input.fq}) && md5sum --status -c $(basename {input.md5})) &> {log}"
