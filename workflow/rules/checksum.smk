# workflow/rules/checksum.smk
#
# This rule validates the checksum of a FASTQ file against a provided checksum
# file. It supports two formats:
# 1. Individual .md5 files (e.g., file.fq.gz.md5) - uses md5sum -c directly
# 2. MD5.txt files (Novogene format) - extracts the relevant line using grep
#    before validation
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
        # We cd into the directory of the fastq file because the checksum file
        # typically contains relative paths. --status suppresses normal output
        # and returns an exit code, which is ideal for scripting.
        # All output (stdout and stderr) is redirected to the log file.
        #
        # If the checksum file is MD5.txt (Novogene format), we need to extract
        # the line matching the FASTQ filename using grep, then pipe to md5sum.
        # Otherwise, use md5sum -c directly on the individual .md5 file.
        """
        (
        cd $(dirname {input.fq})
        md5_file=$(basename {input.md5})
        fq_file=$(basename {input.fq})
        
        if [ "$md5_file" = "MD5.txt" ]; then
            # Novogene format: extract the line matching the FASTQ filename
            # Use -F flag to treat filename as a literal string, not a regex pattern
            grep -F "$fq_file" "$md5_file" | md5sum -c --status -
        else
            # Individual .md5 file format
            md5sum --status -c "$md5_file"
        fi
        ) &> {log}
        """
