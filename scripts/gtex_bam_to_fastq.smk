"""Extract fastq files from bams for a given tissue.

Run with, e.g.:

snakemake -s scripts/gtex_bam_to_fastq.smk --profile slurm --config tissue=ADPSBQ -j200 --retries 2
"""

import os

# Configurable tissue variable (set via --config tissue=XXX)
tissue = config["tissue"]

# Read sample names from samples.txt
with open(f"phenos/{tissue}/input/samples.txt") as f:
    SAMPLES = [line.strip() for line in f if line.strip()]

BAMDIR = "/data/hps/assoc/private/gdml/from_rss/scripps_data/gtex/v8/bams"
FASTQDIR = f"phenos/{tissue}/input/fastq"

rule all:
    input:
        touch_done=f"{FASTQDIR}/touch.done",

rule bam_to_fastq:
    input:
        bam = lambda wildcards: f"{BAMDIR}/{wildcards.sample}.Aligned.sortedByCoord.out.patched.md.bam"
    output:
        fastq1 = f"{FASTQDIR}/{{sample}}_1.fastq.gz",
        fastq2 = f"{FASTQDIR}/{{sample}}_2.fastq.gz"
    params:
        tmpdir = f"{FASTQDIR}/tmp_{{sample}}",
        jar = os.path.expanduser("~/tools/picard.jar")
    resources:
        mem_mb = lambda wildards, attempt: 16000 * attempt,
    shell:
        r"""
        rm -rf {params.tmpdir}
        mkdir {params.tmpdir}
        python3 GTEx/scripts/run_SamToFastq.py {input.bam} \
            --prefix {wildcards.sample} \
            --output_dir {params.tmpdir} \
            --jar {params.jar}
        mv {params.tmpdir}/* {FASTQDIR}/
        rmdir {params.tmpdir}
        """

rule touch_fastqs:
    input:
        fastqs=expand(f"{FASTQDIR}/{{sample}}_{{pair}}.fastq.gz", sample=SAMPLES, pair=[1,2])
    output:
        touch_done=f"{FASTQDIR}/touch.done"
    shell:
        r"""
        touch -t 202501010000 {FASTQDIR}/*.fastq.gz
        touch {output.touch_done}
        """ 
