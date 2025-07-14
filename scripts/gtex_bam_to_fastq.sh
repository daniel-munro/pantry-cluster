#!/bin/bash

tissue=$1

## Make fastqs in temporary directory to allow parallelization
extract_fastq() {
    local tissue=$1
    local sample=$2
    local bamdir=/data/hps/assoc/private/gdml/from_rss/scripps_data/gtex/v8/bams
    local bam=${bamdir}/${sample}.Aligned.sortedByCoord.out.patched.md.bam
    echo $sample
    python3 GTEx/scripts/run_SamToFastq.py $bam \
        --prefix $sample \
        --output_dir phenos/${tissue}/input/fastq/tmp_${sample} \
        --jar ~/tools/picard.jar
    mv phenos/${tissue}/input/fastq/tmp_${sample}/* phenos/${tissue}/input/fastq/
    rmdir phenos/${tissue}/input/fastq/tmp_${sample}
}
export -f extract_fastq

parallel -j 16 extract_fastq $tissue :::: phenos/${tissue}/input/samples.txt
## Set timestamps to be older than bam file symlinks
touch -t 202501010000 phenos/${tissue}/input/fastq/*

