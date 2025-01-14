#!/bin/bash

tissue=$1
# bamdir=/mnt/mohammadi/pejdata/group/mohammadi/data/gtex/v8/bams

# ## NOTE: run_SamToFastq.py creates temporary pipe handles with generic names
# ## inside the output_dir, so samples should be run sequentially, not in parallel.
# cat phenos/$tissue/input/samples.txt | while read sample; do
#     bam=${bamdir}/${sample}.Aligned.sortedByCoord.out.patched.md.bam
#     # This method discards most reads as singletons, I think due to a sorting issue:
#     # fastq1=${tissue}/input/fastq/${sample}.1.fq
#     # fastq2=${tissue}/input/fastq/${sample}.2.fq
#     # samtools fastq -1 $fastq1 -2 $fastq2 -0 /dev/null -s /dev/null -n $bam
#     # bgzip $fastq1
#     # bgzip $fastq2
#     python3 scripts/run_SamToFastq.py $bam \
#         --prefix $sample \
#         --output_dir phenos/$tissue/input/fastq \
#         --jar ~/tools/picard.jar
# done
# ## Set timestamps to be older than bam file symlinks
# touch -t 202301010000 phenos/$tissue/input/fastq/*

## Make fastqs in temporary directory to allow parallelization
extract_fastq() {
    local tissue=$1
    local sample=$2
    local bamdir=/mnt/mohammadi/group/data/gtex/v8/bams
    local bam=${bamdir}/${sample}.Aligned.sortedByCoord.out.patched.md.bam
    echo $sample
    python3 scripts/run_SamToFastq.py $bam \
        --prefix $sample \
        --output_dir phenos/${tissue}/input/fastq/tmp_${sample} \
        --jar ~/tools/picard.jar
    mv phenos/${tissue}/input/fastq/tmp_${sample}/* phenos/${tissue}/input/fastq/
    rmdir phenos/${tissue}/input/fastq/tmp_${sample}
}
export -f extract_fastq

parallel -j 16 extract_fastq $tissue :::: phenos/${tissue}/input/samples.txt
## Set timestamps to be older than bam file symlinks
touch -t 202301010000 phenos/${tissue}/input/fastq/*

