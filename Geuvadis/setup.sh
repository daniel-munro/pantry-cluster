## Use more recent bam file names (touch -h is required to update symlink timestamp)
cat input/samples.txt | while read sample; do
    echo $sample
    ln -s ${sample}.Aligned.sortedByCoord.out.bam intermediate/bam/${sample}.bam
    ln -s ${sample}.Aligned.sortedByCoord.out.bam.bai intermediate/bam/${sample}.bam.bai
    touch -t 202205310000 -h intermediate/bam/${sample}.bam
    touch -t 202205310000 -h intermediate/bam/${sample}.bam.bai
done
