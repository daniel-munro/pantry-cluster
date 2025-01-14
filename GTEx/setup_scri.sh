tissue=$1

## Set up bam files
mkdir -p phenos/$tissue/intermediate/bam
cat phenos/$tissue/input/samples.txt | while read sample; do
    ln -s /active/mohammadi_p/scripps_data/gtex/v8/bams/${sample}.Aligned.sortedByCoord.out.patched.md.bam phenos/$tissue/intermediate/bam/${sample}.bam
done

