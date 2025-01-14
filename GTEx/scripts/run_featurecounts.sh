# Temporary script to run featureCounts for stability phenotypes without Snakemake

TISSUE=$1
TYPE=$2
THREADS=$3

cat phenos/$TISSUE/input/samples.txt | while read sample; do
    bam=phenos/$TISSUE/intermediate/bam/$sample.bam
    out=phenos/$TISSUE/intermediate/stability/$sample.$TYPE.counts.txt
    # If type is constit_exons, extract exons from the constitutive exons GTF
    if [ "$TYPE" == "constit_exons" ]; then
        featureCounts $bam -p -a ref/constit_exons.gtf -t exon --fracOverlap 1 -T $THREADS -o $out
    elif [ "$TYPE" == "introns" ]; then
        featureCounts $bam -p -a ref/introns.gtf -t intron --fracOverlap 0 -T $THREADS -o $out
    fi
    echo "$(date) Wrote $out"    
done
echo "$(date) Done"
