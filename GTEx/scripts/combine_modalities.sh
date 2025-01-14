## Concatenate all phenotype tables produced by Pantry into one, along with a
## groups file that specifies all phenotypes per gene as one group. This allows
## for special downstream analyses such as QTL mapping that reports
## conditionally independent QTLs per gene even across modalities.
##
## Usage: bash scripts/combine_modalities.sh alt_polyA alt_TSS expression isoforms splicing stability
##
## Outputs: output/all.bed.gz, output/all.bed.gz.tbi, output/all.phenotype_groups.txt
##
## The phenotype names are modified by prepending '{modality}:' so they are
## unique and so their modalities can be recovered. After running, add 'all' as
## a modality in the Pheast config and generate QTLs for it.

if [ "$#" -lt 1 ]; then
    echo "Usage: bash scripts/combine_modalities.sh alt_polyA alt_TSS expression isoforms splicing stability"
    exit 1
fi

phenodir=.

## Combine all bed files:
# First confirm that all bed files have the same header:
for modality in "$@"; do
    diff <(zcat $phenodir/$1.bed.gz | head -n 1) <(zcat $phenodir/$modality.bed.gz | head -n 1)
done
rm -f $phenodir/all_tmp.bed
for modality in "$@"; do
    echo $modality
    # Prepend modality name to phenotype names:
    zcat $phenodir/$modality.bed.gz | tail -n+2 | awk -v OFS='\t' '{$4 = "'$modality':"$4; print}' >> $phenodir/all_tmp.bed
done
zcat $phenodir/$1.bed.gz | head -n1 > $phenodir/all.bed
sort -k1,1 -k2,2n $phenodir/all_tmp.bed >> $phenodir/all.bed
rm $phenodir/all_tmp.bed
bgzip $phenodir/all.bed
tabix -p bed $phenodir/all.bed.gz

## Make phenotype_groups file:
zcat $phenodir/all.bed.gz | tail -n+2 | cut -f4 | awk -F'[:.]' '{print $0"\t"$2}' > $phenodir/all.phenotype_groups.txt
