tissue=$1
echo $tissue

mkdir -p pheast/$tissue/input/phenotypes

## Replace sample IDs with individual IDs to match genotypes
for pheno in alt_polyA alt_TSS expression isoforms splicing stability; do
    echo $pheno
    # zcat phenos/$tissue/output/${pheno}.bed.gz \
    #     | head -n1 \
    #     | sed -E 's/(GTEX-[A-Za-z0-9]+)-\S+/\1/g' \
    #     > pheast/$tissue/input/phenotypes/${pheno}.bed
    # zcat phenos/$tissue/output/${pheno}.bed.gz \
    #     | tail -n+2 \
    #     >> pheast/$tissue/input/phenotypes/${pheno}.bed
    # ## Confirm these are the same, i.e. no duplicate individuals
    # head -n1 pheast/$tissue/input/phenotypes/${pheno}.bed \
    #     | cut -f5- \
    #     | sed 's/\t/\n/g' \
    #     | wc -l
    # head -n1 pheast/$tissue/input/phenotypes/${pheno}.bed \
    #     | cut -f5- \
    #     | sed 's/\t/\n/g' \
    #     | sort \
    #     | uniq \
    #     | wc -l
        # --bed pheast/$tissue/input/phenotypes/${pheno}.bed \
    python3 scripts/prepare_gtex_bed.py \
        --bed phenos/$tissue/output/${pheno}.bed.gz \
        --samples geno/ids.txt \
        --out pheast/$tissue/input/phenotypes/${pheno}.bed

    bgzip pheast/$tissue/input/phenotypes/${pheno}.bed
    tabix -p bed pheast/$tissue/input/phenotypes/${pheno}.bed.gz
done
for pheno in alt_polyA alt_TSS isoforms splicing; do
    cp phenos/$tissue/output/${pheno}.phenotype_groups.txt pheast/$tissue/input/phenotypes/
done

## Get samples with genotypes
    # | grep -f geno/ids.txt \
zcat pheast/$tissue/input/phenotypes/expression.bed.gz \
    | head -n1 \
    | cut -f5- \
    | sed 's/\t/\n/g' \
    > pheast/$tissue/input/samples.txt

## Set up Pheast code
cp ~/pantry/Pantry/Pheast/Snakefile pheast/$tissue/
cp -r ~/pantry/Pantry/Pheast/scripts pheast/$tissue/
cp ~/tools/gemma pheast/$tissue/scripts/fusion_twas/
cp -r ~/pantry/Pantry/Pheast/steps pheast/$tissue/
## modified execution parameters:
cp scripts/qtl.smk pheast/$tissue/steps/
cp scripts/twas.smk pheast/$tissue/steps/
cp config_pheast.yaml pheast/$tissue/config.yaml

echo "done"
