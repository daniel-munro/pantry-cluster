# set -e

tissue=$1
echo $tissue

mkdir -p pheast/$tissue/input/phenotypes

## Replace sample IDs with individual IDs to match genotypes
for pheno in alt_polyA alt_TSS expression isoforms splicing stability; do
    echo $pheno
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

## Prepare latent phenotypes
echo "latent"
python3 scripts/prepare_gtex_bed.py \
    --bed phenos/$tissue/output/latent_residual.bed.gz \
    --samples geno/ids.txt \
    --out pheast/$tissue/input/phenotypes/latent_residual_tmp.bed
head -n1 pheast/$tissue/input/phenotypes/latent_residual_tmp.bed > pheast/$tissue/input/phenotypes/latent_residual.bed
grep -E ":PC[1-8]\s" pheast/$tissue/input/phenotypes/latent_residual_tmp.bed >> pheast/$tissue/input/phenotypes/latent_residual.bed
rm pheast/$tissue/input/phenotypes/latent_residual_tmp.bed
bgzip pheast/$tissue/input/phenotypes/latent_residual.bed
tabix -p bed pheast/$tissue/input/phenotypes/latent_residual.bed.gz
grep -E ":PC[1-8]\s" phenos/$tissue/output/latent_residual.phenotype_groups.txt > pheast/$tissue/input/phenotypes/latent_residual.phenotype_groups.txt

## Combine phenotype tables
cd pheast/$tissue/input/phenotypes
bash ../../../../scripts/combine_modalities.sh alt_polyA alt_TSS expression isoforms splicing stability
cd ../../../..

## Get samples with genotypes
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
