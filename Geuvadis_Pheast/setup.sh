# ## Edit copy/pasted code below:
# ## Prepare latent phenotypes
# echo "latent"
# python3 scripts/prepare_gtex_bed.py \
#     --bed phenos/$tissue/output/latent_residual.bed.gz \
#     --samples geno/ids.txt \
#     --out pheast/$tissue/input/phenotypes/latent_residual_tmp.bed
# head -n1 pheast/$tissue/input/phenotypes/latent_residual_tmp.bed > pheast/$tissue/input/phenotypes/latent_residual.bed
# grep -E ":PC[1-8]\s" pheast/$tissue/input/phenotypes/latent_residual_tmp.bed >> pheast/$tissue/input/phenotypes/latent_residual.bed
# rm pheast/$tissue/input/phenotypes/latent_residual_tmp.bed
# bgzip pheast/$tissue/input/phenotypes/latent_residual.bed
# tabix -p bed pheast/$tissue/input/phenotypes/latent_residual.bed.gz
# grep -E ":PC[1-8]\s" phenos/$tissue/output/latent_residual.phenotype_groups.txt > pheast/$tissue/input/phenotypes/latent_residual.phenotype_groups.txt

## Combine phenotype tables
cd ../Geuvadis
bash scripts/combine_modalities.sh alt_polyA alt_TSS expression isoforms splicing stability
cd ../Geuvadis_Pheast

## Prepare genotypes (autosomes only)
plink2 --bfile ../Geuvadis/input/genotype/GEUVADIS.445_samples.GRCh38.20170504.maf01.filtered.nodup \
    --make-bed \
    --max-alleles 2 \
    --chr 1-22 \
    --out ../Geuvadis/input/genotype/GEUVADIS.445_samples.GRCh38.20170504.maf01.filtered.nodup.nochr.chr1-22
