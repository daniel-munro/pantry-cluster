set -e

# cat tissues.phenos.txt | while read tissue; do
#     cp ~/pantry/Pantry/Project/Snakefile phenos/$tissue/
#     cp ~/pantry/Pantry/Project/steps/setup.smk phenos/$tissue/steps/
#     rm -f phenos/$tissue/config.yaml
#     cp config_phenos.yml phenos/$tissue/config.yml
# done

tissue=$1
if [ "$#" -lt 1 ]; then
    echo "Provide tissue name"
    exit 1
fi

# cp ~/pantry/Pantry/Pheast/Snakefile pheast/$tissue/
# cp ~/pantry/Pantry/Pheast/steps/setup.smk pheast/$tissue/steps/
cp pheast/ADPSBQ/Snakefile pheast/$tissue/
cp pheast/ADPSBQ/steps/setup.smk pheast/$tissue/steps/
rm -f pheast/$tissue/config.yaml
cp config_pheast.yml pheast/$tissue/config.yml
# rm -rf pheast/$tissue/intermediate/twas/*stability*
cp scripts/qtl.smk pheast/$tissue/steps/

## Replace sample IDs with individual IDs to match genotypes
python3 scripts/prepare_gtex_bed.py \
    --bed phenos/$tissue/output/stability.bed.gz \
    --samples geno/ids.txt \
    --out pheast/$tissue/input/phenotypes/stability.bed

rm -f pheast/$tissue/input/phenotypes/stability.bed.gz pheast/$tissue/input/phenotypes/stability.bed.gz.tbi
bgzip pheast/$tissue/input/phenotypes/stability.bed
tabix -p bed pheast/$tissue/input/phenotypes/stability.bed.gz

## Combine phenotype tables
cd pheast/$tissue/input/phenotypes
rm -f all.bed all.bed.gz all.bed.gz.tbi all.phenotype_groups.txt
bash ../../../../scripts/combine_modalities.sh alt_polyA alt_TSS expression isoforms splicing stability
cd ../../../..

echo "done"
