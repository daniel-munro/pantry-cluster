set -e

## Set up Pantry Pheast for all tissues
cat GTEx/tissues.pheast.txt | while read tissue; do
    echo $tissue
    mkdir -p pheast/$tissue/input

    ## Set up Pantry code
    cp ~/tools/Pantry/pheast/Snakefile pheast/$tissue/
    cp -r ~/tools/Pantry/pheast/scripts pheast/$tissue/
    cp -r ~/tools/Pantry/pheast/steps pheast/$tissue/
    cp scripts/config_pheast_gtex.yml pheast/$tissue/config.yml

    ## Prepare phenotypes
    cd phenos/$tissue
    bash scripts/combine_modalities.sh alt_polyA alt_TSS expression isoforms splicing stability
    cd -
    cd pheast/$tissue
    python3 scripts/prepare_phenotypes.py \
        --indir ../../phenos/$tissue/output/ \
        --out input/phenotypes/ \
        --map ../../../data/gtex/sample_individual_map.tsv \
        --individuals ../../GTEx/geno/ids.txt
    cd -

    ## Get samples with genotypes
    zcat pheast/$tissue/input/phenotypes/expression.bed.gz \
        | head -n1 \
        | cut -f5- \
        | sed 's/\t/\n/g' \
        > pheast/$tissue/input/samples.txt
done
