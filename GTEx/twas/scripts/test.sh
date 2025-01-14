Rscript ../scripts/FUSION.assoc_test.R \
    --sumstats ../data/sumstats/UKB_50_Standing_height.sumstats \
    --weights weights/expression.pos \
    --weights_dir weights \
    --ref_ld_chr ../data/LDREF_b38ids_chr/1000G.EUR. \
    --chr chr16 \
    --coloc_P 5e-8 \
    --GWASN 337119 \
    --out test.tsv

grep -v 'NA$' test.tsv | less
