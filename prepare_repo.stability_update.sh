set -e

MODALITIES="alt_polyA alt_TSS expression isoforms splicing stability"

echo "=== GTEx Phenotypes ==="
cat GTEx/tissues.phenos.txt | while read TISSUE; do
    echo $TISSUE
    for MODALITY in stability; do
        rsync -a GTEx/phenos/$TISSUE/intermediate/unnorm/$MODALITY.bed.gz repo/RNA_phenotypes/$TISSUE.$MODALITY.unnorm.bed.gz
    done
done

echo "=== Covariates ==="
## Copy all covariate files to make new compressed archives
mv -i repo/covariates repo/covariates_old
mkdir repo/covariates
for MODALITY in $MODALITIES; do
    rsync -a Geuvadis_Pheast/intermediate/covar/$MODALITY.covar.tsv repo/covariates/GEUVADIS.$MODALITY.covar.tsv
    rsync -a Geuvadis_Pheast/intermediate/covar/$MODALITY.covar.plink.tsv repo/covariates/GEUVADIS.$MODALITY.covar.plink.tsv
done
rsync -a Geuvadis_Pheast/intermediate/covar/all.covar.tsv repo/covariates/GEUVADIS.combined.covar.tsv

cat GTEx/tissues.pheast.txt | while read TISSUE; do
    echo $TISSUE
    for MODALITY in $MODALITIES; do
        rsync -a GTEx/pheast/$TISSUE/intermediate/covar/$MODALITY.covar.tsv repo/covariates/$TISSUE.$MODALITY.covar.tsv
        rsync -a GTEx/pheast/$TISSUE/intermediate/covar/$MODALITY.covar.plink.tsv repo/covariates/$TISSUE.$MODALITY.covar.plink.tsv
    done
    rsync -a GTEx/pheast/$TISSUE/intermediate/covar/all.covar.tsv repo/covariates/$TISSUE.combined.covar.tsv
done

## Combine covariate files into compressed archives for easier uploading/downloading
cd repo/covariates
tar -cjf covariates.tar.bz2 *.covar.tsv
tar -cjf covariates.plink.tar.bz2 *.covar.plink.tsv
rm *.covar.tsv
rm *.covar.plink.tsv
cd ../..

echo "=== GTEx QTLs and TWAS ==="
cat GTEx/tissues.pheast.txt | while read TISSUE; do
    echo $TISSUE
    for MODALITY in stability; do
        rsync -a GTEx/pheast/$TISSUE/output/qtl/$MODALITY.cis_qtl.txt.gz repo/QTLs/$TISSUE.$MODALITY.cis_qtl.txt.gz
        rsync -a GTEx/pheast/$TISSUE/output/qtl/$MODALITY.cis_independent_qtl.txt.gz repo/QTLs/$TISSUE.$MODALITY.cis_independent_qtl.txt.gz
        rsync -a GTEx/pheast/$TISSUE/output/twas/$MODALITY.tar.bz2 repo/TWAS_weights/$TISSUE.$MODALITY.twas_weights.tar.bz2
        rsync -a GTEx/pheast/$TISSUE/intermediate/twas/$MODALITY.profile repo/TWAS_weights/$TISSUE.$MODALITY.twas_weights.profile
    done
    rsync -a GTEx/pheast/$TISSUE/output/qtl/all.cis_qtl.txt.gz repo/QTLs/$TISSUE.combined.cis_qtl.txt.gz
    rsync -a GTEx/pheast/$TISSUE/output/qtl/all.cis_independent_qtl.txt.gz repo/QTLs/$TISSUE.combined.cis_independent_qtl.txt.gz
done

## Concatenate combined-mapping phenotype group files to get full list of phenotypes
printf "tissue\tmodality\tgene\tphenotype\n" > repo/info/phenotypes_per_tissue.tsv
sed "s/^/GEUVADIS\t/" Geuvadis/output/all.phenotype_groups.txt \
    | sed "s/:/\t/" \
    | awk '{print $1"\t"$2"\t"$4"\t"$3}' \
    | sort \
    >> repo/info/phenotypes_per_tissue.tsv
cat GTEx/tissues.pheast.txt | while read TISSUE; do
    sed "s/^/$TISSUE\t/" GTEx/pheast/$TISSUE/input/phenotypes/all.phenotype_groups.txt \
        | sed "s/:/\t/" \
        | awk '{print $1"\t"$2"\t"$4"\t"$3}' \
        | sort \
        >> repo/info/phenotypes_per_tissue.tsv
done
gzip repo/info/phenotypes_per_tissue.tsv

echo "=== TWAS associations ==="
mv -i repo/TWAS_associations repo/TWAS_associations_old
mkdir repo/TWAS_associations
cat repo/info/traits.txt | while read TRAIT; do
    echo $TRAIT
    mkdir $TRAIT
    for MODALITY in $MODALITIES; do
        ln -s ../Geuvadis_TWAS/data/output/$MODALITY/fusion.Geuvadis.$MODALITY.$TRAIT.tsv $TRAIT/
        cat GTEx/tissues.pheast.txt | while read TISSUE; do
            ln -s ../GTEx/twas/$TISSUE/output/$MODALITY/$MODALITY.$TRAIT.tsv $TRAIT/fusion.$TISSUE.$MODALITY.$TRAIT.tsv
        done
    done
    ## -h dereferences symlinks to include the actual files
    tar -cjhf repo/TWAS_associations/$TRAIT.tar.bz2 $TRAIT/
    rm -r $TRAIT
done

cd repo && tree > file_tree.txt && cd ..

echo "=== Done ==="