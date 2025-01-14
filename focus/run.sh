set -e

modalities=(alt_TSS alt_polyA expression isoforms splicing stability)

# Append all models to db file (note tissue names etc. get converted to lowercase):
# `convert.py` has been edited to handle Pantry phenotypes.
rm -rf data/geuvadis.db
for modality in ${modalities[@]}; do
    # focus import ../Geuvadis_Pheast/intermediate/twas/$modality.pos fusion --use-ens-id --tissue lcl_$modality --name geuvadis --assay rnaseq --output data/geuvadis
    focus import ../Geuvadis_Pheast/intermediate/twas/$modality.pos fusion --use-ens-id --tissue lcl --name geuvadis --assay rnaseq --output data/geuvadis
done

## View first few lines of tables in geuvadis.db:
# sqlite3 geuvadis.db
# SELECT * FROM molecularfeature LIMIT 5;
