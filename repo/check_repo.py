from pathlib import Path

tissues_gtex = [
    'ADPSBQ', 'ADPVSC', 'ADRNLG', 'ARTAORT', 'ARTCRN', 'ARTTBL', 'BLDDER',
    'BREAST', 'BRNACC', 'BRNAMY', 'BRNCDT', 'BRNCHA', 'BRNCHB', 'BRNCTXA',
    'BRNCTXB', 'BRNHPP', 'BRNHPT', 'BRNNCC', 'BRNPTM', 'BRNSNG', 'BRNSPC',
    'CLNSGM', 'CLNTRN', 'CVSEND', 'CVXECT', 'ESPGEJ', 'ESPMCS', 'ESPMSL',
    'FIBRBLS', 'FLLPNT', 'HRTAA', 'HRTLV', 'KDNCTX', 'KDNMDL', 'LCL', 'LIVER',
    'LUNG', 'MSCLSK', 'NERVET', 'OVARY', 'PNCREAS', 'PRSTTE', 'PTTARY',
    'SKINNS', 'SKINS', 'SLVRYG', 'SNTTRM', 'SPLEEN', 'STMACH', 'TESTIS',
    'THYROID', 'UTERUS', 'VAGINA', 'WHLBLD'
]
## Five GTEx tissues have small sample size and are generally not used for genetic analyses:
excluded = {'BLDDER', 'CVSEND', 'CVXECT', 'FLLPNT', 'KDNMDL'}
tissues_phenos = tissues_gtex + ['GEUVADIS']
tissues_pheast = [t for t in tissues_phenos if t not in excluded]

modalities = ['alt_polyA', 'alt_TSS', 'expression', 'isoforms', 'splicing', 'stability']
## Some modalities can have multiple phenotypes per gene:
modalities_mult = ['alt_polyA', 'alt_TSS', 'isoforms', 'splicing']

with open('info/traits.txt', 'r') as f:
    traits = f.read().splitlines()

files = []

## Metadata and processed results
files_info = ['tissues.phenos.txt', 'tissues.pheast.txt', 'phenotypes_per_tissue.tsv.gz', 'traits.txt', 'gwas_metadata.txt']
files += [f'info/{file}' for file in files_info]
files_processed = [
    'genes_pcg_lncrna.txt', 'geuvadis.sep.assoc.tsv.gz', 'geuvadis.sep.qtls.tsv.gz',
    'geuvadis.comb.qtls.tsv.gz', 'geuvadis.hsq.tsv.gz', 'geuvadis.twas.tsv.gz',
    'gtex.sep.assoc.tsv.gz', 'gtex.sep.qtls.tsv.gz', 'gtex.comb.qtls.tsv.gz',
    'gtex.hsq.tsv.gz', 'gtex.twas_hits.tsv.gz'
]
files += [f'processed/{file}' for file in files_processed]
files.append(f'covariates/covariates.tar.bz2')
files.append(f'covariates/covariates.plink.tar.bz2')

## Phenotype tables
for tissue in tissues_phenos:
    for modality in modalities:
        files.append(f'RNA_phenotypes/{tissue}.{modality}.unnorm.bed.gz')

for tissue in tissues_pheast:
    for modality in modalities + ['combined']:
        ## cis-QTLs
        files.append(f'QTLs/{tissue}.{modality}.cis_qtl.txt.gz')
        files.append(f'QTLs/{tissue}.{modality}.cis_independent_qtl.txt.gz')
    
    ## TWAS weights
    for modality in modalities:
        files.append(f'TWAS_weights/{tissue}.{modality}.twas_weights.tar.bz2')
        files.append(f'TWAS_weights/{tissue}.{modality}.twas_weights.profile')

## TWAS associations
for trait in traits:
    files.append(f'TWAS_associations/{trait}.tar.bz2')

files = [Path(file) for file in files]
for file in files:
    # print(file)
    assert file.exists(), file
print(f'All {len(files)} files present.')

files += [Path(file) for file in ['README.md', 'check_repo.py', 'file_tree.txt']]
for file in Path('.').glob('**/*'):
    if file.is_file() and file not in files:
        print(f'Extra file: {file}')
