#################################
## GWAS loci for TWAS analysis ##
#################################

# Extract genome-wide significant variants from 144 trait GWAS sumstats for GWAS locus-based analyses of TWAS
echo -e "trait\tvariant_id\tpvalue" > gwas_genome_wide_sig.tsv
cat ../Geuvadis_TWAS/data/gwas/gwas_metadata.txt | tail -n+2 | cut -f2 | while read trait; do
    echo $trait
    zcat ../Geuvadis_TWAS/data/gwas/imputed_gwas_hg38_1.1/imputed_${trait}.txt.gz \
        | tail -n+2 \
        | cut -f2,11 \
        | awk '$2 < 5e-8' \
        | sed "s/^/$trait\t/" \
        >> gwas_genome_wide_sig.tsv
done
gzip gwas_genome_wide_sig.tsv

##############################################
## GWAS in cis-window for TWAS example plot ##
##############################################

# TRAIT=GLGC_Mc_LDL
# GENE=ENSG00000142444
# CHROM=chr19
# TSS=10928811

# TRAIT=UKB_50_Standing_height
# GENE=ENSG00000174945
# CHROM=chr7
# TSS=2679522

# TRAIT=Astle_et_al_2016_Eosinophil_counts
# GENE=ENSG00000197536
# CHROM=chr5
# TSS=132410636

TRAIT=UKB_20127_Neuroticism_score
GENE=ENSG00000115947
CHROM=chr2
TSS=148021604

zcat ../Geuvadis_TWAS/data/gwas/imputed_gwas_hg38_1.1/imputed_$TRAIT.txt.gz \
    | awk -v chrom=$CHROM -v tss=$TSS 'NR==1 || ($3 == chrom && $4 >= tss - 1000000 && $4 <= tss + 1000000)' \
    | gzip -c \
    > gwas.$TRAIT.$GENE.$CHROM.$TSS.tsv.gz
