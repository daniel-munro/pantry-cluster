vcf_phased=/mnt/mohammadi/group/data/gtex/v8/ase/variant_calls/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.vcf
gtf=../ref/Homo_sapiens.GRCh38.106.chr.gtf

## VCF must be gzipped for aFC.py
sbatch -c 16 --wrap="bgzip -c -@15 $vcf_phased > geno_phased.vcf.gz"
sbatch --wrap="tabix -p vcf geno_phased.vcf.gz"

## Gene annotation BED for phASER Gene AE
awk '$3 == "gene"' $gtf \
    | cut -f1,4,5,9 \
    | awk -F'\t' '{match($4, /gene_id "([^"]+)";/, arr); print $1"\t"$2"\t"$3+1"\t"arr[1]}' \
    > genes.bed
