## Convert chrom names in Ensembl GTF to match those in BAM files (using Ensembl for compatibility with txrevise)
cat /gpfs/home/dmunro/pantry/Geuvadis/input/human_ref/Homo_sapiens.GRCh38.106.gtf \
    | grep -vP '^(G|K)' \
    | sed 's/^MT/M/' \
    | sed '/^[^#]/s/^/chr/' \
    > ref/Homo_sapiens.GRCh38.106.chr.gtf

## Get available samples per tissue
ls /mnt/mohammadi/group/data/gtex/v8/bams/ \
    | sed 's/.Aligned.sortedByCoord.out.patched.md.bam//g' \
    > samples_with_bam.txt

## Get genotypes and filter to autosomes and MAF > 0.01:
sbatch --wrap="plink2 --make-bed \
    --vcf /mnt/mohammadi/group/data/gtex/v8/ase/variant_calls/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze_MAF001_GTonly.vcf \
    --maf 0.01 \
    --max-alleles 2 \
    --chr 1-22 \
    --output-chr chrM \
    --out geno/gtex"
cut -f2 geno/gtex.fam > geno/ids.txt

## Get tissue name for each sample, then convert to tissue abbreviation
curl 'https://gtexportal.org/rest/v1/dataset/tissueInfo?format=tsv&datasetId=gtex_v8' > tissueInfo.tsv
awk 'BEGIN {FS=OFS="\t"} NR==FNR {h[$1] = $7; next} {print $1, h[$1]}' \
    /mnt/mohammadi/group/data/gtex/v8/annotations/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt \
    samples_with_bam.txt \
    | awk 'BEGIN {FS=OFS="\t"} NR==FNR {h[$3] = $4; next} {print $1, h[$2]}' tissueInfo.tsv - \
    > samples_tissues.txt
rm samples_with_bam.txt

## Set up Pantry for one or more tissues
cut -f2 samples_tissues.txt | sort | uniq | while read tissue; do
    echo $tissue
    ## Set up fastq files
    mkdir -p phenos/$tissue/input/fastq
    awk -v tissue=$tissue '$2 == tissue' samples_tissues.txt | cut -f1 > phenos/$tissue/input/samples.txt
    awk '{ print $1"_1.fastq.gz\t"$1"_2.fastq.gz\t"$1 }' phenos/$tissue/input/samples.txt > phenos/$tissue/input/fastq_map.txt
    ## sbatch --wrap="bash scripts/bam_to_fastq.sh $tissue"
    ## Set up bam files
    mkdir -p phenos/$tissue/intermediate/bam
    cat phenos/$tissue/input/samples.txt | while read sample; do
        ln -s /mnt/mohammadi/group/data/gtex/v8/bams/${sample}.Aligned.sortedByCoord.out.patched.md.bam phenos/$tissue/intermediate/bam/${sample}.bam
    done
    ## Set up Pantry code
    cp ~/pantry/Pantry/Project/Snakefile phenos/$tissue/
    cp -r ~/pantry/Pantry/Project/scripts phenos/$tissue/
    ln -s ~/latent-rna phenos/$tissue/
    cp -r ~/pantry/Pantry/Project/steps phenos/$tissue/
    cp scripts/latent.smk phenos/$tissue/steps/
    cp config_phenos.yaml phenos/$tissue/config.yaml
done

## Get fastq files for one tissue
# sbatch --cpus-per-task=16 --time="2-0" --wrap="bash scripts/bam_to_fastq.sh tissuename" -J tissuename

## Set up Pheast for one tissue (must be done after running phenotyping)
# sbatch --wrap="bash scripts/setup_pheast.sh tissuename"

##########
## TWAS ##
##########

cp -i ../Geuvadis_TWAS/data/gwas/gwas_metadata.txt twas/data/
ln -s ../../../Geuvadis_TWAS/data/intermediate/sumstats twas/data/sumstats

## Add 'chr' to chromosome names in LDREF files
mkdir -p twas/data/LDREF_b38ids_chr
for chrom in {1..22}; do
    echo $chrom
    sed 's/^/chr/' ../Geuvadis_TWAS/data/LDREF_b38ids/1000G.EUR.${chrom}.bim > twas/data/LDREF_b38ids_chr/1000G.EUR.chr${chrom}.bim
    cp -i ../Geuvadis_TWAS/data/LDREF_b38ids/1000G.EUR.${chrom}.bed twas/data/LDREF_b38ids_chr/1000G.EUR.chr${chrom}.bed
    cp -i ../Geuvadis_TWAS/data/LDREF_b38ids/1000G.EUR.${chrom}.fam twas/data/LDREF_b38ids_chr/1000G.EUR.chr${chrom}.fam
done

## Set up TWAS for one or more tissues
cut -f2 samples_tissues.txt | sort | uniq | while read tissue; do
    echo $tissue
    mkdir -p twas/$tissue
    ln -s ../../pheast/$tissue/intermediate/twas twas/$tissue/weights
    cp twas/scripts/Snakefile twas/$tissue/
done
