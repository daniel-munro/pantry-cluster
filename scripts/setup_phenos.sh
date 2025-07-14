set -e

## Set up Pantry phenotyping for all tissues
cut -f2 GTEx/samples_tissues.txt | sort | uniq | while read tissue; do
    echo $tissue
    ## Set up fastq files
    mkdir -p phenos/$tissue/input/fastq
    awk -v tissue=$tissue '$2 == tissue' GTEx/samples_tissues.txt | cut -f1 > phenos/$tissue/input/samples.txt
    awk '{ print $1"_1.fastq.gz\t"$1"_2.fastq.gz\t"$1 }' phenos/$tissue/input/samples.txt > phenos/$tissue/input/fastq_map.txt
    ## Set up bam files
    mkdir -p phenos/$tissue/intermediate/bam
    cat phenos/$tissue/input/samples.txt | while read sample; do
        ln -s /data/hps/assoc/private/gdml/from_rss/scripps_data/gtex/v8/bams/${sample}.Aligned.sortedByCoord.out.patched.md.bam phenos/$tissue/intermediate/bam/${sample}.bam
    done
    ## Set up Pantry code
    cp ~/tools/Pantry/phenotyping/Snakefile phenos/$tissue/
    cp -r ~/tools/Pantry/phenotyping/scripts phenos/$tissue/
    cp -r ~/tools/Pantry/phenotyping/steps phenos/$tissue/
    cp scripts/config_phenos_gtex.yml phenos/$tissue/config.yml
done

## Get fastq files for one tissue
# sbatch --account=cpu-gdml-investor --partition=cpu-gdml-investor -c16 --mem=32000 --time=48:00:00 --wrap="bash scripts/gtex_bam_to_fastq.sh tissuename" -J tissuename
