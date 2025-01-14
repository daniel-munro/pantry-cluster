tissue=$1
# Get sample size from number of lines in samples.txt:
cd pheast/$tissue
samsize=$(wc -l input/samples.txt | cut -d' ' -f1)
for mod in alt_polyA alt_TSS expression isoforms splicing stability; do
    echo 'WGT	ID	CHR	P0	P1	N' > intermediate/twas/$mod.pos
    cut -d'/' -f2 intermediate/twas/$mod.list \
        | sed 's/.wgt.RDat//' \
        | paste intermediate/twas/$mod.list - \
        | join -1 2 -2 4 - <(zcat input/phenotypes/$mod.bed.gz | cut -f1-4 | sort -k4) \
        | awk -v ss=$samsize '{OFS="	"; print $2, $1, $3, $4, $5, ss}' \
        >> intermediate/twas/$mod.pos
    rm -f output/twas/$mod.tar.bz2
    tar -cjf output/twas/$mod.tar.bz2 -C intermediate/twas $mod.list $mod.profile $mod.profile.err $mod.pos $mod/
done


# for tissue in ADPVSC ADRNLG ARTAORT ARTCRN ARTTBL BREAST BRNACC BRNAMY BRNCDT BRNCHA BRNCHB BRNCTXA BRNCTXB BRNHPP BRNHPT BRNNCC BRNPTM BRNSNG BRNSPC CLNSGM CLNTRN; do
#     cp pheast/ADPSBQ/config.yaml pheast/$tissue/
# done

## Run all tissues using alternative script that runs modalities in parallel:
# for tissue in ADRNLG ARTAORT ARTCRN ARTTBL BREAST BRNACC BRNAMY BRNCDT BRNCHA BRNCHB BRNCTXA BRNCTXB BRNHPP BRNHPT BRNNCC BRNPTM BRNSNG BRNSPC CLNSGM CLNTRN; do
#     sbatch -c16 --wrap="bash run2.sh $tissue"
# done
