mod=$1
samsize=$2

echo 'WGT	ID	CHR	P0	P1	N' > intermediate/twas/$mod.pos
cut -d'/' -f2 intermediate/twas/$mod.list \
    | sed 's/.wgt.RDat//' \
    | paste intermediate/twas/$mod.list - \
    | join -1 2 -2 4 - <(zcat input/phenotypes/$mod.bed.gz | cut -f1-4 | sort -k4) \
    | awk -v ss=$samsize '{OFS="	"; print $2, $1, $3, $4, $5, ss}' \
    >> intermediate/twas/$mod.pos
rm -f output/twas/$mod.tar.bz2
tar -cjf output/twas/$mod.tar.bz2 -C intermediate/twas $mod.list $mod.profile $mod.profile.err $mod.pos $mod/
