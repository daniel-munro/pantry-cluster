# A tensorqtl update changed the output format to replace tss_distance (column 8)
# to start_distance and end_distance. Run this to convert to the old format
# (keeping a copy of the old file as {fname}.new_format.txt.gz).

set -e

FILE=$1

mv -i $FILE "$FILE.new_format.txt.gz"
zcat "$FILE.new_format.txt.gz" \
    | cut -f-8,10- \
    | awk 'NR == 1 {gsub("start_distance", "tss_distance", $0)} {print}' \
    | gzip -c \
    > $FILE
