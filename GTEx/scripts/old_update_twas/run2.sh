tissue=$1
cd pheast/$tissue
samsize=$(wc -l input/samples.txt | cut -d' ' -f1)
parallel -j6 bash ../../run3.sh {} $samsize ::: alt_polyA alt_TSS expression isoforms splicing stability
