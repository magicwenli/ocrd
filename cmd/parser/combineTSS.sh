search_dir="/Source/scripts/data/deal"

cat $search_dir/*_TSS.bed | bedtools sort -i | uniq > $search_dir/all_TSS.bed
echo "tss line count: "
wc -l $search_dir/all_TSS.bed

cat $search_dir/*_TSS_low.bed | bedtools sort -i | uniq > $search_dir/all_TSS_low.bed
echo "low tss line count: "
wc -l $search_dir/all_TSS_low.bed
