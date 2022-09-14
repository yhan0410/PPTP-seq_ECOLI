cd ~/work_dir/MPRA_TRN/SAM_files/20210823
# merge all read information from different bins
# promoter read is on P7 end, so if the read is mapped to '-' strand, it corresponds to the promoter on the '+' strand
cat *_result.bed | sort -k2,2n -k4 -k7 -k8,2n | awk '{if($6=="-") print >"for_promoter_read.bed"; else print >"rev_promoter_read.bed"}'
# find the corresponding operon for each promoter read using bedtools closest and 
bedtools closest -D a -S -id -t last -a for_promoter_read.bed  -b ../../External_data/gene_annotation/all_operon.bed >closest_operon_map_for.bed
bedtools closest -D a -S -id -t first -a rev_promoter_read.bed  -b ../../External_data/gene_annotation/all_operon.bed >closest_operon_map_rev.bed
# filter out promoter reads within or outside the operon
# according to Zaslaver et al., 2006, the promoter region include intergenic region plus about 50-150 bp into each flanking coding region
# I found that it is not stricted to 50-150 bp, so I extend this range by extra 50 bp
# the end of a promoter which was sequenced here should be in the first 0-200 bp in the coding region.
awk '{if(($3 >= $10) && ($3 <= $10 + 200)) print $0}' closest_operon_map_for.bed >filtered_closest_operon_map_for.txt
awk '{if(($2 <= $11) && ($2 >= $11 - 200)) print $0}' closest_operon_map_rev.bed >filtered_closest_operon_map_rev.txt
cat filtered_closest_operon_map_for.txt filtered_closest_operon_map_rev.txt >filtered_closest_operon_map_all.txt
