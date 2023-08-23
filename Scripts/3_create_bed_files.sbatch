#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --array=1-16%16
#SBATCH --mem=45G

module load samtools
module load bedtools

ID=${SLURM_ARRAY_TASK_ID}
SAM_FILE_PATH="../SAM_files/Glucose_rep_1"

# extract read id and mapped sgRNA columns
awk '{if($2 == 0)print $1,$3}' $SAM_FILE_PATH/dual_lib_bin${ID}_r1.sam >$SAM_FILE_PATH/dual_lib_bin${ID}_r1.tab
# extract the coordinates of read2 (promoter side) and generate a bed format file
samtools view -bS -F 4 $SAM_FILE_PATH/dual_lib_bin${ID}_r2.sam | bedtools bamtobed -i stdin >$SAM_FILE_PATH/dual_lib_bin${ID}_r2.bed
# replace the read id with the corresponding sgRNA, skip the reads that can not mapped to defined sgRNAs
awk 'NR==FNR {h[$1] = $2; next} {if(h[$4]!="") print $1, $2, $3, h[$4], $5, $6}' $SAM_FILE_PATH/dual_lib_bin${ID}_r1.tab $SAM_FILE_PATH/dual_lib_bin${ID}_r2.bed >$SAM_FILE_PATH/bin${ID}_pair.bed
# count unique pairs of read 2 coordinates and sgRNA
sort $SAM_FILE_PATH/bin${ID}_pair.bed | uniq -c | awk '{OFS="\t"}{print $2,$3,$4,$3":"$4":"$7,$1,$7,$5}' >$SAM_FILE_PATH/bin${ID}_pair_count.bed
awk '{OFS="\t"}{print $0,ARGV[2]}' $SAM_FILE_PATH/bin${ID}_pair_count.bed ${ID} >$SAM_FILE_PATH/bin${ID}_result.bed
