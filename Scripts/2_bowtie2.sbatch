#!/bin/bash
#SBATCH --cpus-per-task=4
#SBATCH --array=1-16%16
#SBATCH --mem=50G

module load bowtie2

ID=${SLURM_ARRAY_TASK_ID}
BARCODE=$( sed -n ${SLURM_ARRAY_TASK_ID}p lookup.txt )
# the pre-fix of raw data file and file path need to be changed for each replicate
RAW_DATA_PATH="../Raw_NGS_data/20210823"
SAM_FILE_PATH="../SAM_files/20210823"
OUT_FILE_PATH="../OUT/20210823"

bowtie2 --very-sensitive-local -p 4 -5 20 -3 100 -x ../Database/IDX_files/sgRNA_library.fa -U ../Raw_NGS_data/20210823/H72VTDSX2_${BARCODE}_L003_R1.fastq.gz -S ../SAM_files/20210823/dual_lib_bin${ID}_r1.sam 2>../OUT/20210823/r1_align_log_bin${ID}.out
bowtie2 --very-sensitive-local -p 4 -5 21 -x ../Database/IDX_files/ECOLI_genome.fa -U ../Raw_NGS_data/20210823/H72VTDSX2_${BARCODE}_L003_R2.fastq.gz -S ../SAM_files/20210823/dual_lib_bin${ID}_r2.sam 2>../OUT/20210823/r2_align_log_bin${ID}.out