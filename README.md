# PPTP-seq
PPTP-seq (pooled promoter responses to transcription factor perturbations sequencing) is a highly scalable method to dissect regulatory genome. 
This repository includes data process scripts, processed data, external data, and jupyter notebooks used in the paper.

## Setup
All the scripts were run on a HPC cluster using Linux system, where SLURM workload manager (version: slurm-wlm 17.11.7) was used. Bowtie2, Samtools, and Bedtools need to be installed in the system. Jupyer notebooks depend on Python 3.8.5 and a few commonly used modules including scipy, numpy, pandas, matplotlib, seaborn, and multiprocess.

## Raw NGS data deposition
NGS data deposited to NCBI GEO as Series GSE213624. 

## Scripts
The scripts #1-#4 were designed to process raw sequencing data for each replicate. It returns a BED file that contains read counts of each sgRNA-promoter pair in each bin. 
```
$ 1_bowtie2-build.sh
$ sbatch 2_bowtie2.sbatch
$ sbatch 3_create_bed_files.sbatch
$ 4_find_closest_operon.sh
```
After obtaining the read counts for all triplicates, the script #5 converted read counts to cell counts and fitted cell counts using log-normal distribution. 
```
$ python 5_data_process.py
```
## Jupyter Notebooks
Exploratory data analysis, differential expression analysis, TF binding site analysis, and validation using a tunable TF library.
