# PPTP-seq
PPTP-seq (pooled promoter responses to transcription factor perturbations sequencing) is a highly scalable method to dissect regulatory genome. 
This repository includes data process scripts, processed data, external data, and jupyter notebooks used in the paper.

## Citation:
Han, Y., Li, W., Filko, A. et al. Genome-wide promoter responses to CRISPR perturbations of regulators reveal regulatory networks in Escherichia coli. Nat Commun 14, 5757 (2023). https://doi.org/10.1038/s41467-023-41572-4

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
$ python 5_calculate_read_counts.py
$ python 6_data_process.py
```
## Jupyter Notebooks
Four files are used for exploratory data analysis, differential expression analysis, TF binding site analysis, and validation using a tunable TF library. Figures were reproduced in these notebooks.

## Processed data
Columns in PPTP-seq_Glu.csv, PPTP-seq_LB.csv, PPTP-seq_Gly.csv:
1. operon: operon expressed by the promoter of interest
1. tf_gene: tf gene targeted by the sgRNA
1. mean: average promoter activity across replicates at natural log scale
1. std: standard deviation of the promoter activity across replicates at natural log scale
1. std_linear: standard deviation of the promoter activity across replicates at linear scale
1. n_rep: number of replicates measured for the variant
1. FC: log2 fold change compared to the promoter activity in negative controls
1. FCOM: log2 fold change compared to the median promoter activity across all CRISPRi perturbations
1. -logP: -log10(adjusted p value)
1. class: 1 represents up-regulated, 0 represents not significant, and -1 represents down-regulated
1. FC_global: average FC of all promoter responses to a TF perturbation
1. FC_specific: FC - FC_global, the gene specific effects that exclude global effects from a TF
