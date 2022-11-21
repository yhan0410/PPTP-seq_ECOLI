# PPTP-seq
PPTP-seq (pooled promoter responses to transcription factor perturbations sequencing) is a highly scalable method to dissect regulatory genome. 
This repository includes data process scripts, processed data, external data, and jupyter notebooks used in the paper.

## Scripts
All the scripts were run on a HPC cluster, where SLURM workload manager was used. The scripts 1-4 were designed to process raw sequencing data for each replicate. It returns a bed file that contains read counts of each sgRNA-promoter pair in each bin. After obtaining the read counts for all triplicates, script 5 converted read counts to cell counts and fitted cell counts using log-normal distribution. 
```
$ 1_bowtie2-build.sh
$ sbatch 2_bowtie2.sbatch
$ sbatch 3_create_bed_files.sbatch
$ 4_find_closest_operon.sh
# after running 1-4 for all triplicates
$ python 5_data_process.py
```
## Jupyter Notebooks
Data processing and vistualization
