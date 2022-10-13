# PPTP-seq
PPTP-seq (pooled promoter responses to transcription factor perturbations sequencing) is a highly scalable method to dissect regulatory genome. 
This repository includes data process scripts, processed data, external data, and jupyter notebooks used in the paper.

## Scripts
All the scripts were run on a HPC cluster, where SLURM workload manager was used. The scripts were designed to process raw sequencing data for each replicate. It returns a bed file that contains read counts of each sgRNA-promoter pair in each bin.
```
$ 1_bowtie2-build.sh
$ sbatch 2_bowtie2.sbatch
$ sbatch 3_create_bed_files.sbatch
$ 4_find_closest_operon.sh
```
## Jupyter Notebooks
Data processing and vistualization
