module load bowtie2
cd ~/work_dir/MPRA_TRN/Database
bowtie2-build ./FASTA_files/genome_MG1655.fasta ./IDX_files/ECOLI_genome.fa
bowtie2-build ./FASTA_files/sgRNA.fasta ./IDX_files/sgRNA_library.fa