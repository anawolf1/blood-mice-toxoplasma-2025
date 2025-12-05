#!/bin/bash
#------------------------------
# Slurm job configuration
#------------------------------
#SBATCH --job-name=hisat2_index
#SBATCH --output=/data/users/awolfruiz/rna_seq/blood-mice-toxoplasma-2025/logs/output/hisat2/hisat2_index.out
#SBATCH --error=/data/users/awolfruiz/rna_seq/blood-mice-toxoplasma-2025/logs/error/hisat2/hisat2_index.err
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --partition=pibu_el8

# REFERENCE is the reference genome fasta file and INDEX_BASENAME is the prefix for the HISAT2 index files created from the reference genome
REFERENCE="/data/users/awolfruiz/rna_seq/blood-mice-toxoplasma-2025/resources/Mus_musculus.GRCm39.dna.primary_assembly.fa"
INDEX_BASENAME="/data/users/awolfruiz/rna_seq/blood-mice-toxoplasma-2025/results/index_files/genome_index"

# Build the HISAT2 index from the reference genome FASTA file
# This creates multiple index files with the prefix specified in $INDEX_BASENAME
# These index files are required for mapping reads with HISAT2
apptainer exec --bind /data/ /containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif hisat2-build $REFERENCE $INDEX_BASENAME