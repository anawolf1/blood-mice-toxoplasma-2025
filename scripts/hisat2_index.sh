#!/bin/bash
#SBATCH --job-name=hisat2_index
#SBATCH --output=logs/output/hisat2_index.out
#SBATCH --error=logs/error/hisat2_index.err
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=08:00:00
#SBATCH --partition=pibu_el8

#check reference and annotation files are intact comparing the to the values in the CHECKSUMS file on the ftp server.
#sum Mus_musculus.GRCm39.dna.primary_assembly.fa.gz
#sum Mus_musculus.GRCm39.115.gtf.gz

#Unizip the fasta file
#gunzip Mus_musculus.GRCm39.dna.primary_assembly.fa.gz

# Paths to reference files
FASTA=/data/users/awolfruiz/rna_seq/blood-mice-toxoplasma-2025/resources/Mus_musculus.GRCm39.dna.primary_assembly.fa
GTF=/data/users/awolfruiz/rna_seq/blood-mice-toxoplasma-2025/resources/Mus_musculus.GRCm39.115.gtf

# Output folder for HISAT2 index
OUTDIR=/data/users/awolfruiz/rna_seq/blood-mice-toxoplasma-2025/data/processed_data/index_files

mkdir -p $OUTDIR
cd $OUTDIR

# Build HISAT2 index using container
apptainer exec /containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif \
  hisat2-build -p $SLURM_CPUS_PER_TASK $FASTA Mus_musculus_index