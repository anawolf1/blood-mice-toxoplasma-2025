#!/bin/bash
#------------------------------
# Slurm job configuration
#------------------------------
#SBATCH --job-name=fastqc
#SBATCH --output=/data/users/awolfruiz/rna_seq/blood-mice-toxoplasma-2025/logs/output/fastqc/fastqc_%A_%a.out
#SBATCH --error=/data/users/awolfruiz/rna_seq/blood-mice-toxoplasma-2025/logs/error/fastqc/fastqc_%A_%a.err 
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:20:00
#SBATCH --partition=pibu_el8

# DIR is the directory where the blood samples FASTQ are located.
DIR="/data/users/awolfruiz/rna_seq/blood-mice-toxoplasma-2025/data/external_data/fastp_output/"

# OUTPUT_DIR is where FastQC will save its reports for this job
OUTPUT_DIR="/data/users/awolfruiz/rna_seq/blood-mice-toxoplasma-2025/data/external_data/fastqc_trimmed_output/"

# Load FASTQC module
module load FastQC/0.11.9-Java-11

# Run FASTQC
#------------------------------
# --threads $SLURM_CPUS_PER_TASK: use the CPU cores allocated in the slurm job configuration
# -o $OUTPUT_DIR: save output files to the output directory
# "$FILE": the FASTQ file to analyze

echo "DB: $FILE"

fastqc --threads $SLURM_CPUS_PER_TASK -o "$OUTPUT_DIR" "$FILE"