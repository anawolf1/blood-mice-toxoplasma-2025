#!/bin/bash
#SBATCH --array=1-12
#SBATCH --job-name= slurm array
#SBATCH --output=array_%J.err
#SBATCH --error=array_%J.err
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=00:10:00
#SBATCH --partition=pibu_el8

# Load FastQC
module load FastQC/0.11.9-Java-11

# Array of all FASTQ files
FILES=(/data/courses/rnaseq_course/toxoplasma_de/reads_Blood/*.fastq.gz)

# Pick the file corresponding to this array task
FILE=${FILES[$SLURM_ARRAY_TASK_ID]}

# Output directory
OUTPUT_DIR=/data/users/awolfruiz/rna_seq/blood-mice-toxoplasma-2025/data/external_data/fastqc_output
mkdir -p $OUTPUT_DIR

# Run FastQC
fastqc --threads $SLURM_CPUS_PER_TASK -o $OUTPUT_DIR "$FILE"


#Goes through all samples and maps the forward and reverse strands in parallel (takes half of the time)
#Each sample has a forward and revers fastqc sample.