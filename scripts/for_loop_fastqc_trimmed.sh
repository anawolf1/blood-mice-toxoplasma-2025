#!/bin/bash

#DIR is the directory where the blood samples FASTQ are located.
DIR="/data/users/awolfruiz/rna_seq/blood-mice-toxoplasma-2025/data/external_data/fastp_output"

# List all FASTQ files in the directory
FILES=`ls $DIR/*.fastq.gz`

# Loop over each FASTQ file and submit a separate Slurm job
for FILE in $FILES; do
    # sbatch submits a job to Slurm
    # --export=FILE=$FILE makes the variable $FILE available to the job script
    # run_fastqc_single_sample_trimmed.sh is the FastQC script
    sbatch --export=FILE=$FILE "run_fastqc_single_sample_trimmed.sh"
done