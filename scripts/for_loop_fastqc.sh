#!/bin/bash

#DIR is the directory where the blood samples FASTQ are located.
DIR="/data/courses/rnaseq_course/toxoplasma_de/reads_Blood/"

# List all FASTQ files in the directory
FILES=`ls $DIR*.fastq.gz`

# Loop over each FASTQ file and submit a separate Slurm job
for FILE in $FILES; do
    # sbatch submits a job to Slurm
    # --export=FILE=$FILE makes the variable $FILE available to the job script
    # run_fastqc_single_sample.sh is your FastQC script
    sbatch --export=FILE=$FILE "run_fastqc_single_sample.sh"
done
