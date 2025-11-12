#!/bin/bash
#SBATCH --job-name=multiqc
#SBATCH --output=logs/output/multiqc_%j.out
#SBATCH --error=logs/error/multiqc_%j.err
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:10:00
#SBATCH --partition=pibu_el8

# Directory with all FastQC outputs
INPUT_DIR=/data/users/awolfruiz/rna_seq/blood-mice-toxoplasma-2025/data/external_data/fastqc_output
OUTPUT_DIR=/data/users/awolfruiz/rna_seq/blood-mice-toxoplasma-2025/data/external_data/multiqc


# Run MultiQC
apptainer exec /containers/path-to/multiqc-1.19.sif multiqc $INPUT_DIR -o $OUTPUT_DIR

