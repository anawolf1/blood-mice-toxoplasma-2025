#!/bin/bash
#------------------------------
# Slurm job configuration
#------------------------------
#SBATCH --job-name=multiqc
#SBATCH --output=/data/users/awolfruiz/rna_seq/blood-mice-toxoplasma-2025/logs/output/multiqc/multiqc_%j.out
#SBATCH --error=/data/users/awolfruiz/rna_seq/blood-mice-toxoplasma-2025/logs/error/multiqc/multiqc_%j.err
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=00:20:00
#SBATCH --partition=pibu_el8

# Load the required module
module load MultiQC/1.11-foss-2021a

WORKDIR="/data/users/awolfruiz/rna_seq/blood-mice-toxoplasma-2025"

# Directory containing the FeatureCounts summary
SUMMARY_FILE="$WORKDIR/results/featureCounts/gene_counts.txt.summary"

# Directory where MultiQC will save the output
OUTPUT_DIR_MULTIQC="$WORKDIR/results/multiqc"

# Create the output directory (if it doesn't exist already)
mkdir -p "/data/users/awolfruiz/rna_seq/blood-mice-toxoplasma-2025/results/multiqc"

# Run multiqc
#------------------------------
# --outdir indicates where the output will be saved
multiqc --outdir "$OUTPUT_DIR_MULTIQC" "$SUMMARY_FILE"