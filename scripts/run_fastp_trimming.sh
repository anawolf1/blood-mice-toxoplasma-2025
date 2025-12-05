#!/bin/bash
#------------------------------
# Slurm job configuration
#------------------------------
#SBATCH --job-name=fastp
#SBATCH --array=0-14
#SBATCH --output=/data/users/awolfruiz/rna_seq/blood-mice-toxoplasma-2025/logs/output/fastp/fastp_trimming_%A_%a.out.out
#SBATCH --error=/data/users/awolfruiz/rna_seq/blood-mice-toxoplasma-2025/logs/error/fastp/fastp_trimming_%A_%a.out.err
#SBATCH --cpus-per-task=2
#SBATCH --mem=2G
#SBATCH --time=02:00:00
#SBATCH --partition=pibu_el8

# Path to the fastp container, to the input directory and to the output directory
CONTAINER=/containers/apptainer/fastp_0.24.1.sif
INPUT_DIR=/data/courses/rnaseq_course/toxoplasma_de/reads_Blood
OUTPUT_DIR=/data/users/awolfruiz/rna_seq/blood-mice-toxoplasma-2025/data/external_data/fastp_output

# Create the output directory (if it doesn't exist already)
mkdir -p "$OUTPUT_DIR"

# Define variables of paired-end reads
R1_FILES=(${INPUT_DIR}/*_1.fastq.gz)          # Create an array of all R1 FASTQ files in the input directory 
R1_FILE=${R1_FILES[$SLURM_ARRAY_TASK_ID]}     # Select the R1 file corresponding to this SLURM array task
R2_FILE=${R1_FILE/_1.fastq.gz/_2.fastq.gz}    # Select the corresponding R2 FASTQ files
SAMPLE=$(basename "$R1_FILE" _1.fastq.gz)     # Obtain sample name

# Run Fastp
apptainer exec \
  --bind "$INPUT_DIR":"$INPUT_DIR" \
  --bind "$OUTPUT_DIR":"$OUTPUT_DIR" \
  "$CONTAINER" \
  fastp \
    -i "$R1_FILE" \
    -I "$R2_FILE" \
    -o "$OUTPUT_DIR/${SAMPLE}_R1_trimmed.fastq.gz" \
    -O "$OUTPUT_DIR/${SAMPLE}_R2_trimmed.fastq.gz" \
    -j "$OUTPUT_DIR/${SAMPLE}.json" \
    -h "$OUTPUT_DIR/${SAMPLE}.html" \
    -w $SLURM_CPUS_PER_TASK \
    --detect_adapter_for_pe \
    --qualified_quality_phred 20 \
    --length_required 30
    
# -i indicates read1 input
# -o indicates read2 input
# -j is used to generate JSON reports
# -h is used to generate html reports
# -w indicates number of threads to use
# --detect_adapter_for_pe automatically detects adapters for paired-end reads
# ---qualified_quality_phred indicates the minimum quality score to keep a base (in this case 20, corresponding to less than 1% error rate)
# --length_required indicates the minimum length of reads to keep after trimming (in this case 30, dicarding reads shorter than 30 bp after trimming)