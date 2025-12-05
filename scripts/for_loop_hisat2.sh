#!/bin/bash
#------------------------------
# Slurm job configuration
#------------------------------
#SBATCH --job-name=for_loop_hisat2
#SBATCH --output=/data/users/awolfruiz/rna_seq/blood-mice-toxoplasma-2025/logs/output/hisat2/for_loop_hisat2_%j.out
#SBATCH --error=/data/users/awolfruiz/rna_seq/blood-mice-toxoplasma-2025/logs/error/hisat2/for_loop_hisat2_%j.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --partition=pibu_el8

# Paths to the hisat2 container, to the genome_index files, to the fastq reads, and to the output directory
CONTAINER="/containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif"
INDEX_BASENAME="/data/users/awolfruiz/rna_seq/blood-mice-toxoplasma-2025/results/index_files/genome_index"
FASTQ_DIR="/data/users/awolfruiz/rna_seq/blood-mice-toxoplasma-2025/data/external_data/fastp_output"
RESULTS_DIR="/data/users/awolfruiz/rna_seq/blood-mice-toxoplasma-2025/results/mapped_trimmed_reads"

# Create the output directory (if it doesn't exist already)
mkdir -p $RESULTS_DIR

# Loop over all paired-end FASTQ files
for R1 in $FASTQ_DIR/*_R1_trimmed.fastq.gz; do     # Loop through all R1 files
    SAMPLE=$(basename $R1 _R1_trimmed.fastq.gz)    # Extract sample name from filename
    R2="$FASTQ_DIR/${SAMPLE}_R2_trimmed.fastq.gz"  # Corresponding R2 file

    # Map reads with HISAT2
    # Produces a SAM file containing all alignments for this sample
    apptainer exec --bind /data/ $CONTAINER hisat2 -x $INDEX_BASENAME -1 $R1 -2 $R2 -S $RESULTS_DIR/${SAMPLE}.sam -p $SLURM_CPUS_PER_TASK

    # Convert SAM to BAM
    # BAM is compressed and easier to work with
    apptainer exec --bind /data/ $CONTAINER samtools view -hbS $RESULTS_DIR/${SAMPLE}.sam > $RESULTS_DIR/${SAMPLE}.bam

    # Sort BAM
    # Sort by genomic coordinates; -m specifies memory per thread
    apptainer exec --bind /data/ $CONTAINER samtools sort -m 2G -@ $SLURM_CPUS_PER_TASK -o $RESULTS_DIR/${SAMPLE}.sorted.bam -T temp $RESULTS_DIR/${SAMPLE}.bam

    # Index BAM
    # Required for visualization and downstream tools
    apptainer exec --bind /data/ $CONTAINER samtools index $RESULTS_DIR/${SAMPLE}.sorted.bam

    # Remove intermediate files
    # Saves space by deleting SAM and unsorted BAM
    rm $RESULTS_DIR/${SAMPLE}.sam $RESULTS_DIR/${SAMPLE}.bam
done