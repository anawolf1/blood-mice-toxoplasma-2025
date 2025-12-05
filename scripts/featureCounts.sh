#!/bin/bash
#------------------------------
# Slurm job configuration
#------------------------------
#SBATCH --job-name=featureCounts
#SBATCH --output=/data/users/awolfruiz/rna_seq/blood-mice-toxoplasma-2025/logs/output/featureCounts/featureCounts_%j.out
#SBATCH --error=/data/users/awolfruiz/rna_seq/blood-mice-toxoplasma-2025/logs/error/featureCounts/featureCounts_%j.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH --partition=pibu_el8


# Paths to container, annotation genome, mapped reads, and the output directory
CONTAINER="/containers/apptainer/subread_2.0.6.sif"
ANNOTATION="/data/users/awolfruiz/rna_seq/blood-mice-toxoplasma-2025/resources/Mus_musculus.GRCm39.115.gtf.gz"
BAM_DIR="/data/users/awolfruiz/rna_seq/blood-mice-toxoplasma-2025/results/mapped_trimmed_reads"
OUT_DIR="/data/users/awolfruiz/rna_seq/blood-mice-toxoplasma-2025/results/featureCounts"

# Create the output directory (if it doesn't exist already)
mkdir -p $OUT_DIR

# List all sorted BAM files
BAM_FILES=$(ls $BAM_DIR/*sorted.bam)

# Run featureCounts command
apptainer exec --bind /data $CONTAINER featureCounts \
    -T $SLURM_CPUS_PER_TASK \
    -p \
    -s 2 \
    -a $ANNOTATION \
    -o $OUT_DIR/gene_counts.txt \
    $BAM_FILES

# -T is the number of threads
# -p indicates our data are paired-end reads
# -s indicates the strandedness (in this case, its 2, meaning reverse-stranded RNA-seq)
# -a indicates the annotation file
# -o indicates the output count matrix
