# blood-mice-toxoplasma-2025

## Overview
This project analyzes RNA-seq data from blood tissue of mice infected with *Toxoplasma gondii*. 
We compare gene expression between different genetic backgrounds (wildtype vs. interferon receptor knockout).

## Folder Structure
rna_seq/blood-mice-toxoplasma-2025  <br>
├─ data/  <br>
│  ├─ external_data/        # raw fastq files or downloaded data (GEO)  <br>
│  ├─ fastqc_output/        # FastQC results  <br>
│  ├─ multiqc/              # MultiQC summary report  <br>
├─ logs/                    # Slurm job outputs and errors  <br>
│  ├─ output/  <br>
│  ├─ error/  <br>
├─ results/                 # Processed results, differential expression tables  <br>
├─ scripts/                 # All analysis scripts  <br>
├─ .gitignore               # Files/folders Git should ignore  <br>
├─ README.md  <br>


## Software / Tools
- FASTQC (quality control)
- HISAT2 (alignment)
- Samtools (processing BAM files)
- featureCounts (counting reads per gene)
- R (differential expression analysis)

## How to Run
### FASTQC

