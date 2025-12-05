# blood-mice-toxoplasma-2025

## Overview
This project analyzes RNA-seq data from blood tissue of mice infected with *Toxoplasma gondii*.  
We compare gene expression between two genotypes: wildtype (WT) and interferon alpha/gamma receptor double knockout (dKO), with animals either infected or uninfected (controls).  

The dataset is a subset of blood samples from Singhania et al. (2019), available on GEO (accession [GSE119855](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE119855)). Libraries are paired-end, strand-specific, and sequenced on Illumina HiSeq 4000.

## Project Structure
rna_seq/blood-mice-toxoplasma-2025
├─ data/
│  ├─ external_data/      # Raw FASTQ files or downloaded data (GEO)
│  ├─ fastqc_output/      # FASTQC results
│  ├─ multiqc/            # MultiQC summary report
├─ logs/                  # Slurm job outputs and errors
│  ├─ output/
│  ├─ error/
├─ results/               # Processed results, differential expression tables
├─ scripts/               # All analysis scripts
├─ .gitignore             # Git ignore rules
├─ README.md

## Tools
### IBU cluster tools

| Tool              | Version / Container / Module                                    | Purpose                        |
|-------------------|-----------------------------------------------------------------|--------------------------------|
| FastQC            | 0.11.9 (Java 11, module) / 0.12.1 container                     | Read quality control           |
| fastp             | 0.24.1 (Apptainer container)                                    | Adapter and quality trimming   |
| HISAT2            | 2.2.1 (Apptainer container)                                     | Read alignment                 |
| Samtools          | 1.15+ (Apptainer container)                                     | SAM/BAM processing             |
| featureCounts     | 2.0.6 (Apptainer container)                                     | Gene-level read quantification |
| MultiQC           | 1.11 (module)                                                   | Aggregation of QC reports      |

### R (local computer) tools and packages

| Package / Tool | Version | Purpose                                                |
|----------------|---------|--------------------------------------------------------|
| DESeq2         | 1.40+   | Differential expression analysis of RNA-seq count data |
| ggplot2        | 3.4+    | Data visualization (PCA, volcano plots)                |
| ggrepel        | 0.9+    | Non-overlapping text labels for plots                  |
| apeglm         | 1.22+   | Log2 fold change shrinkage for DESeq2                  |


## Workflow

### IBU Cluster

1. **Quality Control**  
   - Run FastQC on raw reads: `for_loop_fastqc.sh`  
   - Trim adapters and low-quality bases with `fastp`: `run_fastqc_trimming.sh`  
   - Re-run FastQC on trimmed reads: `for_loop_fastqc_trimmed.sh`  
   - Aggregate FastQC summaries with MultiQC: `multiqc.sh`

2. **Genome Indexing**  
   - Reference genome: `Mus_musculus.GRCm39.dna.primary_assembly.fa`  
   - Annotation: `Mus_musculus.GRCm39.115.gtf`  
   - HISAT2 index generation: `hisat2_index.sh`

3. **Read Mapping**  
   - Map paired-end trimmed reads with HISAT2 (`for_loop_hisat2.sh`)  
   - Convert SAM to BAM, sort, and index using Samtools (inside HISAT2 container)

4. **Gene-level Quantification**  
   - featureCounts (`featureCounts.sh`) 
   - Generate raw gene count matrices and summary reports  
   - Aggregate summary metrics with MultiQC (`multiqc_featureCounts.sh`)

5. **Slurm Resource Settings**

| Step             | CPUs| Memory | Time      |
|------------------|-----|--------|-----------|
| FastQC           | 1   | 4G     | 02:00:00  |
| Trimming         | 2   | 2G     | 00:20:00  |
| HISAT2 indexing  | 1   | 8G     | 04:00:00  |
| HISAT2 mapping   | 4   | 16G    | 12:00:00  |
| featureCounts    | 4   | 4G     | 02:00:00  |
| MultiQC          | 1   | 1G     | 00:20:00  |


### R locally

The main R script is located at 

After generating gene counts, downstream analyses were performed locally in R:

1. **Load count and metadata tables**
2. **Create DESeq2 dataset**  
   - Design: `~ genotype * infection`
3. **Filter low-count genes**  
   - Genes with total counts ≤1 removed
4. **Variance-stabilizing transformation (VST)**
5. **Exploratory analysis**  
   - PCA for sample clustering (by genotype and infection)
6. **Differential expression analysis**  
   - Pairwise contrasts for genotype (dKO vs WT), infection (Case vs Control), and genotype:infection interaction  
   - Shrink log2 fold changes with `apeglm`
   - Volcano plots for all contrasts
7. **Overrepresentaion analysis**
   -


Scripts, not data (data in cloud)
Link  Reference genome LARGE FILE STORAGE (commit files)

Requirements. txt
Environment.jamo?