# blood-mice-toxoplasma-2025

## Project overview
The aim of this study is to investigate how interferon α/γ receptor deficiency alters the systemic transcriptional response to *Toxoplasma gondii* infection in mice. By analyzing blood RNA-seq data from wildtype (WT) and double-knockout (dKO) mice under infected and uninfected conditions, we use a two-factor experimental design to examine the effects of genotype, infection, and their interaction.
 
The dataset used is a subset of blood samples from [Singhania et al. (2019)](https://www.nature.com/articles/s41467-019-10601-6), available on GEO (accession [GSE119855](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE119855)). Libraries are paired-end, strand-specific, and sequenced on Illumina HiSeq 4000.

The first part of the workflow is performed using the IBU cluster of the University of Bern. The scripts ran in the cluster are in the folder [scripts](https://github.com/anawolf1/blood-mice-toxoplasma-2025/tree/main/scripts). The final part of the workflow (differential expression analysis) was executed locally in R, with a single [script](https://github.com/anawolf1/blood-mice-toxoplasma-2025/blob/main/rnaseq.R), found in the main folder of the repository.

### Workflow overview
1. Obtain FASTQ files from GEO (GSE119855) or use the copies provided on the IBU cluster (/data/courses/rnaseq_course/toxoplasma_de).
2. Run the preprocessing and alignment scripts in `scripts/` using Slurm, following the numerical order of the scripts (quality control, read mapping, and gene-level quantification)..
3. Download the resulting gene count matrix and run rnaseq.R locally using the count matrix and the [sample metadata file](https://github.com/anawolf1/blood-mice-toxoplasma-2025/blob/main/metadata.txt)to perform exploratory analysis, differential expression analysis, and GO overrepresentation analysis.

### Setup and usage:
- Update the path of the working directories in all `.sh` and `.R` scripts before running.
- Adjust Slurm `#SBATCH` directives in `.sh` scripts according to your job’s CPU, memory, and time needs. Recommended defaults are provided in each script.
- Optional step: IGV inspection of BAM files.

## Project Structure (IBU cluster)
```
rna_seq/blood-mice-toxoplasma-2025
├─ data/
│  ├─ external_data/            # Raw FASTQ files or downloaded data (GEO)
│     ├─ fastp_output/          # Trimmed reads
│     ├─ fastqc_output/         # FASTQC results
│     ├─ fastqc_trimmed_output/ # Trimmed FASTQC results
├─ logs/                        # Slurm job outputs and errors
│  ├─ output/
│  ├─ error/
├─ resources/                   # Reference genome DNA FASTA file and gtf file.
├─ results/                     # Processed results
├─ scripts/                     # All analysis scripts
├─ .gitignore             
├─ README.md
```

## Tools
### IBU cluster tools

| Tool              | Version                         | Purpose                        |
|-------------------|---------------------------------|--------------------------------|
| FastQC            | 0.11.9 (Module)                 | Read quality control           |
| fastp             | 0.24.1 (Apptainer Container)    | Adapter and quality trimming   |
| HISAT2            | 2.2.1 (Apptainer Container)     | Read alignment                 |
| Samtools          | 1.15  (Apptainer Container)     | SAM/BAM processing             |
| featureCounts     | 2.0.6 (Apptainer Container)     | Gene-level read quantification |
| MultiQC           | 1.11 (Module)                   | Aggregation of QC reports      |

### R (local computer) tools 
R analyses were performed using 2025.09.1+401 on macOS/Linux.

| Package / Library | Version | Purpose                                                |
|-------------------|---------|--------------------------------------------------------|
| DESeq2            | 1.50.2  | Differential expression analysis of RNA-seq count data |
| ggplot2           | 4.0.1   | Data visualization (PCA, normalized counts)            |
| ggrepel           | 0.9.6   | Non-overlapping text labels for plots                  |
| enrichplot        | 1.30.3  | Visualize results of gene enrichment analysis          |
| EnhancedVolcano   | 1.29.1  | Enhanced Volcano Plots                                 |
| clusterprofiler   | 4.18.2  | Functional enrichment analysis                         |
| org.Mm.eg.db      | 3.22.0  | Genome-wide annotation data for *Mus musculus* genome    |


## Workflow

### IBU Cluster workflow

1. **Quality Control**  
   - Run FastQC on raw reads: `for_loop_fastqc.sh`  
   - Trim adapters and low-quality bases with `fastp`: `run_fastqc_trimming.sh`  
   - Run FastQC on trimmed reads: `for_loop_fastqc_trimmed.sh`  
   - Aggregate FastQC summaries with MultiQC: `multiqc.sh`

2. **Genome Indexing**  
   - Reference genome: [Mus_musculus.GRCm39.dna.primary_assembly.fa](https://ftp.ensembl.org/pub/release-115/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz)  
   - Annotation: [Mus_musculus.GRCm39.115.gtf](https://ftp.ensembl.org/pub/release-115/gtf/mus_musculus/Mus_musculus.GRCm39.115.gtf.gz)
   - HISAT2 index generation: `hisat2_index.sh`

3. **Read Mapping**  
   - Map paired-end trimmed reads with HISAT2 (`for_loop_hisat2.sh`)  
   - Convert SAM to BAM, sort, and index using Samtools

4. **Gene-level Quantification**  
   - Generate raw gene count matrices and summary reports with featureCounts (`featureCounts.sh`)
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


### R workflow

1. **Load required libraries**: DESeq2, ggplot2, ggrepel, enrichplot, EnhancedVolcano, clusterProfiler, org.Mm.eg.db.
2. **Set directories**: Input (counts + metadata) and output paths.
3. **Load count and metadata tables**
4. **Data preprocessing**: Trim whitespace, convert metadata columns to factors, set reference levels, ensure counts and metadata match.
5. **Create DESeq2 dataset and run DESeq2 workflow**: The dataset is built using the interaction design, `~ genotype * infection`.
6. **Variance-stabilizing transformation and exploratory analysis**: Perform PCA and visualize sample clustering by genotype and infection.
7. **Differential expression analysis**: Extract results for the genotype × infection interaction term and filter significant genes by FDR < 0.05 and |log2FC| > 1.
8. **Visualize gene level results**: Generate volcano plots and normalized count plots for selected genes.
9. **Gene Ontology (GO) overrepresentation analysis**: Perform enrichment analysis for up and downregulated genes. Produce dotplots and barplots of top enriched biological processes.

