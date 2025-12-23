###############################################################################
# RNA-seq Differential Expression Analysis Pipeline (DESeq2)
#
# Organism: Mus musculus
#
# Loading and preprocessing gene-level count data and sample metadata
# Exploratory data analysis using PCA
# Differential expression analysis using a genotype × infection interaction model
# Visualization of interaction effects (volcano plot, gene-level expression)
# Gene Ontology overrepresentation analysis (GO Biological Process)
#
###############################################################################

# ------------------------------------------------------------------------------
# Load libraries
# ------------------------------------------------------------------------------

library(DESeq2) # Core package for differential expression analysis of RNA-seq count data
library(ggplot2)   # Used for plotting PCA and normalized counts
library(ggrepel)   # Adds non-overlapping text labels to plots
library(enrichplot) # Visualize results of gene enrichment analysis
library(EnhancedVolcano) # Enhanced Volcano Plots
library(clusterProfiler) # Functional enrichment analysis
library(org.Mm.eg.db)    # Genome-wide annotation data for Mus musculus genome

# ------------------------------------------------------------------------------
# Set directories
# ------------------------------------------------------------------------------

# Input directory (counts + metadata)
data_dir <- "/Users/anawolf/Desktop"

# Output directory
base_dir <- "/Users/anawolf/Desktop"

# ------------------------------------------------------------------------------
# Load counts and metadata
# ------------------------------------------------------------------------------

# Read counts. Use read.delim which defaults to tab separator
counts <- read.delim(
  file.path(data_dir, "gene_counts_counts_only.txt"),
  row.names = 1,
  check.names = FALSE
)

# Read metadata
coldata <- read.delim(
  file.path(data_dir, "metadata.txt"),
  row.names = 1,
  check.names = FALSE
)

# ------------------------------------------------------------------------------
#  Exploratory Data Analysis: Prepare Metadata + Counts
# ------------------------------------------------------------------------------

# Trim whitespace from column and row names
colnames(counts) <- trimws(colnames(counts))
rownames(coldata) <- trimws(rownames(coldata))

# Convert metadata columns to categorical factors
# genotype: WT or DKO
# infection: Control or Case
coldata$genotype  <- factor(coldata$genotype)
coldata$infection <- factor(coldata$infection)

# Set the reference level for each factor
# This determines how DESeq2 interprets contrasts
coldata$infection <- relevel(coldata$infection, ref = "Control")
coldata$genotype  <- relevel(coldata$genotype,  ref = "WT")

# Match sample order between counts and metadata
counts <- counts[, rownames(coldata)]

# ------------------------------------------------------------------------------
# Differential Expression Modeling: DESeq2 Setup
# ------------------------------------------------------------------------------

# Create DESeq2 dataset using interaction model
# design = ~ genotype * infection tests:
#   1) main effect of genotype
#   2) main effect of infection
#   3) interaction: specific infection response in DKO vs WT
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData   = coldata,
  design    = ~ genotype * infection
)

# Run DESeq2 workflow (estimation + model fitting)
dds <- DESeq(dds)

# ------------------------------------------------------------------------------
# Sample-Level Visualization: Principal Component Analysis (PCA)
# ------------------------------------------------------------------------------

# Variance-stabilizing transformation
# blind = TRUE ignores experimental design during transformation 
vsd <- vst(dds, blind = TRUE)

# Extract PCA coordinates with metadata
pca_data <- plotPCA(vsd, intgroup = c("genotype", "infection"), returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

# PCA Plot: visualize how samples cluster by genotype and infection
pca_plot <- ggplot(
  pca_data,
  aes(x = PC1, y = PC2, color = genotype, shape = infection, label = name)
) +
  geom_point(size = 3) +
  geom_text_repel(size = 3, max.overlaps = Inf) +
  scale_color_manual(values = c("WT" = "cornflowerblue", "DKO" = "brown")) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid.major = element_line(color = "grey80", size = 0.5),
    panel.grid.minor = element_line(color = "grey90", size = 0.25),
    legend.position = "right"
  ) +
  coord_fixed() +
  ggtitle("Principal Component Analysis of RNA-seq Samples")

# Display PCA plot
print(pca_plot)

# Save PCA plot in base directory
ggsave(
  filename = file.path(base_dir, "PCA_plot.png"),
  plot = pca_plot,
  width = 8,
  height = 6,
  dpi = 300
)

# ------------------------------------------------------------------------------
# Differential Expression Analysis: Interaction Term (Genotype × Infection)
# ------------------------------------------------------------------------------

# Identify interaction term name 
interaction_name <- grep("genotype.*infection.*Case", resultsNames(dds), value = TRUE)

# Extract results for interaction effect
# This tests whether the infection response differs between DKO and WT
# alpha = 0.05 sets the FDR (adjusted p-value) threshold used in summary statistics 
res_interaction <- results(dds, name = interaction_name, alpha = 0.05)

# Count of DE genes based on the default statistical cutoffs
summary(res_interaction)

# Map ENSEMBL IDs to gene symbols
res_interaction$gene_name <- mapIds(
  org.Mm.eg.db,
  keys       = rownames(res_interaction),
  column     = "SYMBOL",
  keytype    = "ENSEMBL",
  multiVals  = "first"
)

# Remove rows with missing gene names
res_interaction <- res_interaction[!is.na(res_interaction$gene_name), ]

# Sort by adjusted p-value (most significant first)
res_interaction <- res_interaction[order(res_interaction$padj), ]

# Save full interaction genes
write.csv(as.data.frame(res_interaction),
          file.path(base_dir, "Interaction_DESeq2_Results.csv"))

# Remove rows with missing p-values or LFC values
res_interaction <- res_interaction[!is.na(res_interaction$padj) &
                                     !is.na(res_interaction$log2FoldChange), ]

# Filter significant interaction genes
# FDR < 0.05 and |log2FC| > 1 threshold
sig_interaction <- res_interaction[
  res_interaction$padj < 0.05 & abs(res_interaction$log2FoldChange) > 1,
]

# Save significant interaction genes
write.csv(as.data.frame(sig_interaction),
          file.path(base_dir, "Interaction_DESeq2_SignificantGenes.csv"))

# ------------------------------------------------------------------------------
# Gene-Level Visualization: Volcano Plot of Differential Expression
# ------------------------------------------------------------------------------

# Split by direction of effect
up_genes   <- sig_interaction[sig_interaction$log2FoldChange > 1, ]
down_genes <- sig_interaction[sig_interaction$log2FoldChange < -1, ]

# Sort each group by significance and effect size
up_sorted   <- up_genes[order(up_genes$padj, -abs(up_genes$log2FoldChange)), ]
down_sorted <- down_genes[order(down_genes$padj, -abs(down_genes$log2FoldChange)), ]

# Select top 10 from each direction (to label them)
top10_up <- head(up_sorted$gene_name, 10)
top10_down <- head(down_sorted$gene_name, 10)

top20_genes <- c(top10_up, top10_down)

# Define custom colors for volcano categories
keyvals <- rep("grey70", nrow(res_interaction))
names(keyvals) <- rep("NS", nrow(res_interaction))

keyvals[res_interaction$padj < 0.05 & res_interaction$log2FoldChange > 1]  <- "brown"
keyvals[res_interaction$padj < 0.05 & res_interaction$log2FoldChange < -1] <- "cornflowerblue"

names(keyvals)[res_interaction$padj < 0.05 & res_interaction$log2FoldChange > 1]  <- "Upregulated"
names(keyvals)[res_interaction$padj < 0.05 & res_interaction$log2FoldChange < -1] <- "Downregulated"

# Create volcano plot
volcano_plot <- EnhancedVolcano(
  res_interaction,
  lab = res_interaction$gene_name,
  selectLab = top20_genes,
  x = "log2FoldChange",
  y = "padj",
  pCutoff = 0.05,
  FCcutoff = 1,
  pointSize = 2,
  labSize = 4.5,
  colCustom = keyvals,
  colAlpha = 0.7,
  drawConnectors = TRUE,
  widthConnectors = 0.6,
  boxedLabels = TRUE,
  title = "Differential Expression for Infection × Genotype Interaction",
  subtitle = "Interaction: (DKO_case – DKO_ctrl) – (WT_case – WT_ctrl)",
  legendPosition = "right",
  legendLabSize = 12,
  legendIconSize = 4.0,
  caption = ""  
)

# Display volcano plot
print(volcano_plot)

# Save volcano plot
ggsave(
  file.path(base_dir, "Volcano_plot.png"),
  plot = volcano_plot,
  width = 10, height = 8, dpi = 300
)

# ------------------------------------------------------------------------------
# Gene-Level Visualization: Normalized Expression Counts
# ------------------------------------------------------------------------------

# Custom function to plot normalized counts for a given gene across groups
# Highlighting the interaction effect (log2 fold change and adjusted p-value)
plotCounts_custom <- function(dds, gene, intgroup, title_gene = NULL) {
  require(ggplot2)
  
  # Extract normalized counts
  d <- plotCounts(dds, gene = gene, intgroup = intgroup, returnData = TRUE)
  
  # Combined group labels for visualization
  d$group <- interaction(d[, intgroup], sep = "_")
  
  # Mean per group
  group_means <- aggregate(count ~ group, d, mean)
  
  # If no title gene provided, show Ensembl ID
  if (is.null(title_gene)) title_gene <- gene
  
  # Retrieve log2FC and padj for this gene from DESeq2 results
  gene_lfc <- res_interaction[gene, "log2FoldChange"]
  gene_padj <- res_interaction[gene, "padj"]
  
  # Legend
  legend_text <- paste0("log2FC = ", round(gene_lfc, 2), "\npadj = ", signif(gene_padj, 3))
  d$Legend <- legend_text
  
  ggplot(d, aes(x = group, y = count)) +
    geom_point(color = "black", size = 3, alpha = 0.85) +
    geom_line(data = group_means, aes(x = group, y = count, group = 1), color = "black", linewidth = 0.7) +
    geom_point(aes(color = Legend), size = 0, alpha = 0) +
    scale_color_manual(values = setNames("black", legend_text)) +
    theme_bw(base_size = 14) +
    theme(
      plot.title = element_text(size = 18, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "right",
      legend.title = element_blank(),
      legend.text = element_text(size = 12)
    ) +
    xlab("Group") +
    ylab("Normalized Counts (log10)") +
    ggtitle(paste0("Expression of ", title_gene)) +
    labs(subtitle = "Interaction Effect (Normalized Counts)")
}


# Choose genes to visualize from "Interaction_DESeq2_Results.csv" or "Interaction_DESeq2_SignificantGenes.csv"

# ----------------------------------------------------------
# IRF1 (Downregulated in interaction contrast)
# ----------------------------------------------------------

plot_irf1 <- plotCounts_custom(
  dds,
  gene = "ENSMUSG00000018899",
  intgroup = c("genotype", "infection"),
  title_gene = "Irf1"
)

print(plot_irf1)

ggsave(
  file.path(base_dir, "Irf1_NormalizedCounts.png"),
  plot_irf1, width = 6, height = 5, dpi = 300
)

# ----------------------------------------------------------
# STAT2 (Downregulated in interaction contrast)
# ----------------------------------------------------------

plot_stat2 <- plotCounts_custom(
  dds,
  gene = "ENSMUSG00000040033",  
  intgroup = c("genotype", "infection"),
  title_gene = "Stat2"
)

print(plot_stat2)

ggsave(
  file.path(base_dir, "Stat2_NormalizedCounts.png"),
  plot_stat2, width = 6, height = 5, dpi = 300
)

# ----------------------------------------------------------
# FECH (Upregulated in interaction contrast)
# ----------------------------------------------------------

plot_fech <- plotCounts_custom(
  dds,
  gene = "ENSMUSG00000024588",
  intgroup = c("genotype", "infection"),
  title_gene = "Fech"
)

print(plot_fech)

ggsave(
  file.path(base_dir, "Fech_NormalizedCounts.png"),
  plot_fech, width = 6, height = 5, dpi = 300
)

# ----------------------------------------------------------
# ISCA1 (Upregulated in interaction contrast)
# ----------------------------------------------------------

plot_isca1 <- plotCounts_custom(
  dds,
  gene = "ENSMUSG00000044792", 
  intgroup = c("genotype", "infection"),
  title_gene = "Isca1"
)

print(plot_isca1)

ggsave(
  file.path(base_dir, "Isca1_NormalizedCounts.png"),
  plot_isca1, width = 6, height = 5, dpi = 300
)

# ------------------------------------------------------------------------------
# Gene Ontology overrepresentation Analysis (Enrichment Analysis)
# ------------------------------------------------------------------------------

# Define DE genes for up- and downregulated separately (Ensembl IDs)
up_genes <- rownames(sig_interaction[sig_interaction$log2FoldChange > 0, ])
down_genes <- rownames(sig_interaction[sig_interaction$log2FoldChange < 0, ])

# Universe: all genes tested in DESeq2
all_genes <- rownames(res_interaction)

# ----------------------------------------------------------
# Upregulated genes
# ----------------------------------------------------------
go_up <- enrichGO(
  gene         = up_genes,
  universe     = all_genes,
  OrgDb        = org.Mm.eg.db,
  ont          = "BP",
  keyType      = "ENSEMBL",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable     = TRUE
)

go_up_df <- as.data.frame(go_up)
write.csv(go_up_df, file.path(base_dir, "GO_Enrichment_Upregulated.csv"))

# Dotplot and Barplot
go_up_dot <- dotplot(go_up, showCategory = 10, title = "Top Enriched GO Terms - Upregulated")
print(go_up_dot)
ggsave(file.path(base_dir, "GO_Upregulated_Dotplot.png"), plot = go_up_dot, width = 10, height = 8, dpi = 300)

go_up_bar <- barplot(go_up, showCategory = 10, title = "Top Enriched GO Terms - Upregulated")
print(go_up_bar)
ggsave(file.path(base_dir, "GO_Upregulated_Barplot.png"), plot = go_up_bar, width = 10, height = 8, dpi = 300)

# ----------------------------------------------------------
# Downregulated genes
# ----------------------------------------------------------
go_down <- enrichGO(
  gene         = down_genes,
  universe     = all_genes,
  OrgDb        = org.Mm.eg.db,
  ont          = "BP",
  keyType      = "ENSEMBL",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable     = TRUE
)

go_down_df <- as.data.frame(go_down)
write.csv(go_down_df, file.path(base_dir, "GO_Enrichment_Downregulated.csv"))

# Dotplot and Barplot
go_down_dot <- dotplot(go_down, showCategory = 10, title = "Top Enriched GO Terms - Downregulated")
print(go_down_dot)
ggsave(file.path(base_dir, "GO_Downregulated_Dotplot.png"), plot = go_down_dot, width = 10, height = 8, dpi = 300)

go_down_bar <- barplot(go_down, showCategory = 10, title = "Top Enriched GO Terms - Downregulated")
print(go_down_bar)
ggsave(file.path(base_dir, "GO_Downregulated_Barplot.png"), plot = go_down_bar, width = 10, height = 8, dpi = 300)
