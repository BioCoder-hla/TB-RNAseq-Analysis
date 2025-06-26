#-------------------------------------------------------------------------------
# Load all necessary libraries for the entire analysis workflow.

# Set working directory (adjust the path as needed)
setwd("/home/wayne/Projects/Data Analytics/Gene Expression")

# Load libraries
library(edgeR)
library(kableExtra)
library(DESeq2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(pheatmap)
library(ggplot2)
library(EnhancedVolcano)
library(readxl)
library(dplyr)
library(tibble)
library(GEOquery)


################################################################################
### RNA-Seq Differential Expression Analysis: TB Progressor Study (GSE107994)
################################################################################


# STAGE 0: DOWNLOADING THE DATASET

# Download GEO data and metadata
gse <- getGEO("GSE107994", GSEMatrix = TRUE, getGPL = TRUE)[[1]]
write.csv(pData(gse), "GSE107994_metadata.csv", row.names = FALSE)
write.csv(exprs(gse), "GSE107994_normalized_expression_data.csv", row.names = TRUE)

# Increase timeout and download raw counts
options(timeout = 300)  # Set timeout to 5 minutes
download.file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE107994&format=file&file=GSE107994_Raw_counts_Leicester_with_progressor_longitudinal.xlsx",
              "GSE107994_Raw_counts_Leicester_with_progressor_longitudinal.xlsx", mode = "wb")
raw_counts_df <- read_excel("GSE107994_Raw_counts_Leicester_with_progressor_longitudinal.xlsx")
write.csv(raw_counts_df, "GSE107994_raw_counts.csv", row.names = FALSE)


# STAGE 1: DATA LOADING & PREPARATION
#-------------------------------------------------------------------------------
# Goal: Load raw counts and metadata, then create clean, aligned objects
#       for analysis: a numeric count matrix and a sample metadata table.

# --- A. DEFINE FILE PATHS AND LOAD RAW DATA ---
counts_file <- "GSE107994_Raw_counts_Leicester_with_progressor_longitudinal.xlsx"
metadata_file <- "GSE107994_metadata.csv"

raw_counts_df <- read_excel(counts_file)
metadata_df <- read.csv(metadata_file)

# --- B. INSPECT RAW DATA ---
print("--- Raw Counts File Head (first 5 columns) ---")
head(raw_counts_df[, 1:5])

print("--- Raw Metadata File Head ---")
head(metadata_df)

# --- C. PROCESS METADATA (colData) ---
# Create the colData table from the metadata file.
colData <- metadata_df %>%
  select(Sample_ID = title, condition = group.ch1) %>%
  column_to_rownames(var = "Sample_ID")

# Recode condition names for clarity and simplicity
colData <- colData %>%
  mutate(condition = recode(condition,
                            "LTBI_Progressor" = "Progressor",
                            "Active_TB" = "Active",
                            "LTBI" = "Latent"))

# Convert 'condition' to a factor and set 'Control' as the reference level
colData$condition <- factor(colData$condition, levels = c("Control", "Latent", "Active", "Progressor"))

# --- D. PROCESS COUNT DATA (counts_matrix) ---
# Convert the raw counts data frame into a numeric matrix with genes as rownames.
counts_matrix <- raw_counts_df %>%
  column_to_rownames(var = "Genes") %>%
  select(starts_with("Leicester_")) %>%
  as.matrix()

# --- E. ALIGN AND VALIDATE ---
# Ensure the columns of the count matrix match the rows of the metadata table.
colData <- colData[colnames(counts_matrix), , drop = FALSE]

# Final validation: This must return TRUE for the analysis to proceed.
print("--- Final Alignment Check ---")
all(colnames(counts_matrix) == rownames(colData))

# Inspect final prepared objects
print("--- Final Processed Metadata (colData) ---")
head(colData)
print("--- Final Condition Counts ---")
table(colData$condition)
print("--- Final Count Matrix Dimensions ---")
dim(counts_matrix)
print("--- Final Metadata Dimensions ---")
dim(colData)


# STAGE 2: QUALITY CONTROL, FILTERING & NORMALIZATION
#-------------------------------------------------------------------------------
# Goal: Remove low-expression genes, normalize for library size differences,
#       and visualize data quality before and after these steps.

# --- A. CREATE DGEList OBJECT AND FILTER LOW-EXPRESSION GENES ---
# Create a DGEList object, which is the main data container for edgeR
dge <- DGEList(counts = counts_matrix, group = colData$condition)

# Calculate CPM (counts per million) to identify low-expressed genes
cpm <- cpm(dge)

# Filter out genes that don't have a CPM > 2 in at least 5 samples
keep <- rowSums(cpm > 2) >= 5
dge <- dge[keep, , keep.lib.sizes = FALSE]

print("--- Dimensions After Filtering ---")
dim(dge$counts)

# --- B. LIBRARY SIZE VISUALIZATION ---
# Define a consistent color palette for all plots
colors <- c("Control" = "#4682B4", "Latent" = "#32CD32", "Active" = "#FF4500", "Progressor" = "#8A2BE2")

# Create simplified sample names for plot axes
simplified_names <- sub("Leicester_with_progressor_longitudinal_", "", colnames(dge))

# Set plot margins to accommodate labels and legend
par(mar = c(8, 5, 4, 8))

barplot(dge$samples$lib.size,
        main = "Library Sizes of Filtered Samples",
        names.arg = simplified_names,
        las = 2,
        cex.names = 0.7,
        col = colors[dge$samples$group],
        border = NA)

title(ylab = "Total Reads", line = 4)

legend("topright",
       legend = names(colors),
       fill = colors,
       title = "Condition",
       inset = c(-0.6, 0),
       xpd = TRUE,
       bty = "o",
       cex = 0.8)

# Reset plot margins to default
par(mar = c(5, 4, 4, 2) + 0.1)

# --- C. TMM NORMALIZATION AND DISTRIBUTION VISUALIZATION ---
# Apply TMM normalization to account for library composition differences
dge_normalized <- calcNormFactors(dge)

# Calculate log-transformed CPM for visualization (before and after normalization)
logcpm_unnormalized <- cpm(dge, log = TRUE)
logcpm_normalized <- cpm(dge_normalized, log = TRUE)

# Generate boxplots to compare sample distributions
par(mfrow = c(1, 2))
boxplot(logcpm_unnormalized, las = 2, col = colors[dge$samples$group],
        main = "Un-normalized Data", ylab = "log-CPM",
        names = simplified_names, cex.axis = 0.7)
boxplot(logcpm_normalized, las = 2, col = colors[dge$samples$group],
        main = "TMM Normalized Data", ylab = "log-CPM",
        names = simplified_names, cex.axis = 0.7)
par(mfrow = c(1, 1))

# --- D. SUMMARY STATISTICS TABLE ---
# Generate a summary table of logCPM values for each condition after normalization
summary_by_condition <- lapply(levels(dge_normalized$samples$group), function(cond) {
  indices <- which(dge_normalized$samples$group == cond)
  apply(logcpm_normalized[, indices], 1, summary)
})
names(summary_by_condition) <- levels(dge_normalized$samples$group)

summary_stats <- lapply(summary_by_condition, function(x) {
  colMeans(t(x), na.rm = TRUE)
})

summary_df <- as.data.frame(do.call(cbind, summary_stats))
rownames(summary_df) <- c("Min.", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max.")

kable(summary_df, format = "html", caption = "Summary Statistics of Normalized logCPM by Condition") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"), full_width = FALSE)


# STAGE 3: EXPLORATORY DATA ANALYSIS (PCA & HEATMAP)
#-------------------------------------------------------------------------------
# Goal: Visualize sample relationships using unsupervised clustering methods
#       to see if they group by biological condition.

# --- A. PRINCIPAL COMPONENT ANALYSIS (PCA) ---
# Perform PCA on the transposed, scaled, normalized log-CPM values
pca <- prcomp(t(logcpm_normalized), scale = TRUE)

# Create a data frame for plotting with ggplot2
pca_data <- data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2], Condition = colData$condition)

# Plot PCA results
# Your original code with the added panel.border
ggplot(pca_data, aes(x = PC1, y = PC2, color = Condition, shape = Condition)) +
  geom_point(size = 4, alpha = 0.8) +
  scale_color_manual(values = colors) +
  theme_minimal(base_size = 14) +
  labs(title = "PCA of TMM-Normalized RNA-seq Data",
       x = paste0("PC1 (", round(summary(pca)$importance[2,1]*100, 1), "%)"),
       y = paste0("PC2 (", round(summary(pca)$importance[2,2]*100, 1), "%)")) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") + 
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  theme(
    legend.position = "bottom",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # This is the new line that adds the border
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1) 
  )

# --- B. HEATMAP OF TOP VARIABLE GENES ---
# Calculate the variance for each gene across all samples
gene_variances <- apply(logcpm_normalized, 1, var)

# Select the top 500 most variable genes
top_var_genes <- order(gene_variances, decreasing = TRUE)[1:200]
top_var_matrix <- logcpm_normalized[top_var_genes, ]

# Create an annotation data frame for the heatmap columns (samples)
annotation_col <- data.frame(Condition = dge_normalized$samples$group)
rownames(annotation_col) <- colnames(dge_normalized)

# Define annotation colors to match other plots
annotation_colors <- list(Condition = colors)

# Define a color scale for the heatmap body
color_scale <- colorRampPalette(c("#4575B4", "#F7F7F7", "#D73027"))(100) # Blue-White-Red

pheatmap(top_var_matrix,
         main = "Heatmap of Top 200 Most Variable Genes",
         color = color_scale,
         scale = "row",                # Scale genes to Z-score
         annotation_col = annotation_col,
         annotation_colors = annotation_colors,
         show_rownames = FALSE,
         show_colnames = FALSE,
         treeheight_row = 30,
         treeheight_col = 30,
         border_color = NA,
         clustering_method = "ward.D2",
         legend_breaks = c(-2, 0, 2),
         legend_labels = c("Low", "Mean", "High"))


# STAGE 4: DIFFERENTIAL EXPRESSION ANALYSIS
#-------------------------------------------------------------------------------
# Goal: Statistically test for gene expression differences between conditions
#       using the edgeR quasi-likelihood pipeline.

# --- A. SETUP MODEL ---
# Create a design matrix for the GLM. `~0+group` fits a mean for each group.
design <- model.matrix(~ 0 + group, data = dge_normalized$samples)
colnames(design) <- levels(dge_normalized$samples$group)

# Estimate gene-wise dispersions, a crucial step for modeling biological variability
dge_normalized <- estimateDisp(dge_normalized, design, robust = TRUE)

# Visualize the dispersion estimates with a BCV (Biological Coefficient of Variation) plot
plotBCV(dge_normalized, main="BCV Plot")

# Fit the Quasi-Likelihood (QL) General-Linear-Model
fit <- glmQLFit(dge_normalized, design, robust = TRUE)

# --- B. DEFINE CONTRASTS AND PERFORM TESTS ---
# Define the biological comparisons of interest
contrasts <- makeContrasts(
  Active_vs_Control = Active - Control,
  Progressor_vs_Latent = Progressor - Latent,
  levels = design
)

# Perform the QL F-test for each contrast
qlf_active_vs_control <- glmQLFTest(fit, contrast = contrasts[, "Active_vs_Control"])
qlf_progressor_vs_latent <- glmQLFTest(fit, contrast = contrasts[, "Progressor_vs_Latent"])

# --- C. EXTRACT AND INSPECT RESULTS ---
# Extract full results tables for each comparison
top_genes_active <- topTags(qlf_active_vs_control, n = Inf)$table
top_genes_progressor <- topTags(qlf_progressor_vs_latent, n = Inf)$table

print("--- Top Genes: Active TB vs Control ---")
head(top_genes_active)

# Filter for statistically significant DEGs (FDR < 0.05 and |logFC| > 1)
sig_genes_active <- top_genes_active[top_genes_active$FDR < 0.05 & abs(top_genes_active$logFC) > 1, ]
sig_genes_progressor <- top_genes_progressor[top_genes_progressor$FDR < 0.05 & abs(top_genes_progressor$logFC) > 1, ]

print(paste("Found", nrow(sig_genes_active), "significant DEGs for Active vs. Control"))
print(paste("Found", nrow(sig_genes_progressor), "significant DEGs for Progressor vs. Latent"))

# Save full results tables
write.csv(top_genes_active, "DEGs_Full_Active_vs_Control.csv")
write.csv(top_genes_progressor, "DEGs_Full_Progressor_vs_Latent.csv")

# --- D. VISUALIZE A SINGLE COMPARISON WITH AN MD PLOT ---
# MD (Mean-Difference) plot for the Active vs Control comparison
plotMD(qlf_active_vs_control, main = "MD Plot: Active TB vs Control")
abline(h = c(-1, 1), col = "blue") # Add lines for logFC cutoffs

# Highlight significant genes on the MD plot
de_genes <- decideTests(qlf_active_vs_control, p.value = 0.05, lfc = 1)
plotMD(qlf_active_vs_control, status = de_genes, values = c(1, -1), col = c("red", "blue"), 
       main = "MD Plot: Active TB vs Control")
legend("topleft", legend = c("Up-regulated", "Down-regulated"), pch = 19, col = c("red", "blue"))

# STAGE 5: VISUALIZATION OF DEGS (VOLCANO PLOTS)
#-------------------------------------------------------------------------------
# Goal: Visualize the results of the differential expression tests, showing
#       both statistical significance (FDR) and biological magnitude (logFC).

# Merge DEGs results with gene names from the original count data for better labels
gene_names_df <- raw_counts_df[, c("Genes", "Gene_name")]
top_genes_active_with_names <- merge(top_genes_active, gene_names_df, by.x = "row.names", by.y = "Genes", all.x = TRUE)
top_genes_progressor_with_names <- merge(top_genes_progressor, gene_names_df, by.x = "row.names", by.y = "Genes", all.x = TRUE)

# Create volcano plots for both comparisons using the user-friendly Gene_name for labels
par(mfrow = c(1, 2), mar = c(5, 5, 5, 2))


# --- Plot 1: Active vs. Control (without grid lines) ---

# 1. Create the volcano plot and store it in a variable (e.g., v_plot1)
v_plot1 <- EnhancedVolcano(top_genes_active_with_names,
                           lab = top_genes_active_with_names$Gene_name,
                           x = 'logFC',
                           y = 'FDR',
                           title = 'Active vs. Control',
                           pCutoff = 0.05,
                           FCcutoff = 1,
                           pointSize = 1.5,
                           labSize = 3.0,
                           legendPosition = 'bottom')

# 2. Add the theme layer to remove grid lines and print the modified plot
v_plot1 + theme(
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank()
)


# --- Plot 2: Progressor vs. Latent (without grid lines) ---

# 1. Create the second volcano plot and store it
v_plot2 <- EnhancedVolcano(top_genes_progressor_with_names,
                           lab = top_genes_progressor_with_names$Gene_name,
                           x = 'logFC',
                           y = 'FDR',
                           title = 'Progressor vs. Latent',
                           pCutoff = 0.05,
                           FCcutoff = 1,
                           pointSize = 1.5,
                           labSize = 3.0,
                           legendPosition = 'bottom')

# 2. Add the theme layer and print the plot
v_plot2 + theme(
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank()
)

par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1) # Reset plotting layout


#-------------------------------------------------------------------------------
# STAGE 6: FUNCTIONAL ENRICHMENT ANALYSIS
#-------------------------------------------------------------------------------
# Goal: Identify biological pathways over-represented in the significant genes
#       from both comparisons.


# --- Best Practice: Define a reusable theme to remove grid lines ---
theme_no_grid <- theme(
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank()
)


#===============================================================================
# Part 1: ANALYSIS FOR 'ACTIVE VS. CONTROL'
#===============================================================================

# --- A. PREPARE GENE LIST (Active vs. Control) ---
deg_ids_active <- rownames(sig_genes_active)

# Convert Ensembl IDs to Entrez IDs
entrez_ids_active <- mapIds(org.Hs.eg.db,
                            keys = deg_ids_active,
                            column = "ENTREZID",
                            keytype = "ENSEMBL",
                            multiVals = "first")
entrez_ids_active <- entrez_ids_active[!is.na(entrez_ids_active)]

# --- B. RUN ENRICHMENT TESTS (Active vs. Control) ---
# GO enrichment
go_results_active <- enrichGO(gene = entrez_ids_active,
                              OrgDb = org.Hs.eg.db, ont = "BP",
                              pAdjustMethod = "BH", qvalueCutoff = 0.05,
                              readable = TRUE)
# KEGG enrichment
kegg_results_active <- enrichKEGG(gene = entrez_ids_active,
                                  organism = "hsa", pAdjustMethod = "BH",
                                  qvalueCutoff = 0.05)

# --- C. VISUALIZE & SAVE RESULTS (Active vs. Control) ---
# GO Plot
if (!is.null(go_results_active) && nrow(go_results_active) > 0) {
  go_plot_active <- dotplot(go_results_active, showCategory = 15) +
    ggtitle("GO Enrichment: Active vs. Control") +
    theme_no_grid
  print(go_plot_active)
} else {
  print("No significant GO terms found for Active vs. Control.")
}

# KEGG Plot
if (!is.null(kegg_results_active) && nrow(kegg_results_active) > 0) {
  kegg_plot_active <- dotplot(kegg_results_active, showCategory = 15) +
    ggtitle("KEGG Enrichment: Active vs. Control") +
    theme_no_grid
  print(kegg_plot_active)
} else {
  print("No significant KEGG pathways found for Active vs. Control.")
}

# Save results to files
write.csv(as.data.frame(go_results_active), "GO_Enrichment_Active_vs_Control.csv", row.names = FALSE)
write.csv(as.data.frame(kegg_results_active), "KEGG_Enrichment_Active_vs_Control.csv", row.names = FALSE)

cat("\n--- Completed analysis for Active vs. Control ---\n\n")


#===============================================================================
# Part 2: ANALYSIS FOR 'PROGRESSOR VS. LATENT'
#===============================================================================
# NOTE: This section assumes you have a dataframe named 'sig_genes_progressor'

# --- A. PREPARE GENE LIST (Progressor vs. Latent) ---
deg_ids_progressor <- rownames(sig_genes_progressor)

# Convert Ensembl IDs to Entrez IDs
entrez_ids_progressor <- mapIds(org.Hs.eg.db,
                                keys = deg_ids_progressor,
                                column = "ENTREZID",
                                keytype = "ENSEMBL",
                                multiVals = "first")
entrez_ids_progressor <- entrez_ids_progressor[!is.na(entrez_ids_progressor)]


# --- B. RUN ENRICHMENT TESTS (Progressor vs. Latent) ---
# GO enrichment
go_results_progressor <- enrichGO(gene = entrez_ids_progressor,
                                  OrgDb = org.Hs.eg.db, ont = "BP",
                                  pAdjustMethod = "BH", qvalueCutoff = 0.05,
                                  readable = TRUE)
# KEGG enrichment
kegg_results_progressor <- enrichKEGG(gene = entrez_ids_progressor,
                                      organism = "hsa", pAdjustMethod = "BH",
                                      qvalueCutoff = 0.05)


# --- C. VISUALIZE & SAVE RESULTS (Progressor vs. Latent) ---
# GO Plot
if (!is.null(go_results_progressor) && nrow(go_results_progressor) > 0) {
  go_plot_progressor <- dotplot(go_results_progressor, showCategory = 15) +
    ggtitle("GO Enrichment: Progressor vs. Latent") +
    theme_no_grid
  print(go_plot_progressor)
} else {
  print("No significant GO terms found for Progressor vs. Latent.")
}

# KEGG Plot
if (!is.null(kegg_results_progressor) && nrow(kegg_results_progressor) > 0) {
  kegg_plot_progressor <- dotplot(kegg_results_progressor, showCategory = 15) +
    ggtitle("KEGG Enrichment: Progressor vs. Latent") +
    theme_no_grid
  print(kegg_plot_progressor)
} else {
  print("No significant KEGG pathways found for Progressor vs. Latent.")
}

# Save results to files
write.csv(as.data.frame(go_results_progressor), "GO_Enrichment_Progressor_vs_Latent.csv", row.names = FALSE)
write.csv(as.data.frame(kegg_results_progressor), "KEGG_Enrichment_Progressor_vs_Latent.csv", row.names = FALSE)


cat("--- Functional Enrichment Analysis Complete ---")
