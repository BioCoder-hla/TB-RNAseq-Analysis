

# 🧬 TB-RNAseq-Analysis

## 📖 Overview
This repository contains the R code, results, and visualizations for a **comprehensive bioinformatics analysis** of RNA-seq data investigating tuberculosis (TB) progression using the Leicester cohort (**GEO: GSE107991**).  
The project explores **gene expression signatures** across four conditions — Active TB, Latent TB, Progressors, and Controls — to identify potential biomarkers and gain deeper insight into immune mechanisms driving disease progression.

---

## 🎯 Key Findings
- **Active TB vs. Control**: Identified 788 significant DEGs (|logFC| ≥ 1, FDR < 0.05), revealing robust up-regulation of interferon-stimulated genes (e.g., `STAT1`, `GBP1`, `GBP6`) and innate immune pathways.
- **Progressor vs. Latent TB**: Found 237 DEGs with similar immune-activation signatures, suggesting early molecular changes predictive of disease progression.
- **Functional Enrichment**: GO and KEGG analysis highlighted key processes, including neutrophil extracellular traps, complement cascade, and interferon signaling, which may serve as targets for further biomarker validation.

---

## 🧪 Analysis Workflow
1. **Data Preparation & QC**  
   - Download raw RNA-seq data and metadata.
   - Filter low-expression genes and normalize using edgeR TMM.

2. **Exploratory Analysis**  
   - Principal Component Analysis (PCA) to visualize group clustering.
   - Heatmaps for top most-variable genes.

3. **Differential Expression**  
   - Fit quasi-likelihood GLM to identify DEGs.
   - Generated MD and volcano plots for Active vs. Control and Progressor vs. Latent comparisons.

4. **Functional Enrichment**  
   - Gene Ontology (GO) and KEGG pathway enrichment using `clusterProfiler`.
   - Visualized top-enriched biological processes with dot plots.

---

## ⚙️ Getting Started

### 📋 Prerequisites:
- R (v4.2.0 or higher)  
- Packages: `edgeR` `clusterProfiler` `ggplot2` `pheatmap` `EnhancedVolcano` `GEOquery`  `kableExtra` `DESeq2` `org.Hs.eg.db` `readxl` `dplyr` `tibble` `GEOquery`

## 📦 Installation:

### From CRAN:
```bash
install.packages(c("tidyverse", "pheatmap", "EnhancedVolcano", "ggrepel"))
```

### From Bioconductor:
```bash
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("edgeR", "clusterProfiler", "org.Hs.eg.db", "GEOquery"))
```

---

## ▶️ Usage

1. Clone this repository:

   ```bash
   git clone https://github.com/YOUR_GITHUB_USERNAME/TB-RNAseq-Analysis.git
   ```
2. Set the working directory.
3. Run the script to reproduce all analyses and figures.

---

## 📜 Citation

Data derived from:  
**Singhania _et al_. (2019)** — *A modular transcriptional signature identifies phenotypic heterogeneity of human tuberculosis infection*  
[DOI: 10.1038/s41467-018-04579-w](https://doi.org/10.1038/s41467-018-04579-w)

---

##  👤 Contact

For questions, contact [hlatwayne@gmail.com] or open an issue on GitHub.

---

## 🤝 Contributing

Contributions, issues, and suggestions are welcome!
Please submit a pull request or open an issue if you’d like to help improve this project.

---

## 📝 License

Distributed under the MIT License — see `LICENSE` for details.

---

**#Bioinformatics #RNAseq #Tuberculosis #GeneExpression #DataAnalysis #RStats #edgeR #clusterProfiler**


---
