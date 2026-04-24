# bulk-transcriptomics-deseq2-analysis
Differential gene expression analysis of bulk RNA-seq.

Data using DESeq2, derived from pleural mesothelioma samples described in Severson et al., Multiomic, Histologic, and scRNA-seq Profiling of Pleural Mesothelioma Reveals a Novel Uncommitted Molecular Phenotype (Journal of Thoracic Oncology, 2025).

## 📌 Bulk RNA-seq Differential Expression Analysis (Colorectal Cancer)
This project performs reproducible bulk RNA-seq differential expression analysis using publicly available colorectal cancer patient data from GEO (GSE327166). The workflow is fully automated using Snakemake and includes preprocessing, statistical modeling, and visualization using DESeq2 in R.

## 🔬 Dataset
**GEO Accession:** GSE327166  
**Samples:** Colorectal cancer tumor vs matched normal tissue  
**Data type:** Bulk RNA-seq (processed counts matrix)  
**Source:** NIH Gene Expression Omnibus (GEO)  

## ⚙️ Workflow Overview
The analysis pipeline includes:
- Input raw count matrix + metadata
- Differential expression analysis (DESeq2)
- Variance stabilizing transformation (VST)
- Exploratory data analysis (PCA)
- Differential expression visualization:
- Volcano plot
- Heatmap of top genes

## 🧪 Tools & Technologies
Snakemake (workflow management)  
R (DESeq2, ggplot2, pheatmap)  
Conda environments  
GEO datasets  
Pandas-style metadata preprocessing  

## 🔁 Workflow Diagram
![Workflow Diagram](workflow/bulk-rnaseq.drawio.svg)

## 📊 Results
Principal Component Analysis (PCA)
Volcano Plot
Heatmap (Top Differential Genes)

## 📈 Key Findings
Distinct clustering of tumor vs normal samples observed in PCA
Identification of significantly differentially expressed genes
Clear transcriptional separation between disease states


