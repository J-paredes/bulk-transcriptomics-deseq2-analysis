# 📌 Bulk RNA-seq Differential Expression Analysis (Colorectal Cancer)
This project performs reproducible bulk RNA-seq differential expression analysis using publicly available colorectal cancer patient data from GEO (GSE327166). The workflow is fully automated using Snakemake and includes preprocessing, statistical modeling, and visualization using DESeq2 in R.

## 🔬 Dataset
**GEO Accession:** GSE327166  
**Samples:** Colorectal cancer tumor vs matched normal tissue  
**Data type:** Bulk RNA-seq (processed counts matrix)  
**Source:** NIH Gene Expression Omnibus (GEO)  
**Contributors:** Severson DT, Freyaldenhoven S, Wadowski B, Richards WG, Bueno R
**GEO Link:** [Link](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE327166)
**Citation:** Severson DT, Freyaldenhoven S, Wadowski B, Hung YP et al. Multiomic, Histologic, and scRNA-seq Profiling of Pleural Mesothelioma Reveals Negative Prognosis Associated With a Novel Uncommitted Molecular Phenotype. J Thorac Oncol 2026 Apr;21(4):103529. [PMID: 41319863](https://www.jto.org/article/S1556-0864(25)02937-5/fulltext)

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


