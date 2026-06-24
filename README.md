# Bulk RNA-seq Analysis Pipeline Using TCGA Data

## Overview

This repository contains a reproducible Snakemake workflow for downloading, processing, and analyzing bulk RNA-seq data from The Cancer Genome Atlas (TCGA). The pipeline performs:

1. TCGA RNA-seq data download using TCGAbiolinks
2. Quality control (QC) and exploratory analysis
3. Differential expression analysis using DESeq2
4. Functional enrichment analysis
5. Publication-ready visualizations

The workflow is designed to be reproducible through Snakemake and Conda environments.

## Data

LUAD molecular subtype annotations were obtained from:
The Cancer Genome Atlas Research Network. [Comprehensive molecular profiling of lung adenocarcinoma.](doi.org/10.1038/nature13385) Nature. 2014. doi:10.1038/nature13385

COAD molecular subtype annotations were obtained from:
The Cancer Genome Atlas Network. [Comprehensive molecular characterization of human colon and rectal cancer.](doi.org/10.1038/nature11252) Nature. 2012. doi:10.1038/nature11252

---

## Workflow

### Step 1: Bioconductor Setup

Installs required Bioconductor and CRAN packages:

- TCGAbiolinks
- DESeq2
- SummarizedExperiment
- clusterProfiler
- org.Hs.eg.db
- ComplexHeatmap
- apeglm
- ggplot2
- readr
- dplyr

Output:

```text
.bioc_installed
```

---

### Step 2: Download TCGA Data

Downloads RNA-seq count data and metadata using TCGAbiolinks.

Projects are specified in:

```text
workflow/env/config.yaml
```

Current implementation supports:

- TCGA-LUAD (Lung Adenocarcinoma)
- TCGA-COAD (Colon Adenocarcinoma)

Outputs:

```text
data/raw/tcga_counts.rds
data/raw/tcga_counts.csv
data/metadata/tcga_metadata.rds
data/metadata/tcga_metadata.csv
```

---

### Step 3: Quality Control

Performs initial quality assessment including:

- DESeq2 object construction
- Sample filtering
- Principal Component Analysis (PCA)
- Association of PC1 with sequencing metrics
- QC reporting

Outputs:

```text
results/qc/qc_deseq2_dds.rds
results/qc/qc_report.txt
results/qc/pca_prelim.png
results/qc/pc1_total_reads.png
results/qc/pc1_multimapped.png
```

---

### Step 4: Differential Expression Analysis

Performs differential gene expression analysis using DESeq2.

Outputs:

```text
results/des/deseq2_results.csv
results/des/normalized_counts.csv
results/des/deseq2_dds.rds
```

---

### Step 5: Functional Enrichment Analysis

Performs Gene Ontology (GO) and KEGG pathway enrichment analysis using clusterProfiler.

Outputs:

```text
results/enrichment/go_enrichment.png
results/enrichment/KEGG_enrichment.png
results/enrichment/go_results.csv
results/enrichment/kegg_results.csv
```

---

### Step 6: Visualization

Generates publication-ready figures including:

- PCA plots
- Volcano plots
- Heatmaps

Outputs:

```text
results/figures/pca.png
results/figures/volcano.png
results/figures/heatmap.png
```

---

## Directory Structure

```text
.
├── workflow/
│   ├── Snakefile
│   ├── env/
│   │   ├── config.yaml
│   │   └── r.yaml
│   └── scripts/
│       ├── setup_bioc.R
│       ├── download_tcga.R
│       ├── qc.R
│       ├── deseq2.R
│       ├── enrichment.R
│       └── plots.R
│
├── data/
│   ├── raw/
│   └── metadata/
│
├── results/
│   ├── qc/
│   ├── des/
│   ├── enrichment/
│   └── figures/
│
└── README.md
```

---

## Installation

### Clone Repository

```bash
git clone <repository-url>
cd bulk-rna-seq
```

### Create Snakemake Environment

```bash
conda create -n bulk-rna-seq snakemake -c conda-forge
conda activate bulk-rna-seq
```

---

## Running the Pipeline

Run the complete workflow:

```bash
snakemake --use-conda --cores 1
```

Run a specific step:

```bash
snakemake --use-conda --cores 1 qc
```

Generate a workflow diagram:

```bash
snakemake --dag | dot -Tpng > dag.png
```

---

## Configuration

Projects analyzed by the pipeline are specified in:

```yaml
projects:
  - TCGA-LUAD
  - TCGA-COAD
```

located in:

```text
workflow/env/config.yaml
```

---

## Methods

### Differential Expression

Differential expression is performed using DESeq2 on raw STAR-generated count matrices downloaded from TCGA.

### Batch Assessment

Metadata fields associated with sequencing center, tissue source site, plate information, and other technical variables are extracted and evaluated for potential batch effects during QC analysis.

### Enrichment Analysis

Gene Ontology and KEGG pathway enrichment analyses are performed using clusterProfiler and organism annotations from org.Hs.eg.db.

---

## Citation

If you use this workflow, please cite:

- Love MI, Huber W, Anders S. DESeq2.
- Colaprico et al. TCGAbiolinks.
- Wu et al. clusterProfiler.
- The Cancer Genome Atlas (TCGA).
