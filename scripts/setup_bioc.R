options(repos = c(CRAN = "https://cloud.r-project.org"))

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages(
    "BiocManager",
    repos = "https://cloud.r-project.org"
  )
}

# Bioconductor packages
BiocManager::install(
  c(
    "TCGAbiolinks",
    "DESeq2",
    "SummarizedExperiment",
    "Biobase",
    "clusterProfiler",
    "org.Hs.eg.db",
    "apeglm",
    "ComplexHeatmap"
  ),
  ask = FALSE,
  update = FALSE
)

# CRAN packages
install.packages(
  c(
    "ggplot2",
    "readr",
    "dplyr",
    "ggrepel"
  ),
  repos = "https://cloud.r-project.org"
)

# Verify critical packages
required_pkgs <- c(
  "TCGAbiolinks",
  "DESeq2",
  "SummarizedExperiment",
  "clusterProfiler",
  "org.Hs.eg.db",
  "apeglm",
  "ComplexHeatmap",
  "ggplot2",
  "readr",
  "dplyr",
  "ggrepel"
)

missing <- required_pkgs[
  !sapply(required_pkgs, requireNamespace, quietly = TRUE)
]

if (length(missing) > 0) {
  stop(
    paste(
      "Failed to install:",
      paste(missing, collapse = ", ")
    )
  )
}

file.create(".bioc_installed")