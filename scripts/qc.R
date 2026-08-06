library(DESeq2)
library(ggplot2)
library(dplyr)
library(vsn)
library(ComplexHeatmap)

qc_bg <- "#F7F9FC"

qc_theme <- theme_classic(base_size = 16) +
  theme(
    panel.background = element_rect(fill = qc_bg, color = NA),
    plot.background = element_rect(fill = qc_bg, color = NA),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold")
  )

# -----------------------------
# INPUTS
# -----------------------------
cancer <- snakemake@wildcards$cancer

counts <- readRDS(snakemake@input[["counts"]])
meta   <- readRDS(snakemake@input[["meta"]])

dir.create(sprintf("results/qc/%s", cancer), recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# OUTPUTS
# -----------------------------
qc_dds_file <- snakemake@output[["dds"]]
pca_plot    <- snakemake@output[["pca"]]
qc_report <- snakemake@output[["report"]]
sample_distance_heatmap <- snakemake@output[["sample_dist"]]
pc1_total   <- snakemake@output[["pc1_total"]]
size_factors <- snakemake@output[["size_factors"]]


dir.create(dirname(qc_dds_file), recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# REPORT
# -----------------------------
report_lines <- c(
  "=== QC REPORT ===",
  paste("Cancer:", cancer),
  paste("Samples:", ncol(counts)),
  paste("Genes:", nrow(counts))
)

# safe column checks
get_safe <- function(df, col) {
  if (col %in% colnames(df)) df[[col]] else NA
}

report_lines <- c(report_lines,
  paste("Tumor samples:",
        sum(get_safe(meta, "sample_type") == "Primary Tumor", na.rm = TRUE)),
  paste("Normal samples:",
        sum(get_safe(meta, "sample_type") == "Solid Tissue Normal", na.rm = TRUE)),
  paste("Median age:",
        median(as.numeric(get_safe(meta, "age_at_diagnosis")), na.rm = TRUE)),
paste("Detected Genes:", sum(rowSums(counts > 0) > 0)),
paste("Zero-count Genes:", sum(rowSums(counts > 0) == 0)),
paste("Top Gene Fraction:", apply(counts, 2, function(x) max(x) / sum(x, na.rm = TRUE)))
)

# -----------------------------
# QC METRICS
# -----------------------------
qc_df <- data.frame(
  sample = colnames(counts),
  total_reads = colSums(counts)
)

qc_df$log_total_reads <- log10(qc_df$total_reads + 1)

# -----------------------------
# DESEQ2 OBJECT
# -----------------------------
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = meta,
  design = ~ 1
)

dds <- estimateSizeFactors(dds)
vsd <- vst(dds, blind = TRUE)

saveRDS(dds, qc_dds_file)

# -----------------------------
# PCA
# -----------------------------
group <- if ("sample_type" %in% colnames(meta)) {
  "sample_type"
} else if ("batch" %in% colnames(meta)) {
  "batch"
} else {
  NULL
}

if (is.null(group)) {
    pca <- plotPCA(vsd, returnData=TRUE)
} else {
    pca <- plotPCA(vsd, intgroup=group, returnData=TRUE)
}

percentVar <- round(100 * attr(pca, "percentVar"), 1)

p <- ggplot(pca, aes(PC1, PC2)) +
  geom_point(aes(color = .data[[group]]), size = 3) +
  labs(
    title = paste(cancer, "PCA"),
    x = paste0("PC1: ", percentVar[1], "%"),
    y = paste0("PC2: ", percentVar[2], "%"),
    color = group
  ) +
  qc_theme

ggsave(pca_plot, p, width = 6, height = 5, bg = qc_bg)

# -----------------------------
# PC1 vs depth
# -----------------------------
pc1_df <- data.frame(
  PC1 = pca$PC1,
  total_reads = qc_df$total_reads
)

p1 <- ggplot(pc1_df, aes(PC1, total_reads)) +
  geom_point() +
  labs(title = "PC1 vs Total Reads") +
  qc_theme

ggsave(pc1_total, p1, bg = qc_bg)

# -----------------------------
# HEATMAP
# -----------------------------
sample_dists <- dist(t(assay(vsd)))

png(filename = sample_distance_heatmap, width = 800, height = 800)
ComplexHeatmap::Heatmap(as.matrix(sample_dists))
dev.off()

# -----------------------------
# Size factor distribution
# -----------------------------

png(filename=size_factors)

boxplot(sizeFactors(dds),
        ylab="Size factor",
        main="Size Factors")

dev.off()

# -----------------------------
# REPORT OUTPUT
# -----------------------------
writeLines(report_lines, qc_report)