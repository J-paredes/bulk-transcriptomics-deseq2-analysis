library(DESeq2)
library(ggplot2)
library(dplyr)
library(ComplexHeatmap)

qc_bg <- "#F7F9FC"

qc_theme <- theme_classic(base_size = 14) +
  theme(
    panel.background = element_rect(fill = qc_bg, color = NA),
    plot.background = element_rect(fill = qc_bg, color = NA),
    legend.background = element_rect(fill = qc_bg, color = NA),
    strip.background = element_rect(fill = qc_bg, color = NA)
  )

# -----------------------------
# INPUTS
# -----------------------------
luad_counts <- readRDS(snakemake@input[["luad_counts"]])
coad_counts <- readRDS(snakemake@input[["coad_counts"]])

luad_meta <- readRDS(snakemake@input[["luad_meta"]])
coad_meta <- readRDS(snakemake@input[["coad_meta"]])

# -----------------------------
# ADD LABELS
# -----------------------------
luad_meta$cancer_type <- "LUAD"
coad_meta$cancer_type <- "COAD"

# -----------------------------
# ALIGN FEATURES (genes)
# -----------------------------
common_genes <- intersect(rownames(luad_counts), rownames(coad_counts))

luad_counts <- luad_counts[common_genes, , drop = FALSE]
coad_counts <- coad_counts[common_genes, , drop = FALSE]

# -----------------------------
# COMBINE MATRICES
# -----------------------------
counts <- cbind(luad_counts, coad_counts)

meta <- bind_rows(luad_meta, coad_meta)

# ensure rownames match sample names
rownames(meta) <- colnames(counts)

# -----------------------------
# SANITY CHECK
# -----------------------------
stopifnot(all(colnames(counts) == rownames(meta)))

# -----------------------------
# DESEQ2 OBJECT
# -----------------------------
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = meta,
  design = ~ cancer_type
)

dds <- estimateSizeFactors(dds)

vsd <- vst(dds, blind = TRUE)

saveRDS(dds, snakemake@output[[1]])

# -----------------------------
# PCA
# -----------------------------
pca <- plotPCA(vsd, intgroup = "cancer_type", returnData = TRUE)

percentVar <- round(100 * attr(pca, "percentVar"))

p_pca <- ggplot(pca, aes(PC1, PC2)) +
  geom_point(aes(color = cancer_type), size = 3, alpha = 0.8) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  labs(title = "Cross-cancer PCA (LUAD vs COAD)") +
  qc_theme

ggsave(
  filename = snakemake@output[[1]],
  plot = p_pca,
  width = 7,
  height = 5,
  bg = qc_bg
)

# -----------------------------
# SAMPLE DISTANCES
# -----------------------------
sample_dists <- dist(t(assay(vsd)))
mat <- as.matrix(sample_dists)

annotation <- data.frame(cancer_type = meta$cancer_type)
rownames(annotation) <- rownames(meta)

png(snakemake@output[[2]], width = 900, height = 900)

ComplexHeatmap::Heatmap(
  mat,
  name = "distance",
  top_annotation = HeatmapAnnotation(df = annotation),
  show_row_names = FALSE,
  show_column_names = FALSE
)

dev.off()

# -----------------------------
# QC METRICS SUMMARY (optional report hook)
# -----------------------------
qc_report <- c(
  "=== CROSS-CANCER QC ===",
  paste("Samples total:", ncol(counts)),
  paste("LUAD samples:", sum(meta$cancer_type == "LUAD")),
  paste("COAD samples:", sum(meta$cancer_type == "COAD")),
  paste("Genes used:", nrow(counts))
)

writeLines(qc_report, snakemake@output[[3]])

