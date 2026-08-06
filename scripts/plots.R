
library(DESeq2)
library(ggplot2)
library(readr)
library(dplyr)
library(ComplexHeatmap)
library(ggrepel)

volcano_plot <- snakemake@output[["volcano"]]
pca_plot <- snakemake@output[["pca"]]
heatmap_plot <- snakemake@output[["heatmap"]]
ma_plot <- snakemake@output[["ma_plot"]]


res <- readr::read_csv(snakemake@input[["results"]], show_col_types = FALSE)
dds <- readRDS(snakemake@input[["dds"]])

# PCA plot ----------------------
vsd <- vst(dds)
pcaData_sample_type <- as.data.frame(
  plotPCA(vsd, intgroup = c("sample_type"), returnData = TRUE)
)
write_csv(pcaData_sample_type, file = "results//pcaData_sample_type.csv")

ggplot(pcaData_sample_type, aes(PC1, PC2, color=sample_type)) +
  geom_point(aes(shape=sample_type), size=10, alpha=0.7) +
  labs(title="PCA by Sample Type", x="PC1 (%)", y="PC2 (%)") +
  stat_ellipse(level=0.95, linetype="dashed", alpha=0.8) +
  theme_classic(base_size=16) +
  # scale_color_manual(values=c("Female" = "pink", "Male" = "blue")) +
  # scale_shape_manual(values=c(16, 17)) +
  theme(legend.title=element_blank(),
        axis.title=element_text(face="bold"))
        
ggsave(pca_plot, width=8, height=6)

# Volcano plot ----------------------
res_df <- as.data.frame(res)
res_df$gene <- res_df$gene_id
res_df$significant <- res_df$padj < 0.05
res_df$label <- ifelse(res_df$padj < 0.01 & abs(res_df$log2FoldChange) > 2,
                        res_df$gene, "")


res_df$significance_label <- dplyr::case_when(
  res_df$padj < 0.01 & abs(res_df$log2FoldChange) > 2 ~ "Significant at q < 0.01 and |Threshold| > 2",
  res_df$padj < 0.05 ~ "Moderately Significant at q < 0.05",
  res_df$log2FoldChange > 0 & res_df$padj > 0.05~ "Upregulated [ns]",
  res_df$log2FoldChange < 0 & res_df$padj > 0.50~ "Downregulated [ns]",
  TRUE ~ "Not Significant"
) 
res_df$significance_label <- factor(res_df$significance_label, levels=c(
  "Significant at q < 0.01 and |Threshold| > 2",
  "Moderately Significant at q < 0.05",
  "Upregulated [ns]",
  "Downregulated [ns]",
  "Not Significant"
))

volcano <- ggplot(res_df, aes(log2FoldChange, -log10(pvalue), color=significance_label)) +
  geom_point(alpha=0.3, size=5) + 
  labs(x=expression(bold(log2(FC))), y=expression(bold(-log10(q))),
  color="Significance") +
  scale_color_manual(values=c("Significant at q < 0.01 and |Threshold| > 2" = "goldenrod",
                      "Moderately Significant at q < 0.05" = "yellow",
                      "Upregulated [ns]" = "darkgreen",
                      "Downregulated [ns]" = "red",
                      "Not Significant" = "grey")) +
geom_hline(yintercept=-log10(0.05), linetype="dashed", color="black") +
  geom_vline(xintercept=c(-2, 2), linetype="dashed", color="black") +
  theme_classic(base_size=16) +  
  theme(legend.title=element_text(face="bold"),
        axis.title=element_text(face="bold")) + 
  geom_text_repel(aes(label=label), color="black", size=3, max.overlaps=10,
   box.padding=0.3, point.padding=0.5, segment.color="grey50")

ggsave(volcano_plot, volcano, width=10, height=6)

# Heatmap ----------------------
top_genes <- res |>
  dplyr::filter(!is.na(padj), !is.na(gene_id)) |>
  dplyr::arrange(padj) |>
  dplyr::slice_head(n = 20) |>
  dplyr::pull(gene_id)

top_genes <- intersect(top_genes, rownames(assay(vsd)))
mat <- assay(vsd)[top_genes, , drop = FALSE]

anno <- as.data.frame(colData(vsd))[, "sample_type", drop = FALSE]
anno <- anno[colnames(mat), , drop = FALSE]

library(ComplexHeatmap)

mat_scaled <- t(scale(t(mat)))

ht <- Heatmap(
mat_scaled,
name = "Z-score",
show_row_names = FALSE,
show_column_names = FALSE,
top_annotation = HeatmapAnnotation(df = anno)
)

png(heatmap_plot, width = 800, height = 600)
draw(ht, merge_legends = TRUE)
dev.off()


# MA plot ----------------------
resLFC <- lfcShrink(dds, coef=2, type="apeglm")
resLFC_df <- as.data.frame(resLFC)

# Customizable MA plot settings
ma_padj_cutoff <- 0.05
ma_lfc_cutoff <- 1
ma_point_size <- 1.8
ma_alpha <- 0.7
ma_label_top_n <- 10

resLFC_df$gene <- rownames(resLFC_df)
resLFC_df <- subset(
  resLFC_df,
  !is.na(baseMean) & baseMean > 0 & !is.na(log2FoldChange)
)

resLFC_df$significance <- ifelse(
  !is.na(resLFC_df$padj) & resLFC_df$padj < ma_padj_cutoff & abs(resLFC_df$log2FoldChange) >= ma_lfc_cutoff,
  "Significant",
  "Not significant"
)

sig_for_labels <- resLFC_df[
  resLFC_df$significance == "Significant",
]
sig_for_labels <- sig_for_labels[order(sig_for_labels$padj), ]
label_genes <- head(sig_for_labels$gene, ma_label_top_n)
resLFC_df$label <- ifelse(resLFC_df$gene %in% label_genes, resLFC_df$gene, "")

ma <- ggplot(resLFC_df, aes(x = log10(baseMean), y = log2FoldChange, color = significance)) +
  geom_hline(yintercept = c(-ma_lfc_cutoff, ma_lfc_cutoff), linetype = "dashed", color = "grey40") +
  geom_point(size = ma_point_size, alpha = ma_alpha) +
  geom_text_repel(aes(label = label), color = "black", size = 3, max.overlaps = Inf) +
  # scale_color_manual(values = c("Not significant" = "grey70", "Significant" = "red")) +
  theme_classic(base_size = 16) +
  labs(
    title = "MA Plot",
    subtitle = paste0("Significant: padj < ", ma_padj_cutoff, " and |log2FC| >= ", ma_lfc_cutoff),
    x = expression(log[10](baseMean)),
    y = expression(log[2](Fold~Change)),
    color = "Gene set"
  ) +
  theme(
    legend.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold")
  )

ggsave(snakemake@output[["ma_plot"]], ma, width = 8, height = 6, dpi = 300)
