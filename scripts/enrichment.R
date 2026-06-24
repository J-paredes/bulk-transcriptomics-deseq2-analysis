

library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(readr)
library(dplyr)

go_plot <- snakemake@output[[1]]
kegg_plot <- snakemake@output[[2]]
go_csv <- snakemake@output[[3]]
kegg_csv <- snakemake@output[[4]]

for (output_path in as.character(snakemake@output)) {
  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
}

res <- readr::read_csv(snakemake@input[[1]], show_col_types = FALSE)
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

save_empty_result <- function(plot_file, csv_file, title_text) {
  empty_plot <- ggplot() +
    annotate("text", x = 0, y = 0, label = title_text, size = 5) +
    xlim(-1, 1) + ylim(-1, 1) +
    theme_void()
  ggsave(plot_file, empty_plot, width = 8, height = 6)
  readr::write_csv(data.frame(), csv_file)
}

# GO Enrichment ---------------------
sig <- res[which(res$padj < 0.05), ]
genes <- unique(gsub("\\..*", "", sig$gene_id))
genes <- genes[!is.na(genes) & genes != ""]

ego <- tryCatch(
  enrichGO(gene          = genes,
                OrgDb         = org.Hs.eg.db,
                keyType       = "ENSEMBL",
                ont           = "BP",
                pAdjustMethod = "BH",
                readable      = TRUE),
  error = function(e) NULL
)

if (is.null(ego) || nrow(as.data.frame(ego)) == 0) {
  save_empty_result(go_plot, go_csv, "No GO terms could be mapped")
} else {
  dotplot(ego, showCategory = 15) + 
    labs(title = "GO Enrichment of Significant DE Genes") +
    theme_classic(base_size = 16) +
    theme(legend.title = element_text(face = "bold"),
          axis.title = element_text(face = "bold")) 
  ggsave(go_plot)
  readr::write_csv(as.data.frame(ego), go_csv)
}
# KEGG Pathway ---------------------
gene_df <- bitr(genes,
                fromType = "ENSEMBL",   # or "SYMBOL"
                toType = "ENTREZID",
                OrgDb = org.Hs.eg.db)

if (is.null(gene_df) || nrow(gene_df) == 0) {
  save_empty_result(kegg_plot, kegg_csv, "No KEGG terms could be mapped")
} else {
  gene_df <- gene_df[!is.na(gene_df$ENTREZID), ]
  ekegg <- tryCatch(
    enrichKEGG(gene = gene_df$ENTREZID, organism = "hsa"),
    error = function(e) NULL
  )

  if (is.null(ekegg) || nrow(as.data.frame(ekegg)) == 0) {
    save_empty_result(kegg_plot, kegg_csv, "No KEGG terms were enriched")
  } else {
    dotplot(ekegg) + 
      labs(title = "KEGG Enrichment of Significant DE Genes") +
      theme_classic(base_size = 16) +
      theme(legend.title = element_text(face = "bold"),
            axis.title = element_text(face = "bold")) 
    ggsave(kegg_plot)
    readr::write_csv(as.data.frame(ekegg), kegg_csv)
  }
}
