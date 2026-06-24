
library(DESeq2)


# Create output subdirectories before writing files.
for (output_path in as.character(snakemake@output)) {
  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
}

counts <- readRDS(snakemake@input[["counts"]])
meta <- readRDS(snakemake@input[["meta"]])

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = meta,
  design = ~ cancer_type
)

dds <- DESeq(dds)
saveRDS(dds, file = snakemake@output[["dds"]])

normalized_counts <- counts(dds, normalized = TRUE)
normalized_counts_df <- as.data.frame(normalized_counts)
normalized_counts_df$gene_id <- rownames(normalized_counts_df)
normalized_counts_df <- normalized_counts_df[, c("gene_id", colnames(normalized_counts))]
write.csv(normalized_counts_df,
          file = snakemake@output[["normalized"]],
          row.names = FALSE)

res <- results(dds)
res_df <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)
res_df <- res_df[, c("gene_id", setdiff(colnames(res_df), "gene_id"))]
write.csv(res_df,
          file = snakemake@output[["results"]],
          row.names = FALSE)