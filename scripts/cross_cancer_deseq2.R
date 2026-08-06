library(DESeq2)


# Create output subdirectories before writing files.
for (output_path in as.character(snakemake@output)) {
  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
}

luad_counts <- readRDS(snakemake@input[["luad_counts"]])
coad_counts <- readRDS(snakemake@input[["coad_counts"]])

luad_meta <- readRDS(snakemake@input[["luad_meta"]])
coad_meta <- readRDS(snakemake@input[["coad_meta"]])

# Add cancer labels
luad_meta$cancer_type <- "LUAD"
coad_meta$cancer_type <- "COAD"

rownames(luad_meta) <- colnames(luad_counts)
rownames(coad_meta) <- colnames(coad_counts)

# Keep only columns needed for DESeq2
luad_meta <- luad_meta[, "cancer_type", drop = FALSE]
coad_meta <- coad_meta[, "cancer_type", drop = FALSE]

# Preserve sample IDs
rownames(luad_meta) <- rownames(readRDS(snakemake@input[["luad_meta"]]))
rownames(coad_meta) <- rownames(readRDS(snakemake@input[["coad_meta"]]))

# Combine metadata
meta <- rbind(luad_meta, coad_meta)

# Make sure order matches count matrix
meta <- meta[colnames(counts), , drop = FALSE]

# Keep only shared genes
common_genes <- intersect(
  rownames(luad_counts),
  rownames(coad_counts)
)

luad_counts <- luad_counts[common_genes, , drop = FALSE]
coad_counts <- coad_counts[common_genes, , drop = FALSE]

# Combine counts
counts <- cbind(luad_counts, coad_counts)

# Combine metadata
meta <- rbind(luad_meta, coad_meta)

# Ensure metadata order matches counts
meta <- meta[colnames(counts), , drop = FALSE]

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

res <- results(dds, contrast=c("cancer_type", "LUAD", "COAD"))
res_df <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)
res_df <- res_df[, c("gene_id", setdiff(colnames(res_df), "gene_id"))]
write.csv(res_df,
          file = snakemake@output[["results"]],
          row.names = FALSE)

