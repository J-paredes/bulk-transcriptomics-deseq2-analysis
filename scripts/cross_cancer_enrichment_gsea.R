library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(readr)
library(ggplot2)
library(dplyr)


res <- read_csv(
    snakemake@input[["results"]],
    show_col_types = FALSE
)


# Remove missing values
gsea_df <- res %>%
    filter(
        !is.na(log2FoldChange),
        !is.na(gene_id),
        !is.na(stat)
    )


# Remove Ensembl versions
gsea_df$gene_id <- gsub(
    "\\..*",
    "",
    gsea_df$gene_id
)


gene_rank <- gsea_df$stat

# Safety check
gene_rank <- gene_rank[is.finite(gene_rank)]


names(gene_rank) <- gsea_df$gene_id

gene_rank <- sort(
    gene_rank,
    decreasing = TRUE
)

# GO GSEA
ego <- gseGO(
    geneList = gene_rank,
    OrgDb = org.Hs.eg.db,
    keyType = "ENSEMBL",
    ont = "BP",
    minGSSize = 10,
    maxGSSize = 500,
    pAdjustMethod = "BH"
)


write_csv(
    as.data.frame(ego),
    snakemake@output[["go_csv"]]
)


ggsave(
    snakemake@output[["go_plot"]],
    dotplot(ego, showCategory=15),
    width=8,
    height=6
)


# KEGG GSEA
gene_df <- bitr(
    names(gene_rank),
    fromType="ENSEMBL",
    toType="ENTREZID",
    OrgDb=org.Hs.eg.db
)

gene_rank_kegg <- gene_rank[
    names(gene_rank) %in% gene_df$ENSEMBL
]

# Convert Ensembl names to Entrez
names(gene_rank_kegg) <- gene_df$ENTREZID[
    match(names(gene_rank_kegg), gene_df$ENSEMBL)
]

gene_rank_kegg <- sort(
    gene_rank_kegg,
    decreasing = TRUE
)

gene_rank_kegg <- gene_rank_kegg[
    !duplicated(names(gene_rank_kegg))
]

ekk <- gseKEGG(
    geneList=gene_rank_kegg,
    organism="hsa"
)

write_csv(
    as.data.frame(ekk),
    snakemake@output[["kegg_csv"]]
)


ggsave(
    snakemake@output[["kegg_plot"]],
    dotplot(ekk, showCategory=15),
    width=8,
    height=6
)