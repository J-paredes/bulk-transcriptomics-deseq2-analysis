library(DESeq2)
library(ggplot2)
library(dplyr)

qc_bg <- "#F7F9FC"
qc_theme <- theme_classic(base_size = 16) +
    theme(
        panel.background = element_rect(fill = qc_bg, color = NA),
        plot.background = element_rect(fill = qc_bg, color = NA),
        legend.background = element_rect(fill = qc_bg, color = NA),
        strip.background = element_rect(fill = qc_bg, color = NA),
        axis.title = element_text(face = "bold"),
        plot.title = element_text(face = "bold")
    )

print("Starting QC script...")
# -----------------------------
# INPUTS (Snakemake style)
# -----------------------------
counts <- readRDS(snakemake@input[["counts"]])
meta   <- readRDS(snakemake@input[["meta"]])

# -----------------------------
# OUTPUT FILES
# -----------------------------
qc_dds_file <- snakemake@output[[1]]
qc_report <- snakemake@output[[2]]
pca_plot <- snakemake@output[[3]]
pc1_total_plot <- snakemake@output[[4]]
pc1_multimapped_plot <- snakemake@output[[5]]

for (output_path in as.character(snakemake@output)) {
    dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
}

print(qc_report)
print(pca_plot)

# -----------------------------
# 1. BASIC SANITY CHECKS
# -----------------------------
report_lines <- c()

report_lines <- c(report_lines, "=== QC REPORT ===")

# check dimensions
report_lines <- c(report_lines,
                   paste("Samples:", ncol(counts)),
                   paste("Genes:", nrow(counts)))

# check sample matching
common <- intersect(colnames(counts), rownames(meta))
report_lines <- c(report_lines,
                   paste("Matched samples:", length(common)))

if (length(common) == 0) stop("No overlap between counts and metadata!")

counts <- counts[, common]
meta <- meta[common, , drop = FALSE]

# -----------------------------
# 2. CHECK IF AGE EXISTS
# -----------------------------
if ("age_at_diagnosis" %in% colnames(meta)) {
    report_lines <- c(report_lines, "Age column FOUND: age_at_diagnosis")
    age_available <- TRUE
} else if ("age" %in% colnames(meta)) {
    report_lines <- c(report_lines, "Age column FOUND: age")
    age_available <- TRUE
} else {
    report_lines <- c(report_lines, "WARNING: No age column found")
    age_available <- FALSE
}

# -----------------------------
# 3. QC METRICS EXTRACTION
# -----------------------------
qc_df <- meta %>%
    mutate(
        total_reads = colSums(counts),
        log_total_reads = log10(total_reads + 1)
    )

# check if multimapped exists
if ("multimapped_reads" %in% colnames(meta)) {
    qc_df$multimapped_reads <- meta$multimapped_reads
} else {
    qc_df$multimapped_reads <- NA
    report_lines <- c(report_lines, "WARNING: multimapped_reads not found")
}

# batch check
if ("batch" %in% colnames(meta)) {
    report_lines <- c(report_lines, "Batch column FOUND")
    if ("batch_source" %in% colnames(meta)) {
        report_lines <- c(report_lines,
                          paste("Batch source:", unique(na.omit(meta$batch_source))))
    }
} else {
    report_lines <- c(report_lines, "WARNING: No batch column found")
}

# -----------------------------
# 4. DESEQ2 PCA (CORE QC STEP)
# -----------------------------
dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData = meta,
    design = ~ 1
)

saveRDS(dds, file = qc_dds_file)


vsd <- vst(dds, blind = TRUE)
pca <- plotPCA(vsd, intgroup = c("cancer_type"), returnData = TRUE)

percentVar <- round(100 * attr(pca, "percentVar"))

if ("batch" %in% colnames(meta)) {
    pca$batch <- meta$batch
    p <- ggplot(pca, aes(PC1, PC2)) +
        geom_point(size = 3, aes(color = batch)) +
        xlab(paste0("PC1: ", percentVar[1], "% variance")) +
        ylab(paste0("PC2: ", percentVar[2], "% variance")) +
        qc_theme +
        facet_wrap(~ meta$cancer_type)

    ggsave(pca_plot, p, width = 6, height = 5, bg = qc_bg)
} else {
    p <- ggplot(pca, aes(PC1, PC2)) +
        geom_point(size = 3) +
        xlab(paste0("PC1: ", percentVar[1], "% variance")) +
        ylab(paste0("PC2: ", percentVar[2], "% variance")) +
        qc_theme +
        facet_wrap(~ meta$cancer_type)
    ggsave(pca_plot, p, width = 6, height = 5, bg = qc_bg)
}

# -----------------------------
# 5. PC1 vs QC METRICS PLOTS
# -----------------------------

pc1_df <- data.frame(
    PC1 = pca$PC1,
    total_reads = qc_df$total_reads,
    multimapped_reads = qc_df$multimapped_reads
)

p1 <- ggplot(pc1_df, aes(PC1, total_reads)) +
    geom_point(color = "#2C7FB8") +
    labs(title = "PC1 vs Total Reads", x = "PC1", y = "Total Reads") +
    qc_theme

p2 <- ggplot(pc1_df, aes(PC1, multimapped_reads)) +
    geom_point(color = "#D95F0E") +
    labs(title = "PC1 vs Multimapped Reads", x = "PC1", y = "Multimapped Reads") +
    qc_theme

ggsave(pc1_total_plot, p1, bg = qc_bg)
ggsave(pc1_multimapped_plot, p2, bg = qc_bg)

# -----------------------------
# 6. WRITE REPORT
# -----------------------------
writeLines(report_lines, qc_report)