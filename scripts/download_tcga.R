options(repos = c(CRAN = "https://cloud.r-project.org"))

dir.create("data", recursive = TRUE, showWarnings = FALSE)

ensure_bioc_package <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager")
    }

    BiocManager::install(pkg, ask = FALSE, update = FALSE)

    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(sprintf("Failed to install required package: %s", pkg))
    }
  }
}

ensure_bioc_package("TCGAbiolinks")
ensure_bioc_package("SummarizedExperiment")

suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(SummarizedExperiment)
})

with_retry <- function(expr_fun, label, max_attempts = 6, initial_wait = 10) {
  last_err <- NULL

  for (attempt in seq_len(max_attempts)) {
    call_res <- tryCatch(
      list(ok = TRUE, value = expr_fun()),
      error = function(e) {
        last_err <<- e
        list(ok = FALSE, value = NULL)
      }
    )

    if (isTRUE(call_res$ok)) {
      if (attempt > 1) {
        message(sprintf("%s succeeded on attempt %d.", label, attempt))
      }
      return(call_res$value)
    }

    msg <- conditionMessage(last_err)
    transient <- grepl(
      "Could not resolve host|Could not resolve hostname|cannot open the connection|getURL\\(\\) failed|timeout",
      msg,
      ignore.case = TRUE
    )

    if (!transient || attempt == max_attempts) {
      stop(last_err)
    }

    wait_s <- min(120, initial_wait * (2 ^ (attempt - 1)))
    message(sprintf(
      "%s failed on attempt %d/%d with transient network error. Retrying in %ds...",
      label, attempt, max_attempts, wait_s
    ))
    Sys.sleep(wait_s)
  }
}


extract_batch <- function(meta) {
  batch_fields <- c(
    "batch",
    "center",
    "plate",
    "tissue_source_site",
    "tss_code",
    "portion",
    "vial",
    "analyte",
    "is_ffpe"
  )

  available_fields <- intersect(batch_fields, colnames(meta))
  available_fields <- available_fields[available_fields != "batch"]

  if (length(available_fields) == 0) {
    meta$batch <- NA_character_
    meta$batch_source <- NA_character_
    return(meta)
  }

  batch_values <- do.call(
    interaction,
    c(
      meta[available_fields],
      list(drop = TRUE, sep = "_")
    )
  )

  meta$batch <- as.character(batch_values)
  meta$batch_source <- paste(available_fields, collapse = ",")
  meta
}

out_dir <- "data/raw"
gdc_dir <- file.path(out_dir, "GDCdata")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(gdc_dir, recursive = TRUE, showWarnings = FALSE)

query <- GDCquery(
  project = c("TCGA-LUAD", "TCGA-COAD"),
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

with_retry(
  function() GDCdownload(query, method = "api", files.per.chunk = 20, directory = gdc_dir),
  "GDCdownload"
)

data <- with_retry(
  function() GDCprepare(query, directory = gdc_dir),
  "GDCprepare"
)

available_assays <- SummarizedExperiment::assayNames(data)


# Available assays in SummarizedExperiment : 
#   => unstranded
## TCGA standard pipeline output
## raw counts (integer-like)
## compatible with DESeq2 assumptions
## most commonly used in publications
#   => stranded_first
## TCGA standard pipeline output
## raw counts (integer-like)
#   => stranded_second

#   => tpm_unstrand
## TPM and FPKM are normalized expression values that account for gene length and sequencing depth, but they are not raw counts and may not be suitable for DESeq2 analysis.
#   => fpkm_unstrand
## TPM and FPKM are normalized expression values that account for gene length and sequencing depth, but they are not raw counts and may not be suitable for DESeq2 analysis.
#   => fpkm_uq_unstrand
## upper-quartile normalized FPKM values (not raw counts)
## not suitable for DESeq2

assay_to_use <- if ("unstranded" %in% available_assays) {
  "unstranded"
} else {
  available_assays[[1]]
}

counts <- SummarizedExperiment::assay(data, assay_to_use)

meta <- as.data.frame(SummarizedExperiment::colData(data))
meta <- extract_batch(meta)

meta$cancer_type <- ifelse(meta$project_id == "TCGA-LUAD",
                           "LUAD", "COAD")

dir.create("data/raw", recursive = TRUE, showWarnings = FALSE)
dir.create("data/metadata", recursive = TRUE, showWarnings = FALSE)

saveRDS(counts, "data/raw/tcga_counts.rds")
saveRDS(meta, "data/metadata/tcga_metadata.rds")

library(readr)
write_csv(
  as.data.frame(counts, check.names = FALSE),
  "data/raw/tcga_counts.csv"
)

write_csv(
  meta,
  "data/metadata/tcga_metadata.csv"
)
