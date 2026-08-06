options(repos = c(CRAN = "https://cloud.r-project.org"))

if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

print(R.version.string)
library(BiocManager)

pkgs <- c(
    "TCGAbiolinks",
    "DESeq2",
    "SummarizedExperiment",
    "Biobase",
        "vsn",
    # "clusterProfiler",
    # "org.Hs.eg.db",
    "apeglm",
    "ComplexHeatmap"
)

for (pkg in pkgs) {

    if (!requireNamespace(pkg, quietly = TRUE)) {

        message("Installing ", pkg)

        BiocManager::install(
            pkg,
            ask = FALSE,
            update = FALSE,
            force = TRUE
        )
    }
}


missing <- pkgs[
    !sapply(pkgs, requireNamespace, quietly = TRUE)
]

if (length(missing) > 0) {
    stop(
        "Missing packages: ",
        paste(missing, collapse = ", ")
    )
}


cat("LIBPATHS\n")
print(.libPaths())

cat("R HOME\n")
print(R.home())

cat("DESeq2\n")
print(find.package("DESeq2"))

message(.libPaths())
installed.packages()[, "Package"]
print(.libPaths())

print(find.package("DESeq2"))

print(requireNamespace("TCGAbiolinks", quietly = TRUE))

dir.create("results/setup", recursive = TRUE, showWarnings = FALSE)
write.csv(installed.packages(),
          "results/setup/installed_packages.csv")

writeLines(
    capture.output(sessionInfo()),
    "results/setup/sessionInfo.txt"
)

stopifnot(requireNamespace("DESeq2", quietly = TRUE))
stopifnot(requireNamespace("TCGAbiolinks", quietly = TRUE))
stopifnot(requireNamespace("ComplexHeatmap", quietly = TRUE))

file.create("results/setup/.bioc_installed")