options(repos = c(CRAN = "https://cloud.r-project.org"))

dir.create(
  "results/setup/enrichment",
  recursive = TRUE,
  showWarnings = FALSE
)

bioc_pkgs <- c(
  "clusterProfiler",
  "org.Hs.eg.db",
  "enrichplot",
  "DOSE"
)

cran_pkgs <- c(
  "shadowtext"
)

for(pkg in bioc_pkgs){
    if(!requireNamespace(pkg, quietly=TRUE)){
        BiocManager::install(
            pkg,
            ask=FALSE,
            update=FALSE
        )
    }
}

for(pkg in cran_pkgs){
    if(!requireNamespace(pkg, quietly=TRUE)){
        install.packages(
            pkg,
            repos="https://cloud.r-project.org"
        )
    }
}

file.create(
  "results/setup/enrichment/.packages_installed"
)