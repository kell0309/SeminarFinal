#  Setup script: install + load required packages
#  Purpose:
#   *Ensures required CRAN and Bioconductor packages exist
#   *Loads libraries used in downstream analysis


#  CRAN packages download
#  download openxlsx which is used to write out excel files
install.packages("openxlsx")


#  Bioconductor packages
#  Download biomanager
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("edgeR", "clusterProfiler", "enrichplot", "org.Hs.eg.db"))


#  edgeR: differential expression analysis for count data
#  clusterProfiler: functional enrichment analysis
#  enrichplot: visualisation of enrichment results
#  org.Hs.eg.db: human gene annotation database


# Load libraries
library(edgeR)
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
library(openxlsx)

