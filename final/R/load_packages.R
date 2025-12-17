install.packages("openxlsx")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("edgeR", "clusterProfiler", "enrichplot", "org.Hs.eg.db"))

library(edgeR)
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
library(openxlsx)

