# RNA-seq Differential Expression Analysis Pipeline
## Overview
The repository contains an R-based RNA-seq analysis pipeline. The pipeline performs data loading, filtering of low-expressed genes, normalisation and differential gene expression analysis using the edgeR package. The main aim of this project is code clearity, modular design, and proper documentation of functions.

## Repository Structure

1. load_counts_new.R -Load raw count matrix 
2. load_samples_new.R -Load sample metadata
3. filter_low.R -Filter low-expressed genes and normalise counts
4. run_deg_new.R -Run differential expression analysis
5. Enrichment_kegg.R -Run GO and KEGG enrichment analysis
6. Enrichment_plot.R -Visualize enrichment results
7. excel_table.R -Export DEG results to Excel
8. load.packages.R -Load required R packages 

## Usage
A typical analysis workflow is:

```r
counts <- load_counts("E-MTAB-2523.counts.txt")
samples <- load_samples("E-MTAB-2523_sample table.txt")

flt <- filter_low_expression(counts, samples)

deg <- run_deg(flt$dge, flt$samples)
```

## GO and KEGG Enrichment Analysis

```r
enr <- run_enrichment_go_kegg(deg, use_only_pass = TRUE)
```
### This function:

- Selects genes from the differential expression table
- Runs GO Biological process enrichment 
- Runs KEGG pathway enrichment
- Returns enrichment results and gene ID mappings

## Visualisation

```r
plot_enrichment(enr$go, type = "dotplot", showCategory = 10)
plot_enrichment(enr$kegg, type = "dotplot", showCategory = 10)
```

## Results of Export

```r
export_deg_excel(deg, out_xlsx = "DEG_results.xlsx")
```

## Documentation
All functions are documented using roxygen2 style comments, containing:
- Function descriptions
- Parameter explanations
- Return values

## Version Control
Development was carried out on a separate Git branch and merged via pull request to ensure safe and traceable changes, following collaborative software development best practices.

