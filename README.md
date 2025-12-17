# RNA-seq Differential Expression Analysis Pipeline
## Overview
The repository contains an R-based RNA-seq analysis pipeline. The pipeline performs data loading, filtering of low-expressed genes, normalization and differential gene expression analysis using the edgeR package. The main aim of this project is code clearity, modular design, and proper documentation of functions.

## Repository Structure
1.load_counts:new.R #Load raw count matrix 
2.load_samples_new.R #Load sample metadata
3.filter_low.R #Filter low-expressed genes and normalize counts
4.run_deg_new.R #Run differential expression analysis

## Usage
A typical analysis workflow is:

```r
counts <- load_counts("E-MTAB-2523.counts.txt")
samples <- load_samples("E-MTAB-2523_sample table.txt")

flt <- filter_low_expression(counts, samples)

deg <- run_deg(flt$dge, flt$samples)
```
## Documentation
All functions are documented using roxygen2 style comments, containing:
- Function descriptions
- Parameter explanations
- Return values

## Version Control
Development was carried out on a separate Git branch and merged via pull request to ensure safe and traceable changes, following collaborative software develepment best practices.
