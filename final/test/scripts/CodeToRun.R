## This the pipeline that runs the full RNA-seq

# 1. Load the input data
counts  <- load_counts("E-MTAB-2523.counts.txt")
samples <- load_samples("E-MTAB-2523_sample table.txt")


# 2. Filter the lowley expressed genes and normalise the counts
filt <- filter_low_expression(counts, samples, group_col = "disease")

# 3. Differential expression analysis
deg <- run_deg(
  dge = filt$dge,
  samples = filt$samples,
  group_col = "disease",
  ref_level = "normal",
  fdr_cutoff = 0.05,
  logfc_cutoff = 1
)


# 4.export DEG results to excel
export_deg_excel(deg, "DEG_results.xlsx")

# 5. Functional enrichment using significant genes
enr <- run_enrichment_go_kegg(deg, use_only_pass = TRUE)

# 6. Plots to visualise the enrichment results
plot_enrichment(enr$go,   type = "dotplot", showCategory = 10)
plot_enrichment(enr$kegg, type = "dotplot", showCategory = 10)

