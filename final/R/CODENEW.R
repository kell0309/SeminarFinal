counts  <- load_counts("E-MTAB-2523.counts.txt")

samples <- load_samples("E-MTAB-2523_sample table.txt")

filt <- filter_low_expression(counts, samples, group_col = "disease")

deg <- run_deg(
  dge = filt$dge,
  samples = filt$samples,
  group_col = "disease",
  ref_level = "normal",
  fdr_cutoff = 0.05,
  logfc_cutoff = 1
)

export_deg_excel(deg, "DEG_results.xlsx")

enr <- run_enrichment_go_kegg(deg, use_only_pass = TRUE)

# Plots
plot_enrichment(enr$go,   type = "dotplot", showCategory = 10)
plot_enrichment(enr$kegg, type = "dotplot", showCategory = 10)

