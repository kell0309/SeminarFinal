filter_low_expression <- function(counts, samples, group_col = "disease") {
  samples <- samples[match(colnames(counts), samples$sample), ]
  group <- samples[[group_col]]
  
  y <- DGEList(counts = counts, group = group)
  keep <- filterByExpr(y, group = group)
  
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y)
  
  list(dge = y, samples = samples)
}
