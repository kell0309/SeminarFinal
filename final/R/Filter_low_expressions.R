#'Filters lowly expressed genes and creates a normalised list using edgeR and DGElist

#' #´@description
#' This code aligns the sample metadata and the count matrix, also removes lowly expressed genes.
#'
#' @param counts A data frame of raw counts
#' @samples a data frame with sample table
#' @group col a column names with in samples

#' @return
filter_low_expression <- function(counts, samples, group_col = "disease") {
  samples <- samples[match(colnames(counts), samples$sample), ]
  group <- samples[[group_col]]

  y <- DGEList(counts = counts, group = group)
  keep <- filterByExpr(y, group = group)

  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y)

  list(dge = y, samples = samples)
}
