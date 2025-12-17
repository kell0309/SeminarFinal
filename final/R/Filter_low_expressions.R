#'Filters lowly expressed genes and creates a normalised list using edgeR and DGElist

#' @description
#' This code aligns the sample metadata and the count matrix, also removes lowly expressed genes.
#'
#' @param counts A data frame of raw counts
#' @param samples a data frame with sample table
#' @param group_col a column names with in samples

#' @return a list of filtered expression from samples
#' @export

filter_low_expression <- function(counts, samples, group_col = "disease") {

  # Align sample metadata rows to the count matrix columns
  samples <- samples[match(colnames(counts), samples$sample), ]
  group <- samples[[group_col]]


  # identify genes to keep, removes genes with too-low expression across group
  y <- DGEList(counts = counts, group = group)
  keep <- filterByExpr(y, group = group)


  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y)



  # Return the filtered + normalized DGEList and  the aligned sample table
  list(dge = y, samples = samples)
}

