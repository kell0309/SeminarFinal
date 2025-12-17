#' Filter low-expressed genes and create a normalized edgeR DGEList
#'
#' @description
#' It matches sample metadata to the count matrix column, then filter low expressed gene.
#
#' @param counts A numeric matrix including raw gene expression counts, with genes as rows and samples
#' @param samples A data.frame containing sample metadata.
#' @param group_col A character string specifying the column name.
#'
#' @return A named list containing an edgeR \code{DGEList} object after filtering and normalization,
#' and a reordered sample metadata data.frame matched to the count matrix columns.

filter_low_expression <- function(counts, samples, group_col = "disease") {
  samples <- samples[match(colnames(counts), samples$sample), ]
  group <- samples[[group_col]]

  y <- DGEList(counts = counts, group = group)
  keep <- filterByExpr(y, group = group)

  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y)

  list(dge = y, samples = samples)
}
