#' create plots for go and kegg
#'
#' @description
#' Creates a visualisation of an enrichment result using either a dotplot or a cnetplot.
#'
#' @param enrich_result An enrichment result object .
#' @param type Character. Plot type: `"dotplot"` or `"cnetplot"` (default: `"dotplot"`).
#' @param showCategory Integer. Number of categories to display (default: 10).
#'
#' @return A dotplot  or a plot produced by cnetplot()
#'
#' @export
plot_enrichment <- function(enrich_result,
                            type = c("dotplot", "cnetplot"),
                            showCategory = 10) {

  type <- match.arg(type)

  if (type == "dotplot") {
    dotplot(enrich_result, showCategory = showCategory)
  } else {
    cnetplot(enrich_result, showCategory = showCategory)
  }
}
