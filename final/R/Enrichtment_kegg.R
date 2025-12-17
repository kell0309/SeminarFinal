#' Enrichemt script for kegg
#'
#' @description
#' Takes a differential expression table  selects the  genes
#' runs GO  enrichment and KEGG pathway enrichment, and returns the results
#' @param deg_table A data frame containing at least `gene  and, if `use_only_pass = TRUE`, a logical `pass` column.
#' @param use_only_pass Logical. If TRUE, only genes with `deg_table$pass == TRUE` are used.
#' @param kegg_organism  KEGG organism code.
#' @param show_top Number of top categories to show in downstream plotting
#'
#' @return A list with Go enrichment results, kegg enrichment results, id maps from SYMBOL to ENTREZID, Show_top for top 10 value from kegg
#' @export
run_enrichment_go_kegg <- function(deg_table,
                                   use_only_pass = TRUE,
                                   kegg_organism = "hsa",
                                   show_top = 10) {

  # Get the gene list from the DEG table
  genes <- deg_table$gene

  # Restrict to genes that are signficant
  if (use_only_pass) genes <- deg_table$gene[deg_table$pass]

  # Run GO enrich analysis process using gene symbols
  go <- enrichGO(
    gene = genes,
    OrgDb = org.Hs.eg.db,
    keyType = "SYMBOL",
    ont = "BP",
    pvalueCutoff = 0.05,
    readable = TRUE
  )

  # Convert gene symbols to Entrez IDs needed for KEGG enrichmment
  conv <- bitr(genes,
               fromType = "SYMBOL",
               toType = "ENTREZID",
               OrgDb = org.Hs.eg.db)


  # Run KEGG pathway enrichment using Entrez IDs
  kegg <- enrichKEGG(
    gene = conv$ENTREZID,
    organism = kegg_organism,
    pvalueCutoff = 0.05
  )

  # Return both the enrichments results for GO and KEGG and the ID mappings
  list(go = go, kegg = kegg, id_map = conv, show_top = show_top)
}
