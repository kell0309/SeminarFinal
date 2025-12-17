run_enrichment_go_kegg <- function(deg_table,
                                   use_only_pass = TRUE,
                                   kegg_organism = "hsa",
                                   show_top = 10) {

  genes <- deg_table$gene
  if (use_only_pass) genes <- deg_table$gene[deg_table$pass]

  go <- enrichGO(
    gene = genes,
    OrgDb = org.Hs.eg.db,
    keyType = "SYMBOL",
    ont = "BP",
    pvalueCutoff = 0.05,
    readable = TRUE
  )

  conv <- bitr(genes,
               fromType = "SYMBOL",
               toType = "ENTREZID",
               OrgDb = org.Hs.eg.db)

  kegg <- enrichKEGG(
    gene = conv$ENTREZID,
    organism = kegg_organism,
    pvalueCutoff = 0.05
  )

  list(go = go, kegg = kegg, id_map = conv, show_top = show_top)
}
