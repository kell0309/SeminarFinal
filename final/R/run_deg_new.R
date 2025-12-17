#' Differential expression analysis using edgeR
#'
#' @description
#' It performs differential expression analysis between two groups.
#'
#' @param dge An edgeR \code{DGEList} object containing filtered and normalized counts.
#' @param samples A data.frame containing sample metadata matched to the count data.
#' @param group_col A character string specifying the grouping column.
#' @param ref_level A character string specifying the reference group.
#' @param fdr_cutoff A numeric value specifying the FDR significance threshold by 0.05.
#' @param logfc_cutoff A numeric value specifying the absolute log2 fold-change.
#
#' @return A data.frame containing differential expression results and a logical column
#' indicating whether each gene passes the specified thresholds.



run_deg <- function(dge, samples,
                    group_col = "disease",
                    ref_level = "normal",
                    fdr_cutoff = 0.05,
                    logfc_cutoff = 1) {

  grp <- factor(samples[[group_col]])
  grp <- relevel(grp, ref = ref_level)

  design <- model.matrix(~ grp)
  dge <- estimateDisp(dge, design)
  fit <- glmQLFit(dge, design)
  qlf <- glmQLFTest(fit, coef = 2)

  tab <- topTags(qlf, n = Inf)$table
  tab$gene <- rownames(tab)
  tab$pass <- (tab$FDR <= fdr_cutoff) & (abs(tab$logFC) >= logfc_cutoff)

  tab[, c("gene", "logFC", "logCPM", "PValue", "FDR", "pass")]
}
