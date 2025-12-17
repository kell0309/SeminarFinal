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
