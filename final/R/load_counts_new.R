load_counts <- function(counts_path) {
  x <- read.table(counts_path, header = TRUE, sep = "\t", check.names = FALSE)
  rownames(x) <- x$Gene
  x$Gene <- NULL
  as.matrix(x)
}
