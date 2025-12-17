#' Load RNA-seq count matrix from a tab-delimited file
#'
#'  @description
#'  It reads a tab-separated counts table where the first column is gene identifiers.
#'
#'  @param counts_path It paths to the counts file.
#'
#'  @return A numeric matrix of counts with genes as rownames and samples IDs as a colnames.
#'  @export


load_counts <- function(counts_path) {

  # header = TRUE: first row contains column names
  # check.names = FALSE: keep sample names exactly as they are
  x <- read.table(counts_path, header = TRUE, sep = "\t", check.names = FALSE)


  # Set the row names to the gene identifiers from the "Gene" column
  rownames(x) <- x$Gene
  x$Gene <- NULL
  as.matrix(x) #Convert the data frame to a matrix
}
