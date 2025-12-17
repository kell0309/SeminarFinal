#' Load samples metadata samples
#'
#' @description
#' It reads a tab-separated sample table and converts.
#'
#' @param sample_path It pats to the sample table
#'
#' @return A data.frame contains sample annotations.
#'
#' @export

load_samples <- function(samples_path) {

  # header = TRUE: first row contains column names
  s <- read.table(samples_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  s$disease <- factor(s$disease)


  # Return the processed sample metadata
  s
}
