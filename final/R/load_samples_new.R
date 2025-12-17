load_samples <- function(samples_path) {
  s <- read.table(samples_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  s$disease <- factor(s$disease)
  s
}
