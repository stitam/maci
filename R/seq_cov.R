#' Calculate sequencing coverage
#'
#' This function takes an assembly fasta file and attempts to calculate
#' sequencing coverage from contig headers.
#' @param file character; file path to the fasta file
#' @return A list with two elements. The first element is a tibble and contains
#' coverage data for each contig. The second element is the overall coverage for
#' the whole assembly. Overall coverage is calculated as the weighted average of
#' contig coverage, weighted by contig length.
#' @export
seq_cov <- function(file) {
  contigs <- Biostrings::readDNAStringSet(file)
  meta <- strsplit(names(contigs), "_")
  meta <- lapply(meta, function(x) {
    var_names <- x[which(seq_along(x) %% 2 == 1)]
    var_values <- x[which(seq_along(x) %% 2 == 0)]
    df <- tibble::as_tibble(t(var_values))
    names(df) <- var_names
    return(df)
  })
  meta <- dplyr::bind_rows(meta)
  meta$length <- as.numeric(meta$length)
  meta$cov <- as.numeric(meta$cov)
  seq_cov <- signif(sum(meta$length*meta$cov)/sum(meta$length), 4)
  return(list("cov_by_contig" = meta, "cov_overall" = seq_cov))
}