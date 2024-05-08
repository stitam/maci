#' Add coverage metrics to BLAST results
#' 
#' This function calculates coverage metrics based on BLAST results and adds
#' them to the output table.
#' @param sample_id character; sample ID
#' @param blastres data.frame; BLAST results
#' @return data.frame; BLAST results with new columns
#' @importFrom dplyr bind_cols
#' @export
add_coverage <- function(sample_id, blastres) {
  cov_vars <- c("qstart", "qend", "sstart", "send", "slen")
  if (all(cov_vars %in% names(blastres))) {
    if (any(is.na(cov_vars))) {
      blastres$covs <- NA
      blastres$covq <- NA
      blastres$qstrand <- NA
    } else {
      # calculate coverage 1 - fraction of subject in 'subject' match
      blastres$covs <- with(blastres, ifelse(
        send > sstart, (send - sstart + 1)/slen, (sstart-send +1)/slen
      ))
      blastres$covs <- round(blastres$covs, 5)
      # calculate coverage 2 - fraction of subject in 'query' match
      blastres$covq <- with(blastres, ifelse(
        qend > qstart, (qend - qstart + 1)/slen, (qstart-qend +1)/slen
      ))
      blastres$covq <- round(blastres$covq, 5)
    }
  } else {
    stop("Variable missing for custom blast coverage calculation.")
  }
  blastres <- dplyr::bind_cols(assembly = sample_id, blastres)
  return(blastres)
}
