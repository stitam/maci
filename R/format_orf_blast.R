#' Format ORF BLAST results
#' 
#' This function uses ORF prediction and subsequent BLAST results to find BLAST
#' hit locations within genomes.
#' @param sample_id character; sample ID
#' @param orfs data.frame; predicted ORFs parsed from a GFF file
#' @param blastres data.frame; BLAST results
#' @return A formatted data frame.
#' @details Formatting ORF BLAST results is a bit tricky because BLAST is called
#' on the Open Reading Frames (ORFs) and not the genomes. Because of this, query
#' locations will be locations within an ORF and not within a contig, and we
#' will never find hits on negative strands because those have already been
#' reverse complemented in the ORF prediction process. This function joins
#' results from ORF prediction and BLAST to trace results to locations within
#' the original genome (before ORF prediction).
#' @importFrom dplyr rename left_join bind_cols
#' @export
format_orf_blast <- function(sample_id, orfs, blastres) {
  if ("qseqid" %in% names(blastres) == FALSE) {
    stop("BLAST results must include the qseqid column.")
  }
  if (nrow(orfs) == 1 && all(is.na(orfs))) {
  } else {
    ID <- gsub("^[0-9]+_", "_", get_attribute(orfs$attributes, "ID"))
    orfs$seqid <- paste0(orfs$seqid, ID)
    orfs$strand <- ifelse(orfs$strand == "+", 1, -1)
  }
  orfs <- dplyr::rename(
    orfs, qseqid = seqid, qstart = start, qend = end, qstrand = strand)
  orfs <- orfs[, c("qseqid", "qstart", "qend", "qstrand")]
  index <- which(names(blastres) %in% c("qstart", "qend"))
  if (length(index) > 0) {
    blastres <- blastres[, -index]
  }
  df <- dplyr::left_join(blastres, orfs, by = "qseqid")
  df <- dplyr::bind_cols(assembly = sample_id, df)
  for (i in 1:nrow(df)) {
    if (!is.na(df$qstrand[i])) {
      if (df$qstrand[i] == -1) {
        aux <- df$send[i]
        df$send[i] <- df$sstart[i]
        df$sstart[i] <- aux
      }
    }
  }
  return(df)
}
