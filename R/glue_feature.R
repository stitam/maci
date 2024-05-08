glue_feature <- function(assembly, features, feature){
  features <- features[which(features$sseqid == feature),]
  features$forward <- features$send >= features$sstart
  features$sequence <- NA
  for (i in 1:nrow(features)) {
    index <- features$qstart[i]:features$qend[i]
    s <- paste0(assembly[[features$qseqid[i]]][index], collapse = "")
    if(features$forward[i] == FALSE) {
      s <- Biostrings::DNAString(s)
      s <- Biostrings::reverseComplement(s)
      s <- as.character(s)
    }
    features$sequence[i] <- toupper(s)
  }
  for (i in 1:nrow(features)) {
    if (features$forward[i] == FALSE) {
      a <- features$send[i]
      features$send[i] <- features$sstart[i]
      features$sstart[i] <- a
    }
  }
  merge_order <- order(features$sstart)
  post_N_count <- features$send[merge_order[length(merge_order)]] %% 3
  sequence <- c(
    rep("N", times = features$sstart[merge_order[1]]-1),
    strsplit(features$sequence[merge_order[1]], "")[[1]]
  )
  for (i in 2:length(merge_order)) {
    chunk <- strsplit(features$sequence[merge_order[i]], "")[[1]]
    d <- features$sstart[merge_order[i]]-features$send[merge_order[i-1]]
    if (d == 1) {
      # chunk follows previous chunk without gaps
      sequence <- c(sequence, chunk)
    }
    if (d > 1) {
      # chunk follows previous chunk but there are gaps
      sequence <- c(
        sequence,
        rep("N", times = d - 1),
        chunk
      )
    }
    if (d < 1) {
      # chunk follows previous chunk but there is overlap
      pre_overlap <-  sequence[(length(sequence) + d):length(sequence)]
      post_overlap <- chunk[1:(1 - d)]
      overlap <- ifelse(pre_overlap == post_overlap, pre_overlap, "N")
      sequence <- c(
        sequence[1:(length(sequence) + d - 1)],
        overlap,
        chunk[(2-d):length(chunk)]
      )
      print(sequence)
    }
  }
  sequence <- c(sequence, rep("N", times = post_N_count))
  sequence <- paste(sequence, collapse = "")
  return(sequence)
}
