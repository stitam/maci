rm(list=ls())

args <- commandArgs(trailingOnly = TRUE)

# import R functions from the R directory
lapply(list.files(args[1], full.names = TRUE), source)

sample_id <- args[2]
blast_iseqs <- try(read.table(args[3]), silent = TRUE)
coverage <- as.numeric(args[4])

blast_names <- c(
  "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart",
  "qend", "sstart", "send", "slen", "evalue", "bitscore", "qseq", "sseq"
)

if (inherits(blast_iseqs, "try-error")) {
  blast_iseqs <- data.frame(matrix(NA, nrow = 1, ncol = length(blast_names)))
}

names(blast_iseqs) <- blast_names

blast_iseqs$sseqid <- sapply(blast_iseqs$sseqid, function(x) {
  if (is.na(x)) return(x)
  strsplit(x, "_")[[1]][1]
})

blast_iseqs <- add_coverage(sample_id, blast_iseqs)

blast_iseqs$qstrand <- sign(blast_iseqs$send-blast_iseqs$sstart)

# filter results by coverage
# if blast_iseqs inherits try-error, this line will produce and empty data.frame
blast_iseqs <- blast_iseqs[which(blast_iseqs$covs >= coverage),]

write.table(
  blast_iseqs,
  file = paste0(sample_id, "_format_iseqs.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
