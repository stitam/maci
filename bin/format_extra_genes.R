rm(list=ls())

args <- commandArgs(trailingOnly = TRUE)

# import R functions from the R directory
lapply(list.files(args[1], full.names = TRUE), source)

sample_id <- args[2]
blast_extra_genes <- try(read.table(args[3]), silent = TRUE)
coverage <- as.numeric(args[4])

blast_names <- c(
  "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart",
  "qend", "sstart", "send", "slen", "evalue", "bitscore", "qseq", "sseq"
)

if (inherits(blast_extra_genes, "try-error")) {
  blast_extra_genes <- data.frame(matrix(NA, nrow = 1, ncol = length(blast_names)))
}

names(blast_extra_genes) <- blast_names

blast_extra_genes$sseqid <- sapply(blast_extra_genes$sseqid, function(x) {
  if (is.na(x)) return(x)
  strsplit(x, "_")[[1]][1]
})

blast_extra_genes <- add_coverage(sample_id, blast_extra_genes)

blast_extra_genes$qstrand <- sign(blast_extra_genes$send-blast_extra_genes$sstart)

# filter results by coverage
index <- which(blast_extra_genes$covs >= coverage)
blast_extra_genes <- blast_extra_genes[index, ]

write.table(
  blast_extra_genes,
  file = paste0(sample_id, "_format_extra_genes.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
