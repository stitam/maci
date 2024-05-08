library(Biostrings)
library(dplyr)
rm(list = ls())

args <- commandArgs(trailingOnly = TRUE)

# import R functions from the R directory
lapply(list.files(args[1], full.names = TRUE), source)
# sample_id
sample_id <- args[2]
# genome assembly
assembly <- seqinr::read.fasta(args[3])
# blast results for fq resgenes
blast_hit <- try(read.csv(args[4], sep = "\t", header = TRUE), silent = TRUE)
# fq resgenes
ref_faa <- seqinr::read.fasta(
  args[5], as.string = TRUE, forceDNAtolower = FALSE)

feature <-  c("gyrA", "gyrB", "parC")
blast_hit <- blast_hit[which(blast_hit$sseqid %in% feature),]

seqs <- list()
for (i in seq_along(blast_hit$sseqid)) {
  index <- which(blast_hit$sseqid == blast_hit$sseqid[i])
  # If there are multiple hits, keep the one with the highest bitscore
  # If there are still multiple hits, keep the first
  if (length(index) > 1) {
    seq_dna <- glue_feature(assembly, blast_hit, blast_hit$sseqid[i])
  } else {
    index2 <- blast_hit$qstart[index]:blast_hit$qend[index]
    seq_dna <- assembly[[blast_hit$qseqid[index]]][index2]
    seq_dna <- paste0(seq_dna, collapse = "")
    # If blast hit is in reversed order, calculate reverse complement
    if (blast_hit$sstart[index] > blast_hit$send[index]) {
      seq_dna <- Biostrings::reverseComplement(Biostrings::DNAString(seq_dna))
      seq_dna <- as.character(seq_dna)
    }
    seq_dna <- toupper(seq_dna)
  }
  sample_feature <- paste0(sample_id, "_", blast_hit$sseqid[i])
  seqs[[sample_feature]] <- seq_dna
}

# translate nucleotide sequence to amino acid sequence
sample_faa <- seqs %>%
  as.character() %>%
  Biostrings::DNAStringSet() %>%
  Biostrings::translate(., if.fuzzy.codon = "solve") %>% 
  as.character() %>%
  gsub("\\*$", "", .) %>% 
  as.list()

names(sample_faa) <- gsub(paste0("^", sample_id, "_"), "", names(seqs))

# These are the point mutations we are looking for
fqdet <- data.frame(
  name = c("gyrAS81L", "gyrBA414T", "parCS84L"),
  # The name of the feature must be exactly this for the script to find it
  feature = c("gyrA", "gyrB", "parC"),
  position = c(81, 414, 84),
  reference = c("S", "A", "S"),
  determinant = c("L", "T", "L")
)

# TODO put this function into R folder and add tests
# NOTE throws and error if Biostrings is not imported
find_determinants <- function(sample_seq, ref_seq, determinants) {
  
  if (mean(determinants$feature %in% names(ref_seq)) < 1) {
    index <- which(determinant$feature %in% names(ref_seq) == FALSE)
    missing <- paste(determinants$feature[index], collapse = ", ")
    stop("Features missing from reference sequences: ", missing, ".")
  }
  
  # return the character
  foo <- function(feature, position) {
    if (feature %in% names(sample_seq)) {
      aln <- Biostrings::pairwiseAlignment(
        sample_seq[[feature]], ref_seq[[feature]])
      mdf <- summary(aln)@mismatchSummary$subject
      index <- which(mdf$SubjectPosition == position)
      if (length(index) == 0) {
        return(strsplit(as.character(sample_seq[[feature]]), "")[[1]][position])
      }
      if (length(index) == 1) {
        return(mdf$Pattern[index])
      }
    }
    return(NA)
  }
  
  res <- mapply(foo, determinants$feature, determinants$position)
  
  res_df <- data.frame(
    feature = names(res),
    sample = unname(res)
  )
  
  out <- dplyr::left_join(determinants, res_df, by = "feature")
  out$hit <- out$sample == out$determinant
  
  return(out)
}

fqres <- find_determinants(sample_faa, ref_faa, fqdet)

df <- data.frame (
  assembly = sample_id
)

index <- which(fqres$hit == TRUE)

if (length(index) == 0) {
  df$fluoroquinolone <- NA
}

if (length(index) >= 1) {
  df$fluoroquinolone = paste(fqres$name[index], collapse = "|")
}

write.table(
  df,
  file = paste0(sample_id, "_predict_fq.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
