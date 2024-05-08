rm(list=ls())
library(Biostrings)
library(dplyr)
library(reshape2)
library(seqinr)

args <- commandArgs(trailingOnly = TRUE)

resdir <- paste0("./jobs/", args[1], "/output/resfinder/resfinder_db")

filenames <- dir(resdir)[grep(".fsa$", dir(resdir))]

resdb <- vector()
resclasses <- data.frame()
for (i in seq_along(filenames)) {
  compclass <- strsplit(filenames[i], split = ".fsa$")[[1]][1]
  f <- read.fasta(paste0(resdir, "/", filenames[i]))
  f_aa <- lapply(f, function(x) {
    dna <- x %>% as.character() %>% paste(collapse = "") %>% DNAString()
    aa <- try(dna %>% Biostrings::translate(), silent = TRUE)
    if (inherits(aa, "try-error")) {
      return(NA)
    } else {
      aa <- gsub("\\*$", "", aa)
      return(aa)
    }
  })
  names(f_aa) <- sapply(f, function(x) attr(x, "name"))
  resdb <- c(resdb, f_aa)
  
  df <- data.frame(Resfinder_Gene_Accession = names(f_aa),
                   Resfinder_Class = compclass)
  resclasses <- rbind(resclasses, df)
}

write.fasta(sequences = resdb,
            names = names(resdb),
            as.string = TRUE,
            file.out = paste0("./jobs/", args[1], "/output/resfinder/resdb.fasta"))

resclasses_wide <- dcast(resclasses,
                         Resfinder_Gene_Accession ~ Resfinder_Class,
                         fun.aggregate = length)

write.csv(resclasses,
          file = paste0("./jobs/", args[1], "/output/resfinder/resclasses_long.csv"),
          row.names = FALSE)

write.csv(resclasses_wide,
          file = paste0("./jobs/", args[1], "/output/resfinder/resclasses_wide.csv"),
          row.names = FALSE)
