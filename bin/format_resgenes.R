rm(list=ls())
library(dplyr)
library(tidyr)

args <- commandArgs(trailingOnly = TRUE)

# import R functions from the R directory
lapply(list.files(args[1], full.names = TRUE), source)

dbdir <- paste0(args[2], "/resfinder/2.0.0/")
sample_id <- args[3]
diamond <- try(read.table(args[4]), silent = TRUE)
coverage <- as.numeric(args[5])

myPaths <- .libPaths()
myPaths <- c(myPaths, args[6])

orfs <- try(ape::read.gff(args[7]), silent = TRUE)

# format blast results

blast_names <- c(
  "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart",
  "qend", "sstart", "send", "slen", "evalue", "bitscore"
)

if (inherits(diamond, "try-error")) {
  diamond <- data.frame(matrix(NA, nrow = 1, ncol = length(blast_names)))
}

names(diamond) <- blast_names

if (inherits(orfs, "try-error")) {
  orfs <- data.frame(matrix(NA, nrow = 1, ncol = 9))
  names(orfs) <- c(
    "seqid", "source", "type", "start", "end", "score", "strand", "phase",
    "attributes"
  )
}

blast_resfinder <- format_orf_blast(sample_id, orfs, diamond)

if (any(!is.na(blast_resfinder$sseqid))) {
  blast_resfinder$sseqid <- sapply(blast_resfinder$sseqid, function(x) {
    if (is.na(x)) return(x)
    strsplit(x, "_")[[1]][1]
  })
}

write.table(
  blast_resfinder,
  file = paste0(sample_id, "_format_resgenes.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# format phenotypes

phenotypes <- read.delim(paste0(dbdir, "resfinder_db/phenotypes.txt"))
phenotypes$Gene_accession.no.<- sapply(
  phenotypes$Gene_accession.no., function(x) strsplit(x, "_")[[1]][1])

abs <- phenotypes$Phenotype %>% 
  strsplit(split = ", +") %>% 
  unlist() %>% 
  tolower() %>% 
  gsub("-", "_", .) %>%
  gsub(" +", "_", .) %>%
  gsub("\\+", "_", .) %>%
  gsub(",", "_", .) %>%
  gsub("_+", "_", .) %>%
  unique() %>% 
  sort()

resclasses <- readRDS(paste0(dbdir, "resclasses.rds"))
resclasses$Resfinder_Gene_Accession <- sapply(
  resclasses$Resfinder_Gene_Accession, function(x) strsplit(x, "_")[[1]][1])
resclasses$Resfinder_Class <- gsub("-", "_", resclasses$Resfinder_Class)
resclasses <- dplyr::distinct(resclasses)

classes <- resclasses$Resfinder_Class %>% unique() %>% sort()

genes <- sapply(phenotypes$Gene_accession.no., function(x) {
  strsplit(x, "_")[[1]][1]
}) %>% unique() %>% sort()
  
blast_resfinder <- dplyr::rename(blast_resfinder, gene = sseqid)

# EACH COLUMN IS AN ANTIBIOTIC, VALUES ARE RESISTANCE GENES

if (nrow(diamond) == 1 && all(is.na(diamond))) {
  df_wider <- as.data.frame(matrix(NA, nrow = 1, ncol = length(abs) + 1))
  names(df_wider) <- c("assembly", abs)
  df_wider$assembly <- sample_id
} else {
  df <- blast_resfinder
  # if a gene is found more than once only include it once
  df <- df[,c("assembly", "gene")]
  df <- dplyr::distinct(df)
  # merge diamond results with phenotypes
  df <- dplyr::left_join(df, phenotypes, by = c("gene" = "Gene_accession.no." ))
  # keep relevant columns
  df <- df[, which(names(df) %in% c("assembly", "gene", "Phenotype"))]
  # remove any rows that contain NAs. Note that Phenotype for gene formA is NA.
  df <- df[complete.cases(df),]

  longer <- function(row) {
    if (grepl(", +", row$Phenotype)){
      values <- strsplit(row$Phenotype, ", +")[[1]]
      newrows <- data.frame(
        assembly = row$assembly,
        gene = row$gene,
        Phenotype = values
      )
    } else {
      newrows = row
    }
    return(newrows)
  }

  df_longer <- data.frame()
  for (i in 1:nrow(df)){
    df_longer <- dplyr::bind_rows(df_longer, longer(df[i,]))
  }

  df_longer$Phenotype <- df_longer$Phenotype %>%
    tolower() %>% 
    gsub("-", "_", .) %>%
    gsub(" +", "_", .) %>%
    gsub("\\+", "_", .) %>%
    gsub(",", "_", .) %>%
    gsub("_+", "_", .)

  df_wider <- tidyr::pivot_wider(
    df_longer,
    id_cols = assembly,
    names_from = Phenotype,
    names_sort = TRUE,
    values_from = gene,
    values_fn = function(x) paste(sort(unique(x)), collapse = "|"),
  )

  tmp <- as.data.frame(matrix(NA, nrow = 1, ncol = length(abs) + 1))
  names(tmp) <- c("assembly", abs)

  df_wider <- dplyr::bind_rows(tmp, df_wider)[-1,]
}

write.table(
  df_wider,
  file = paste0(sample_id, "_abs.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# EACH COLUMN IS A CLASS, VALUES ARE RESISTANCE GENES

if (nrow(diamond) == 1 && all(is.na(diamond))) {
  df_wider <- as.data.frame(matrix(NA, nrow = 1, ncol = length(classes) + 1))
  names(df_wider) <- c("assembly", classes)
  df_wider$assembly <- sample_id
} else {
  df <- blast_resfinder
  # if a gene is found more than once only include it once
  df <- df[,c("assembly", "gene")]
  df <- dplyr::distinct(df)
  # merge diamond results with resclasses
  df <- dplyr::left_join(df, resclasses, by = c("gene" = "Resfinder_Gene_Accession" ))
  # keep relevant columns
  df <- df[, which(names(df) %in% c("assembly", "gene", "Resfinder_Class"))]
  
  df_wider <- tidyr::pivot_wider(
    df,
    id_cols = assembly,
    names_from = Resfinder_Class,
    names_sort = TRUE,
    values_from = gene,
    values_fn = function(x) paste(sort(unique(x)), collapse = "|"),
  )
  
  tmp <- as.data.frame(matrix(NA, nrow = 1, ncol = length(classes) + 1))
  names(tmp) <- c("assembly", classes)
  
  df_wider <- dplyr::bind_rows(tmp, df_wider)[-1,]
}

write.table(
  df_wider,
  file = paste0(sample_id, "_classes.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# EACH COLUMN IS AN ANTIBIOTIC RESISTANCE GENE, VALUES ARE RESISTANCE GENES

if (nrow(diamond) == 1 && all(is.na(diamond))) {
  df_wider <- as.data.frame(matrix(NA, nrow = 1, ncol = length(genes) + 1))
  names(df_wider) <- c("assembly", genes)
  df_wider$assembly <- sample_id
} else {
  df <- blast_resfinder
  # if a gene is found more than once only include it once
  df <- df[,c("assembly", "gene")]
  df <- dplyr::distinct(df)
  # merge diamond results with phenotypes
  df <- dplyr::left_join(df, phenotypes, by = c("gene" = "Gene_accession.no." ))
  # keep relevant columns
  df <- df[, which(names(df) %in% c("assembly", "gene"))]
  # because of shortened gene names some genes may occur multiple times
  df <- dplyr::distinct(df)
  
  df_wider <- tidyr::pivot_wider(
    df,
    id_cols = assembly,
    names_from = gene,
    names_sort = TRUE,
    values_from = gene,
    values_fn = function(x) !is.na(x),
  )
  
  tmp <- as.data.frame(matrix(NA, nrow = 1, ncol = length(genes) + 1))
  names(tmp) <- c("assembly", genes)
  
  df_wider <- dplyr::bind_rows(tmp, df_wider)[-1,]
}

write.table(
  df_wider,
  file = paste0(sample_id, "_genes.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
