library(dplyr)
library(optparse)
rm(list = ls())

args_list <- list(
  make_option(
    c("-s", "--sample_id"),
    type = "character",
    help = "Sample ID."
  ),
  make_option(
    c("-r", "--resgenes"),
    type = "character",
    help = "Path to formated resfinder results."
  ),
  make_option(
    c("-x", "--extra_genes"),
    type = "character",
    help = "Path to formated extra genes blast results."
  ),
  make_option(
    c("-i", "--insertion_sequences"),
    type = "character",
    help = "Path to formated insertion sequences blast results."
  ),
  make_option(
    c("-g", "--resgene_groups"),
    type = "character",
    help = "Path to a data frame that maps resgenes to resgene groups."
  )
)

args_parser  <- OptionParser(option_list = args_list)

if (!interactive()) {
  args  <- parse_args(args_parser)
} else {
  args <- list(
    sample_id = "Aci1",
    resgenes = "Aci1_format_resgenes.tsv",
    extra_genes = "Aci1_format_extra_genes.tsv",
    insertion_sequences = "Aci1_format_iseqs.tsv",
    resgene_groups = "~/Methods/aci/data/resgene_groups.csv"
  )
}

# define sample_id
sample_id <- args$sample_id
# read blast results for resgenes
blast_resgenes <- read.csv(args$resgenes, sep = "\t")
# read blast results for extra genes
blast_extra_genes <- read.csv(args$extra_genes, sep = "\t")
# read blast results for extra features
blast_iseqs <- read.csv(args$insertion_sequences, sep = "\t")
# import resgene groups dictionary
resgroups <- read.csv(args$resgene_groups, sep = "\t")

blast_resgenes$sseqid <- sapply(blast_resgenes$sseqid, function(x) {
  index <- which(resgroups$gene == x)
  # if gene name is not in db return gene name.
  if (length(index) == 0) {
    return(x)
  }
  # if gene name is in db return group name
  if (length(index) == 1) {
    return(resgroups$group[index])
  }
  # if gene name is in db multiple times, return gene name and a warning
  if (length(index) >= 1) {
    msg <- paste0(
      x, " belongs to multiple groups: ",
      paste(resgroups$group[index], collapse = ", "),
      ". Returning ", x, "."
    )
    warning(msg)
    return(x)
  }
})

# Currently the only redundant entry in the resgene group database is blaOXA-199:
# blaOXA-199 is categorised as blaOXA-48-like in one publication, blaOXA-51-like
# in another. Here we decide to categorise blaOXA-199 as blaOXA-51-like because
# blaOXA-51-like genes may confer resistance against carbapenems.
# TODO find a general solution for clashes.
for (i in 1:nrow(blast_resgenes)){
  if (!is.na(blast_resgenes$sseqid[i]) && blast_resgenes$sseqid[i] == "blaOXA-199") {
    blast_resgenes$sseqid[i] == "blaOXA-51-like"
  }
}

# Contains at least one of the following carbapenem resistance genes:

# 2023-01-01 updates to previous version:
# blaOXA-23: renamed, blaOXA-23-like
# blaOXA-24: removed, blaOXA-24 was renamed to blaOXA-40
# blaOXA-25: removed, blaOXA-40-like
# blaOXA-26: removed, blaOXA-40-like
# blaOXA-40: renamed, blaOXA-40-like
# blaOXA-58: renamed, blaOXA-58-like
# blaOXA-143: renamed, blaOXA-143-like
# blaOXA-235: renamed, blaOXA-235-like
# blaOXA-6: removed, not in ResFinder
# blaOXA-66, removed, blaOXA-51-like but that group is handled separately
# blaOXA-72: removed, blaOXA-40-like

# 2023-07-20
# After analysing measured carbapenem resistance using the BV-BRC database the
# following changes were made to the resistance definition:
# blaOXA-312: added, 16 True Positive, 0 False Positive

carbapenem_db <- c(
  "blaOXA-23-like",
  "blaOXA-40-like",
  "blaOXA-58-like",
  "blaOXA-143-like",
  "blaOXA-235-like",
  "blaOXA-312",
  "NDM", # MULTIPLE PRESENT in resfinder_db, beta-lactam.fsa
  "VIM", # MULTIPLE PRESENT in resfinder_db, beta-lactam.fsa
  "IMP", # MULTIPLE PRESENT in resfinder_db, beta-lactam.fsa
  "KPC", # MULTIPLE PRESENT in resfinder_db, beta-lactam.fsa
  "GES" # MULTIPLE PRESENT in resfinder_db, beta-lactam.fsa
)

carbapenem_resfinder <- blast_resgenes$sseqid[which(
  blast_resgenes$sseqid %in% carbapenem_db | # present in carbapenem_db object
    grepl("^blaNDM-[0-9]+", blast_resgenes$sseqid) | # any NDM gene
    grepl("^blaVIM-[0-9]+", blast_resgenes$sseqid) | # any VIM gene
    grepl("^blaIMP-[0-9]+", blast_resgenes$sseqid) | # any IMP gene
    grepl("^blaKPC-[0-9]+", blast_resgenes$sseqid) | # any KPC gene
    grepl("^blaGES-[0-9]+", blast_resgenes$sseqid) # any GES gene
)] %>% unique() %>% sort()

# carO gene is missing
# https://doi.org/10.1371/journal.pone.0122793

# 2023-07-20
# After analysing measured carbapenem resistance using the BV-BRC database
# presence/absence of carO was not found to be a good predictor of carbapenem
# 663 True Positives, 611 False Positives

if (nrow(blast_extra_genes) == 1 && is.na(blast_extra_genes$qseqid)) {
  carO <- TRUE
  #carbapenem_extra_genes <- vector()
} else {
  if ("carO" %in% blast_extra_genes$sseqid) {
    carO <- TRUE
    #carbapenem_extra_genes <- vector()
  } else {
    carO <- FALSE
    #carbapenem_extra_genes <- "missing-carO"
  }
}

# Contains ISAba1 upstream of a blaOXA-51-like gene, in the opposite direction
# https://doi.org/10.1128/CMR.00117-13

# If no insertion sequences can be found return empty vector
if (nrow(blast_iseqs) == 0) {
  ISAba1 <- FALSE
  # create an empty data frame
  blast_iseqs_range <- as.data.frame(matrix(NA, nrow = 1, ncol = 5))
  names(blast_iseqs_range) <- c(
    "assembly","contig","iseq_qstart", "iseq_qend", "iseq_dir")
  blast_iseqs_range$assembly <- sample_id
} else {
  ISAba1 <- TRUE
  # filter and rename blast_iseqs
  blast_iseqs_range <- blast_iseqs[, c("qseqid","qstart", "qend", "qstrand")]
  blast_iseqs_range <- dplyr::distinct(blast_iseqs_range)
  blast_iseqs_range <- dplyr::rename(
    blast_iseqs_range,
    contig = qseqid,
    iseq_qstart = qstart,
    iseq_qend = qend,
    iseq_dir = qstrand
  )
  blast_iseqs_range <- dplyr::bind_cols(
   data.frame("assembly" = sample_id), blast_iseqs_range)
}

# If no blaOXA-51-like gene can be found return empty vector
if ("blaOXA-51-like" %in% blast_resgenes$sseqid == FALSE) {
  blast_combos <- as.data.frame(matrix(NA, nrow = 1, ncol = 6))
  names(blast_combos) <- c(
    "assembly", "contig", "sseqid", "gene_qstart", "gene_qend", "gene_dir")
  blast_combos$assembly <- sample_id
} else {
  # filter resgenes to blaOXA-51
  blast_combos <- blast_resgenes[which(blast_resgenes$sseqid %in% c(
    "blaOXA-51-like"
  )),]
  blast_combos <- blast_combos[,c("qseqid", "sseqid", "qstart", "qend", "qstrand")]
  blast_combos$qseqid <- gsub("_[0-9]+$", "", blast_combos$qseqid)
  blast_combos <- dplyr::rename(
    blast_combos,
    contig = qseqid,
    gene_qstart = qstart,
    gene_qend = qend,
    gene_dir = qstrand
  )
  blast_combos <- dplyr::bind_cols(
    data.frame("assembly" = sample_id), blast_combos)
}

# join the two tables
df <- dplyr::full_join(
  blast_iseqs_range, blast_combos, by = c("assembly", "contig"))

# look for valid combinations
# source: https://doi.org/10.1128/CMR.00117-13
# article mentions 25db for blaOXA-23 and 7bp for blaOXA-51
# let's be permissive
threshold <- 25

df$threshold <- NA
df$combo <- FALSE
for (i in 1:nrow(df)) {
  if (complete.cases(df[i,-which(names(df) == "threshold")])) {
    # if gene and iseq are on the same contig
    if (df$gene_dir[i] != df$iseq_dir[i]) {
      # if gene and iseq are in opposite directions
      if (df$gene_dir[i] == 1) {
        df$threshold[i] <- df$gene_qstart[i] - df$iseq_qend[i]
      } else {
        df$threshold[i] <- df$iseq_qstart[i] - df$gene_qend[i]
      }
      if (df$threshold[i] > 0 & df$threshold[i] < threshold) {
        df$combo[i] <- TRUE
      }
    }
  }
}
    
write.table(
  df,
  file = paste0(sample_id, "_iseq_debug.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

if (any(df$combo == TRUE)) {
  carbapenem_combos = "ISAba1_blaOXA-51-like"
} else {
  carbapenem_combos = vector()
}

# Collect resistance determinants

# 2023-07-20
# After analysing measured carbapenem resistance using the BV-BRC database
# presence/absence of carO was not found to be a good predictor of carbapenem
# 663 True Positives, 611 False Positives

carbapenem_sample <- c(
  carbapenem_resfinder,
  #carbapenem_extra_genes,
  carbapenem_combos
) %>% unique() %>% sort() %>% paste(collapse = "|")

carbapenem <- data.frame(
  sample_id = sample_id,
  carO = carO,
  ISAba1 = ISAba1,
  blaOXA51like = "blaOXA-51-like" %in% blast_resgenes$sseqid,
  carbapenem = carbapenem_sample
)

write.table(
  carbapenem,
  file = paste0(sample_id, "_predict_carbapenem.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
