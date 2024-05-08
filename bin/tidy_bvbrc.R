# This script is used to prepare an AMR data set for subsequently training a
# prediction model. The script uses data downloaded from BV-BRC
# (https://www.bv-brc.org/).

library(dplyr)
library(optparse)
rm(list = ls())

args_list <- list(
  make_option(
    c("-g", "--genome"),
    type = "character",
    help = "BVBRC_genome.csv downloaded from https://www.bv-brc.org/"
  ),
  make_option(
    c("-a", "--amr"),
    type = "character",
    help = "BVBRC_genome_amr.csv downloaded from https://www.bv-brc.org/"
  )
)

args_parser  <- OptionParser(option_list = args_list)

if (!interactive()) {
  args  <- parse_args(args_parser)
} else {
  args <- list(
    genome = "BVBRC_genome.csv",
    amr = "BVBRC_genome_amr.csv"
  )
}

genome <- read.csv(args$genome, na.strings = "")
genome <- genome[,which(names(genome) %in% c(
  "Genome.ID", "BioSample.Accession" 
))]
genome <- genome[complete.cases(genome),]

amr <- read.csv(args$amr, na.strings = "") %>%
  filter(Evidence == "Laboratory Method" | 
         Laboratory.Typing.Method == "Broth dilution"
  ) %>% select("Genome.ID", "Antibiotic", "Resistant.Phenotype")

df <- dplyr::left_join(
  genome,
  amr,
  by = "Genome.ID",
  relationship = "many-to-many"
)

df <- df[-which(names(df) == "Genome.ID")]
df <- df[complete.cases(df),]

testthat::expect_true(all(df$Resistant.Phenotype %in% c(
  "Resistant", "Susceptible", "Intermediate"
)))

saveRDS(df, "amr_bvbrc.rds")
