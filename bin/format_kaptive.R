rm(list=ls())

args <- commandArgs(trailingOnly = TRUE)
sample_id <- args[1]
k <- args[2]
o <- args[3]

myPaths <- .libPaths()
myPaths <- c(myPaths, args[4])

k <- read.table(k, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
names(k) <- c(
  "assembly", "k_serotype", "k_type", "k_confidence", "k_problems", "k_coverage",
  "k_identity", "k_length_discrepancy", "k_expected_genes_in_locus",
  "k_expected_genes_in_locus_details", "k_missing_expected_genes", 
  "k_other_genes_in_locus", "k_other_genes_in_locus_details",
  "k_expected_genes_outside_locus", "k_expected_genes_outside_locus_details",
  "k_other_genes_outside_locus", "k_other_genes_outside_locus_details")

o <- read.table(o, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
names(o) <- c(
  "assembly", "oc_serotype", "oc_type","oc_confidence", "oc_problems", "oc_coverage",
  "oc_identity", "oc_length_discrepancy", "oc_expected_genes_in_locus",
  "oc_expected_genes_in_locus_details", "oc_missing_expected_genes",
  "oc_other_genes_in_locus", "oc_other_genes_in_locus_details",
  "oc_expected_genes_outside_locus", "oc_expected_genes_outside_locus_details",
  "oc_other_genes_outside_locus", "oc_other_genes_outside_locus_details")

kaptive <- plyr::join(k, o, by = "assembly", type = "full")

write.table(
  kaptive,
  file = paste0(sample_id, "_formatted.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
