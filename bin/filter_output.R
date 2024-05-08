rm(list=ls())

args <- commandArgs(trailingOnly = TRUE)
pipeline_version <- args[1]
results <- read.csv(args[2], sep = "\t", stringsAsFactors = FALSE)

myPaths <- .libPaths()
myPaths <- c(myPaths, args[3])

results <- results[order(results$assembly), ]
results <- tibble::as_tibble(results)

exclude <- c(
  "k_type",
  "k_problems",
  "k_coverage",
  "k_identity",
  "k_length_discrepancy",
  "k_expected_genes_in_locus",
  "k_expected_genes_in_locus_details", 
  "k_missing_expected_genes",
  "k_other_genes_in_locus",
  "k_other_genes_in_locus_details",  
  "k_expected_genes_outside_locus",
  "k_expected_genes_outside_locus_details",
  "k_other_genes_outside_locus",
  "k_other_genes_outside_locus_details",
  "oc_type",
  "oc_problems",
  "oc_coverage",
  "oc_identity",
  "oc_length_discrepancy",
  "oc_expected_genes_in_locus",
  "oc_expected_genes_in_locus_details",
  "oc_missing_expected_genes",
  "oc_other_genes_in_locus",
  "oc_other_genes_in_locus_details",
  "oc_expected_genes_outside_locus",
  "oc_expected_genes_outside_locus_details",
  "oc_other_genes_outside_locus",
  "oc_other_genes_outside_locus_details"
)

results_short <- results[, which(names(results) %in% exclude == FALSE)]

attr(results_short, "pipeline_version") <- pipeline_version
saveRDS(results_short, file = "results_short.rds")

output_file = "results_short.tsv"
header_message <- paste0("# Pipeline version ", pipeline_version, ". \n")
cat(header_message, file = output_file)
suppressWarnings(write.table(results_short,
                             file = output_file,
                             append = TRUE,
                             quote = FALSE,
                             sep = "\t",
                             row.names = FALSE))

attr(results, "pipeline_version") <- pipeline_version
saveRDS(results, file = "results_long.rds")

output_file = "results_long.tsv"
header_message <- paste0("# Pipeline version: ", pipeline_version, " \n")
cat(header_message, file = output_file)
suppressWarnings(write.table(results,
                             file = output_file,
                             append = TRUE,
                             quote = FALSE,
                             sep = "\t",
                             row.names = FALSE))