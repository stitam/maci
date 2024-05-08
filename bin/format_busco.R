# This script is used to format the BUSCO output into a table. BUSCO returns a
# JSON file with the results, which is parsed and converted into a table. The
# table is then written to a file. The script is currently not used in the
# pipeline.

rm(list = ls())

if (!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  busco_dir <- args[1]
  busco_db_ver <- args[2]
} else {
  test_dir <- "~/Methods/aci-tests/test-format_busco"
  sample_id <- "Aci370"
  busco_results_path <- paste0(test_dir, "/Aci370")
  busco_db_ver <- "pseudomonadales_odb10"
}

resfile <- paste0(
  busco_dir,
  "/short_summary.specific.",
  busco_db_ver,
  ".",
  busco_dir,
  ".json"
)

cont <- jsonlite::fromJSON(resfile)

res <- data.frame(
  busco_complete = cont$results$Complete,
  busco_single_copy = cont$results$`Single copy`,
  busco_multi_copy = cont$results$`Multi copy`,
  busco_fragmented = cont$results$Fragmented,
  busco_missing = cont$results$Missing,
  busco_n_markers = cont$results$n_markers
)

write.table(
  res,
  file = paste0(sample_id, "_busco.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
