# This script runs a number of QC checks, adds a "qc_pass" variable to the 
# input data frame and exports other files as well. See ?check_qc() for more
# information.

library(optparse)
rm(list = ls())

args_list <- list(
  make_option(
    c("-p", "--project_dir"),
    type = "character",
    help = "Path to project directory which contains the pipeline scripts."
  ),
  make_option(
    c("-f", "--file"),
    type = "character",
    help = "Path to aci prediction results."
  ),
  make_option(
    c("-t", "--taxid"),
    type = "character",
    help = "The NCBI Taxonomy ID of the organism."
  )
)

args_parser  <- OptionParser(option_list = args_list)

if (!interactive()) {
  args  <- parse_args(args_parser)
} else {
  args <- list(
    project_dir = "aci",
    file = "data/typing_summary_tables/TableS1_aci_study_redacted_renamed2.rds",
    taxid= 573
  )
}

# create log file and start logging
con <- file("log.txt", open = "wt")
sink(con, type = "output", split = TRUE)
sink(con, type = "message") # comment when debugging

library(dplyr)
library(devtools)
load_all(args$project_dir)

# import data set
aci <- read_df(args$file)

# run qc checks

qc <- check_qc(aci, taxid = args$taxid, new_na = "ignore", verbose = TRUE)

# export data frame with QC results

saveRDS(qc$data, file = "aci_with_qc.rds")

write_tsv(qc$data, file = "aci_with_qc.tsv")

# export data frame with QC failures

saveRDS(qc$failed, file = "aci_qc_failed.rds")

write_tsv(qc$failed, file = "aci_qc_failed.tsv")

sink(type = "output")
sink(type = "message")
close(con)
