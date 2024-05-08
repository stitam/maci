# This script filters the assemblies:
# - only include rows where QC passes
# - only include rows where continent is known
# - only include rows where region is known
# - only include rows where country is known
# - only include rows where collection year is known
#
# This script also adds new temporal variables:
# - yearweek
# - yearmonth
#
# This script also adds a new variable:
# - serotype

library(optparse)
rm(list = ls())

args_list <- list(
  make_option(
    c("-p", "--project_dir"),
    type = "character",
    help = "Path to the project directory."
  ),
  make_option(
    c("-f", "--file"),
    type = "character",
    help = "Path to aci prediction results."
  ),
  make_option(
    c("-s", "--shortlist"),
    type = "character",
    help = "Path to shortlisted regions and countries."
  )
)

args_parser  <- OptionParser(option_list = args_list)

if (!interactive()) {
  args  <- parse_args(args_parser)
} else {
  args <- list(
    project_dir = "~/Methods/aci",
    file = "aci_with_qc.rds",
    shortlist = "geographic_locations_in_study.tsv"
  )
}

# create log file and start logging
con <- file(paste0("log.txt"))
sink(con, split = TRUE)

library(dplyr)
library(devtools)
load_all(args$project_dir)

# import data set
aci <- readRDS(args$file)

# filter by qc
aci <- aci[which(aci$qc_pass == TRUE), ]

### EDIT
# Do not filter by resistance here, because some analysis downstream may require
# both crab and non crab isolates. Instead filter in each downstream step.

# filter by resistance
# aci <- aci[which(aci$crab == TRUE), ]
###

# filter by environment, only include clinical isolates

print("Only include clinical isolates.")

if (any(aci$human_related == FALSE)) {
  aci <- aci[-which(aci$human_related == FALSE), ]
}

print("Number of rows:")
print(nrow(aci))

# filter by geographic locations

print("Only include rows where continent is known.")

aci <- aci[which(!is.na(aci$continent)), ]

print("Number of rows:")
print(nrow(aci))

print("Only include rows where region is known.")

aci <- aci[which(!is.na(aci$region23)), ]

print("Number of rows:")
print(nrow(aci))

print("Only include rows where country is known.")

aci <- aci[which(!is.na(aci$country)), ]

print("Number of rows:")
print(nrow(aci))

# update temporal variables
set.seed(0)
# if date is complete keep it
# if date is incomplete, simulate complete date between boundaries
aci$collection_day <- aci$collection_date %>%
  date_runif() %>%
  lubridate::date_decimal() %>%
  as.Date()

# collection year is the year of the (sometimes simulated) collection day
aci$collection_year <- as.integer(lubridate::year(aci$collection_day))

# filter by collection_year

print("Only include rows where collecion year is known.")

aci <- aci[which(!is.na(aci$collection_year)), ]

print("Number of rows:")
print(nrow(aci))

# add combined serotype
aci$serotype <- paste(aci$mlst, aci$k_serotype, sep = "_")

# eliminate VUB duplicates
# It seems some VUB samples were sequenced more twice. Both sequences were 
# uploaded to NCBI as separate assemblies with separate biosample metadata and
# the metadata are sometimes not identical which produces errors downstream.
# Approach to handle this:
# - When merging tables using scripts/merge_prediction_results.R, look at
# `collection_date` for each VUB strain; if one is NA and the other is not,
# replace NA with the other.
# - Here, if both samples would be kept after all the filtering, keep the one
# with the lowest number of contigs.

vub_strains <- unique(aci$strain[grep("AB[0-9]+-VUB", aci$strain)])

for (i in vub_strains) {
  index <- which(aci$strain == i)
  if (length(index) == 1) {
    next()
  } else {
    index_keep <- which(
      aci$contig_count[index] == min(aci$contig_count[index])
    )[1]
    index_remove <- index[-index_keep]
    aci <- aci[-index_remove, ]
  }
}

# check that all VUB duplicates were removed
for (i in vub_strains) {
  testthat::expect_equal(length(which(aci$strain == i)), 1)
}

shortlist <- read.csv(args$shortlist, sep = "\t", na = "")
shortlist_country <- shortlist %>% dplyr::filter(variable == "country")
shortlist_region <- shortlist %>% dplyr::filter(variable == "region23")

# check that all countries and regions are in the shortlist
testthat::expect_true(all(aci$country %in% shortlist_country$term))
testthat::expect_true(all(aci$region23 %in% shortlist_region$term))

# check that all countries and regions in the shortlist are in the db.
testthat::expect_true(all(shortlist_country$term %in% aci$country))
testthat::expect_true(all(shortlist_region$term %in% aci$region23))

# export data set
saveRDS(aci, "aci_filtered.rds")

write.table(
  aci,
  file = "aci_filtered.tsv",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# end logging
sink(con)
