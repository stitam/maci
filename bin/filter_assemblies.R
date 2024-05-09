# This script adds new variables:
# - collection_day
# - collection_year
# - serotype
# - country_iso2c
# - region23
# - continent
#
# This script filters assemblies:
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
    project_dir = "aci",
    file = "results_redacted/run_qc_checks/aci_with_qc.rds",
    shortlist = "data/geographic_locations_in_study.tsv"
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

# REMOVE SOME VARIABLES IF THEY EXIST

if ("filtered" %in% names(aci)) {
  aci <- dplyr::select(aci, -filtered)
}

if ("continent" %in% names(aci)) {
  aci <- dplyr::select(aci, -continent)
}

if ("region23" %in% names(aci)) {
  aci <- dplyr::select(aci, -region23)
}

if ("country_iso2c" %in% names(aci)) {
  aci <- dplyr::select(aci, -country_iso2c)
}

if ("collection_day" %in% names(aci)) {
  aci <- dplyr::select(aci, -collection_day)
}

if ("collection_year" %in% names(aci)) {
  aci <- dplyr::select(aci, -collection_year)
}

# ADD SOME NEW VARIABLES

aci <- validate_var(
  df = aci,
  varname = "country",
  varclass = "character",
  coerce_fun = as.character,
  verbose = TRUE
)

# Query continents
aci$continent <- countrycode::countrycode(
  aci$country, 
  origin = "country.name",
  destination = "continent"
)

# Query regions
aci$region23 <- countrycode::countrycode(
  aci$country, 
  origin = "country.name",
  destination = "region23"
)

aci$country_iso2c <- countrycode::countrycode(
  aci$country, 
  origin = "country.name",
  destination = c("iso2c", "ecb")
)

# Add non ISO standard countries
if (any(aci$country == "Kosovo", na.rm = TRUE)) {
  index <- which(aci$country == "Kosovo")
  aci$country_iso2c[index] <- "XK"
  aci$continent[index] <- "Europe"
  aci$region23[index] <- "Southern Europe"
}

# add combined serotype

aci <- validate_var(
  df = aci,
  varname = "mlst",
  varclass = "character",
  coerce_fun = as.character,
  verbose = TRUE
)

aci <- validate_var(
  df = aci,
  varname = "k_serotype",
  varclass = "character",
  coerce_fun = as.character,
  verbose = TRUE
)

aci$serotype <- paste(aci$mlst, aci$k_serotype, sep = "-")

# FILTER

aci_all <- aci

# filter by qc
aci <- aci[which(aci$qc_pass == TRUE), ]

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

# ADD NEW TEMPORAL VARIABLE

# Note: this could be moved to the beginning but that would break previous
# calculations because of the random nature of the process. Keep this here
# for now and maybe tidy it for a new project.

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

aci <- validate_var(
  df = aci,
  varname = "qc_pass",
  varclass = "logical",
  verbose = TRUE
)

acif <- aci %>%
  dplyr::filter(!is.na(qc_pass) & qc_pass) %>%
  dplyr::filter(is.na(human_related) | human_related) %>% 
  dplyr::filter(!is.na(continent)) %>%
  dplyr::filter(!is.na(region23)) %>%
  dplyr::filter(!is.na(country)) %>%
  dplyr::filter(!is.na(collection_year))

# TODO move this assertion dowstream it should not break filtering
if (file.exists(args$shortlist)) {
  shortlist <- read.csv(args$shortlist, sep = "\t", na = "")
  shortlist_country <- shortlist %>% dplyr::filter(variable == "country")
  shortlist_region <- shortlist %>% dplyr::filter(variable == "region23")
  
  # check that all countries and regions are in the shortlist
  testthat::expect_true(all(acif$country %in% shortlist_country$term_pretty))
  testthat::expect_true(all(acif$region23 %in% shortlist_region$term_pretty))
  
  # TODO is this necessary? Will it break anything? 
  # check that all countries and regions in the shortlist are in the db.
  # testthat::expect_true(all(shortlist_country$term %in% acif$country))
  # testthat::expect_true(all(shortlist_region$term %in% acif$region23))
} else {
  print("No table for shortlisted geographic locations found.")
}

set.seed(0)

aci_all$collection_day <- as.Date("1500-01-01")
aci_all$collection_year <- NA_integer_
aci_all$filtered <- NA

for (i in 1:nrow(aci_all)) {
  index <- which(acif$assembly == aci_all$assembly[i])
  if (length(index) == 0) {
    aci_all$collection_day[i] <- date_runif(aci_all$collection_date[i]) %>%
      lubridate::date_decimal() %>%
      as.Date()
    aci_all$collection_year[i] <- as.integer(
      lubridate::year(aci_all$collection_day[i]))
    aci_all$filtered[i] <- FALSE
  } else if (length(index) > 1) {
    stop()
  } else {
    aci_all$collection_day[i] <- acif$collection_day[index]
    aci_all$collection_year[i] <- acif$collection_year[index]
    aci_all$filtered[i] <- TRUE
  }
}

testthat::expect_false(min(aci_all$collection_day, na.rm = TRUE) == as.Date("1500-01-01"))

# export data set
saveRDS(aci_all, "aci_filtered.rds")

write.table(
  aci_all,
  file = "aci_filtered.tsv",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# end logging
sink(con)
