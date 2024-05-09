library(optparse)
rm(list = ls())

args_list <- list(
  make_option(
    "--project_dir",
    type = "character",
    help = "Path to project directory."
  ),
  make_option(
    "--file",
    type = "character",
    help = "Path to a file which contains typing results for each sample."
  ),
  make_option(
    "--strategy",
    type = "character",
    help = "Downsampling strategy"
  ),
  make_option(
    "--minyear",
    type = "integer",
    help = "Lowest year to include."
  )
)

args_parser  <- OptionParser(option_list = args_list)

if (!interactive()) {
  args  <- parse_args(args_parser)
} else {
  args <- list(
    project_dir = "aci",
    file = "results_redacted/downsample_isolates/aci_crab_ds_geodate.tsv",
    strategy = "geodate",
    minyear = 2016
  )
}

if (!interactive()) {
  # create log file and start logging
  con <- file(paste0("log_", args$strategy, ".txt"))
  sink(con, split = TRUE)
}

library(devtools)
library(dplyr)
library(ggplot2)

load_all(args$project_dir)

# import data
aci <- read_df(args$file) %>% 
  dplyr::filter(filtered & crab & downsampled & downsampled_by_pop)

# filter to minimum year
aci <- aci[which(aci$collection_year >= args$minyear), ]

print("Collection year for earliest sample:")
print(min(aci$collection_year))

print("Are all samples CRAB?")
print(all(aci$crab))

continent <- prep_serotop_tbl(df = aci, type = "serotype", group_by = "continent")
region23 <- prep_serotop_tbl(df = aci, type = "serotype", group_by = "region23")
country <- prep_serotop_tbl(df = aci, type = "serotype", group_by = "country")

# set random seed for jitter
set.seed(0)

index <- which(continent$continent != "all")
continent$quantile[index] <- jitter(continent$quantile[index])

index <- which(region23$region23 != "all")
region23$quantile[index] <- jitter(region23$quantile[index])

index <- which(country$country != "all")
country$quantile[index] <- jitter(country$quantile[index])

if (!interactive()) {
  # export continent summary
  write.table(
    continent,
    file = paste0("serotop_continent_ds_", args$strategy, ".tsv"),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
  # export region summary
  write.table(
    region23,
    file = paste0("serotop_region23_ds_", args$strategy, ".tsv"),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
  # export country summary
  write.table(
    country,
    file = paste0("serotop_country_ds_", args$strategy, ".tsv"),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
}

if (!interactive()) {
  # end logging
  sink(con)
}
