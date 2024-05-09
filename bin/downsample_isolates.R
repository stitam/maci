library(devtools)
library(dplyr)
library(optparse)
rm(list = ls())

args_list <- list(
  make_option(
    c("-p", "--project_dir"),
    type = "character",
    help = "Path to pipeline directory."
  ),
  make_option(
    c("-a", "--aci_path"),
    type = "character",
    help = "Path to aci prediction results."
  ),
  make_option(
   c("-o", "--pp_path"),
   type = "character",
   help = "Path to PopPUNK prediction results."
  ),
  make_option(
    c("-c", "--downsampling_strategy"),
    type = "character",
    help = "Downsampling strategy, e.g. geodate."
  ),
  make_option(
    "--geographic_locations",
    type = "character",
    help = "Path to geographic locations."
  ),
  make_option(
    "--population_sampling_rate",
    type = "numeric",
    help = "Number of samples to draw per million inhabitants in each country."
  )
)

args_parser  <- OptionParser(option_list = args_list)

if (!interactive()) {
  args  <- parse_args(args_parser)
} else {
  args <- list(
    project_dir = "aci",
    aci_path = "results_redacted/filter_assemblies/aci_filtered.rds",
    pp_path = "ppdb_clusters.csv",
    downsampling_strategy = "geodate2",
    geographic_locations = "data/geographic_locations_in_study.tsv",
    population_sampling_rate = 1
  )
}

# import custom functions
load_all(args$project_dir)

# import data set
aci_all <- readRDS(args$aci_path)

# consistency checks

testthat::expect_equal(class(aci_all$filtered), "logical")
testthat::expect_equal(class(aci_all$crab), "logical")

# remove column if it already exists
if ("downsampled" %in% names(aci_all)) {
  aci_all <- dplyr::select(aci_all, -downsampled)
}

# downsample
if (args$downsampling_strategy == "none") {
  aci_all$downsampled <- TRUE
}

# filter
acif <- aci_all %>% dplyr::filter(filtered) %>% dplyr::filter(crab)

# set random seed
set.seed(0)

# downsample
if (args$downsampling_strategy %in% c("geodate", "geodate2", "geodate3")) {
  acid <- downsample(acif, strategy = args$downsampling_strategy)
}

if (args$downsampling_strategy == "poppunk") {
  pp <- read.csv(args$pp_path)
  pp <- dplyr::rename(pp, assembly = Taxon, pp = Cluster)
  acid <- downsample(acif, strategy = "poppunk", pp = pp)
}

# create "downsampled" variable

aci_all$downsampled <- sapply(aci_all$assembly, function(x) {
  if (x %in% acid$assembly) {
    return(TRUE)
  } else {
    return(FALSE)
  }
})

# SUBSEQUENTLY DOWNSAMPLE BY POPULATION

# prep/inport input data
acid <- aci_all %>% dplyr::filter(filtered & crab & downsampled)

population <- read_df(args$geographic_locations) %>%
  dplyr::filter(variable == "country") %>%
  dplyr::rename(country = term_pretty) %>%
  dplyr::select(country, population)

# reset seed
set.seed(0)

# downsample
acid_pop <- downsample_population(
  acid,
  population = population,
  rate = as.numeric(args$population_sampling_rate)
)

# create "downsampled_by_pop" variable
aci_all$downsampled_by_pop <- sapply(aci_all$assembly, function(x) {
  if (x %in% acid_pop$assembly) {
    return(TRUE)
  } else {
    return(FALSE)
  }
})

# export
if (!interactive()) {
  write_tsv(
    aci_all,
    file = paste0("aci_crab_ds_", args$downsampling_strategy, ".tsv")
  )
}

# prepare table of country county by filter type

count_all <- aci_all %>%
  dplyr::group_by(country) %>%
  dplyr::summarise(raw = n())

count_filtered <- aci_all %>%
  dplyr::filter(filtered) %>%
  dplyr::group_by(country) %>%
  dplyr::summarise(filtered = n())

count_crab <- aci_all %>%
  dplyr::filter(filtered & crab) %>%
  dplyr::group_by(country) %>%
  dplyr::summarise(crab = n())

count_ds <- aci_all %>%
  dplyr::filter(filtered & crab & downsampled) %>%
  dplyr::group_by(country) %>%
  dplyr::summarise(downsampled = n())

count_ds_pop <- aci_all %>%
  dplyr::filter(filtered & crab & downsampled & downsampled_by_pop) %>%
  dplyr::group_by(country) %>%
  dplyr::summarise(downsampled_by_pop = n())

counts <- count_all %>%
  dplyr::left_join(count_filtered, by = "country") %>%
  dplyr::left_join(count_crab, by = "country") %>%
  dplyr::left_join(count_ds, by = "country") %>%
  dplyr::left_join(count_ds_pop, by = "country") %>%
  dplyr::arrange(country)

for (i in 1:nrow(counts)) {
  for (j in 2:ncol(counts)) {
    if (is.na(counts[i,j])) {
      counts[i,j] <- 0
    }
  }
}

# export
if (!interactive()) {
  write.table(
    counts,
    file = paste0("country_counts_", args$downsampling_strategy, ".tsv"),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
}
