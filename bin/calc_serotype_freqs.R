library(optparse)
rm(list = ls())

args_list <- list(
  make_option(
    "--project_dir",
    type = "character",
    help = "Path to project directory."
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
  ),
  make_option(
    c("-c", "--downsampling_strategy"),
    type = "character",
    help = "Downsampling strategy, e.g. geodate."
  ),
  make_option(
    c("-m", "--minyear"),
    type = "integer",
    help = "Lowest year to include."
  ),
  make_option(
    c("-r", "--minyear_recent"),
    type = "integer",
    help = "Lowest year to include for recent time period."
  )
)

args_parser  <- OptionParser(option_list = args_list)

if (!interactive()) {
  args  <- parse_args(args_parser)
} else {
  args <- list(
    project_dir = "aci",
    file = "results_redacted/downsample_isolates/aci_crab_ds_geodate2.tsv",
    shortlist = "data/geographic_locations_in_study.tsv",
    downsampling_strategy = "geodate2",
    minyear = 2009,
    minyear_recent = 2016
  )
}

# create log file and start logging
con <- file(paste0("log.txt"))
sink(con, split = TRUE)

library(dplyr)
library(ggplot2)

library(devtools)
load_all(args$project_dir)

# import aci data
aci <- read_df(args$file) %>% 
  dplyr::filter(filtered & crab & downsampled & downsampled_by_pop)

aci$serotype <- gsub("-", "_", aci$serotype)

if (file.exists(args$shortlist)) {
# import short list of regions and countries that were manually selected after
# rarefaction
  shortlist <- read.csv(args$shortlist, sep = "\t") %>% 
    filter(shortlisted == TRUE)
  if (nrow(shortlist) == 0) rm("shortlist")
}

# filter to minimum year
aci <- aci[which(aci$collection_year >= args$minyear_recent), ]

saveRDS(
  list(minyear = args$minyear_recent, maxyear = max(aci$collection_year)),
  file = "meta.rds"
)

# check that all samples are recent, stop if not
testthat::expect_true(all(aci$collection_year >= args$minyear_recent))

# check that all samples are crab, stop if not
testthat::expect_true(all(aci$crab))

# REGIONS

region_continent_dict <- aci %>%
  select(region23, continent) %>%
  distinct()

serotop_region <- serotype_freqs(aci, group_by = "region23")

if (exists("shortlist")) {
# filter to regions that were manually shortlisted after rarefaction
shortlisted_regions <- shortlist$term_pretty[which(shortlist$variable == "region23")]
index <- which(serotop_region$region23 %in% shortlisted_regions)
serotop_region <- serotop_region[index,]
}

##### TODO The summary table is no longer used. Remove? #####

# WARNING: if serotop_region is not filtered to most prevalent serotypes then
# serotop_region_summary is not about whether serotypes are "prevalent" in a 
# region but whether they are "present".
serotop_region_summary <- serotop_region[, c("region23", "serotype")] %>%
  tidyr::pivot_wider(names_from = serotype, values_from = serotype)


for (j in 2:ncol(serotop_region_summary)) {
  serotop_region_summary[,j] <- as.numeric(!is.na(serotop_region_summary[,j]))
}

sums <- apply(serotop_region_summary[,2:ncol(serotop_region_summary)], 2, sum)

serotop_region_summary <- serotop_region_summary[,c(1, (order(sums, decreasing = TRUE)+1))]

#############################################################

serotop_region_exp <- serotop_region
serotop_region_exp$serotype <- gsub("_", "-", serotop_region_exp$serotype)

# export region23 top serotypes
write.table(
  serotop_region_exp,
  file = paste0("TableS2A_top_serotypes_region23_ds_", args$downsampling_strategy, ".tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
# export region23 top serotype summaries
write.table(
  serotop_region_summary,
  file = paste0("top_serotypes_summary_region23_ds_", args$downsampling_strategy, ".tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# DEFINE PREVALENT AND GLOBAL SEROTYPES AND EXPORT COLOR CODES

if (exists("shortlist")) {
  # filter to regions that were manually shortlisted after rarefaction
  shortlisted_regions <- shortlist$term_pretty[which(shortlist$variable == "region23")]
  index <- which(aci$region23 %in% shortlisted_regions)
  aci_shortlisted <- aci[index,]
}

global_prevalent <- global_prevalent(aci_shortlisted)

global_prevalent_exp <- global_prevalent
global_prevalent_exp$serotype <- gsub("_", "-", global_prevalent_exp$serotype)

write.table(
  global_prevalent_exp,
  file = paste0(
    "global_or_prevalent_serotypes_region23_ds_", args$downsampling_strategy, ".tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# DEFINE PREVALENT AND GLOBAL SEROTYPES FOR THE PREVIOUS TIME PERIOD

# beginning of the previous time period
minyear2 <- args$minyear

# import aci data
aci_prev <- read_df(args$file) %>% 
  dplyr::filter(filtered & crab & downsampled & downsampled_by_pop)

aci_prev$serotype <- gsub("-", "_", aci_prev$serotype)

# filter to time period
aci_prev <- aci_prev[which(
  aci_prev$collection_year >= minyear2 & 
  aci_prev$collection_year < args$minyear_recent
), ]
  
serotop_region_prev <- serotype_freqs(aci_prev, group_by = "region23")

if (exists("shortlist")) {
# filter to regions that were manually shortlisted after rarefaction
shortlisted_regions <- shortlist$term_pretty[which(shortlist$variable == "region23")]
index <- which(serotop_region_prev$region23 %in% shortlisted_regions)
serotop_region_prev <- serotop_region_prev[index,]
}

# DEFINE PREVALENT AND GLOBAL SEROTYPES FOR PREVIOUS TIME PERIOD

# Here the question is if we define global and prevalent serotypes based on the
# previous time period, how many of them are still global and prevalent in the
# current time period?

# prevalent_clones
# a serotype is considered prevalent if it has a frequency of at least 5% in at
# least 1 region
prevalent_clones_prev <- serotop_region_prev %>%
  filter(ratio > 0.05) %>%
  summarise(serotype = unique(serotype)) %>% 
  mutate(color = "orange3")

# global clones
# a serotype is considered a global if it is prevalent and has a frequency of at
# least 2% in at least 3 regions

min2perc_min3reg_prev <- serotop_region_prev %>% 
  filter(ratio > 0.02) %>% 
  group_by(serotype) %>% 
  summarise(n_region = length(unique(region23))) %>% 
  filter(n_region >= 3) %>% 
  subset(select = 1) %>% 
  mutate(color = "#DC143C")

global_index <- which(
  min2perc_min3reg_prev$serotype %in% prevalent_clones_prev$serotype)
global_clones_prev <- min2perc_min3reg_prev[global_index,]

prevalent_clones_prev <- anti_join(prevalent_clones_prev, global_clones_prev, by="serotype")

global_prevalent_prev <- dplyr::bind_rows(global_clones_prev, prevalent_clones_prev)

period1 <- global_prevalent %>% dplyr::rename(period = color)
period1_name <- paste0(
  args$minyear_recent, "-", max(aci$collection_year)
)
period1$period <- period1_name

period2 <- global_prevalent_prev %>% dplyr::rename(period = color)
period2_name <- paste0(
  minyear2, "-", args$minyear_recent-1
)
period2$period <- period2_name

period <- dplyr::bind_rows(period1, period2)
period <- tidyr::pivot_wider(
  period,
  id_cols = serotype,
  names_from = period,
  values_from = period
)
period[[period1_name]] <- !is.na(period[[period1_name]])
period[[period2_name]] <- !is.na(period[[period2_name]])
period$overlap <- period[[period1_name]] & period[[period2_name]]

period$serotype <- gsub("_", "-", period$serotype)

write.table(
  period,
  file = paste0(
    "global_or_prevalent_overlap_region23_ds_", args$downsampling_strategy, ".tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# Another question when we look at prevalent clones are they on the same
# continents in the two time periods?

prevalent_period1 <- serotop_region %>%
  filter(ratio > 0.05) %>%
  tidyr::pivot_wider(
    id_cols = serotype,
    names_from = region23,
    values_from = region23
  )
prevalent_period1[, 2: ncol(prevalent_period1)] <- apply(
  prevalent_period1[, 2: ncol(prevalent_period1)], 2, function(x) !is.na(x)
)

prevalent_period2 <- serotop_region_prev %>%
  filter(ratio > 0.05) %>%
  tidyr::pivot_wider(
    id_cols = serotype,
    names_from = region23,
    values_from = region23
  )
prevalent_period2[, 2: ncol(prevalent_period2)] <- apply(
  prevalent_period2[, 2: ncol(prevalent_period2)], 2, function(x) !is.na(x)
)

for (i in seq_along(names(prevalent_period1))) {
  p1name <- names(prevalent_period1)[i]
  if (p1name %in% names(prevalent_period2) == FALSE) {
    # TODO CHECK IF THIS IS MEANS WHAT IT SHOULD
    prevalent_period2[[p1name]] <- FALSE
  }
  index <- which(names(prevalent_period2) == p1name)
  if (i == 1) {
    prevalent_period2 <- dplyr::relocate(prevalent_period2, p1name)
  } else {
    prevalent_period2 <- dplyr::relocate(
      prevalent_period2,
      p1name,
      .after = names(prevalent_period1)[i-1])
  }
}

for (i in seq_along(names(prevalent_period2))) {
  p2name <- names(prevalent_period2)[i]
  if (p2name %in% names(prevalent_period1) == FALSE) {
    # TODO CHECK IF THIS IS MEANS WHAT IT SHOULD
    prevalent_period1[[p2name]] <- FALSE
  }
  index <- which(names(prevalent_period1) == p2name)
  if (i == 1) {
    prevalent_period1 <- dplyr::relocate(prevalent_period1, p2name)
  } else {
    prevalent_period1 <- dplyr::relocate(
      prevalent_period1,
      p2name,
      .after = names(prevalent_period2)[i-1])
  }
}

# test that variable names are in the same order
testthat::expect_true(all(names(prevalent_period1) == names(prevalent_period2)))

prevres <- matrix(
  NA_character_, nrow = nrow(prevalent_period1), ncol = ncol(prevalent_period1))
prevres[,1] <- prevalent_period1$serotype
prevres <- as.data.frame(prevres)
names(prevres) <- names(prevalent_period1)

get_freq_change <- function() {
  # get ratio of serotype in region for previous time period
  index_ratio_old <- which(
    serotop_region_prev$serotype == serotype &
      serotop_region_prev$region23 == names(prevres)[j])
  if (length(index_ratio_old) == 0) {
    ratio_old <- 0
  } else {
    ratio_old <- round(serotop_region_prev$ratio[index_ratio_old], 3)
  }
  # get ratio of serotype in region for latest time period
  index_ratio_new <- which(
    serotop_region$serotype == serotype &
      serotop_region$region23 == names(prevres)[j])
  if (length(index_ratio_new) == 0) {
    ratio_new <- 0
  } else {
    ratio_new <- round(serotop_region$ratio[index_ratio_new], 3)
  }
  ratio_diff <- round(ratio_new - ratio_old, 3)
  ratio_diff_sign <- ifelse(sign(ratio_diff) == 1, "+", "")
  freq_change <- paste(
    ratio_old, "/", ratio_new, "(", ratio_diff_sign, ratio_diff, ")",
    sep = "",
    collapse = ""
  )
  return(freq_change)
}

for (i in 1:nrow(prevres)) {
  serotype <- prevres$serotype[i]
  index_x <- which(prevalent_period2$serotype == serotype)
  # if serotype cannot be found in previous period, each TRUE entry is emergence
  if (length(index_x) == 0) {
    for (j in 2:ncol(prevres)) {
      p1true <- as.logical(prevalent_period1[i,j])
      if (p1true) {
        prevres[i,j] <- get_freq_change()
      } else {
        prevres[i,j] <- NA_character_ 
      }
    }
  } else {
    for (j in 2:ncol(prevres)) {
      p1true <- as.logical(prevalent_period1[i,j])
      p2true <- as.logical(prevalent_period2[index_x, j])
      if (p1true & p2true) {
        prevres[i,j] <- get_freq_change()
      }
      if (p1true == FALSE & p2true == FALSE) {
        prevres[i,j] <- NA_character_
      }
      if (p1true == FALSE & p2true) {
        prevres[i,j] <- get_freq_change()
      }
      if (p1true & p2true == FALSE) {
        prevres[i,j] <- get_freq_change()
      }
    }
  }
}

prevres$serotype = gsub("_", "-", prevres$serotype)

write.table(
  prevres,
  file = paste0(
    "TableS2_prevalent_overlap_region23_ds_", args$downsampling_strategy, ".tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE,
  na = ""
)

# COUNTRIES

serotop_country <- serotype_freqs(aci, group_by = "country")

if (exists("shortlist")) {
# filter to countries that were manually shortlisted after rarefaction
shortlisted_countries <- shortlist$term_pretty[which(shortlist$variable == "country")]
index <- which(serotop_country$country %in% shortlisted_countries)
serotop_country <- serotop_country[index,]
}

# export country top serotypes
write.table(
  serotop_country,
  file = paste0("top_serotypes_country_ds_", args$downsampling_strategy, ".tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# pairwise comparisons - bray curtis dissimilarity and overlap prevalence

if (exists("shortlisted_countries")) {
pairs <- combn(shortlisted_countries, 2) %>% t() %>% as.data.frame()
names(pairs) <- c("country1", "country2")
pairs$country1_continent <- sapply(pairs$country1, function(x) {
  countrycode::countrycode(
    x,
    origin = "country.name",
    destination = "continent"
  ) 
}) %>% tolower() %>% gsub(" +", "_", .)
pairs$country2_continent <- sapply(pairs$country2, function(x) {
  countrycode::countrycode(
    x,
    origin = "country.name",
    destination = "continent"
  ) 
}) %>% tolower() %>% gsub(" +", "_", .)

# calculate bray-curtis dissimilarity

bray <- serotop_country  %>% 
  group_by(country, serotype) %>% 
  summarise(count=sum(count))

diversity_matrix <- tidyr::pivot_wider(
  bray, 
  names_from=serotype, 
  values_from=count
  ) %>%
  tibble::column_to_rownames("country") %>%
  replace(is.na(.), 0)

beta_diversity <- vegan::vegdist(diversity_matrix, "bray")

beta_diversity_long <- reshape2::melt(as.matrix(beta_diversity))

# merge bray-curtis results with pairwise comparison table

pairs$bray_curtis <- NA
for (i in 1:nrow(pairs)) {
  index <- which(
    beta_diversity_long$Var1 == pairs$country1[i] &
    beta_diversity_long$Var2 == pairs$country2[i]
  )
  if (length(index) == 1) {
    pairs$bray_curtis[i] <- beta_diversity_long$value[index]
  }
}
# Stop with an error message if any of the values are still NA.
if (any(is.na(pairs$bray_curtis))) stop("Check Bray-Curtis.")

pairs$morisita <- NA
for (i in 1:nrow(pairs)) {
  index1 <- which(rownames(diversity_matrix) == pairs$country1[i])
  index2 <- which(rownames(diversity_matrix) == pairs$country2[i])
  pairs$morisita[i] <- abdiv::horn_morisita(
    diversity_matrix[index1, ],
    diversity_matrix[index2, ]
  )
}
# Stop with an error message if any of the values are still NA.
if (any(is.na(pairs$morisita))) stop("Check Morisita.")

# calculate prevalence of overlapping serotypes
# method:
# for each pair, find serotypes that are present in both
# calculate overall prevalence of overlapping serotypes in country 1
# calculare overall prevalence of overlapping serotypes in country 2
# calculate the mean of the two values

pairs$country1_overlap_prevalence <- NA
pairs$country2_overlap_prevalence <- NA
pairs$mean_overlap_prevalence <- NA
for (i in 1:nrow(pairs)) {
  country1 <- pairs$country1[i]
  country2 <- pairs$country2[i]
  # find serotypes that are present in both countries
  serotypes1 <- serotop_country$serotype[which(serotop_country$country == country1)]
  serotypes2 <- serotop_country$serotype[which(serotop_country$country == country2)]
  overlap <- intersect(serotypes1, serotypes2)
  # calculate overall prevalence of overlapping serotypes in country 1
  index_country1_overlap <- which(
    serotop_country$country == country1 & serotop_country$serotype %in% overlap)
  pairs$country1_overlap_prevalence[i] <- sum(serotop_country$ratio[index_country1_overlap])
  # calculare overall prevalence of overlapping serotypes in country 2
  index_country2_overlap <- which(
    serotop_country$country == country2 & serotop_country$serotype %in% overlap)
  pairs$country2_overlap_prevalence[i] <- sum(serotop_country$ratio[index_country2_overlap])
}
# calculate the mean of the two values
pairs$mean_overlap_prevalence <- apply(pairs[,c("country1_overlap_prevalence", "country2_overlap_prevalence")], 1, mean)
} else {
  pairs <- data.frame()
}

write.table(
  pairs,
  file = paste0("country_comparisons_ds_", args$downsampling_strategy, ".tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# REGIONS BY YEAR

# import aci data
aci <- read_df(args$file) %>% 
  dplyr::filter(filtered & crab & downsampled & downsampled_by_pop)

serotop_region_year <- data.frame()
for (k in sort(unique(aci$collection_year))) {
  for (i in unique(aci$region23)) {
  aci_small <- aci[which(aci$region23 == i & aci$collection_year == k),]
  if (nrow(aci_small) > 0) {
    df <- aci_small %>%
      group_by(serotype) %>%
      summarise(count = length(assembly))
    df$ratio <- signif(df$count/sum(df$count), 4)
    df <- df[order(df$ratio, decreasing = TRUE),]
    df$region23 <- i
    df$collection_year <- k
    df <- dplyr::relocate(df, region23, collection_year)
    serotop_region_year <- dplyr::bind_rows(
      serotop_region_year,
      df
    )
  }
  
  }
}

if (exists("shortlist")) {
# filter to regions that were manually shortlisted after rarefaction
shortlisted_regions <- shortlist$term_pretty[which(shortlist$variable == "region23")]
index <- which(serotop_region_year$region23 %in% shortlisted_regions)
serotop_region_year <- serotop_region_year[index,]
}

write.csv(
  serotop_region_year,
  file = paste0(
    "serotypes_region23_year_ds_", args$downsampling_strategy, ".tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# COUNTRIES BY YEAR

# import aci data
aci <- read_df(args$file) %>% 
  dplyr::filter(filtered & crab & downsampled & downsampled_by_pop)

serotop_country_year <- data.frame()
for (k in sort(unique(aci$collection_year))) {
  for (i in unique(aci$country)) {
    aci_small <- aci[which(aci$country == i & aci$collection_year == k),]
    if (nrow(aci_small) > 0) {
      df <- aci_small %>%
        group_by(serotype) %>%
        summarise(count = length(assembly))
      df$ratio <- signif(df$count/sum(df$count), 4)
      df <- df[order(df$ratio, decreasing = TRUE),]
      df$country <- i
      df$collection_year <- k
      df <- dplyr::relocate(df, country, collection_year)
      serotop_country_year <- dplyr::bind_rows(
        serotop_country_year,
        df
      )
    }
  }
}

if( exists("shortlist")) {
# filter to countries that were manually shortlisted after rarefaction
shortlisted_countries <- shortlist$term_pretty[which(shortlist$variable == "country")]
index <- which(serotop_country_year$country %in% shortlisted_countries)
serotop_country_year <- serotop_country_year[index,]
}

write.csv(
  serotop_country_year,
  file = paste0(
    "serotypes_country_year_ds_", args$downsampling_strategy, ".tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
                  
# end logging
sink(con)
