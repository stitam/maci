library(optparse)
rm(list = ls())

args_list <- list(
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
    c("-c", "--collapse_strategy"),
    type = "character",
    help = "Collapse strategy, e.g. geodate."
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
    file = "results/filter_crab/aci_collapse_geodate_crab.rds",
    shortlist = "aci/data/geographic_locations_in_study.tsv",
    collapse_strategy = "geodate",
    minyear = 2009,
    minyear_recent = 2016
  )
}

# create log file and start logging
con <- file(paste0("log.txt"))
sink(con, split = TRUE)

library(dplyr)
library(ggplot2)

# import aci data
aci <- readRDS(args$file)

# import short list of regions and countries that were manually selected after
# rarefaction
shortlist <- read.csv(args$shortlist, sep = "\t") %>% filter(shortlisted == TRUE)

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

serotop_region <- data.frame()
for (i in unique(aci$region23)) {
  aci_small <- aci[which(aci$region23 == i),]
  df <- aci_small %>%
    group_by(serotype) %>%
    summarise(count = length(assembly))
  df$ratio = signif(df$count/sum(df$count), 4)
  df <- df[order(df$ratio, decreasing = TRUE),]
  df$ratio_cumsum <- signif(cumsum(df$ratio), 4)
  df$order <- 1:nrow(df)
  df$region23 <- i
  df <- dplyr::relocate(df, region23)
  serotop_region = dplyr::bind_rows(
    serotop_region,
    df
  )
}

# filter to regions that were manually shortlisted after rarefaction
shortlisted_regions <- shortlist$term[which(shortlist$variable == "region23")]
index <- which(serotop_region$region23 %in% shortlisted_regions)
serotop_region <- serotop_region[index,]

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

# export region23 top serotypes
write.table(
  serotop_region,
  file = paste0("TableS2A_top_serotypes_region23_collapse_", args$collapse_strategy, ".tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
# export region23 top serotype summaries
write.table(
  serotop_region_summary,
  file = paste0("top_serotypes_summary_region23_collapse_", args$collapse_strategy, ".tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# DEFINE PREVALENT AND GLOBAL SEROTYPES AND EXPORT COLOR CODES

serotop_region2 <- dplyr::left_join(
  serotop_region,
  region_continent_dict,
  by = "region23"
)

# prevalent_clones
# a serotype is considered prevalent if it has a frequency of at least 5% in at
# least 1 region
prevalent_clones <- serotop_region2 %>%
  filter(ratio > 0.05) %>%
  summarise(serotype = unique(serotype)) %>% 
  mutate(color = "orange3")

# global clones
# a serotype is considered a global if it is prevalent and has a frequency of at
# least 2% in at least 3 regions and at least 2 continents

min2perc_min3reg <- serotop_region2 %>% 
  filter(ratio > 0.02) %>% 
  group_by(serotype) %>% 
  summarise(
    n_region = length(unique(region23)),
    n_continent = length(unique(continent))
  ) %>% 
  filter(n_region >= 3 & n_continent >= 2) %>% 
  subset(select = 1) %>% 
  mutate(color = "#DC143C")

global_index <- which(min2perc_min3reg$serotype %in% prevalent_clones$serotype)
global_clones <- min2perc_min3reg[global_index,]

# combine the two tables
prevalent_clones <- anti_join(prevalent_clones, global_clones, by="serotype")
global_prevalent <- dplyr::bind_rows(global_clones, prevalent_clones)

write.table(
  global_prevalent,
  file = paste0(
    "global_or_prevalent_serotypes_region23_collapse_", args$collapse_strategy, ".tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# DEFINE PREVALENT AND GLOBAL SEROTYPES FOR THE PREVIOUS TIME PERIOD

# beginning of the previous time period
minyear2 <- args$minyear

# import aci data
aci_prev <- readRDS(args$file)

# filter to time period
aci_prev <- aci_prev[which(
  aci_prev$collection_year >= minyear2 & 
  aci_prev$collection_year < args$minyear_recent
), ]

serotop_region_prev <- data.frame()
for (i in unique(aci_prev$region23)) {
  aci_small <- aci_prev[which(aci_prev$region23 == i),]
  df <- aci_small %>%
    group_by(serotype) %>%
    summarise(count = length(assembly))
  df$ratio = signif(df$count/sum(df$count), 4)
  df <- df[order(df$ratio, decreasing = TRUE),]
  df$ratio_cumsum <- signif(cumsum(df$ratio), 4)
  df$order <- 1:nrow(df)
  df$region23 <- i
  df <- dplyr::relocate(df, region23)
  serotop_region_prev <- dplyr::bind_rows(
    serotop_region_prev,
    df
  )
}

# filter to regions that were manually shortlisted after rarefaction
shortlisted_regions <- shortlist$term[which(shortlist$variable == "region23")]
index <- which(serotop_region_prev$region23 %in% shortlisted_regions)
serotop_region_prev <- serotop_region_prev[index,]

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

write.table(
  period,
  file = paste0(
    "global_or_prevalent_overlap_region23_collapse_", args$collapse_strategy, ".tsv"),
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

index <- sapply(names(prevalent_period1), function(x) which(names(prevalent_period2) == x))
prevalent_period2 <- prevalent_period2[,index]

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

write.table(
  prevres,
  file = paste0(
    "TableS2_prevalent_overlap_region23_collapse_", args$collapse_strategy, ".tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE,
  na = ""
)

# COUNTRIES

serotop_country <- data.frame()
for (i in unique(aci$country)) {
  aci_small <- aci[which(aci$country == i),]
  df <- aci_small %>%
    group_by(serotype) %>%
    summarise(count = length(assembly))
  df$ratio = signif(df$count/sum(df$count), 4)
  df <- df[order(df$ratio, decreasing = TRUE),]
  df$ratio_cumsum <- signif(cumsum(df$ratio), 4)
  df$order <- 1:nrow(df)
  df$country <- i
  df <- dplyr::relocate(df, country)
  serotop_country = dplyr::bind_rows(
    serotop_country,
    df
  )
}

# filter to countries that were manually shortlisted after rarefaction
shortlisted_countries <- shortlist$term[which(shortlist$variable == "country")]
index <- which(serotop_country$country %in% shortlisted_countries)
serotop_country <- serotop_country[index,]

# export country top serotypes
write.table(
  serotop_country,
  file = paste0("top_serotypes_country_collapse_", args$collapse_strategy, ".tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# pairwise comparisons - bray curtis dissimilarity and overlap prevalence

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
  pairs$morisita[i] <- abdiv::morisita(
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

write.table(
  pairs,
  file = paste0("country_comparisons_collapse_", args$collapse_strategy, ".tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# REGIONS BY YEAR

# import aci data
aci <- readRDS(args$file)

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
# filter to regions that were manually shortlisted after rarefaction
shortlisted_regions <- shortlist$term[which(shortlist$variable == "region23")]
index <- which(serotop_region_year$region23 %in% shortlisted_regions)
serotop_region_year <- serotop_region_year[index,]

write.csv(
  serotop_region_year,
  file = paste0(
    "serotypes_region23_year_collapse_", args$collapse_strategy, ".tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# COUNTRIES BY YEAR

# import aci data
aci <- readRDS(args$file)

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

# filter to countries that were manually shortlisted after rarefaction
shortlisted_countries <- shortlist$term[which(shortlist$variable == "country")]
index <- which(serotop_country_year$country %in% shortlisted_countries)
serotop_country_year <- serotop_country_year[index,]

write.csv(
  serotop_country_year,
  file = paste0(
    "serotypes_country_year_collapse_", args$collapse_strategy, ".tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
                  
# end logging
sink(con)
