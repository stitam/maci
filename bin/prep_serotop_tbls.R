rm(list = ls())

if (!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  aci_path <- args[1]
  collapse_strategy <- args[2]
  minyear <- args[3]
} else {
  test_dir <- "~/Methods/aci-tests/test-prep_serotop_tbls"
  aci_path <- paste0(test_dir, "/aci_filtered.rds")
  collapse_strategy <- "none"
  minyear <- 2016
}

if (!interactive()) {
  # create log file and start logging
  con <- file(paste0("log_", collapse_strategy, ".txt"))
  sink(con, split = TRUE)
}

library(dplyr)
library(ggplot2)

# import data
aci <- readRDS(aci_path)

# filter to minimum year
aci <- aci[which(aci$collection_year >= minyear), ]

print("Collection year for earliest sample:")
print(min(aci$collection_year))

print("Are all samples CRAB?")
print(all(aci$crab))

# prepare summary tables
foo <- function(data, quantile) {
  tbl <- data %>% table() %>% sort(., decreasing = TRUE)
  tbl <- tbl/length(data)
  tbl <- cumsum(tbl)
  res <- min(which(tbl >= quantile))
  return(res)
}

bar <- function(variable) {
  df <- aci %>% group_by(get(variable)) %>% summarise(
    q10 = foo(serotype, 0.1),
    q20 = foo(serotype, 0.2),
    q30 = foo(serotype, 0.3),
    q40 = foo(serotype, 0.4),
    q50 = foo(serotype, 0.5),
    q60 = foo(serotype, 0.6),
    q70 = foo(serotype, 0.7),
    q80 = foo(serotype, 0.8),
    q90 = foo(serotype, 0.9),
    q100 = foo(serotype, 1)
  )
  df_all <- aci %>% summarise(
    q10 = foo(serotype, 0.1),
    q20 = foo(serotype, 0.2),
    q30 = foo(serotype, 0.3),
    q40 = foo(serotype, 0.4),
    q50 = foo(serotype, 0.5),
    q60 = foo(serotype, 0.6),
    q70 = foo(serotype, 0.7),
    q80 = foo(serotype, 0.8),
    q90 = foo(serotype, 0.9),
    q100 = foo(serotype, 1)
  )
  df <- dplyr::bind_rows(df, df_all)
  df[nrow(df), 1] <- "all"
  df2 <- tidyr::pivot_longer(
    df,
    cols = 2:ncol(df),
    names_to = "quantile",
    names_prefix = "q",
    values_to = "serotop"
  )
  df2$quantile <- as.numeric(df2$quantile) / 100
  names(df2)[1] <- variable
  index <- which(df2[[variable]] != "all")
  df2$quantile[index] = jitter(df2$quantile[index])
  return(df2)
}

# set random seed for jitter
set.seed(0)

continent <- bar("continent")
region23 <- bar("region23")
country <- bar("country")

if (!interactive()) {
  # export continent summary
  write.table(
    continent,
    file = paste0("serotop_continent_collapse_", collapse_strategy, ".tsv"),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
  # export region summary
  write.table(
    region23,
    file = paste0("serotop_region23_collapse_", collapse_strategy, ".tsv"),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
  # export country summary
  write.table(
    country,
    file = paste0("serotop_country_collapse_", collapse_strategy, ".tsv"),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
}

if (!interactive()) {
  # end logging
  sink(con)
}
