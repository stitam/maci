library(optparse)
rm(list = ls())

args_list <- list(
  make_option(
    c("-f", "--file"),
    type = "character",
    help = "Path to aci prediction results."
  ),
  make_option(
    c("-p", "--project_dir"),
    type = "character",
    help = "Path to the project directory."
  ),
  make_option(
    c("-s", "--downsampling_strategy"),
    type = "character",
    help = "Downsampling strategy."
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
    file = "aci_filtered.rds",
    project_dir = "~/Methods/aci",
    downsampling_strategy = "none",
    minyear = 2009,
    minyear_recent = 2016
    
  )
}

if (!interactive()) {
  # create log file and start logging
  con <- file(paste0("log_", args$downsampling_strategy, ".txt"))
  sink(con, split = TRUE)
}

library(devtools)
library(dplyr)
library(tidyr)

#produces an error
load_all(args$project_dir)

# import data set
aci <- read_df(args$file) %>% 
  dplyr::filter(filtered & crab & downsampled & downsampled_by_pop)

scopes <- c("continent", "region23", "country")
types <- c("all", "crab", "nocrab")

index <- which(aci$crab == TRUE)

if (length(index) == 0) {
  warning("Input table for diversity analysis does not contain crab isolates.")
}

if (length(index) == nrow(aci)) {
  warning("Input table for diversity analysis does not contain non-crab isolates.")
}

df <- data.frame()
for (i in scopes) {
  for (j in types) {
    # do not add new df if new df generation would fail
    add_new_df <- TRUE
    if (length(index) == 0 & j == "crab") {
      add_new_df <- FALSE
    }
    if (length(index) == nrow(aci) & j == "nocrab") {
      add_new_df <- FALSE
    }
    # add new df
    if (add_new_df) {
      newdf <- data.frame(
        scope = i,
        type = j,
        switch(
          j,
          all = div(
            aci, 
            by = i, 
            minyear = args$minyear, 
            minyear_recent = args$minyear_recent
          ),
          crab = div(
            aci[index, ], 
            by = i, 
            minyear = args$minyear, 
            minyear_recent = args$minyear_recent
          ),
          nocrab = div(
            aci[-index, ], 
            by = i,
            minyear = args$minyear, 
            minyear_recent = args$minyear_recent
          )
        )
      )
      # TODO include this in div()
      names(newdf)[3] <- "location"
      df <- dplyr::bind_rows(df, newdf)
    }
  }
}

df <- tibble::as_tibble(df)

write.table(
  df,
  file = paste0("diversity_stats_ds_", args$downsampling_strategy, ".tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

if (!interactive()) {
  # end logging
  sink(con)
}
