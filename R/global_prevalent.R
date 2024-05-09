#' Identify prevalent and global types
#' 
#' This function takes a data frame which contains one row for each sample and
#' return a data frame of prevalent and global serotypes.
#' @param df data.frame, a table with one row for each sample. The table must
#' have a `serotype` variable.
#' @return A data frame.
global_prevalent <- function(df) {
  if ("continent" %in% names(df) == FALSE) {
    stop("Input data frame must have a variable called 'continent'.")
  }
  region_continent_dict <- df %>%
    dplyr::select(region23, continent) %>%
    dplyr::distinct()
  serotop_region <- serotype_freqs(df, group_by = "region23")
  serotop_region2 <- dplyr::left_join(
    serotop_region,
    region_continent_dict,
    by = "region23"
  )
  # prevalent_clones: a serotype is considered prevalent if it has a frequency
  # > 5% in at  least 1 region
  prevalent_clones <- serotop_region2 %>%
    filter(ratio > 0.05) %>%
    summarise(serotype = unique(serotype)) %>% 
    mutate(color = "orange3")
  # global clones: a serotype is considered a global if it is prevalent and has
  # a frequency > 2% in at least 3 regions and at least 2 continents
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
  return(global_prevalent)
}
