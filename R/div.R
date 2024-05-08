#' Calculate diversity metrics
#' 
#' This function takes a data frame and a variable of the data frame and
#' calculates various diversity metrics e.g. Minor Allele Frequency (MAF) or
#' Shannon diversity.
#' @param df data.frame; a data frame which contains a row for each sample
#' @param by character; a variable of the data frame to be investigated
#' @param minyear numeric; the lowest year to include
#' @param minyear_recent numeric; the lowest year to include for calculations
#' on recent samples
#' @return a data frame.
div <- function(df, by, minyear, minyear_recent) {
  # TODO the data frame must have obligatory columns. Include these in doc.
  # TODO functions used by this function have extra arguments, include these
  # arguments here as well?
  df_summary <- df %>% 
    group_by(get(by)) %>%
    summarise(
      
      samples = n(),
      divAll_count = length(unique(serotype)),
      div0.5_count = div_count(serotype, collection_year, year_th = NULL, 0.5),
      div0.8_count = div_count(serotype, collection_year, year_th = NULL, 0.8),
      div3_ratio = div_ratio(serotype, collection_year, year_th = NULL, 3),
      div5_ratio = div_ratio(serotype, collection_year, year_th = NULL, 5),
      maf_ratio = maf_ratio(serotype, collection_year, year_th = NULL),
      shannon = shannon(serotype, collection_year, year_th = NULL),
      
      samples_ALL = length(which(collection_year >= minyear)),
      divAll_ALL_count = length(
        unique(serotype[which(collection_year >= minyear)])),
      div0.5_ALL_count = div_count(serotype, collection_year, minyear, 0.5),
      div0.8_ALL_count = div_count(serotype, collection_year, minyear, 0.8),
      div3_ALL_ratio = div_ratio(serotype, collection_year, minyear, 3),
      div5_ALL_ratio = div_ratio(serotype, collection_year, minyear, 5),
      maf_ALL_ratio = maf_ratio(serotype, collection_year, minyear),
      shannon_ALL = shannon(serotype, collection_year, year_th = minyear),
      
      samples_RECENT = length(which(collection_year >= minyear_recent)),
      divAll_RECENT_count = length(unique(serotype[which(collection_year >= minyear_recent)])),
      div0.5_RECENT_count = div_count(serotype, collection_year, minyear_recent, 0.5),
      div0.8_RECENT_count = div_count(serotype, collection_year, minyear_recent, 0.8),
      div3_RECENT_ratio = div_ratio(serotype, collection_year, minyear_recent, 3),
      div5_RECENT_ratio = div_ratio(serotype, collection_year, minyear_recent, 5),
      maf_RECENT_ratio = maf_ratio(serotype, collection_year, minyear_recent),
      shannon_RECENT = shannon(serotype, collection_year, year_th = minyear_recent),
    )
  names(df_summary)[1] = by
  
  df_summary$divShr_count = NA
  df_summary$divShr_ratio = NA
  for (i in 1:nrow(df_summary)) {
    if (is.na(df_summary[[by]][i])) {
      sero <- unique(df$serotype[which(is.na(df[[by]]))])
      sero_others <- unique(df$serotype[which(!is.na(df[[by]]))])
    } else {
      sero <- unique(df$serotype[which(df[[by]] == df_summary[[by]][i])])
      sero_others <- unique(df$serotype[which(df[[by]] != df_summary[[by]][i])])
    }
    df_summary$divShr_count[i] <- sum(sero %in% sero_others)
    df_summary$divShr_ratio[i] <- round(mean(sero %in% sero_others),3)
  }
  
  df_summary$divShr_ALL_count = NA
  df_summary$divShr_ALL_ratio = NA
  for (i in 1:nrow(df_summary)) {
    if (is.na(df_summary[[by]][i])) {
      sero <- unique(df$serotype[which(is.na(df[[by]]) & df$collection_year >= minyear)])
      sero_others <- unique(df$serotype[which(!is.na(df[[by]]) & df$collection_year >= minyear)])
    } else {
      sero <- unique(df$serotype[which(df[[by]] == df_summary[[by]][i] & df$collection_year >= minyear)])
      sero_others <- unique(df$serotype[which(df[[by]] != df_summary[[by]][i] & df$collection_year >= minyear)])
    }
    df_summary$divShr_ALL_count[i] <- sum(sero %in% sero_others)
    df_summary$divShr_ALL_ratio[i] <- round(mean(sero %in% sero_others),3)
  }
  
  df_summary$divShr_RECENT_count = NA
  df_summary$divShr_RECENT_ratio = NA
  for (i in 1:nrow(df_summary)) {
    if (is.na(df_summary[[by]][i])) {
      sero <- unique(df$serotype[which(is.na(df[[by]]) & df$collection_year >= minyear_recent)])
      sero_others <- unique(df$serotype[which(!is.na(df[[by]]) & df$collection_year >= minyear_recent)])
    } else {
      sero <- unique(df$serotype[which(df[[by]] == df_summary[[by]][i] & df$collection_year >= minyear_recent)])
      sero_others <- unique(df$serotype[which(df[[by]] != df_summary[[by]][i] & df$collection_year >= minyear_recent)])
    }
    df_summary$divShr_RECENT_count[i] <- sum(sero %in% sero_others)
    df_summary$divShr_RECENT_ratio[i] <- round(mean(sero %in% sero_others),3)
  }
  names(df_summary)[c(10:17, 28:29)] <- gsub("ALL", as.character(minyear), names(df_summary)[c(10:17, 28:29)])
  names(df_summary)[c(18:25, 30:31)] <- gsub("RECENT", as.character(minyear_recent), names(df_summary)[c(18:25, 30:31)])
  return(df_summary) 
}
