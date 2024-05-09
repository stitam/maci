#' Cumulative relative prevalence of most prevalent types
#' 
#' This function prepares a table about the number of most prevalent types which
#' account for a certain ratio of all isolates.
#' @param df data.frame; a table which contains one row for each sample and at
#' at least two columns: the `type` of the sample and a grouping variable.
#' @param type character; the name of a variable within `df` we are interested
#' in.
#' @param group_by character; the name of the grouping variable.
#' @return a data.frame
#' @importFrom magrittr %>%
prep_serotop_tbl <- function(df, type, group_by) {
  
  foo <- function(data, quantile) {
    tbl <- data %>% table() %>% sort(., decreasing = TRUE)
    tbl <- tbl/length(data)
    tbl <- cumsum(tbl)
    res <- min(which(tbl >= quantile))
    return(res)
  }
  
  top <- df %>% group_by(get(group_by)) %>% summarise(
    q10 = foo(get(type), 0.1),
    q20 = foo(get(type), 0.2),
    q30 = foo(get(type), 0.3),
    q40 = foo(get(type), 0.4),
    q50 = foo(get(type), 0.5),
    q60 = foo(get(type), 0.6),
    q70 = foo(get(type), 0.7),
    q80 = foo(get(type), 0.8),
    q90 = foo(get(type), 0.9),
    q100 = foo(get(type), 1)
  )
  
  top_all <- df %>% summarise(
    q10 = foo(get(type), 0.1),
    q20 = foo(get(type), 0.2),
    q30 = foo(get(type), 0.3),
    q40 = foo(get(type), 0.4),
    q50 = foo(get(type), 0.5),
    q60 = foo(get(type), 0.6),
    q70 = foo(get(type), 0.7),
    q80 = foo(get(type), 0.8),
    q90 = foo(get(type), 0.9),
    q100 = foo(get(type), 1)
  )
  
  out <- dplyr::bind_rows(top, top_all)
  names(out)[1] <- "group"
  out[nrow(out), 1] <- "All"
  out <- tidyr::pivot_longer(
    out,
    cols = 2:ncol(out),
    names_to = "quantile",
    names_prefix = "q",
    values_to = "top"
  )
  out$quantile <- as.numeric(out$quantile) / 100
  
  return(out)
}
