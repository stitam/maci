#' Calculate serotype frequences
#'
#' This function takes a data frame with a single grouping variable, e.g.
#' `region23` or `country` and prepares a table in which for each level of the
#' grouping variable serotypes are listed in decreasing frequence order.
#' @param df data.frame, a table with one row for each sample. The table must
#' have a `serotype` variable and an optional grouping variable.
#' @param group_by character; the name of the grouping variable. If not `NULL`
#' the input data frame must contain a variable with this name.
#' @return A data frame.
#' @examples
#' serotype_freqs(df, group_by = "region23")
serotype_freqs <- function(df, group_by = NULL) {
  if (length(group_by) > 1) stop()
  if ("serotype" %in% names(df) == FALSE) {
    stop("Input data frame must have a variable called 'serotype'.")
  }
  if (group_by %in% names(df) == FALSE) {
    stop(paste0(
      "Input data frame must have a variable called '",
      group_by,
      "'."
    ))
  }
  foo <- function(tbl) {
    tbl_sum <- tbl %>%
      group_by(serotype) %>%
      summarise(count = n())
    tbl_sum$ratio = signif(tbl_sum$count/sum(tbl_sum$count), 4)
    tbl_sum <- tbl_sum[order(tbl_sum$ratio, decreasing = TRUE),]
    tbl_sum$ratio_cumsum <- signif(cumsum(tbl_sum$ratio), 4)
    tbl_sum$order <- 1:nrow(tbl_sum)
    return(tbl_sum)
  }
  if (is.null(group_by)) {
    freqs <- foo(df)
  } else {
    freqs <- data.frame()
    for (i in unique(df[[group_by]])) {
      df_small <- df[which(df[[group_by]] == i),]
      df_sum <- foo(df_small)
      df_sum[[group_by]] <- i
      df_sum <- dplyr::relocate(df_sum, group_by)
      freqs = dplyr::bind_rows(
        freqs,
        df_sum
      )
    }
  }
  return(freqs)
}
