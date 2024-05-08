get_ts_df <- function(df, var, timevar, window) {
  if ("data.frame" %in% class(df) == FALSE) {
    stop ("Argument 'df' must be a data.frame.")
  }
  if (var %in% names(df) == FALSE) {
    stop(paste0("Variable '", var, "' not found."))
  }
  if (timevar %in% names(df) == FALSE) {
    stop(paste0("Variable '", timevar, "' not found."))
  }
  if (!is.numeric(df[[timevar]])) {
    stop(paste0("Variable '", timevar, "' must be numeric."))
  }
  if (!is.numeric(window)) {
    stop(paste("Variable '", window, "' must be numeric."))
  }
  df <- tibble::as_tibble(as.data.frame(
    table(df[[var]], df[[timevar]], useNA = "ifany")))
  names(df) <- c(var, timevar, "count")
  df[[timevar]] <- as.numeric(as.character(df[[timevar]]))
  df$yearcount <- NA
  for (i in 1:nrow(df)) {
    df$yearcount[i] <- sum(df$count[which(df[[timevar]] == df[[timevar]][i])])
  }
  df$count_window <- NA
  df$yearcount_window <- NA
  for (i in 1:nrow(df)) {
    if (is.na(df[[var]][i])) {
      index <- which(
        is.na(df[[var]]) &
        df[[timevar]] <= df[[timevar]][i] &  
        df[[timevar]] > df[[timevar]][i]-window
      )
    } else {
      index <- which(
        df[[var]] == df[[var]][i] & 
        df[[timevar]] <= df[[timevar]][i] &  
        df[[timevar]] > df[[timevar]][i]-window
      )
    }
    df$count_window[i] <- sum(df$count[index])
    index_year <- which(
      df[[timevar]] <= df[[timevar]][i] & 
      df[[timevar]] > df[[timevar]][i]-window
    )
    df$yearcount_window[i] <- sum(df$count[index_year])
  }
  df$ratio_window <- df$count_window/df$yearcount_window
  return(df)
}

