#' Downsample isolates
#' 
#' Downsample isolates using various downsampling strategies.
#' @param df data.frame; a table with one row for each isolate
#' @param strategy character; the downsampling strategy, can be `geodate`,
#' `geodate2`, `geodate3`, `poppunk`, or `population`. 
#' @return The function returns a data frame which contains the same variables
#' as the input data frame but does not contain all isolates only the
#' downsampled isolates.
downsample <- function(df, strategy, ...) {
  if (!"assembly" %in% names(df)) {
    stop("'df' must have a variable called 'assembly'.")
  }
  df <- switch(
    strategy,
    geodate = downsample_geodate(df),
    geodate2 = downsample_geodate2(df),
    geodate3 = downsample_geodate3(df),
    poppunk = downsample_poppunk(df, pp),
    population = downsample_population(df, population)
  )
  return(df)
}