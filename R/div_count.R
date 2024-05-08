div_count <- function(var, year, year_th, ratio_th) {
  if (is.null(year_th)) {
    index <- seq_along(var)
  } else {
    index <- which(year >= year_th)
  }
  if (length(index) > 0) {
    cs_count <- cumsum(sort(table(var[index]), decreasing = TRUE))
    cs_ratio <- cs_count/length(index)
    m <- min(which(cs_ratio >= ratio_th))
    return(m)
  }
  return(NA)
}