div_ratio <- function(var, year, year_th, top) {
  if (is.null(year_th)) {
    index <- seq_along(var)
  } else {
    index <- which(year >= year_th)
  }
  if (length(index) > 0) {
    tbl <- sort(table(var[index]), decreasing = TRUE)
    if (length(tbl) >= top) {
      count <- sum(tbl[1:top])
      ratio <- round(count/length(index), 3)
      return(ratio)
    } else {
      return(1)
    }
  } 
  return(NA)
}