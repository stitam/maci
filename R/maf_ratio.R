maf_ratio <- function(var, year, year_th) {
  if (is.null(year_th)) {
    index <- seq_along(var)
  } else {
    index <- which(year >= year_th)
  }
  if (length(index) > 0) {
    tbl <- sort(table(var[index]), decreasing = TRUE)
    if (length(tbl) > 1) {
      count <- unname(tbl[2])
      ratio <- round(count/length(index), 3)
      return(ratio)
    }
  }
  return(NA)
}
