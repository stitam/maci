shannon <- function(var, year, year_th) {
  if (is.null(year_th)) {
    index <- seq_along(var)
  } else {
    index <- which(year >= year_th)
  }
  if (length(index) > 0) {
    p <- table(var[index])/length(index)
    sh <- round(-sum(p*log(p)), 3)
    return(sh)
  }
  return(NA)
}
