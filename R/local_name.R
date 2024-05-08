local_name <- function(city, context = NULL) {
  data(citydict, envir = rlang::current_env())
  if (is.null(context)) {
    return(city) 
  } else {
    res <- unname(sapply(city, function(x) {
      index <- which(citydict$international == x & citydict$context == context)
      if (length(index) > 0) {
        out <- citydict$local[index]
      } else {
        out <- x
      }
      return(out)
    }))
    return(res)
  }
}