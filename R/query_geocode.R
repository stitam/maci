#' Query a geocode from Google API
#' 
#' This function is a wrapper around ggmap::geocode(). The Google Geocoding API
#' is not a free service, therefore any queries towards the API should make use
#' of caching. This function implements such caching approach.
#' @param location character; a vector of locations.
#' @param cache_file character; the name of the cache file without the file
#' extension.
#' @param verbose logical; should a verbose output be printed on the console?
#' @examples
#' \donttest{
#' query_geocode("Budapest")
#' }
#' @import ggmap
#' @export
query_geocode <- function(query,
                          cache_file = "geocodes",
                          verbose = getOption("verbose")) {
  foo <- function(query, verbose = verbose) {
    res <- suppressMessages(ggmap::geocode(query, output = "latlona"))
    res <- cbind(data.frame(query = query), res)
    if (inherits(res, "try-error")) {
      return(NA)
    }
    if (verbose) message("Done.")
    return(res)
  }
  cfpath <- paste0(cache_file, ".rds")
  if (file.exists(cfpath)) {
    geocodes <- readRDS(file = cfpath)
  } else {
    geocodes <- data.frame()
  }
  out <- lapply(query, function(x) {
    if (verbose) message("Querying ", x, ". ", appendLF = FALSE)
    
    if (x %in% geocodes$query) {
      if (verbose) message("Already retrieved.")
      return(geocodes[which(geocodes$query == x),])
    } else {
      new <- foo(x, verbose = verbose)
      if (!is.na(x)) {
        geocodes <<- rbind(geocodes, new)
        geocodes <<- geocodes[order(geocodes$query, decreasing = FALSE),]
        geocodes <<- tibble::as_tibble(geocodes)
        saveRDS(geocodes, file = cfpath)
      }
      return(new)
    }
  })
  out <- dplyr::bind_rows(out)
  return(out)
}
