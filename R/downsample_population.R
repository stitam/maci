#' Downsample by population
#' 
#' This function randomly selects isolates in proportion with the number of
#' inhabitants in each country.
#' @param aci data.frame; a table with one row for each sample.
#' @param population data.frame; a table with country names and populations.
#' @param rate numeric; sampling rate, samples per million inhabitants.
#' @return A data frame which contains the same variables as the input data 
#' frame but does not contain all isolates only the downsampled isolates.
downsample_population <- function(aci, population, rate) {
  if (any(is.na(aci$country))) stop()
  testthat::expect_true(all(aci$country %in% population$country))
  countries <- unique(aci$country)
  aci_pop <- data.frame()
  for (i in countries) {
    index <- which(aci$country == i)
    country_population <- population$population[which(population$country == i)]
    # sample size depends on population and rate but is at least one
    sample_size <- as.integer(max(round(country_population * rate), 0), 1)
    if (sample_size >= length(index)) {
      index_sample <- index
    } else {
      index_sample <- index[sample(1:length(index), sample_size)]
    }
    aci_pop <- dplyr::bind_rows(
      aci_pop,
      aci[index_sample, ]
    )
  }
  out <- aci[sort(which(aci$assembly %in% aci_pop$assembly)), ]
  return(out)
}
