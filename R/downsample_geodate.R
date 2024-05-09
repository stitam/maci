downsample_geodate <- function(aci) {

  # Strategy for collapsing outbreaks:
  
  # If isolates from the same city and the same year and month have the same ST-K
  # combined serotype, then they are considered to be samples from the same
  # outbreak and only one of them is kept per month, at random.
  
  # For isolates where the city is not known, if isolates from the same country 
  # and from the same year and same week of the year have the same ST-K combined
  # serotype, then they are considered to be samples from the same outbreak and
  # only one of them is kept per week, at random.
  
  # Isolates with unknown collection dates will not be eliminated from the data
  # set.
  
  # add new temporal variables
  
  aci <- validate_var(
    df = aci,
    varname = "collection_day",
    varclass = "Date",
    coerce_fun = function(x) as.Date(x, origin = "1970-01-01"),
    verbose = TRUE
  )
  
  aci$yearweek <- paste(
    aci$collection_year, lubridate::week(aci$collection_day), sep = "_")
  aci$yearweek <- ifelse(grepl("NA", aci$yearweek), NA, aci$yearweek)
  
  aci$yearmonth <- paste(
    aci$collection_year, lubridate::month(aci$collection_day), sep = "_")
  aci$yearmonth <- ifelse(grepl("NA", aci$yearmonth), NA, aci$yearmonth)
  
  # collapse samples where city is known
  
  aci <- validate_var(
    df = aci,
    varname = "city",
    varclass = "character",
    coerce_fun = as.character,
    verbose = TRUE
  )
  
  aci <- validate_var(
    df = aci,
    varname = "serotype",
    varclass = "character",
    coerce_fun = as.character,
    verbose = TRUE
  )
  
  cities <- unique(aci$city)[which(!is.na(unique(aci$city)))]
  keep_city <- vector()
  for (i in cities) {
    geo <- aci[which(aci$city == i),]
    serotypes <- unique(geo$serotype)
    for (j in serotypes) {
      geo_sero <- geo[which(geo$serotype == j),]
      for (k in unique(geo_sero$yearmonth)) {
        if (is.na(k)) {
          # if either year or month or both are not known, keep assembly
          keep_city <- c(
            keep_city, geo_sero$assembly[which(is.na(geo_sero$yearmonth))])
        } else {
          # otherwise keep one assembly per year per month
          geo_sero_yearmonth <- geo_sero[which(geo_sero$yearmonth == k),]
          keep_city <- c(keep_city, sample(geo_sero_yearmonth$assembly, 1))
        }
      }
    }
  }
  
  aci <- validate_var(
    df = aci,
    varname = "country",
    varclass = "character",
    coerce_fun = as.character,
    verbose = TRUE
  )
  
  # collapse samples where city is not known but country is known
  countries <- unique(aci$country)[which(!is.na(unique(aci$country)))]
  keep_country <- vector()
  for (i in countries) {
    # only consider samples where the country is known but the city is not
    geo <- aci[which(aci$country == i & is.na(aci$city)),]
    serotypes <- unique(geo$serotype)
    for (j in serotypes) {
      geo_sero <- geo[which(geo$serotype == j),]
      for (k in unique(geo_sero$yearweek)) {
        if (is.na(k)) {
          # if either year or week or both are not known, keep assembly
          keep_country <- c(
            keep_country, geo_sero$assembly[which(is.na(geo_sero$yearweek))])
        } else {
          # otherwise keep one assembly per year per week
          geo_sero_yearweek <- geo_sero[which(geo_sero$yearweek == k),]
          keep_country <- c(keep_country, sample(geo_sero_yearweek$assembly, 1))
        }
      }
    }
  }
  
  index  <- which(
    aci$assembly %in% keep_city | 
      aci$assembly %in% keep_country |
      is.na(aci$country)
  )
  
  aci_collapsed <- aci[index, ]
  
  return(aci_collapsed)
}