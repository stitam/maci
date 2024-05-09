downsample_geodate3 <- function(aci) {
  # geodate3 is different from geodate2 in that it only looks at countries to 
  # avoid oversampling from countries where many isolates have cities.
  
  # add new temporal variables
  
  aci <- validate_var(
    df = aci,
    varname = "collection_day",
    varclass = "Date",
    coerce_fun = function(x) as.Date(x, origin = "1970-01-01"),
    verbose = TRUE
  )
  
  aci$yearweek <- paste(
    aci$collection_year, 
    lubridate::week(aci$collection_day),
    sep = "_"
  )
  
  aci$yearweek <- ifelse(grepl("NA", aci$yearweek), NA, aci$yearweek)
  
  if (any(is.na(aci$yearweek))) stop()
  
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
    geo <- aci[which(aci$country == i),]
    for (k in unique(geo$yearweek)) {
      geo_yearweek <- geo[which(geo$yearweek == k),]
      if (nrow(geo_yearweek) ==1) {
        keep_country <- c(keep_country, geo_yearweek$assembly)
      } else {
        keep_country <- c(keep_country, sample(geo_yearweek$assembly, 1))
      }
      
    }
  }
  
  index  <- which(aci$assembly %in% keep_country)
  
  aci_collapsed <- aci[index, ]
  
  return(aci_collapsed)
}