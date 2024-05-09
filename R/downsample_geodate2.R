downsample_geodate2 <- function(aci) {
  # geodate2 is different from geodate in that it always keeps one sample per
  # yearweek (both when cities are know and when they are not). Also geodate2
  # does not look at serotypes, only at yearweek.

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
  
  # Downsample where city is known
  
  aci <- validate_var(
    df = aci,
    varname = "city",
    varclass = "character",
    coerce_fun = as.character,
    verbose = TRUE
  )
  
  cities <- unique(aci$city)[which(!is.na(unique(aci$city)))]
  keep_city <- vector()
  for (i in cities) {
    geo <- aci[which(aci$city == i),]
    for (k in unique(geo$yearweek)) {
      geo_yearweek <- geo[which(geo$yearweek == k),]
      keep_city <- c(keep_city, sample(geo_yearweek$assembly, 1))
    }
  }
  
  # Downsample where city is not known but country is known
  
  aci <- validate_var(
    df = aci,
    varname = "country",
    varclass = "character",
    coerce_fun = as.character,
    verbose = TRUE
  )
   
  countries <- unique(aci$country)[which(!is.na(unique(aci$country)))]
  keep_country <- vector()
  for (i in countries) {
    # only consider samples where the country is known but the city is not
    geo <- aci[which(aci$country == i & is.na(aci$city)),]
    for (k in unique(geo$yearweek)) {
      geo_yearweek <- geo[which(geo$yearweek == k),]
      if (nrow(geo_yearweek) == 1) {
        keep_country <- c(keep_country, geo_yearweek$assembly)
      } else {
        # Keep max 2 samples per country per yearweek
        keep_country <- c(keep_country, sample(geo_yearweek$assembly, 2))
      }
    }
  }
  
  index  <- which(aci$assembly %in% c(keep_city, keep_country))
  
  aci_collapsed <- aci[index, ]
  
  return(aci_collapsed)
}