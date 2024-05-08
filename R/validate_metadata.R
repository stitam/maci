#' Validate metadata
#' 
#' This function checks whether the data frame which contains the sample
#' metadata contains certain columns and whether these columns have the right
#' classes.
#' @param df data.frame; a data.frame of sample metadata
#' @param keep_extra logical; if the data.frame contains extra variables which
#' will not be checked by this function, should these columns be kept in the
#' returned object? Note that this may result in merge errors.
#' @param reorder logical; should the function return variables in a specific
#' order? If extra variables are kept, these variables will be put to the end of
#' the data frame.
#' @return  If the data frame meets all the requirements, it will be converted
#' to an object of class \code{meta_df}. Subsequent elements of the data
#' analytical pipeline will only accept metadata as a \code{meta_df} object.
#' @details Some variables are required and the function will stop with an
#' error if any of these variables are missing or not provided appropriately.
#' Other variables are optional and the function will print a warning message
#' when any of these variables are missing or stop with an error if any of them
#' are not provided appropriately. The function will look for the following
#' variables and classes:
#' \itemize{
#' \item \code{assembly} (required) character
#' \item \code{coverage} (optional) character
#' \item \code{organism} (required) character
#' \item \code{collection_date} (optional) character
#' \item \code{collection_year} (optional) integer
#' \item \code{collection_day} (optional) Date
#' \item \code{continent} (optional) character
#' \item \code{region23} (optional) character
#' \item \code{country} (optional)  character
#' \item \code{country_iso2c} (optional) character
#' \item \code{city} (optional) character
#' \item \code{ward} (optional) character
#' \item \code{covid} (optional) logical
#' \item \code{assembly_source} (required) character
#' \item \code{in_lab} (optional) logical
#' \item \code{lon} (optional) numeric
#' \item \code{lat} (optional) numeric
#' \item \code{jobname} (required) character
#' \item \code{relative_path} (required) character
 #' }
#' @note The purpose of this function is to ensure that results from analysing
#' multiple batches of sequences are compatible and can be merged easily.
#' @export
validate_metadata <- function(df, keep_extra = FALSE, reorder = TRUE){
  if ("assembly" %in% names(df) == FALSE) {
    stop("Required column 'assembly' is missing.")
  }
  if (class(df$assembly) != "character"){
    stop("Required column 'assembly' must be a 'character'.")
  }
  if (sum(grepl("root", df$assembly, ignore.case = TRUE)) > 0) {
    stop("Assembly name cannot be 'root'. Specify another name.")
  }
  if (sum(!is.na(suppressWarnings(as.numeric(df$assembly)))) > 0) {
    stop("Assembly name cannot be a number. Specify another name.")
  }
  if (length(unique(df$assembly)) != length(df$assembly)) {
    stop("Assembly names must be unique. Check.")
  }
  if ("coverage" %in% names(df) == FALSE) {
    warning("Optional column 'coverage' is missing. \n")
  } else if (class(df$coverage) != "numeric"){
    stop("Optional column 'coverage' must be a 'numeric'.")
  }
  if ("organism" %in% names(df) == FALSE) {
    stop("Required column 'organism' is missing.")
  }
  if (class(df$organism) != "character"){
    stop("Required column 'organism' must be a 'character'.")
  }
  if (length(unique(df$organism)) > 1) {
    stop("'organism' must be 'Acinetobacter baumannii'.")
  }
  if (unique(df$organism) != "Acinetobacter baumannii") {
    stop("'organism' must be 'Acinetobacter baumannii'.")
  }
  if ("strain" %in% names(df) == FALSE) {
    warning("Optional column 'strain' is missing. \n")
  } else if (class(df$strain) != "character"){
    stop("Optional column 'strain' must be a 'character'.")
  }
  if ("collection_date" %in% names(df) == FALSE) {
    warning("Optional column 'collection_date' is missing. \n")
  } else if (class(df$collection_date) != "character"){
    stop("Optional column 'collection_date' must be a 'character'.")
  }
  if ("collection_year" %in% names(df) == FALSE) {
    warning("Optional column 'collection_year' is missing. \n")
  } else if (class(df$collection_year) != "integer"){
    stop("Optional column 'collection_year' must be an 'integer'.")
  }
  if ("collection_day" %in% names(df) == FALSE) {
    warning("Optional column 'collection_day' is missing. \n")
  } else if (class(df$collection_day) != "Date"){
    stop("Optional column 'collection_day' must be a 'Date'.")
  }
  if ("continent" %in% names(df) == FALSE) {
    warning("Optional column 'continent' is missing. \n")
  } else if (class(df$continent) != "character"){
    stop("Optional column 'continent' must be a 'character'.")
  }
  if ("region23" %in% names(df) == FALSE) {
    warning("Optional column 'region23' is missing. \n")
  } else if (class(df$region23) != "character") {
    stop("Optional column 'region23' must be a 'character'.")
  }
  if ("country" %in% names(df) == FALSE) {
    warning("Optional column 'country' is missing. \n")
  } else if (class(df$country) != "character"){
    stop("Optional column 'country' must be a 'character'.")
  }
  if ("country_iso2c" %in% names(df) == FALSE) {
    warning("Optional column 'country_iso2c' is missing. \n")
  } else if (class(df$country_iso2c) != "character"){
    stop("Optional column 'country_iso2c' must be a 'character'.")
  }
  if ("city" %in% names(df) == FALSE) {
    warning("Optional column 'city' is missing. \n")
  } else if (class(df$city) != "character"){
    stop("Optional column 'city' must be a 'character'.")
  }
  if ("ward" %in% names(df) == FALSE) {
    warning("Optional column 'ward' is missing. \n")
  } else if (class(df$ward) != "character"){
    stop("Optional column 'ward' must be a 'character'.")
  }
  if ("covid" %in% names(df) == FALSE) {
    warning("Optional column 'covid' is missing. \n")
  } else if (class(df$covid) != "logical"){
    stop("Optional column 'covid' must be a 'logical'.")
  }
  if ("assembly_source" %in% names(df) == FALSE) {
    stop("Required column 'assembly_source' is missing.")
  }
  if (class(df$assembly_source) != "character"){
    stop("Required column 'assembly_source' must be a 'character'.")
  }
  if ("in_lab" %in% names(df) == FALSE) {
    warning("Optional column 'in_lab' is missing. \n")
  } else if (class(df$in_lab) != "logical"){
    stop("Optional column 'in_lab' must be a 'logical'.")
  }
  if ("lon" %in% names(df) == FALSE) {
    warning("Optional column 'lon' is missing. \n")
  } else if (class(df$lon) != "numeric"){
    stop("Optional column 'lon' must be a 'numeric'.")
  }
  if ("lat" %in% names(df) == FALSE) {
    warning("Optional column 'lat' is missing. \n")
  } else if (class(df$lat) != "numeric"){
    stop("Optional column 'lat' must be a 'numeric'.")
  }
  if ("jobname" %in% names(df) == FALSE) {
    stop("Required column 'jobname' is missing.")
  }
  if (class(df$jobname) != "character"){
    stop("Required column 'jobname' must be a 'character'.")
  }
  if (sum(is.na(df$jobname)) > 0) {
    stop("Required column 'jobname' cannot contain NA values.")
  }
  if ("relative_path" %in% names(df) == FALSE) {
    stop("Required column 'relative_path' is missing.")
  }
  if (class(df$relative_path) != "character") {
    stop("Required column 'relative_path' must be a 'character'.")
  }
  if (any(is.na(df$relative_path))) {
    stop("Required column 'relative_path' cannot contain NA values.")
  }
  meta_vars <- c(
    "assembly",
    "coverage",
    "organism",
    "strain",
    "collection_date",
    "collection_year",
    "collection_day",
    "continent",
    "region23",
    "country",
    "country_iso2c",
    "city",
    "ward",
    "covid",
    "assembly_source",
    "in_lab",
    "lon",
    "lat",
    "jobname",
    "relative_path"
  )
  if (keep_extra == TRUE) {
    extra_vars <- names(df)[which(names(df) %in% meta_vars == FALSE)]
    extra_vars_collapsed <- paste(extra_vars, collapse = ", ")
    warning(
      "The following extra columns are included in the output: ", extra_vars_collapsed,  "\n"
    )
    df <- df[, which(names(df) %in% c(meta_vars, extra_vars))]
  } else {
    df <- df[,which(names(df) %in% meta_vars)]
  }
  if (reorder == TRUE) {
    if (keep_extra == TRUE) {
      df1 <- df[, meta_vars[which(meta_vars %in% names(df))]]
      df2 <- df[, which(names(df) %in% extra_vars)]
      df_out <- dplyr::bind_cols(df1, df2)
    } else {
      df_out <- df
    }
  }
  class(df_out) = c("meta_df", "data.frame")
  return(df_out)
}
