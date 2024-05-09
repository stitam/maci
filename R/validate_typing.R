#' Validate a typing summary table
#' 
#' There are functions and scripts in the package that take a typing summary
#' table as input. In order for these to work properly, the summary table must 
#' contain certain variables and these must be of the right type. This function
#' validates such a table, checks for variables names and classes, attempts to
#' resolve any issues and stops with an error if necessary. 
#' @param df data.frame; an assembly summary table.
#' @param taxid character; The NCBI Taxonomy ID for the organism.
#' @param new_na character; what should the function do if validation introduces
#' new NAs? Either \code{"ignore"} or \code{"stop"}.
#' @param verbose logical; should verbose messages be printed to the console?
#' @return a data frame
validate_typing <- function(
    df,
    taxid,
    new_na = "ignore",
    verbose = getOption("verbose")
    ) {
  taxid <- as.numeric(taxid)
  new_na <- match.arg(new_na, choices = c("ignore", "stop"))
  character_variables <- c(
    "assembly", "collection_date", "country", "city", "mlst", "kraken2_taxid"
  )
  if (taxid == 470) {
    character_variables <- c(character_variables, "mlst_genes_detected")
  }
  if (taxid %in% c(470, 573)) {
    character_variables = c(character_variables, "k_serotype", "k_confidence")
  }
  for (i in character_variables) {
    df <- validate_var(
      df, 
      varname = i,
      varclass = "character",
      coerce_fun = as.character,
      new_na = new_na,
      verbose = verbose
    )
  }
  if (any(is.na(df[["assembly"]]))) {
    stop("'assembly' variable within 'df' cannot be NA.")
  }
  if (nrow(df) != length(unique(df[["assembly"]]))) {
    stop("'assembly' variable within 'df' must be unique.")
  }
  # TODO this is a metadata issue that should be fixed in the metadata
  # This is a temporary fix. Include in webseq e.g. validate parsed
  # NCBI BioSample metadata?
  if (any(df[["country"]] == "Korea", na.rm = TRUE)) {
    df[["country"]][which(df[["country"]] == "Korea")] <- "South Korea"
  }
  numeric_variables <- c(
    "coverage", "busco_complete_perc", "kraken2_taxid_freq", "genome_size",
    "gc_content", "contig_count", "longest_contig", "N50", "N95", "n_count"
  )
  for (i in numeric_variables) {
    df <- validate_var(
      df,
      varname = i,
      varclass = "numeric",
      new_na = new_na,
      coerce_fun = function(x) suppressWarnings(as.numeric(as.character(x))),
      verbose = verbose
    )
  }
  logical_variables <- c("human_related", "crab")
  for (i in logical_variables) {
    df <- validate_var(
      df,
      varname = i,
      varclass = "logical",
      coerce_fun = function(x) sapply(x, function(A) {
        if (is.na(A)) {
          return(NA) 
        }
        if (A %in% c("TRUE", "True", "true", "YES", "Yes", "yes")) {
          return(TRUE)
        } else if (A %in% c("FALSE", "False", "false", "NO", "No", "no")) {
          return(FALSE)
        } else {
          return(NA)
        }
      }),
      new_na = new_na,
      verbose = verbose
    )
  }
  opt_logical_variables <- c("qc_pass")
  for (i in opt_logical_variables) {
    if (i %in% names(df)) {
      df <- validate_var(
        df,
        varname = i,
        varclass = "logical",
        coerce_fun = function(x) sapply(x, function(A) {
          if (is.na(A)) {
            return(NA) 
          }
          if (A %in% c("TRUE", "True", "true", "YES", "Yes", "yes")) {
            return(TRUE)
          } else if (A %in% c("FALSE", "False", "false", "NO", "No", "no")) {
            return(FALSE)
          } else {
            return(NA)
          }
        }),
        new_na = new_na,
        verbose = verbose
      )
    }
  }
  return(df)
}