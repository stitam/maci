#' Group samples using text data
#' 
#' This function looks at inclusion and exclusion terms in a character vector to
#' determine whether samples belong to a group or not.
#' @param x character; a vector of strings used for grouping samples.
#' @param include character; a vector of inclusion terms which clearly identify
#' samples that belong to the group.
#' @param exclude character; a vector of exclusion terms which clearly identify
#' samples that do not belong to the group.
#' @return A logical vector. Elements of the vector can be \code{TRUE} if an 
#' inclusion term was matched, \code{FALSE} if an exclusion term was matched,
#' or \code{NA} if neither inclusion nor exclusion terms were matched.
#' @examples
#' \dontrun{
#' # Identify human associated samples using metadata terms from multiple
#' variables.
#' 
#' # Step 1. Look at variable called "host"
#' include <- "human"
#' exclude <- "dog"
#' meta$human_related <- in_group(
#'   meta$host,
#'   include = include,
#'   exclude = exclude
#' )
#' 
#' # Step 2. Look at variable "isolation_source"
#' # Only group samples that were not grouped in the previous step.
#' include <- "human blood"
#' exclude <- "wastewater"
#' meta$human_related <- in_group(
#'   meta$isolation_source,
#'   include = include,
#'   exclude = exclude,
#'   cache = meta$human_related)
#' )
#' }
#' @export
in_group <- function(x,
                     include = character(),
                     exclude = character(),
                     cache = NULL
                     ) {
  if (!is.null(cache)) {
    in_group <- cache
  } else {
    in_group <- rep(NA, times = length(x))
  }
  for (i in 1:length(x)) {
    if (!is.na(in_group[i])) {
      next()
    } else {
      if (x[i] %in% include & x[i] %in% exclude) {
        msg <- paste0(
          "The following string is both included and excluded: ",
          x[i]
        )
        stop(msg)
      }
      if (x[i] %in% include) {
        in_group[i] <- TRUE
      }
      if (x[i] %in% exclude) {
        in_group[i] <- FALSE
      }
    }
  }
  return(in_group)
}