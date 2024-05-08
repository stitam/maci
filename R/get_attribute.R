#' Get any attribute from a parsed GFF file
#' 
#' GFF files include a field which contains a number of attributes concatenated
#' into a single variable. This function can be used to recover individual
#' attributes from this concatenated variable.
#' @param x character; attributes variable of the parsed gff file
#' @param field character; the attribute to be recovered
#' @param attrsep character; the separator which separates the attributes
#' @return character; attribute values
#' @note https://stat.ethz.ch/pipermail/bioconductor/2008-October/024669.html
#' @examples 
#' \dontrun {
#' # TODO include a working example
#' gff <- ape::read.gff("path_to_gff_file")
#' ID <- get_attribute(gff$attributes, "ID")
#' }
get_attribute <- function (x, field, attrsep = ";") {
  s = strsplit(x, split = attrsep, fixed = TRUE)
  sapply(s, function(atts) {
    a = strsplit(atts, split = "=", fixed = TRUE)
    m = match(field, sapply(a, "[", 1))
    if (!is.na(m)) {
      rv = a[[m]][2]
    }
    else {
      rv = as.character(NA)
    }
    return(rv)
  })
}