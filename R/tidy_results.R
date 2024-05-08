#' Combine prediction results with metadata
#' 
#' This function is the final step in an analysis job. It takes the metadata
#' file, validates its content, takes the prediction results file, merges the
#' two and exports the combined data set in two formats with proper names.
#' @param dir character; path to the directory which contains the metadata and
#' the prediction results in standardized format.
#' @return combined data set in .rds and .tsv format.
#' @examples 
#' \dontrun{
#' tidy_results(<path to the directory where analysis data is stored>)
#' }
#' @export
tidy_results <- function(dir) {
  jobname <- basename(dir)
  meta <- readRDS(paste0(jobname, "/metadata/metadata.rds"))
  meta <- validate_metadata(meta)
  results <- readRDS(paste0(jobname, "/results/results_short.rds"))
  out <- dplyr::left_join(results, meta, by = "assembly")
  saveRDS(out, file = paste0(dir, "/", jobname, ".rds"))
  write.table(
    out,
    file = paste0(dir, "/", jobname, ".tsv"),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
}
