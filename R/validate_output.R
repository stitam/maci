# NOT YET SURE THIS IS NEEDED
validate_output <- function(predictions, metadata) {
  pred_ids <- unique(predictions$assembly)
  meta_ids <- unique(metadata$assembly)
  if (sum(pred_ids %in% meta_ids == FALSE) > 0) {
    index <- which(pred_ids %in% meta_ids == FALSE)
    msg <- paste(
      "Metadata missing for the following entries: ",
      pred_ids[index],
      collapse = " "
    )
    warning(msg)
  }
  if (sum(meta_ids %in% pred_ids == FALSE) > 0) {
    index <- which(meta_ids %in% pred_ids == FALSE)
    msg <- paste(
      "Predictions missing for the following entries: ",
      pred_ids[index],
      collapse = " "
    )
    warning(msg)
  }
}