rm(list=ls())

args <- commandArgs(trailingOnly = TRUE)
jobname <- args[1]
pipeline_version <- args[2]

# load mlst results
mlst <- read.table(
  file = paste0("./jobs/", jobname, "/output/mlst_summary.tsv"),
  stringsAsFactors = FALSE)
mlst <- dplyr::rename(mlst, assembly = V1)
mlst <- dplyr::rename(mlst, mlst = V3)
mlst <- mlst[,c("assembly", "mlst")]

# load kaptive results
load(paste0("./jobs/", jobname, "/output/kaptive/kaptive_results_short.rda"))
kaptive_short <- kaptive_short[,c(
  "assembly",
  "k_serotype",
  "k_confidence",
  "oc_serotype",
  "oc_confidence"
)]

# load diamond results
diamond <- read.csv(
  paste0("./jobs/", jobname, "/output/diamond_names.csv"),
  stringsAsFactors = FALSE)

# merge results
aci <- dplyr::full_join(mlst, kaptive_short, by = "assembly")
aci <- dplyr::full_join(aci, diamond, by = "assembly")
aci <- tibble::as_tibble(aci)

# if metadata exists, add it

metafile <- paste0("./jobs/", jobname, "/metadata/metadata.tsv")

if (file.exists(metafile)) {
  meta <- read.csv(metafile, sep = "\t", stringsAsFactors = FALSE)
  for (i in 1:ncol(meta)) {
    metacol <- meta[, i]
    metacol <- metacol[which(!is.na(metacol))]
    is_date <- mean(grepl("[0-9]{4}-[0-9]{2}-[0-9]{2}", metacol), na.rm = TRUE)
    if (!is.nan(is_date) && is_date == 1) meta[, i] <- as.Date(meta[, i])
  }
  aci <- dplyr::left_join(aci, meta, by = "assembly")
  index1 <- which(names(aci) == "assembly")
  index2 <- which(names(aci) != "assembly")
  aci <- aci[ ,c(index1, index2)]
  if (!is.null(aci$collection_date)) {
    aci$collection_date <- as.character(aci$collection_date)
  }
} else {
  aci$jobname <- jobname
}
aci <- tibble::as_tibble(aci)
attr(aci, "pipeline.version") <- pipeline_version

# export results
save(aci, file = paste0("./jobs/",jobname, "/output/", jobname, "_results.rda"))

output_file <- paste0("./jobs/",jobname, "/output/", jobname, "_results.tsv")
header_message <- paste0("# Pipeline version ", pipeline_version, ". \n")
cat(header_message, file = output_file)
suppressWarnings(write.table(aci,
                             file = output_file,
                             append = TRUE,
                             quote = FALSE,
                             sep = "\t",
                             row.names = FALSE))
