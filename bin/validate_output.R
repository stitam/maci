rm(list = ls())

args <- commandArgs(trailingOnly = TRUE)

if (grepl("rds$", args[1])) {
  assemblies <- readRDS(args[1])
}
if (grepl("csv$", args[1])) {
  assemblies <- read.csv(args[1], sep = ",")
}
if (grepl("tsv$", args[1])) {
  assemblies <- read.csv(args[1], sep = "\t")
}

aci <- read.csv(args[2], sep = "\t")

pipeline_version <- args[3]

# Validate output
if (nrow(aci) != nrow(assemblies)) {
  stop("Number of rows in output is different from number of rows in input.")
}

write.table(
  aci,
  file = "validated_results.tsv",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)