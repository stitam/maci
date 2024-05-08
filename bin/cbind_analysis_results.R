rm(list=ls())

args <- commandArgs(trailingOnly = TRUE)
assembly_stats <- args[1]
busco <- args[2]
kraken2 <- args[3]
mlst <- args[4]
kaptive <- args[5]
diamond <- args[6]
crab <- args[7]
metafile <- args[8]

myPaths <- .libPaths()
myPaths <- c(myPaths, args[9])

# load tables
assembly_stats <- read.csv(assembly_stats, sep = "\t")
busco <- read.csv(busco, sep = "\t")
kraken2 <- read.csv(kraken2, sep = "\t")
mlst <- read.csv(mlst, sep = "\t")
kaptive <- read.csv(kaptive, sep = "\t")
diamond <- read.csv(diamond, sep = "\t")
crab <- read.csv(crab, sep = "\t")

if (grepl("rds$", metafile)) {
  meta <- readRDS(metafile)
}
if (grepl("csv$", metafile)) {
  meta <- read.csv(metafile, sep = ",")
}
if (grepl("tsv$", metafile)) {
  meta <- read.csv(metafile, sep = "\t")
}

# merge tables
aci <- dplyr::left_join(meta, assembly_stats)

# remove columns
busco <- busco[,which(names(busco) %in% c(
  "assembly", "busco_complete_perc"
))]

aci <- dplyr::left_join(aci, busco)
aci <- dplyr::left_join(aci, kraken2)
aci <- dplyr::left_join(aci, mlst)
aci <- dplyr::left_join(aci, kaptive)
aci <- dplyr::left_join(aci, diamond)
aci <- dplyr::left_join(aci, crab)
aci <- tibble::as_tibble(aci)

write.table(
  aci,
  file = paste0("analysis_results.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
