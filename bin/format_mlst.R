rm(list=ls())

args <- commandArgs(trailingOnly = TRUE)
sample_id <- args[1]
mlst <- args[2]

mlst <- read.table(file = mlst)

mlst_small <- mlst[,c(1,3)]
names(mlst_small) <- c("assembly", "mlst")
mlst_small$mlst <- ifelse(mlst_small$mlst == "-", NA, paste0("ST", mlst_small$mlst))

mlst_small$mlst_genes_detected <- ncol(mlst) - 3

write.table(
    mlst_small,
    file = paste0(sample_id, "_formatted.tsv"),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
)
