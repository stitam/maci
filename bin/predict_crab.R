library(dplyr)
rm(list=ls())

args <- commandArgs(trailingOnly = TRUE)
sample_id <- args[1]
resfinder <- read.csv(args[2], sep = "\t", header = TRUE)
carbapenem <- read.csv(args[3], sep = "\t", header = TRUE)
fq <- read.csv(args[4], sep = "\t", header = TRUE)

## INDICATORS OF CARBAPENEM RESISTANCE MUST BE IMPORTED HERE

# Local definition of CRAB:

df <- data.frame(
  assembly = sample_id,
  # Contains at least one aminoglycoside resistance gene
  aminoglycoside = resfinder$aminoglycoside,
  # Contains at least one carbapenem resistance indicator
  carO = carbapenem$carO,
  ISAba1 = carbapenem$ISAba1,
  blaOXA51like = carbapenem$blaOXA51like,
  carbapenem = carbapenem$carbapenem,
  # Contains at least one fluoroquinolone resitance indicator
  fluoroquinolone = fq$fluoroquinolone
)
df$crab = complete.cases(df[,c("aminoglycoside", "carbapenem", "fluoroquinolone")])

write.table(
  df,
  file = paste0(sample_id, "_crab.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
