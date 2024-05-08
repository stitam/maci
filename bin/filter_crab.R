rm(list = ls())

if (!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  aci_path <- args[1]
  collapse_strategy <- args[2]
} else {
  test_dir <- "~/Methods/aci-tests/test-filter_crab"
  aci_path <- paste0(test_dir, "/aci_collapse_geodate.rds")
  collapse_strategy <- "geodate"
}

# import data set
aci <- readRDS(aci_path)

# filter to crab
aci_filtered <- aci[which(aci$crab == TRUE), ]

# export
saveRDS(
  aci_filtered,
  paste0("aci_collapse_", collapse_strategy, "_crab.rds")
)