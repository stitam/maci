# This script checks assemblies based on a number of QC metrics:

# DETERMINISTIC QUALITY FILTERS
# - mlst is not NA and exactly 7 genes detected
# - k_serotype is not NA
# - k_confidence is at least 'Good'
# - kraken2_taxid contains '470'
# - coverage is at least 25
# - BUSCO complete score is at least 95% 
# STOCHASTIC QUALITY FILTERS
# - GC content - use IQR extreme outliers, keep middle
# - contig count - use IQR extreme outliers, keep lower
# - longest contig - use IQR extreme outliers, keep higher
# - N50 - use IQR extreme outliers, keep higher
# - N95 - use IQR extreme outliers, keep higher
# - n_count - use IQR extreme outliers, keep lower, exclude 0
# - kraken2_taxid_freq - use IQR extreme outliers, keep higher
#
# Deterministic quality filters are applied once. Stochastic quality filters
# are applied iteratively until no more rows are removed.
#
# The script does not remove any rows from the data set, but adds a new variable
# called 'qc_pass' that is TRUE if the assembly passes all QC checks and FALSE
# otherwise. The script exports the data set and also a data frame that contains
# the assemblies that failed QC checks and the reason why they failed.

library(optparse)
rm(list = ls())

args_list <- list(
  make_option(
    c("-p", "--project_dir"),
    type = "character",
    help = "Path to project directory which contains the pipeline scripts."
  ),
  make_option(
    c("-f", "--file"),
    type = "character",
    help = "Path to aci prediction results."
  )
)

args_parser  <- OptionParser(option_list = args_list)

if (!interactive()) {
  args  <- parse_args(args_parser)
} else {
  args <- list(
    project_dir = "~/Methods/aci",
    file = "aci_all.rds"
  )
}

# create log file and start logging
con <- file(paste0("log.txt"))
sink(con, split = TRUE)

library(dplyr)
library(devtools)
load_all(args$project_dir)

# import data set
aci <- readRDS(args$file)
aci_failed <- data.frame()

# add new variable for qc
aci$qc_pass <- TRUE

print("Number of rows:")
print(nrow(aci))

# DETERMINISTIC QUALITY FILTERS

# mlst
print("Keep rows where MLST is not NA.")

if (length(which(is.na(aci$mlst))) > 0) {
  aci_rm <- data.frame(
    assembly = aci$assembly[which(is.na(aci$mlst))],
    reason = "MLST is NA"
  )
  aci_failed <- dplyr::bind_rows(
    aci_failed,
    aci_rm
  )
  aci$qc_pass[which(aci$assembly %in% aci_rm$assembly)] <- FALSE
}

print("QC passed overall:")
print(sum(aci$qc_pass))

print("Keep rows where exactly 7 MLST genes are detected.")
if (length(which(aci$mlst_genes_detected != 7)) > 0) {
  aci_rm <- data.frame(
    assembly = aci$assembly[which(aci$mlst_genes_detected != 7)],
    reason = "More or less than 7 MLST genes detected"
  )
  aci_failed <- dplyr::bind_rows(
    aci_failed,
    aci_rm
  )
  aci$qc_pass[which(aci$assembly %in% aci_rm$assembly)] <- FALSE
}

print("QC passed overall:")
print(sum(aci$qc_pass))

# k_serotype

print("Keep rows where K serotype is not NA.")
if (length(which(is.na(aci$k_serotype))) > 0) {
  aci_rm <- data.frame(
    assembly = aci$assembly[which(is.na(aci$k_serotype))],
    reason = "K serotype is NA"
  )
  aci_failed <- dplyr::bind_rows(
    aci_failed,
    aci_rm
  )
  aci$qc_pass[which(aci$assembly %in% aci_rm$assembly)] <- FALSE
}

print("QC passed overall:")
print(sum(aci$qc_pass))

# k_confidence

# print("Keep rows where K confidence is at least 'Low'.")
# if (length(which(aci$k_confidence %in% "None")) > 0) {
#   aci_rm <- data.frame(
#     assembly = aci$assembly[which(aci$k_confidence %in% "None")],
#     reason = "K confidence is less than 'Low'"
#   )
#   aci_failed <- dplyr::bind_rows(
#     aci_failed,
#     aci_rm
#   )
#   aci$qc_pass[which(aci$assembly %in% aci_rm$assembly)] <- FALSE
# }

# print("QC passed overall:")
# print(sum(aci$qc_pass))

print("Keep rows where K confidence is at least 'Good'.")
if (length(which(aci$k_confidence %in% c("None", "Low"))) > 0) {
  aci_rm <- data.frame(
    assembly = aci$assembly[which(aci$k_confidence %in% c("None", "Low"))],
    reason = "K confidence is less than 'Good'"
  )
  aci_failed <- dplyr::bind_rows(
    aci_failed,
    aci_rm
  )
  aci$qc_pass[which(aci$assembly %in% aci_rm$assembly)] <- FALSE
}

print("QC passed overall:")
print(sum(aci$qc_pass))

# kraken2_taxid

print("Keep rows where kraken2_taxid contains '470'.")
index_keep <- strsplit(aci$kraken2_taxid, "\\|") %>%
  sapply(function(x) "470" %in% x) %>%
  which()

if (length(which(!index_keep)) > 0) {
  aci_rm <- data.frame(
    assembly = aci$assembly[which(!index_keep)],
    reason = "kraken2_taxid does not contain '470'"
  )
  aci_failed <- dplyr::bind_rows(
    aci_failed,
    aci_rm
  )
  aci$qc_pass[which(aci$assembly %in% aci_rm$assembly)] <- FALSE
}

print("QC passed overall:")
print(sum(aci$qc_pass))

# coverage

print("Keep rows where coverage is at least 25")

if (length(which(aci$coverage < 25)) > 0) {
  aci_rm <- data.frame(
    assembly = aci$assembly[which(aci$coverage < 25)],
    reason = "Coverage is too low"
  )
  aci_failed <- dplyr::bind_rows(
    aci_failed,
    aci_rm
  )
  aci$qc_pass[which(aci$assembly %in% aci_rm$assembly)] <- FALSE
}

print("QC passed overall:")
print(sum(aci$qc_pass))

# busco

print("Keep rows where BUSCO Complete score is at least 95")

if (length(which(aci$busco_complete < 95)) > 0) {
  aci_rm <- data.frame(
    assembly = aci$assembly[which(aci$busco_complete < 95)],
    reason = "BUSCO Complete score is too low"
  )
  aci_failed <- dplyr::bind_rows(
    aci_failed,
    aci_rm
  )
  aci$qc_pass[which(aci$assembly %in% aci_rm$assembly)] <- FALSE
}

print("QC passed overall:")
print(sum(aci$qc_pass))

# STOCHASTIC QUALITY FILTERS

# IMPORTANT: THE THRESHOLDS BELOW ARE CALCULATED USING ISOLATES THAT HAVE PASSED
# THE DETERMINISTIC QUALITY FILTERS ABOVE. THIS IS TO ENSURE THAT THE THRESHOLDS
# ARE NOT BIASED BY ISOLATES THAT ARE CLEARLY OF POOR QUALITY.

ip0 <- 1
ip1 <- 0
counter <- 0
while (ip1 < ip0 & counter < 10) {
  
  print(paste0("Iteration: ", counter))
  
  ip0 <- sum(aci$qc_pass)
  
  aci_qc <- aci[which(aci$qc_pass), ]

  gc_content <- flag_outliers(
    aci_qc$gc_content,
    qlow = 0.25,
    qhigh = 0.75,
    m = 3,
    keep = "middle"
  )
  contig_count <- flag_outliers(
    aci_qc$contig_count,
    qlow = 0.25,
    qhigh = 0.75,
    m = 3,
    keep = "lower"
  )

longest_contig <- flag_outliers(
  aci_qc$longest_contig,
  qlow = 0.25,
  qhigh = 0.75,
  m = 3,
  keep = "higher"
)

N50 <- flag_outliers(
  aci_qc$N50,
  qlow = 0.25,
  qhigh = 0.75,
  m = 3,
  keep = "higher"
)

N95 <- flag_outliers(
  aci_qc$N95,
  qlow = 0.25,
  qhigh = 0.75,
  m = 3,
  keep = "higher"
)

n_count <- flag_outliers(
  aci_qc$n_count,
  qlow = 0.25,
  qhigh = 0.75,
  m = 3,
  keep = "lower",
  exclude = 0
)

kraken2_taxid_freq <- flag_outliers(
  aci_qc$kraken2_taxid_freq,
  qlow = 0.25,
  qhigh = 0.75,
  m = 3,
  keep = "higher"
)

qc <- data.frame(
  assembly = aci_qc$assembly,
  gc_content = gc_content$keep,
  contig_count = contig_count$keep,
  longest_contig = longest_contig$keep,
  N50 = N50$keep,
  N95 = N95$keep,
  n_count = n_count$keep,
  kraken2_taxid_freq = kraken2_taxid_freq$keep
)

print(apply(qc[,-1], 2, table))

thresholds <- dplyr::bind_rows(
  gc_content = data.frame(
    name = "gc_content",
    qlow = gc_content$range[1],
    qhigh = gc_content$range[2]
  ),
  contig_count = data.frame(
    name = "contig_count",
    qlow = contig_count$range[1],
    qhigh = contig_count$range[2]
  ),
  longest_contig = data.frame(
    name = "longest_contig",
    qlow = longest_contig$range[1],
    qhigh = longest_contig$range[2]
  ),
  N50 = data.frame(
    name = "N50",
    qlow = N50$range[1],
    qhigh = N50$range[2]
  ),
  N95 = data.frame(
    name = "N95",
    qlow = N95$range[1],
    qhigh = N95$range[2]
  ),
  n_count = data.frame(
    name = "n_count",
    qlow = n_count$range[1],
    qhigh = n_count$range[2]
  ),
  kraken2_taxid_freq = data.frame(
    name = "kraken2_taxid_freq",
    qlow = kraken2_taxid_freq$range[1],
    qhigh = kraken2_taxid_freq$range[2]
  )
)

print(thresholds)

# gc_content

print("Keep rows where GC content is not too high or too low")

if (length(which(qc$gc_content == FALSE)) > 0) {
  aci_rm <- data.frame(
    assembly = qc$assembly[which(qc$gc_content == FALSE)],
    reason = "GC content is too high or too low"
  )
  aci_failed <- dplyr::bind_rows(
    aci_failed,
    aci_rm
  )
  aci$qc_pass[which(aci$assembly %in% aci_rm$assembly)] <- FALSE
}

print("QC passed overall:")
print(sum(aci$qc_pass))

# contig_count

print("Keep rows where contig count is not too high")

if (length(which(qc$contig_count == FALSE)) > 0) {
  aci_rm <- data.frame(
    assembly = qc$assembly[which(qc$contig_count == FALSE)],
    reason = "Contig count is too high"
  )
  aci_failed <- dplyr::bind_rows(
    aci_failed,
    aci_rm
  )
  aci$qc_pass[which(aci$assembly %in% aci_rm$assembly)] <- FALSE
}

print("QC passed overall:")
print(sum(aci$qc_pass))

# longest_contig

print("Keep rows where longest contig is not too short")

if (length(which(qc$longest_contig == FALSE)) > 0) {
  aci_rm <- data.frame(
    assembly = qc$assembly[which(qc$longest_contig == FALSE)],
    reason = "Longest contig is too short"
  )
  aci_failed <- dplyr::bind_rows(
    aci_failed,
    aci_rm
  )
  aci$qc_pass[which(aci$assembly %in% aci_rm$assembly)] <- FALSE
}

print("QC passed overall:")
print(sum(aci$qc_pass))

# N50

print("Keep rows where N50 is not too short")

if (length(which(qc$N50 == FALSE)) > 0) {
  aci_rm <- data.frame(
    assembly = qc$assembly[which(qc$N50 == FALSE)],
    reason = "N50 is too short"
  )
  aci_failed <- dplyr::bind_rows(
    aci_failed,
    aci_rm
  )
  aci$qc_pass[which(aci$assembly %in% aci_rm$assembly)] <- FALSE
}

print("QC passed overall:")
print(sum(aci$qc_pass))

# N95

print("Keep rows where N95 is not too short")

if (length(which(qc$N95 == FALSE)) > 0) {
  aci_rm <- data.frame(
    assembly = qc$assembly[which(qc$N95 == FALSE)],
    reason = "N95 is too short"
  )
  aci_failed <- dplyr::bind_rows(
    aci_failed,
    aci_rm
  )
  aci$qc_pass[which(aci$assembly %in% aci_rm$assembly)] <- FALSE
}

print("QC passed overall:")
print(sum(aci$qc_pass))

# n_count

print("Keep rows where N count is not too high")

if (length(which(qc$n_count == FALSE)) > 0) {
  aci_rm <- data.frame(
    assembly = qc$assembly[which(qc$n_count == FALSE)],
    reason = "N count is too high"
  )
  aci_failed <- dplyr::bind_rows(
    aci_failed,
    aci_rm
  )
  aci$qc_pass[which(aci$assembly %in% aci_rm$assembly)] <- FALSE
}

print("QC passed overall:")
print(sum(aci$qc_pass))

# kraken2_taxid_freq

print("Keep rows where kraken2_taxid_freq is not too low")

if (length(which(qc$kraken2_taxid_freq == FALSE)) > 0) {
  aci_rm <- data.frame(
    assembly = qc$assembly[which(qc$kraken2_taxid_freq == FALSE)],
    reason = "kraken2_taxid_freq is too low"
  )
  aci_failed <- dplyr::bind_rows(
    aci_failed,
    aci_rm
  )
  aci$qc_pass[which(aci$assembly %in% aci_rm$assembly)] <- FALSE
}

print("QC passed overall:")
print(sum(aci$qc_pass))

ip1 <- sum(aci$qc_pass)
counter <- counter + 1

}




# export data frame with QC results

saveRDS(aci, file = "aci_with_qc.rds")

write.table(
  aci,
  file = "aci_with_qc.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# export data frame with QC failures

aci_failed <- aci_failed[order(aci_failed$assembly), ]

saveRDS(aci_failed, file = "aci_qc_failed.rds")

write.table(
  aci_failed,
  file = "aci_qc_failed.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# stop logging
sink(con)
