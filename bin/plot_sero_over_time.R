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
  ),
  make_option(
    c("-s", "--serotype_file"),
    type = "character",
    help = "Path to a tbl which lists global or prevalent serotypes."
  ),
  make_option(
    c("-m", "--minyear"),
    type = "character",
    help = "Minimum year to plot."
  ),
  make_option(
    "--maxyear",
    type = "character",
    help = "Maximum year to plot."
  )
)

args_parser  <- OptionParser(option_list = args_list)

if (!interactive()) {
  args  <- parse_args(args_parser)
} else {
  args <- list(
    project_dir = "aci",
    file = "results_redacted/downsample_isolates/aci_crab_ds_geodate.tsv",
    serotype_file = "results_redacted/calc_serotype_freqs/geodate/global_or_prevalent_serotypes_region23_ds_geodate.tsv",
    minyear = "2009",
    maxyear = "2020"
  )
}

library(dplyr)

downsampling_strategy <- args$file %>% basename() %>% strsplit(split = "_") %>% unlist()
downsampling_strategy <- gsub("\\.tsv", "", downsampling_strategy[4])

if (!interactive()) {
  # create log file and start logging
  con <- file(paste0("log_", downsampling_strategy, ".txt"))
  sink(con, split = TRUE)
}

library(devtools)
load_all(args$project_dir)

# import 
aci <- read_df(args$file) %>% 
  dplyr::filter(filtered & crab & downsampled & downsampled_by_pop) %>%
  dplyr::filter(collection_year <= as.numeric(args$maxyear))

aci$serotype <- gsub("-", " ", aci$serotype)

global_prevalent <- read.csv(args$serotype_file, sep = "\t")

global_prevalent$serotype <- gsub("-", " ", global_prevalent$serotype)

keep <- sort(unique(aci$serotype[which(aci$serotype %in% global_prevalent$serotype)]))

aci$serotype <- ifelse(aci$serotype %in% keep, aci$serotype, "Other")
aci$serotype <- factor(aci$serotype, levels = c(keep, "Other"))

# check that plotted serotypes are either global or prevalent or pooled as "Other"
testthat::expect_true(all(aci$serotype %in% c(global_prevalent$serotype, "Other")))

g <- plot_ts_area(
  aci,
  var = "serotype",
  timevar = "collection_year",
  minyear = as.numeric(args$minyear),
  window = 2,
  show_count = FALSE
) + 
guides(fill=guide_legend(title="serotype", ncol = 1)) + 
ylab("Relative prevalence")

ggsave(
  filename = paste0("sero_over_time_ds_", downsampling_strategy, ".pdf"),
  plot = g,
  width = 10,
  height= 8
)
ggsave(
  filename = paste0("sero_over_time_ds_", downsampling_strategy, ".png"),
  plot = g,
  width = 10,
  height= 8
)

g <- g +
  guides(fill=guide_legend(title = "serotype", nrow = 3)) + 
  theme(
    panel.background = element_rect(
      colour = "#000000",
      linewidth = 0.1
    ),
    axis.title = element_text(size = 5),
    axis.text = element_text(size = 5),
    axis.ticks = element_line(colour = "#000000", linewidth = 0.1),
    legend.title = element_text(size = 5),
    legend.text = element_text(size = 5),
    legend.key.size = unit(2, "mm"),
    legend.margin = margin(t = 0, r = 0, b = 0, l = 0)
  )

saveRDS(g, file = paste0("sero_over_time_ds_", downsampling_strategy, ".rds"))

# NO LONGER NEEDED

# KEEP THE TOP N FROM EACH YEAR
# aci <- read_df(args$file) %>% dplyr::filter(filtered & crab & downsampled & downsampled_by_pop)

# keep2 <- vector()
# k <- 3
# for (i in unique(aci$collection_year)) {
#   aci_crop <- aci[which(aci$collection_year == i),]
#   if (length(unique(aci_crop$serotype)) >= k) {
#     keep2 <- c(
#       keep2,
#       names(sort(table(aci_crop$serotype), decreasing = TRUE))[1:k]
#     )
#   } else {
#     keep2 <- c(keep2, unique(aci_crop$serotype2))
#   }
# }
# keep2 <- sort(unique(keep2))

# keep3 <- sort(unique(c(keep, keep2)))

# aci$serotype <- ifelse(aci$serotype %in% keep3, aci$serotype, "Other")
# aci$serotype <- factor(aci$serotype, levels = c(keep3, "Other"))

# g <- plot_ts(
#   aci,
#   var = "serotype",
#   timevar = "collection_year",
#   minyear = as.numeric(args$minyear),
#   type = "ratio",
#   geom = "geom_area",
#   show_count = FALSE
#   ) + 
#   theme(
#     legend.key.size = unit(0.4, "cm")
#   )

# ggsave(
#   filename = paste0(
#     "sero_over_time_ds_", downsampling_strategy, "_top_", k, "_per_year.pdf"),
#   plot = g,
#   width = 10,
#   height= 8
# )
# ggsave(
#   filename = paste0(
#     "sero_over_time_ds_", downsampling_strategy, "_top_", k, "_per_year.png"),
#   plot = g,
#   width = 10,
#   height= 8
# )
# saveRDS(g, file = paste0(
#   "sero_over_time_ds_", downsampling_strategy, "_top_", k, "_per_year.rds"))