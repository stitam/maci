library(optparse)
rm(list = ls())

args_list <- list(
  make_option(
    c("-f", "--region_file"),
    type = "character",
    help = "Path to a tbl prepared for serotop by regions, e.g. top_serotypes_region23_ds_geodate.tsv."
  ),
  make_option(
    c("-e", "--sero_over_time"),
    type = "character",
    help = "Path to an rds file for the sero over time plot."
  ),
  make_option(
    c("-r", "--minyear_recent"),
    type = "integer",
    help = "Lowest year to include for recent time period."
  )
)

args_parser  <- OptionParser(option_list = args_list)

if (!interactive()) {
  args  <- parse_args(args_parser)
} else {
  args <- list(
    region_file = "results_redacted/calc_serotype_freqs/geodate/top_serotypes_region23_ds_geodate.tsv",
    sero_over_time = "results_redacted/plot_sero_over_time/geodate/sero_over_time_ds_geodate.rds",
    minyear_recent = 2016
  )
}

library(dplyr)
library(ggplot2)
library(ggrepel)

# import data set
regions <- read.csv(args$region_file, sep = "\t")

regions$serotype <- gsub("-", " ", regions$serotype)

# get overall number of serotypes
all_sero_count <- length(unique(regions$serotype))

# get list of prevalent serotype and their colors
sero_over_time <- readRDS(args$sero_over_time)
data <- sero_over_time$data %>%
  filter(collection_year >= as.numeric(args$minyear_recent))

# get number of prevalent serotypes
prev_sero_count <- length(unique(data$serotype))
if ("Other" %in% data$serotype) {
  prev_sero_count <- prev_sero_count - 1
}

# get number of serotypes covered by the "Other" category
non_prev_sero_count <- all_sero_count - prev_sero_count

# rename non-prevalent serotypes to "Other"
df <- regions
df$serotype <- ifelse(
  df$serotype %in% as.character(data$serotype),
  df$serotype,
  "Other"
)

df <- df %>%
  group_by(serotype) %>%
  summarise(count = sum(count))

plotdf <- dplyr::left_join(
  df,
  sero_over_time$plot_env$cls_df,
  by = "serotype"
)

# Rename "Other" to include the number of serotypes
if ("Other" %in% plotdf$serotype) {
  plotdf$serotype <- as.character(plotdf$serotype)
  index <- which(plotdf$serotype == "Other")
  plotdf$serotype[index] <- paste0(
    "Other (",
    non_prev_sero_count,
    ")"
  )
  plotdf$serotype <- factor(plotdf$serotype, levels = c(
    plotdf$serotype[-index], plotdf$serotype[index]
  ))
  plotdf <- dplyr::bind_rows(
    plotdf[-index, ],
    plotdf[index, ]
  )
}

plotdf$percentage <- round(plotdf$count/sum(plotdf$count)*100, 1)

sero_colors <- plotdf$color
names(sero_colors) <- plotdf$serotype

# Get the positions for printing percentage
df2 <- plotdf %>% 
  mutate(
    csum = rev(cumsum(rev(percentage))), 
    pos = percentage/2 + lead(csum, 1),
    pos = if_else(is.na(pos), percentage/2, pos)
  )

# Only keep percentage for "Other" category
index <- grep("Other", df2$serotype)
df2$pos[-index] <- NA

# Prepary plot
g <- ggplot(plotdf, aes(x="", y=percentage, fill=serotype)) +
  geom_bar(stat="identity", width=1, color = "white") +
  coord_polar("y", start=0) +
  geom_label_repel(
    data = df2,
    aes(y = pos, label = paste0(percentage, "%")),
    size = 4.5, nudge_x = 1, show.legend = FALSE
  ) +
  scale_fill_manual(values = sero_colors) +
  guides(fill = guide_legend(title = "MLST-CPS")) +
  theme_void() +
  theme(
    plot.margin = unit(c(0,0,0,0),"cm")
  )
  
ggsave(
  filename = "prevalence_piechart.pdf",
  plot = g
)

ggsave(
  filename = "prevalence_piechart.png",
  plot = g
)

saveRDS(g, "prevalence_piechart.rds")
