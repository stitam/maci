library(optparse)
rm(list = ls())

args_list <- list(
  make_option(
    c("-R", "--serotop_region23"),
    type = "character",
    help = "Path to a tbl prepared for serotop, e.g. serotop_region23_ds_geodate.tsv."
  ),
  make_option(
    c("-C", "--serotop_country"),
    type = "character",
    help = "Path to a tbl prepared for serotop, e.g. serotop_country_ds_geodate.tsv."
  ),
  make_option(
    c("-r", "--regions"),
    type = "character",
    help = "Path to a table with region names, coordinates, colors."
  )
)

args_parser  <- OptionParser(option_list = args_list)

if (!interactive()) {
  args  <- parse_args(args_parser)
} else {
  args <- list(
    serotop_region23 = "results_redacted_new/prep_serotop_tbls/serotop_region23_ds_geodate.tsv",
    serotop_country = "results_redacted_new/prep_serotop_tbls/serotop_country_ds_geodate.tsv",
    regions = "data/geographic_locations_in_study.tsv"
  )
}

library(tidyverse)

# extract strategy for filename
strategy <- strsplit(basename(args$serotop_region23), split = "_")[[1]][4]
strategy <- gsub("\\.tsv$", "", strategy)

# Read data frame for geographic regions, shortlist, etc.
geoloc <- read.csv(args$regions, sep = "\t", na = "")

# Country
country_geodate <- read_delim(args$serotop_country)

names(country_geodate) <- c("country", "quantile", "serotop")

# Only include shortlisted countries

shortlist <- geoloc %>% filter(variable == "country" & shortlisted == TRUE)

# data frame "with" the "All" category
# df <- data.frame(
#   country = c(shortlist$term, "all"),
#   country_pretty = c(shortlist$term_pretty, "All"),
#   colour = c(shortlist$colour, "black"),
#   linetype = c(rep("solid", times = nrow(shortlist)), "dashed")
# )

# data frame "without" to "All" category
df <- data.frame(
  country = shortlist$term,
  country_pretty = shortlist$term_pretty,
  colour = shortlist$colour,
  linetype = rep("solid", times = nrow(shortlist))
)

testthat::expect_true(all(df$country_pretty %in% country_geodate$country))

country_geodate <- country_geodate[which(country_geodate$country %in% df$country_pretty),]
country_geodate <- dplyr::left_join(
  country_geodate,
  df,
  by = "country"
)

# country_geodate %>% group_by(country) %>% summarise(x=sum(serotop)) %>% filter(x>40)

country_colors <- df$colour
names(country_colors) <- df$country_pretty

country_linetypes <- df$linetype
names(country_linetypes) <- df$country_pretty

country_geodate_fig <- ggplot(country_geodate, aes(quantile, serotop)) +
  geom_line(aes(col=country_pretty, linetype = country))+
  geom_point(aes(col=country))+
  theme_classic()+
  ylab("Number of MLST-CPS types")+
  xlab("Proportion of genomes")+
  scale_y_log10()+
  scale_x_continuous(breaks=c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1))+
  scale_colour_manual(values = country_colors, name="Country")+
  scale_linetype_manual(values = country_linetypes, name = "Country")+
  theme(panel.grid.major.y=element_line(color="gray", linetype = 3, linewidth = 0.3),
        axis.line=element_line(color="black", linewidth=0.1),
        axis.ticks = element_line(linewidth=0.1),
        legend.key.size = unit(1,"point"))+
  annotation_logticks(sides="l", size=0.1)

ggsave(
  filename = paste0("serotop_country_ds_", strategy, ".pdf"),
  width = 7,
  height = 4
)
ggsave(
  filename = paste0("serotop_country_ds_", strategy, ".png"),
  width = 7,
  height = 4
)
saveRDS(country_geodate_fig, file = paste0("serotop_country_ds_", strategy, ".rds"))

# Region

region23_geodate <- read_delim(args$serotop_region23)

names(region23_geodate) <- c("region23", "quantile", "serotop")

# Only include shortlisted countries

shortlist <- geoloc %>% filter(variable == "region23" & shortlisted == TRUE)

# data frame "with" the "All" category
# df <- data.frame(
#   region23 = c(shortlist$term, "all"),
#   region23_pretty = c(shortlist$term_pretty, "All"),
#   colour = c(shortlist$colour, "black"),
#   linetype = c(rep("solid", times = nrow(shortlist)), "dashed")
# )

# data frame "without" to "All" category
df <- data.frame(
  region23 = shortlist$term,
  region23_pretty = shortlist$term_pretty,
  colour = shortlist$colour,
  linetype = rep("solid", times = nrow(shortlist))
)

testthat::expect_true(all(df$region23_pretty %in% region23_geodate$region23))

region23_geodate <- region23_geodate[which(region23_geodate$region23 %in% df$region23_pretty),]
region23_geodate <- dplyr::left_join(
  region23_geodate,
  df,
  by = "region23"
)

region23_colors <- df$colour
names(region23_colors) <- df$region23_pretty

region23_linetypes <- df$linetype
names(region23_linetypes) <- df$region23_pretty

region23_geodate_fig <- ggplot(region23_geodate, aes(quantile, serotop)) +
  geom_line(aes(col=region23, linetype = region23), linewidth = 0.1)+
  geom_point(aes(col=region23), size = 0.5)+
  theme_classic()+
  ylab("Number of MLST-CPS types")+
  xlab("Proportion of genomes")+
  scale_y_log10()+
  scale_x_continuous(breaks=c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1))+
  scale_colour_manual(values = region23_colors, name="")+
  scale_linetype_manual(values = region23_linetypes, name = "")+
  guides(
    colour = guide_legend(keyheight = unit(2, "mm")),
    linetype = guide_legend(keyheight = unit(2, "mm"))
  )+
  theme(
    axis.title = element_text(size = 5),
    axis.title.x = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0)),
    axis.text = element_text(size = 5),
    axis.line=element_line(color="black", linewidth=0.1),
    legend.title = element_text(size = 5),
    legend.text = element_text(size = 5),
    panel.grid.major.y=element_line(color="gray", linetype = 3, linewidth = 0.1),
    axis.ticks = element_line(linewidth=0.1),
    plot.margin = unit(c(0,0,0,0), "cm"),
    legend.margin = margin(t = 0, r = 0, b = 0, l = 0),
    panel.background = element_rect(fill='transparent'), #transparent panel bg
    plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
    legend.background = element_rect(fill='transparent', color=NA)
  )+
  annotation_logticks(sides="l", size=0.1)

ggsave(
  filename = paste0("serotop_region23_ds_", strategy, ".pdf"),
  units = "cm",
  width = 10,
  height = 6
)
ggsave(
  filename = paste0("serotop_region23_ds_", strategy, ".png"),
  units = "cm",
  width = 10,
  height = 6,
  bg = 
)
saveRDS(region23_geodate_fig, file = paste0("serotop_region23_ds_", strategy, ".rds"))
