library(optparse)
rm(list = ls())

args_list <- list(
  make_option(
    c("-f", "--file"),
    type = "character",
    help = "Path to a tbl prepared for serotop by regions, e.g. top_serotypes_region23_collapse_geodate.tsv."
  ),
  make_option(
    c("-r", "--regions"),
    type = "character",
    help = "Path to a table with region names, coordinates, colors."
  ),
  make_option(
    c("-s", "--serotype_file"),
    type = "character",
    help = "Path to a tbl which lists global or prevalent serotypes."
  )
)
args_parser  <- OptionParser(option_list = args_list)

if (!interactive()) {
  args  <- parse_args(args_parser)
} else {
  args <- list(
    file = "results/calc_serotype_freqs/geodate/TableS2A_top_serotypes_region23_collapse_geodate.tsv",
    regions = "aci/data/geographic_locations_in_study.tsv",
    serotype_file = "results/calc_serotype_freqs/geodate/global_or_prevalent_serotypes_region23_collapse_geodate.tsv"
  )
}

library(tidyverse)

strategy <- args$file %>% strsplit(., "_") %>% unlist()
strategy <- gsub(".tsv", "", strategy[6])

regions_serotypes_geodate <- read_delim(args$file)
regions_serotypes_geodate$serotype <- gsub("_"," ", regions_serotypes_geodate$serotype)

regions_serotypes_geodate_hm <- regions_serotypes_geodate  %>%  
  mutate(regions_serotypes_geodate, serotype2=if_else(ratio_cumsum>0.9, "Other", serotype)) %>% 
  group_by(region23, serotype2) %>% 
  summarise(percent=sum(ratio)) # remove *100 to get ratio instead of percent

global_prevalent <- read.csv(args$serotype_file, sep = "\t")
global_prevalent$serotype <- gsub("_", " ", global_prevalent$serotype)

global_clones <- global_prevalent %>%
  filter(color == "#DC143C") %>%
  select(serotype)

regions_serotypes_geodate_boxplot <- regions_serotypes_geodate_hm[regions_serotypes_geodate_hm$serotype2 %in% global_clones$serotype,]
regions_serotypes_geodate_boxplot <- regions_serotypes_geodate_boxplot %>% pivot_wider(names_from=serotype2, values_from = percent)
regions_serotypes_geodate_boxplot[is.na(regions_serotypes_geodate_boxplot)] <- 0
regions_serotypes_geodate_boxplot <- pivot_longer(regions_serotypes_geodate_boxplot, cols = c(2:ncol(regions_serotypes_geodate_boxplot)), values_to = "percent")

# import region23 tbl with region colors
region_colors <- read.csv(args$regions, sep = "\t", na = "")
region_colors <- region_colors %>% 
  filter(variable == "region23") %>%
  select(term, term_pretty, colour)

region_colors_vector <- region_colors$colour
names(region_colors_vector) <- region_colors$term_pretty

# merge tables to add colors
regions_serotypes_geodate_boxplot <- dplyr::left_join(
  regions_serotypes_geodate_boxplot,
  region_colors,
  by = c("region23" = "term")
)

# check that all levels have colors (missing colors would be excluded!)
testthat::expect_true(
  all(regions_serotypes_geodate_boxplot$term_pretty %in% names(region_colors_vector)))

g <- ggplot(regions_serotypes_geodate_boxplot, aes(reorder(name, -percent), percent))+ 
  geom_boxplot(outlier.color = "transparent")+
  geom_point(position=position_jitter(0.2), aes(color=term_pretty))+
  scale_color_manual(values = region_colors_vector, name = "")+
  geom_hline(yintercept=0.02, colour="red", linetype=2)+
  ylab("Relative prevalence")+
  xlab("")+
  theme_classic()+
  theme(axis.text.x = element_text(angle=60, hjust=1, vjust=1),
        legend.position = "right")

ggsave(
  g,
  file = paste0("FigS3_boxplot_global_serotypes_collapse_", strategy, ".pdf"),
  width = 6.5,
  height = 4.5
)

ggsave(
  g,
  file = paste0("FigS3_boxplot_global_serotypes_collapse_", strategy, ".png"),
  width = 6.5,
  height = 4.5
)

saveRDS(g, file = paste0("FigS3_boxplot_global_serotypes_collapse_", strategy, ".rds"))
