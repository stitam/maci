# What this script does:
# - Plots a phylogenetic tree using the fan layout
# - The tree may be the whole tree or a random subsample
# - Adds three heatmaps: mlst, k_serotype, continent 

library(devtools)
library(dplyr)
library(ggnewscale)
library(ggimage)
library(ggplot2)
library(ggtree)
library(ggtreeExtra)
library(optparse)
library(treeio)

rm(list = ls())

args_list <- list(
  make_option(
    c("-p", "--project_dir"),
    type = "character",
    help = "Path to project directory which contains the pipeline scripts."
  ),
  make_option(
    c("-t", "--tree"),
    type = "character",
    help = "Path to a dated phylogenetic tree in rds format."
  ),
  make_option(
    c("-m", "--metadata"),
    type = "character",
    help = "Path to a metadata file with .rds extension e.g. aci_all.rds."
  ),
  make_option(
    c("-s", "--sample_size"),
    type = "character",
    help = "Sample size, for testing purposes. Will only do sampling if sample size > 0."
  ),
  make_option(
    c("-e", "--sero_over_time"),
    type = "character",
    help = "Path to an rds file for the sero over time plot."
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
    project_dir = "aci",
    tree = "data/global_ST2_tree/dated_tree.rds",
    metadata = "results_redacted/filter_assemblies/aci_filtered.rds",
    sample_size = 50,
    sero_over_time = "results_redacted/plot_sero_over_time/geodate/sero_over_time_collapse_geodate.rds",
    regions = "aci/data/geographic_locations_in_study.tsv"
  )
}

project_dir <- args$project_dir
intree_path <- args$tree
meta_path <- args$metadata
sample_size <- args$sample_size

# load functions
load_all(project_dir)
# load funcions from submodule
load_all(paste0(project_dir, "/prophyl"))

tree <- readRDS(intree_path)
meta <- readRDS(meta_path)

# drop tips which do not have related metadata
if (all(tree$tip.label %in% meta$assembly) == FALSE) {
  labels <- tree$tip.label[which(tree$tip.label %in% meta$assembly == FALSE)]
  msg <- paste(labels, collapse = ", ")
  warning("One or more tip labels cannot be found in metadata table: ", msg)
  # drop tips that are not in the table
  tree <- ape::drop.tip(tree, tip = labels)
}

# only keep metadata for assemblies that are on the tree
meta <- meta[which(meta$assembly %in% tree$tip.label),]

# Capitalise the first letters of continent names
meta$continent <- stringr::str_to_title(meta$continent)

# collapse serotype and define colors

sero_over_time <- readRDS(args$sero_over_time)

seros <- sero_over_time$data %>% 
  group_by(serotype) %>% 
  summarise(count = sum(count)) %>%
  arrange(desc(count))
  
sero_colors <- sero_over_time$plot_env$cls_df

seros <- dplyr::left_join(
  seros,
  sero_colors,
  by = "serotype"
)
seros$serotype <- as.character(seros$serotype)

top_sero_count <- length(unique(meta$serotype))
if (length(unique(meta$serotype)) >= top_sero_count) {
  
  if ("Other" %in% seros$serotype) {
    top_s <- seros$serotype[-which(seros$serotype == "Other")][1:top_sero_count]
  } else {
    top_s <- seros$serotype[1:top_sero_count]
  }
  
  meta$serotype <- gsub("_", " ", meta$serotype)

  if (any(top_s %in% meta$serotype) == FALSE) {
    stop("Serotype ", top_s[which(!top_s %in% meta$serotype)]," not found in metadata table.")
  }

  meta$serotype <- ifelse(meta$serotype %in% top_s, meta$serotype, "Other")
  
  sero_colors <- data.frame(
    serotype = c(sort(top_s), "Other")
  )
  sero_colors <- dplyr::left_join(
    sero_colors,
    seros,
    by = "serotype"
  )

} else {

  sero_colors <- data.frame(
    serotype = c(sort(unique(meta$serotype)), "Other")
  )
  sero_colors <- dplyr::left_join(
    sero_colors,
    seros,
    by = "serotype"
  )
}

# define colors for continents

# read geoloc database
geoloc <- read.csv(args$regions, sep = "\t", na = "")
geoloc_continents <- geoloc %>% 
  filter(variable == "continent" & colour != "#000000") %>%
  dplyr::rename(continent = term, continent_pretty = term_pretty)

continent_colors <- data.frame(
  continent = geoloc_continents$continent_pretty,
  color = geoloc_continents$colour
)

# define colors for countries

# read geoloc database
geoloc <- read.csv(args$regions, sep = "\t", na = "")
geoloc_countries <- geoloc %>% 
  filter(variable == "country") %>%
  dplyr::rename(country = term, country_pretty = term_pretty)

# make country names in meta "pretty"
meta$country <- sapply(meta$country, function(x) {
  index <- which(geoloc_countries$country == x)
  if (length(index) == 0) {
    return(NA)
  } else {
    geoloc_countries$country_pretty[index]
  }
})

geoloc_countries <- geoloc_countries %>% filter(colour != "#000000")

country_colors <- data.frame(
  country = geoloc_countries$country_pretty,
  color = geoloc_countries$colour
)

# subsample tree
if (sample_size > 0) {
  set.seed(0)
  sample_size <- round(as.numeric(sample_size),0)                    
  tree <- ape::keep.tip(tree, sample(tree$tip.label, sample_size))
}

tree_tbl <- as_tibble(tree)

tree_tbl <- dplyr::left_join(
  tree_tbl,
  meta[,c(
    "assembly", "continent", "country", "serotype", "crab"
  )],
  by = c("label" = "assembly")
)

# test that none of the countries are missing
testthat::expect_false(any(is.na(tree_tbl$country[1:ape::Ntip(tree)])))

# update colors
sero_colors <- get_colors(tree_tbl, "serotype", colors = sero_colors, verbose = TRUE)
continent_colors <- get_colors(tree_tbl, "continent", colors = continent_colors, verbose = TRUE)
country_colors <- get_colors(tree_tbl, "country", colors = country_colors, verbose = TRUE)

testthat::expect_true(all(tree_tbl$serotype %in% sero_colors$serotype))
testthat::expect_true(all(tree_tbl$continent %in% continent_colors$continent))
testthat::expect_true(all(tree_tbl$country %in% country_colors$country))

# get a list of countries to show on legend

country_counts <- tree_tbl[1:ape::Ntip(tree),] %>% 
  group_by(country) %>% 
  summarise(count = n()) %>% 
  arrange(desc(count))

index <- which(country_counts$count >= 100)
country_breaks <- country_counts$country[index]

g <- plot_tree_fan(
  tree_tbl,
  linewidth = 0.1,
  open_angle = 25 ,
  heatmap_var = c("continent", "country", "serotype"),
  heatmap_colors = list(
    "serotype" = sero_colors,
    "continent" = continent_colors,
    "country" = country_colors
  ),
  heatmap_colnames_hjust = 1.05,
  legend_show = c("continent", "country"),
  legend_breaks = list(
    "country" = country_breaks
  ),
  legend_guides = list(
    "country" = list(title = "", nrow = 4),
    "continent" = list(title = "", nrow = 4)
  ),
  heatmap_colnames_font_size = 1.8,
  verbose = TRUE,
  heatmap_width = 7
)

ggsave(
  filename = "global_tree.pdf",
  plot = g
)

ggsave(
  filename = "global_tree.png",
  plot = g
)

saveRDS(g, "global_tree.rds")
